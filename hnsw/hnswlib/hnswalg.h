#pragma once

#include "visited_list_pool.h"
#include "hnswlib.h"
#include <random>
#include <stdlib.h>
#include <unordered_set>
#include <list>
#include <cmath>

//#define max_level_ 2
#define size_n 10000
#define L 256
#define OFF 3 //from 128bits
#define min_book 4
#define cen 16
#define nnum 4 //from 128 to 32
#define KK 101
#define fan 10
//#define delta 0.5

struct elem{
	unsigned int id;
	float dist;
	bool flag;
};

namespace hnswlib {
    typedef unsigned int tableint;
    typedef unsigned int linklistsizeint;
	//typedef unsigned int labeltype;

    template<typename dist_t>
    class HierarchicalNSW : public AlgorithmInterface<dist_t> {
    public:
        HierarchicalNSW(SpaceInterface<dist_t> *s) {

        }

        HierarchicalNSW(SpaceInterface<dist_t> *s, const std::string &location, const std::string &location2,  bool nmslib = false, size_t max_elements=0) {
            loadIndex(location, location2, s, max_elements);
        }

        HierarchicalNSW(int level_, SpaceInterface<dist_t> *s, size_t max_elements, int vecdim_, size_t M = 16, size_t ef_construction = 500, size_t random_seed = 100) :
                link_list_locks_(max_elements), element_levels_(max_elements) {
            max_elements_ = max_elements;

			vecdim = vecdim_;
			max_level_ = level_;
            has_deletions_=false;
            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();
			dist_func_param2_ = s->get_dist_func_param2();
			
            M_ = M;
            maxM_ = M_;
            maxM0_ = M_ * 2;
			maxQ0_ = 32;
			
            ef_construction_ = std::max(ef_construction,M_);
            ef_ = 10;

            level_generator_.seed(random_seed);

            size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
			size_quan_level0_ = maxQ0_ * sizeof(tableint) + sizeof(linklistsizeint);
			
            size_data_per_element_ = size_links_level0_ + sizeof(float) + data_size_ + sizeof(labeltype);
            //offsetData_ = size_links_level0_;
            //label_offset_ = size_links_level0_ + data_size_;
            offsetLevel0_ = 0;
			
			offsetNorm_= size_links_level0_;
            offsetData_ = sizeof(float) + size_links_level0_;
            label_offset_ = size_links_level0_ + sizeof(float) + data_size_;
			
			
			
			//------------------------------------------------------------------------------
			length = new int[max_level_];
			quan_size_ = new size_t[max_level_];
            size_quan_per_element_ = new size_t[max_level_];
			quan_offset_ = new size_t[max_level_];
			quan_level0_memory_ = new char* [max_level_];
			
			for (int i = 0; i < max_level_; i++)
			    length[i] = pow(2, max_level_ -i + OFF);
			

            data_level0_memory_ = (char *) malloc(max_elements_ * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            cur_element_count = 0;

            visited_list_pool_ = new VisitedListPool(1, max_elements);



            //initializations for special treatment of the first node
            enterpoint_node_ = -1;
            maxlevel_ = -1;

            linkLists_ = (char **) malloc(sizeof(void *) * max_elements_);
            if (linkLists_ == nullptr)
                throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate linklists");
            size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
            mult_ = 1 / log(1.0 * M_);
            revSize_ = 1.0 / mult_;
			count_tmp = new int[max_level_];
			for(int i = 0; i < max_level_; i++) {count_tmp[i] = 0;}
			count0 = 0;
        }

    struct Neighbor {
        unsigned id;
        float distance;
        bool flag;
        bool flag2;
	
        Neighbor() = default;
        Neighbor(unsigned id, float distance, bool f, bool f2) : id{id}, distance{distance}, flag(f), flag2(f2) {}

        inline bool operator<(const Neighbor &other) const {
            return distance < other.distance;
        }
    };

    static inline int InsertIntoPool (Neighbor *addr, unsigned K, Neighbor nn) {
        // find the location to insert
        int left=0,right=K-1;
        if(addr[left].distance>nn.distance){
            memmove((char *)&addr[left+1], &addr[left],K * sizeof(Neighbor));
            addr[left] = nn;
            return left;
        }
        if(addr[right].distance<nn.distance){
            addr[K] = nn;
            return K;
        }
        while(left<right-1){
            int mid=(left+right)/2;
            if(addr[mid].distance>nn.distance)right=mid;
            else left=mid;
        }
        //check equal ID

        while (left > 0){
            if (addr[left].distance < nn.distance) break;
            if (addr[left].id == nn.id) return K + 1;
            left--;
        }
        if(addr[left].id == nn.id||addr[right].id==nn.id)return K+1;
            memmove((char *)&addr[right+1], &addr[right],(K-right) * sizeof(Neighbor));
            addr[right]=nn;
            return right;
    }

        struct CompareByFirst {
            constexpr bool operator()(std::pair<dist_t, tableint> const &a,
                                      std::pair<dist_t, tableint> const &b) const noexcept {
                return a.first < b.first;
            }
        };

        ~HierarchicalNSW() {    

            free(data_level0_memory_);
            for (tableint i = 0; i < cur_element_count; i++) {
                if (element_levels_[i] > 0)
                    free(linkLists_[i]);
            }
            free(linkLists_);
            delete visited_list_pool_;
        }

        size_t max_elements_;
        size_t cur_element_count;
        size_t size_data_per_element_;
        size_t size_links_per_element_;

		int max_level_;
        size_t M_;
        size_t maxM_;
        size_t maxM0_;
		size_t maxQ0_;
		
        size_t ef_construction_;

        double mult_, revSize_;
        int maxlevel_;

		//--------------------------------------------
		int* length;
		int vecdim;
		size_t*	quan_size_;
        size_t* size_quan_per_element_;
		size_t*	quan_offset_ ;
		char**	quan_level0_memory_;
        int* count_tmp;	
        int count0;		
		//-------------------------------------------

        VisitedListPool *visited_list_pool_;
        std::mutex cur_element_count_guard_;

        std::vector<std::mutex> link_list_locks_;
        tableint enterpoint_node_;


        size_t size_links_level0_;
		size_t size_quan_level0_;
        size_t offsetNorm_, offsetData_, offsetLevel0_;


        char *data_level0_memory_;
        char **linkLists_;
        std::vector<int> element_levels_;

        size_t data_size_;

        bool has_deletions_;


        size_t label_offset_;
        DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_;
		void *dist_func_param2_;
		
        std::unordered_map<labeltype, tableint> label_lookup_;

        std::default_random_engine level_generator_;

        inline labeltype getExternalLabel(tableint internal_id, int level) const {
            labeltype return_label;
			if(level == -1)
                memcpy(&return_label,(data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), sizeof(labeltype));
			else{
				memcpy(&return_label,(quan_level0_memory_[level] + internal_id * size_quan_per_element_[level] + quan_offset_[level]), sizeof(labeltype));
			}
            return return_label;
        }

        inline void setExternalLabel(tableint internal_id, labeltype label, int level) const {
			if(level == -1)
                memcpy((data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), &label, sizeof(labeltype));
			else{
				memcpy((quan_level0_memory_[level] + internal_id * size_quan_per_element_[level] + quan_offset_[level]), &label, sizeof(labeltype));
			}
        }

        inline labeltype *getExternalLabeLp(tableint internal_id, int level) const {
			if(level == -1)
                return (labeltype *) (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_);
			else {return (labeltype *) (quan_level0_memory_[level] + internal_id * size_quan_per_element_[level] + quan_offset_[level]);}
        }

        inline char *getDataByInternalId(tableint internal_id, int level) const {
			if(level == -1)
                return (data_level0_memory_ + internal_id * size_data_per_element_ + offsetData_);
			else {return (quan_level0_memory_[level] + internal_id * size_quan_per_element_[level]);}
        }
		
        inline float getNormByInternalId(tableint internal_id) const {
		    return  *( (float *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetNorm_) );
        }	

        inline char *getNormByInternalId2(tableint internal_id) const {
		    return (data_level0_memory_ + internal_id * size_data_per_element_ + offsetNorm_);
        }		

        int getRandomLevel(double reverse_size) {
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            double r = -log(distribution(level_generator_)) * reverse_size;
            return (int) r;
        }

        float compare_(const float* a, const float* b, unsigned size) const {
			float result = 0;
			float diff0, diff1, diff2, diff3;
            const float* last = a + size;
            const float* unroll_group = last - 3;

            /* Process 4 items with each loop for efficiency. */
            while (a < unroll_group) {
                diff0 = a[0] - b[0];
                diff1 = a[1] - b[1];
                diff2 = a[2] - b[2];
                diff3 = a[3] - b[3];
                result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
                a += 4;
                b += 4;
            }
            /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
            while (a < last) {
                diff0 = *a++ - *b++;
                result += diff0 * diff0;
            }
            return result;
        }
		
		void rotation_(int vecsize, int dim_, float** data, float* data2, float* R){
			
                    #pragma omp parallel for		
			for(int i = 0; i < vecsize; i++){
                            float* data3 = new float[dim_]; 
		        for(int j = 0; j < dim_; j++){
			        data3[j] = 0;
			        //for(int l = 0; l < dim_; l++){
			            //data2[j] += R[j * dim_ + l] * data[i][l];	
						data3[j] = fstdistfunc_( (const void*) (R + j * dim_), (const void*) (data[i]), dist_func_param2_);	
			        //}
		        }
		
		        for(int l = 0; l < dim_; l++){
			        data[i][l] = data3[l];
		        }
				delete[] data3;
	        }			
		}
			

        void restore_index(float *query_data, float* array0, float* book, 
            float** start_book, unsigned char*** merge, unsigned char** merge0, int* length, int* sdim, int tol_dim, float* temp_arr){
	        
		//printf("restore1\n");	
            float* temp_arr2 = new float[sdim[0]];

            int base_dim = sdim[0];
		
	        float* ind2 = array0 ;
			
	        for(int j = 0; j < base_dim; j++){
		        temp_arr2[j] = fstdistfunc_( (const void*) ind2, (const void*) query_data, dist_func_param_);
                ind2 += tol_dim;			
	        }   
 
	        for(int j = 0; j < L; j++){
		        _mm_prefetch(ind2, _MM_HINT_T0);

		        book[j] = compare_(ind2, temp_arr2, base_dim);

	            ind2 += base_dim;		
	        }
           	
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>		
        searchBaseLayer(tableint ep_id, const void *data_point, float temp_norm, int layer, int level, float*** book) {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

            dist_t lowerBound;
			dist_t dist;
			dist_t dist1;
            //if (!isMarkedDeleted(ep_id)) {
				if(level == -1)
                    dist = temp_norm + getNormByInternalId(ep_id) - 2 * fstdistfunc_(data_point, getDataByInternalId(ep_id, -1), dist_func_param_);
				else{dist = quandistfunc_((unsigned char*) data_point, (unsigned char*) getDataByInternalId(ep_id, level), length[level], book);}
				
                top_candidates.emplace(dist, ep_id);
                lowerBound = dist;
                candidateSet.emplace(-dist, ep_id);

            visited_array[ep_id] = visited_array_tag;

            while (!candidateSet.empty()) {
                std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
                if ((-curr_el_pair.first) > lowerBound) {
                    break;
                }
                candidateSet.pop();

                tableint curNodeNum = curr_el_pair.second;

                std::unique_lock <std::mutex> lock(link_list_locks_[curNodeNum]);

                int *data;// = (int *)(linkList0_ + curNodeNum * size_links_per_element0_);
                if (layer == 0) {
                    data = (int*)get_linklist0(curNodeNum, level);
                } else {
                    data = (int*)get_linklist(curNodeNum, layer);
                }
                size_t size = getListCount((linklistsizeint*)data);
                tableint *datal = (tableint *) (data + 1);
				
                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = *(datal + j);

                    if (visited_array[candidate_id] == visited_array_tag) continue;
                    visited_array[candidate_id] = visited_array_tag;
                    //char *currObj1 = (getDataByInternalId(candidate_id, level));

					if(level == -1)
                        dist1 = temp_norm + getNormByInternalId(candidate_id) - 2 * fstdistfunc_(data_point, getDataByInternalId(candidate_id, level), dist_func_param_);
					else{dist1 = quandistfunc_((unsigned char*) data_point, (unsigned char*) getDataByInternalId(candidate_id, level), length[level], book);}
					
                    if (top_candidates.size() < ef_construction_ || lowerBound > dist1) {
                        candidateSet.emplace(-dist1, candidate_id);

                            top_candidates.emplace(dist1, candidate_id);

                        if (top_candidates.size() > ef_construction_)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }
            visited_list_pool_->releaseVisitedList(vl);

            return top_candidates;
        }

        template <bool has_deletions>
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerST(tableint ep_id, const void *data_point, size_t ef) const {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
			dist_t dist;
            dist = fstdistfunc_(data_point, getDataByInternalId(ep_id, -1), dist_func_param_);

            lowerBound = dist;
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);

            visited_array[ep_id] = visited_array_tag;

            while (!candidate_set.empty()) {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound) {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;
                int *data = (int *) get_linklist0(current_node_id, -1);
                size_t size = getListCount((linklistsizeint*)data);

                for (size_t j = 1; j <= size; j++) {
                    int candidate_id = *(data + j);

                    if (!(visited_array[candidate_id] == visited_array_tag)) {

                        visited_array[candidate_id] = visited_array_tag;

                        char *currObj1 = (getDataByInternalId(candidate_id, -1));
						
                        dist = fstdistfunc_(data_point, currObj1, dist_func_param_);

                        if (top_candidates.size() < ef || lowerBound > dist) {
                            candidate_set.emplace(-dist, candidate_id);

                                top_candidates.emplace(dist, candidate_id);

                            if (top_candidates.size() > ef)
                                top_candidates.pop();

                            if (!top_candidates.empty())
                                lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            return top_candidates;
        }

        void getNeighborsByHeuristic2(
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates,
                const size_t M, int level, float*** book) {
            if (top_candidates.size() < M) {
                return;
            }
            std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
            std::vector<std::pair<dist_t, tableint>> return_list;
            while (top_candidates.size() > 0) {
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
                top_candidates.pop();
            }
            dist_t curdist;
            while (queue_closest.size()) {
                if (return_list.size() >= M)
                    break;
                std::pair<dist_t, tableint> curent_pair = queue_closest.top();
                dist_t dist_to_query = -curent_pair.first;
                queue_closest.pop();
                bool good = true;
                for (std::pair<dist_t, tableint> second_pair : return_list) {
					if(level == -1)
                    curdist =
                            getNormByInternalId(second_pair.second) + getNormByInternalId(curent_pair.second) - 2 *
							             fstdistfunc_(getDataByInternalId(second_pair.second, level),
                                         getDataByInternalId(curent_pair.second, level),
                                         dist_func_param_); //two 
					else{
					curdist =
                            quandistfunc_((unsigned char*) getDataByInternalId(second_pair.second, level),
                                         (unsigned char*) getDataByInternalId(curent_pair.second, level),
                                         length[level], book);
					}					 
										 
										 
                    if (curdist < dist_to_query) {
                        good = false;
                        break;
                    }
                }
                if (good) {
                    return_list.push_back(curent_pair);
                }


            }

            for (std::pair<dist_t, tableint> curent_pair : return_list) {

                top_candidates.emplace(-curent_pair.first, curent_pair.second);
            }
        }


        linklistsizeint *get_linklist0(tableint internal_id, int level) const {
			if(level == -1)
                return (linklistsizeint *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetLevel0_);
			else{
				return (linklistsizeint *) (quan_level0_memory_[level] + internal_id * size_quan_per_element_[level] + quan_size_[level]);
			}
        };

        linklistsizeint *get_linklist(tableint internal_id, int level) const {
            return (linklistsizeint *) (linkLists_[internal_id] + (level - 1) * size_links_per_element_);
        };

        void mutuallyConnectNewElement(const void *data_point, tableint cur_c,
                                       std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates,
                                       int level, int level2, float*** book) {

            size_t Mcurmax;	
            if(level2 == -1)			
			    Mcurmax = level ? maxM_ : maxM0_;
			else{
				Mcurmax = level ? maxM_ : maxQ0_;
			}
			
            getNeighborsByHeuristic2(top_candidates, M_, level2, book);
            if (top_candidates.size() > M_)
                throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

            std::vector<tableint> selectedNeighbors;
            selectedNeighbors.reserve(M_);
            while (top_candidates.size() > 0) {
                selectedNeighbors.push_back(top_candidates.top().second);
                top_candidates.pop();
            }

            {
                linklistsizeint *ll_cur;
                if (level == 0)
                    ll_cur = get_linklist0(cur_c, level2);
                else
                    ll_cur = get_linklist(cur_c, level);

                if (*ll_cur) {
                    throw std::runtime_error("The newly inserted element should have blank link list");
                }
                setListCount(ll_cur,selectedNeighbors.size());
                tableint *data = (tableint *) (ll_cur + 1);


                for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                    if (data[idx])
                        throw std::runtime_error("Possible memory corruption");
                    if (level > element_levels_[selectedNeighbors[idx]])
                        throw std::runtime_error("Trying to make a link on a non-existent level");

                    data[idx] = selectedNeighbors[idx];

                }
            }
			dist_t d_max;
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {

                std::unique_lock <std::mutex> lock(link_list_locks_[selectedNeighbors[idx]]);
                
                linklistsizeint *ll_other;
                if (level == 0)
                    ll_other = get_linklist0(selectedNeighbors[idx], level2);
                else
                    ll_other = get_linklist(selectedNeighbors[idx], level);

                size_t sz_link_list_other = getListCount(ll_other);

                if (sz_link_list_other > Mcurmax)
                    throw std::runtime_error("Bad value of sz_link_list_other");
                if (selectedNeighbors[idx] == cur_c)
                    throw std::runtime_error("Trying to connect an element to itself");
                if (level > element_levels_[selectedNeighbors[idx]])
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                tableint *data = (tableint *) (ll_other + 1);
                if (sz_link_list_other < Mcurmax) {
                    data[sz_link_list_other] = cur_c;
                    setListCount(ll_other, sz_link_list_other + 1);
                } else {
                    // finding the "weakest" element to replace it with the new one
					if(level2 == -1)
                        d_max =  getNormByInternalId(cur_c) + getNormByInternalId(selectedNeighbors[idx]) - 2 * 
					            fstdistfunc_(getDataByInternalId(cur_c, level2), getDataByInternalId(selectedNeighbors[idx], level2),
                                                dist_func_param_);
					else{
                        d_max = quandistfunc_((unsigned char*) getDataByInternalId(cur_c, level2), (unsigned char*) getDataByInternalId(selectedNeighbors[idx], level2),
                                                length[level2], book);						
					}							
												
                    // Heuristic:
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidates;
                    candidates.emplace(d_max, cur_c);

                    for (size_t j = 0; j < sz_link_list_other; j++) {
						if(level2 == -1)
                        candidates.emplace(
					            getNormByInternalId(data[j]) + getNormByInternalId(selectedNeighbors[idx]) - 2 *
                                fstdistfunc_(getDataByInternalId(data[j], level2), getDataByInternalId(selectedNeighbors[idx], level2),
                                             dist_func_param_), data[j]);
						else{
                        candidates.emplace(
                                quandistfunc_((unsigned char*) getDataByInternalId(data[j], level2), (unsigned char*) getDataByInternalId(selectedNeighbors[idx], level2),
                                             length[level2], book), data[j]);							
						}					 
                    }

                    getNeighborsByHeuristic2(candidates, Mcurmax, level2, book);

                    int indx = 0;
                    while (candidates.size() > 0) {
                        data[indx] = candidates.top().second;
                        candidates.pop();
                        indx++;
                    }
                    setListCount(ll_other, indx);
                    // Nearest K:
                    /*int indx = -1;
                    for (int j = 0; j < sz_link_list_other; j++) {
                        dist_t d = fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(rez[idx]), dist_func_param_);
                        if (d > d_max) {
                            indx = j;
                            d_max = d;
                        }
                    }
                    if (indx >= 0) {
                        data[indx] = cur_c;
                    } */
                }

            }
        }

        std::mutex global;
        size_t ef_;

        void setEf(size_t ef) {
            ef_ = ef;
        }

        void permutation(int** link_, int vecsize_, int max_level, int* num_points, int level2){  
		    //--------------------------------------------------------
		    //labeltype *getExternalLabeLp(tableint internal_id, int level);
		    
			if(level2 == -1){
			    for(int i = 0; i < vecsize_; i++){
				    link_[0][ *(getExternalLabeLp(i, -1)) ] = i;
			    }			
			}
			else{
			    for(int i = 0; i < max_level; i++){
				    for(int j = 0; j < num_points[i]; j++){
				        link_[i+1][ *(getExternalLabeLp(j, i)) ] = j;
                    }				
			    }
			}
                                   			            			
		    //--------------------------------------------------------
        }
		
		void est_density(float* density, int vecsize_, int vecdim0_){  
		    //--------------------------------------------------------
		    //labeltype *getExternalLabeLp(tableint internal_id, int level);
		    linklistsizeint *ll_other;
			float sigma_ = -1 * vecdim0_ / 2 / log(2);
			for(int i = 0; i < vecsize_; i++){
			    ll_other = get_linklist0(i, -1);
				size_t sz_link_list_other = getListCount(ll_other);
				tableint *data = (tableint *) (ll_other + 1);
				tableint a= data[sz_link_list_other - 1];
				labeltype b = * ( getExternalLabeLp(i, -1) );
				
				dist_t dist1 = -2 * fstdistfunc_(getDataByInternalId(i, -1), getDataByInternalId(a, -1), dist_func_param_) + getNormByInternalId(i) + getNormByInternalId(a);
				if(dist1 < 0) {continue;}
			    else{
			        density[b] = sigma_ * log(dist1);
				}
			}                         
			             
        }
		

        void resizeIndex(size_t new_max_elements){   // not use
            if (new_max_elements<cur_element_count)
                throw std::runtime_error("Cannot resize, max element is less than the current number of elements");


            delete visited_list_pool_;
            visited_list_pool_ = new VisitedListPool(1, new_max_elements);



            element_levels_.resize(new_max_elements);

            std::vector<std::mutex>(new_max_elements).swap(link_list_locks_);


            // Reallocate base layer
            char * data_level0_memory_new = (char *) malloc(new_max_elements * size_data_per_element_);
            if (data_level0_memory_new == nullptr)
                throw std::runtime_error("Not enough memory: resizeIndex failed to allocate base layer");
            memcpy(data_level0_memory_new, data_level0_memory_,cur_element_count * size_data_per_element_);
            free(data_level0_memory_);
            data_level0_memory_=data_level0_memory_new;

            // Reallocate all other layers
            char ** linkLists_new = (char **) malloc(sizeof(void *) * new_max_elements);
            if (linkLists_new == nullptr)
                throw std::runtime_error("Not enough memory: resizeIndex failed to allocate other layers");
            memcpy(linkLists_new, linkLists_,cur_element_count * sizeof(void *));
            free(linkLists_);
            linkLists_=linkLists_new;

            max_elements_=new_max_elements;

        }

        void saveIndex(const std::string &location, float** data) {
            std::ofstream output(location, std::ios::binary);
            std::streampos position;

            writeBinaryPOD(output, max_level_);
            writeBinaryPOD(output, offsetLevel0_);
          //  writeBinaryPOD(output, max_elements_);
			
			for(int i = 0; i < max_level_; i++)
                writeBinaryPOD(output, count_tmp[i]);
			
          //  writeBinaryPOD(output, size_data_per_element_);
            writeBinaryPOD(output, label_offset_);
			writeBinaryPOD(output, offsetNorm_);
            writeBinaryPOD(output, offsetData_);
            writeBinaryPOD(output, maxM_);

            writeBinaryPOD(output, maxM0_);
			writeBinaryPOD(output, maxQ0_);
            writeBinaryPOD(output, M_);
            writeBinaryPOD(output, mult_);
            writeBinaryPOD(output, ef_construction_);

	    //printf("ef_construction = %d elem = %d, size = %d\n", ef_construction_, max_elements_, size_data_per_element_);
			for(int i = 0; i < max_level_; i++){
                            
	                      // printf("count = %d, size = %d\n", count_tmp[i], size_quan_per_element_[i]);
				output.write(quan_level0_memory_[i], count_tmp[i] * size_quan_per_element_[i]);

			}			             
            output.close();
        }

        void loadIndex(const std::string &location1, const std::string &location2, SpaceInterface<dist_t> *s, size_t max_elements_i=0) {

            std::ifstream input(location1, std::ios::binary);
          
            std::ifstream input2(location2, std::ios::binary);
            readBinaryPOD(input2, max_elements_);
            readBinaryPOD(input2, size_data_per_element_);

	    
            data_level0_memory_ = (char *) malloc(max_elements_ * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
            input2.read(data_level0_memory_, max_elements_ * size_data_per_element_);

            if (!input.is_open())
                throw std::runtime_error("Cannot open file");

            // get file size:
            input.seekg(0,input.end);
            std::streampos total_filesize=input.tellg();
            input.seekg(0,input.beg);


            readBinaryPOD(input, max_level_);
            readBinaryPOD(input, offsetLevel0_);
//            readBinaryPOD(input, max_elements_);

	        count_tmp = new int[max_level_];
	    
			for(int i = 0; i < max_level_; i++)
                readBinaryPOD(input, count_tmp[i]);
			
			 size_t max_elements=max_elements_;
  //          readBinaryPOD(input, size_data_per_element_);
            readBinaryPOD(input, label_offset_);
			
			readBinaryPOD(input, offsetNorm_);
            readBinaryPOD(input, offsetData_);
            readBinaryPOD(input, maxM_);
            readBinaryPOD(input, maxM0_);
			readBinaryPOD(input, maxQ0_);
            readBinaryPOD(input, M_);
            readBinaryPOD(input, mult_);
            readBinaryPOD(input, ef_construction_);
      
			//---------------------------------------------------------
			size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
            size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
			size_quan_level0_ = maxQ0_ * sizeof(tableint) + sizeof(linklistsizeint);
			
			length = new int[max_level_];
			quan_size_ = new size_t[max_level_];
            size_quan_per_element_ = new size_t[max_level_];
			quan_offset_ = new size_t[max_level_];
			quan_level0_memory_ = new char* [max_level_];
			
			for (int i = 0; i < max_level_; i++)
			    length[i] = pow(2, max_level_ -i + OFF);
			
			for(int i = 0; i < max_level_; i++){
			    quan_size_[i] = length[i] * sizeof(unsigned char);	
				size_quan_per_element_[i] = size_quan_level0_ + quan_size_[i] + sizeof(labeltype);
				quan_offset_[i] = size_quan_level0_ + quan_size_[i];

	                     //  printf("count = %d, size = %d\n", count_tmp[i], size_quan_per_element_[i]);
				quan_level0_memory_[i] = (char *) malloc(count_tmp[i] * size_quan_per_element_[i]);
                if (quan_level0_memory_[i] == nullptr)
                throw std::runtime_error("Not enough memory");				
			}
			//---------------------------------------------------------------------
				
            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();     


	    for(int i = 0; i < max_level_; i++){
                input.read(quan_level0_memory_[i], count_tmp[i] * size_quan_per_element_[i]);
            }

		
            std::vector<std::mutex>(max_elements).swap(link_list_locks_);


            visited_list_pool_ = new VisitedListPool(1, max_elements);

			
            element_levels_ = std::vector<int>(max_elements);
            revSize_ = 1.0 / mult_;
            ef_ = 10;
					
            has_deletions_=false;
            input.close();
            return;
        }

        template<typename data_t>
        std::vector<data_t> getDataByLabel(labeltype label)  //not use
        {
            tableint label_c;
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end() ) {  //|| isMarkedDeleted(search->second)
                throw std::runtime_error("Label not found");
            }
            label_c = search->second;

            char* data_ptrv = getDataByInternalId(label_c);
            size_t dim = *((size_t *) dist_func_param_);
            std::vector<data_t> data;
            data_t* data_ptr = (data_t*) data_ptrv;
            for (int i = 0; i < dim; i++) {
                data.push_back(*data_ptr);
                data_ptr += 1;
            }
            return data;
        }


 //       static const unsigned char DELETE_MARK = 0x01;
//        static const unsigned char REUSE_MARK = 0x10;
        /**
         * Marks an element with the given label deleted, does NOT really change the current graph.
         * @param label
         */
        void markDelete(labeltype label)
        {

        }

        /**
         * Uses the first 8 bits of the memory for the linked list to store the mark,
         * whereas maxM0_ has to be limited to the lower 24 bits, however, still large enough in almost all cases.
         * @param internalId
         */
        void markDeletedInternal(tableint internalId) {

        }

        /**
         * Remove the deleted mark of the node.
         * @param internalId
         */
        void unmarkDeletedInternal(tableint internalId) {

        }

        /**
         * Checks the first 8 bits of the memory to see if the element is marked deleted.
         * @param internalId
         * @return
         */
        bool isMarkedDeleted(tableint internalId) const {

        }

        unsigned short int getListCount(linklistsizeint * ptr) const {
            return *((unsigned short int *)ptr);
        }

        void setListCount(linklistsizeint * ptr, unsigned short int size) const {
            *((unsigned short int*)(ptr))=*((unsigned short int *)&size);
        }

        void addPoint(const void *data_point, labeltype label, float*** book, int level2, bool flag) {
		addPoint(data_point, label, -1, book, level2, flag);
        }
			
		float quandistfunc_(unsigned char* obj1, unsigned char* obj2, int num, float*** dist_book) const{
			float dist = 0;
			for(int i = 0; i < num; i++)
			    dist += dist_book[i][ obj1[i] ][ obj2[i] ];
			
			return dist;
		}

		float querydistfunc_(unsigned char* id, int num, float** book) const{
			float dist = 0;
            //_mm_prefetch( id, _MM_HINT_T0);
			 
			for(int i = 0; i < num; i++){
			    dist += book[i][ id[i] ];
			}                   
			return dist;
		}

		
		void deleteLinklist(int ii, std::string location){
			for(int i = 0; i < max_elements_; i++){
                if (element_levels_[i] > 0){
                    free(linkLists_[i]);
					element_levels_[i] = 0;
				}					
			}
			
                        if(ii == -1){
                               std::ofstream output(location, std::ios::binary);

                                   writeBinaryPOD(output, max_elements_);
                                 writeBinaryPOD(output, size_data_per_element_);
                                output.write(data_level0_memory_, max_elements_ * size_data_per_element_);
                                output.close();

				free(data_level0_memory_);
			    data_level0_memory_ = NULL;
			}
			
			int k = ii + 1;
            if(k <= max_level_ - 1){
			    quan_size_[k] = length[k] * sizeof(unsigned char);	
				size_quan_per_element_[k] = size_quan_level0_ + quan_size_[k] + sizeof(labeltype);
				quan_offset_[k] = size_quan_level0_ + quan_size_[k];
				quan_level0_memory_[k] = (char *) malloc(max_elements_ * size_quan_per_element_[k]);
                if (quan_level0_memory_[k] == nullptr)
                throw std::runtime_error("Not enough memory");				
          	}		
		}
		
		//-----------------------------------------------------------------
        tableint addPoint(const void *data_point, labeltype label, int level, float*** book, int level2, bool flag) {	
		tableint cur_c = 0;
			dist_t curdist;
			if(flag == true){    //the first point for every level
                enterpoint_node_ = -1;
                maxlevel_ = -1;
			    count0 = 0;
				if(visited_list_pool_!= NULL){
				    delete visited_list_pool_;
				}	
                visited_list_pool_ = new VisitedListPool(1, max_elements_);				
			}
            
            {
		        std::unique_lock <std::mutex> lock(cur_element_count_guard_);
		
		        if ( (level2 >= 0 && count_tmp[level2] >= max_elements_) || count0 >= max_elements_) {
                    throw std::runtime_error("The number of elements exceeds the specified limit");
                };
		
				if(level2 >= 0){
                    cur_c = count_tmp[level2];  //set to 0
                    count_tmp[level2]++;
		        }
		    	else{
		    		cur_c = count0;
		    		count0++;
		    	}
            }

            std::unique_lock <std::mutex> lock_el(link_list_locks_[cur_c]);
            int curlevel = getRandomLevel(mult_);
            if (level > 0)
                curlevel = level;

            element_levels_[cur_c] = curlevel;


            std::unique_lock <std::mutex> templock(global);
            int maxlevelcopy = maxlevel_;
            if (curlevel <= maxlevelcopy)
                templock.unlock();
            tableint currObj = enterpoint_node_;
            tableint enterpoint_copy = enterpoint_node_;

               // printf("ccheck2\n");	
            if(level2 == -1)
                memset(data_level0_memory_ + cur_c * size_data_per_element_ + offsetLevel0_, 0, size_data_per_element_);
            else{
				memset(quan_level0_memory_[level2] + cur_c * size_quan_per_element_[level2] + offsetLevel0_, 0, size_quan_per_element_[level2]);
			}
			
            // Initialisation of the data and label
			
			float temp_norm = 0;
			float* vec_;
			
			if(level2 == -1){
				vec_ = (float*) data_point;
				for(int i = 0; i < vecdim; i++)
				    temp_norm += vec_[i] * vec_[i];
                memcpy(getExternalLabeLp(cur_c, level2), &label, sizeof(labeltype));
				memcpy(getNormByInternalId2(cur_c), &temp_norm, sizeof(float));				
                memcpy(getDataByInternalId(cur_c, level2), data_point, data_size_);	
                              
			}
			else{
                memcpy(getExternalLabeLp(cur_c, level2), &label, sizeof(labeltype));
                memcpy(getDataByInternalId(cur_c, level2), data_point, quan_size_[level2]);				
			}
			
            if (curlevel) {
                linkLists_[cur_c] = (char *) malloc(size_links_per_element_ * curlevel + 1);
                if (linkLists_[cur_c] == nullptr)
                    throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
                memset(linkLists_[cur_c], 0, size_links_per_element_ * curlevel + 1);
            }
			dist_t d;
            if ((signed)currObj != -1) {

                if (curlevel < maxlevelcopy) {
                    if(level2 == -1)
                        curdist = temp_norm +  getNormByInternalId(currObj) - 2 * fstdistfunc_(data_point, getDataByInternalId(currObj, level2), dist_func_param_);
					else{
						curdist = quandistfunc_((unsigned char*) data_point, (unsigned char*) getDataByInternalId(currObj, level2), length[level2], book);
					}
					
                    for (int level = maxlevelcopy; level > curlevel; level--) {


                        bool changed = true;
                        while (changed) {
                            changed = false;
                            unsigned int *data;
                            std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                            data = get_linklist(currObj,level);
                            int size = getListCount(data);

                            tableint *datal = (tableint *) (data + 1);
                            for (int i = 0; i < size; i++) {
                                tableint cand = datal[i];
                                if (cand < 0 || cand > max_elements_)
                                    throw std::runtime_error("cand error");
								if(level2 == -1)
                                    d = temp_norm +  getNormByInternalId(cand) - 2 * fstdistfunc_(data_point, getDataByInternalId(cand, level2), dist_func_param_);
								else{
									d = quandistfunc_((unsigned char*) data_point, (unsigned char*) getDataByInternalId(cand, level2), length[level2], book);
								}
								
                                if (d < curdist) {
                                    curdist = d;
                                    currObj = cand;
                                    changed = true;
                                }
                            }
                        }
                    }
                }

				bool epDeleted = false;
				
                for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
                    if (level > maxlevelcopy || level < 0)  // possible?
                        throw std::runtime_error("Level error");

                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                            currObj, data_point, temp_norm, level, level2, book);

                    mutuallyConnectNewElement(data_point, cur_c, top_candidates, level, level2, book);

                    currObj = top_candidates.top().second;
                }


            } else {
                // Do nothing for the first element
                enterpoint_node_ = 0;
                maxlevel_ = curlevel;

            }

            //Releasing lock for the maximum level
            if (curlevel > maxlevelcopy) {
                enterpoint_node_ = cur_c;
                maxlevel_ = curlevel;
            }
            return cur_c;
        };

        std::priority_queue<std::pair<dist_t, labeltype >>
        searchKnn(const void *query_data, size_t k, float*** book, int** start_book, int** init_obj, int* ord) const {
			
            std::priority_queue<std::pair<dist_t, labeltype >> result;
  
			float min_sum;
			int min_id, id1, id2;
			for( int i = 0; i < L; i++){
				if(i == 0) {min_id = 0; min_sum = book[max_level_-1][0][i];}
			    else{
					if(book[max_level_-1][0][i] < min_sum) { min_id = i; min_sum = book[max_level_-1][0][i];}
				}
			}
			id1 = min_id;

			for( int i = 0; i < L; i++){
				if(i == 0) {min_id = 0; min_sum = book[max_level_-1][1][i];}
			    else{
					if(book[max_level_-1][1][i] < min_sum) { min_id = i; min_sum = book[max_level_-1][1][i];}
				}
			}
			id2 = min_id;

            tableint currObj = start_book[id1][id2]; //revise
			
            dist_t curdist = querydistfunc_((unsigned char*) getDataByInternalId(currObj, max_level_-2), length[max_level_-2], book[max_level_-2]);
			
			int length_;
			float** book_;
            for (int level = max_level_-2; level >= 0; level--) { 
                bool changed = true;
				length_ = length[level];
				book_ = book[level];
				
                while (changed) {
                    changed = false;
                    unsigned int *data;

					data = (unsigned int *) get_linklist0(currObj, level);
					
                    int size = getListCount(data);
                    tableint *datal = (tableint *) (data + 1);
                    for (int i = 0; i < size; i++) {
                        tableint cand = datal[i];
                        if (cand < 0 || cand > max_elements_){
			            printf("cand = %d\n",cand);
                            throw std::runtime_error("cand error");
			        }
                        //dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_);
						dist_t d = querydistfunc_((unsigned char*) getDataByInternalId(cand, level), length_, book_);

                        if (d < curdist) {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
					if(changed == false ){
						if(level > 0){
						    currObj = init_obj[level][ *(getExternalLabeLp(currObj, level)) ];  //currobj is internal id , should be replaced first!!!!
						    curdist = querydistfunc_((unsigned char*) getDataByInternalId(currObj, level-1), length[level-1], book[level-1]);
						}
						else{
							currObj = *(getExternalLabeLp(currObj, 0));
							curdist = fstdistfunc_(query_data, getDataByInternalId(currObj, -1), dist_func_param_);
						}
					}
                }
            }

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates1=searchBaseLayerST<false>(
                        currObj, query_data, std::max(ef_, k));
                top_candidates.swap(top_candidates1);
            
            while (top_candidates.size() > k) {
                top_candidates.pop();
            }
            while (top_candidates.size() > 0) {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second, -1)));
                top_candidates.pop();
            }
            return result;
        };
		
        void connect(const void *query_data, float*** book, unsigned int* init_obj, int ii, bool flag){
						
            std::priority_queue<std::pair<dist_t, labeltype >> result;
           
            tableint currObj = enterpoint_node_;
          
			dist_t curdist;
			if(ii > -1)
  			    curdist = quandistfunc_((unsigned char*) query_data, 
			    (unsigned char*) getDataByInternalId(enterpoint_node_, ii), length[ii], book);
            else{
				curdist = fstdistfunc_(query_data, getDataByInternalId(enterpoint_node_, ii), dist_func_param_);
			}
				
		    dist_t d;
            for (int level = maxlevel_; level > 0; level--) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    unsigned int *data;

                    data = (unsigned int *) get_linklist(currObj, level);
                    int size = getListCount(data);
                    tableint *datal = (tableint *) (data + 1);
                    for (int i = 0; i < size; i++) {
                        tableint cand = datal[i];
                        if (cand < 0 || cand > max_elements_)
                            throw std::runtime_error("cand error");
						if(ii > -1)
                            d = quandistfunc_((unsigned char*) query_data, 
						    (unsigned char*) getDataByInternalId(cand, ii), length[ii], book);
						else{	
							d = fstdistfunc_(query_data, getDataByInternalId(cand, ii), dist_func_param_);
						}
						
                        if (d < curdist) {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
                }
            }

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates1=searchQuanST(
                        currObj, query_data, book, ii);   //ef_ has been set to 500
            
			top_candidates.swap(top_candidates1);
            int rec;
			
			if(flag == true){
				rec = KK - 1;
                while (top_candidates.size() > KK) {
                    top_candidates.pop();
                }
			}
			else{
				rec = fan -1;
                while (top_candidates.size() > fan) {
                    top_candidates.pop();
                }				
			}
			
            while (top_candidates.size() > 0) {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                init_obj[rec] = rez.second;  //internal ID
                top_candidates.pop();
				rec--;
            }
			
			
			if(rec >= 0){
				printf("rec =%d -- connect_error\n",rec);
				exit(-1);
		    }
						
        };	

		std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> 
	    searchQuanST(tableint ep_id, const void *data_point, float*** book, int ii) const {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
			dist_t dist;
 		
		    if(ii > -1)
                dist = quandistfunc_((unsigned char*) data_point, (unsigned char*) getDataByInternalId(ep_id, ii), length[ii], book);		
			else{
				dist = fstdistfunc_(data_point, getDataByInternalId(ep_id, ii), dist_func_param_);
			}
			
			lowerBound = dist;
			
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);
		
            visited_array[ ep_id] = visited_array_tag;

            while (!candidate_set.empty()) {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound) {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;
                int *data = (int *) get_linklist0(current_node_id, ii);
                size_t size = getListCount((linklistsizeint*)data);

				
                for (size_t j = 1; j <= size; j++) {
                    int candidate_id = *(data + j);

                    if (!(visited_array[candidate_id] == visited_array_tag)) {

                        visited_array[candidate_id] = visited_array_tag;

						if(ii > -1)
                        dist = quandistfunc_((unsigned char*) data_point, 
						        (unsigned char*) getDataByInternalId(candidate_id, ii), length[ii], book);
                        else{
							dist = fstdistfunc_(data_point, getDataByInternalId(candidate_id, ii), dist_func_param_);
						}
								
                        if (top_candidates.size() < ef_ || lowerBound > dist) {
                            candidate_set.emplace(-dist, candidate_id);

                                top_candidates.emplace(dist, candidate_id);

                            if (top_candidates.size() > ef_)
                                top_candidates.pop();

                            if (!top_candidates.empty())
                                lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }
			
			if (top_candidates.size() > ef_)
                top_candidates.pop();

            visited_list_pool_->releaseVisitedList(vl);
            return top_candidates;
        };		

        template <typename Comp>
        std::vector<std::pair<dist_t, labeltype>>
        searchKnn(const void* query_data, size_t k, Comp comp) {
            std::vector<std::pair<dist_t, labeltype>> result;
            //auto ret = searchKnn(query_data, k);
			/*
            auto ret;
            while (!ret.empty()) {
                result.push_back(ret.top());
                ret.pop();
            }

            std::sort(result.begin(), result.end(), comp);
            */
            return result;
        }
		
        float querydistfunc_(unsigned char* id, int num, float** book){
	        float dist = 0;
	        for(int i = 0; i < num; i++){
	            dist += book[i][ id[i] ];
	        }
        	return dist;
        }


void SearchWithsingleGraph(float** book_, int X, unsigned* points, 
    unsigned* enter_point, unsigned int* trans, size_t ef, VisitedListPool *visited_list_pool_) {
  //unsigned L = parameters.Get<unsigned>("L_search");
  unsigned LL = ef;
   VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
  
  
  float sum, dist;
  unsigned int currObj, min_id, fan_new;
  if(LL < fan){
	fan_new = LL;
  }
  else{
	fan_new = fan;
  }
  std::vector<Neighbor> retset(LL + 1);
  
  size_t quan_len_ = quan_size_[0];
  size_t size_quan_per_element = size_quan_per_element_[0];
  char* quan_graph_ = quan_level0_memory_[0];
  	
    for(int i = 0; i< fan_new; i++){
       sum = querydistfunc_( (unsigned char*) quan_graph_ + size_quan_per_element * points[i], X, book_);
       retset[i] = Neighbor(points[i], sum, true, true);
	   visited_array[points[i]] = visited_array_tag;
    }
	
    std::sort(retset.begin(), retset.begin() + fan_new);

	int L2 = fan_new;
    int k = 0;

    while (k < (int)LL) {
    int nk = LL;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

       _mm_prefetch(quan_graph_ + size_quan_per_element * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + size_quan_per_element * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;

      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + size_quan_per_element * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
        if (visited_array[id] == visited_array_tag) continue;
        visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + size_quan_per_element * id, X, book_);
		int r;
        if (L2 >= LL && dist >= retset[LL - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < LL && dist >= retset[L2-1].distance){
                        retset[L2] = nn; 
                        r = L2;
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), L2, nn); }
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }

   for (size_t i = 0; i < LL; i++) {		
	    enter_point[i] = trans[retset[i].id];
   } 
   visited_list_pool_->releaseVisitedList(vl);
}




void SearchWithquanGraph(float** book_, int X, unsigned* points, elem* init_point, 
        unsigned int* trans, size_t ef, char* fflag, int level, VisitedListPool *visited_list_pool_) {
 // unsigned L = parameters.Get<unsigned>("L_search");
     unsigned LL = ef;
    VisitedList *vl = visited_list_pool_->getFreeVisitedList();
    vl_type *visited_array = vl->mass;
    vl_type visited_array_tag = vl->curV;  
  

  float sum, dist;
  unsigned int currObj, min_id, fan_new;
  if(LL < fan){
	fan_new = LL;
  }
  else{
	fan_new = fan;
  }
  std::vector<Neighbor> retset(LL + 1);

  size_t quan_len_ = quan_size_[level];
  size_t size_quan_per_element = size_quan_per_element_[level];
  char* quan_graph_ = quan_level0_memory_[level];	
    	
    for(int i = 0; i< fan_new; i++){
       sum = querydistfunc_( (unsigned char*) quan_graph_ + size_quan_per_element * points[i], X, book_);
       retset[i] = Neighbor(points[i], sum, true, true);
	   visited_array[points[i]] = visited_array_tag;
    }
	
    std::sort(retset.begin(), retset.begin() + fan_new);

	int L2 = fan_new;
    int k = 0;

    while (k < (int)LL) {

    int nk = LL;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + size_quan_per_element * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + size_quan_per_element * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;

      for (unsigned m = 0; m < MaxM; ++m)
	   _mm_prefetch(quan_graph_ + size_quan_per_element * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + size_quan_per_element * id, X, book_);
		int r;
        if (L2 >= LL && dist >= retset[LL - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < LL && dist >= retset[L2-1].distance){
            retset[L2] = nn; 
            r = L2;
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), L2, nn); }
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }

   for (size_t i = 0; i < LL; i++) {	
	if(fflag[ retset[i].id ] == 1){
		init_point[i].flag = true;
    }
	
    else{ 
	    init_point[i].flag = false;
    }
	
	init_point[i].id = trans[retset[i].id];
	init_point[i].dist = retset[i].distance;
 }
  visited_list_pool_->releaseVisitedList(vl); 
}

void SearchWithquanGraph2(float** book_, int X, elem* points, 
    unsigned* init_point, unsigned int* trans, size_t ef, VisitedListPool *visited_list_pool_) {
    unsigned LL = ef;

   VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
   
  float dist, sum;
   std::vector<Neighbor> retset(LL + 1);
 
	size_t quan_len_ = quan_size_[0];
    size_t size_quan_per_element = size_quan_per_element_[0];
    char* quan_graph_ = quan_level0_memory_[0];	
	
	for(int i = 0; i < LL; i++){
		if(points[i].flag == true){
		    sum = querydistfunc_( (unsigned char*) quan_graph_ + size_quan_per_element * points[i].id, X, book_);
            retset[i] = Neighbor(points[i].id, sum, true, true);
			visited_array[retset[i].id] = visited_array_tag;
		}
	    else{
			retset[i] = Neighbor(points[i].id, points[i].dist, false, false);
		}
	}
	
	std::sort(retset.begin(), retset.begin() + LL);
	int L2= LL;
	
    int k = 0;
    while (k < (int)LL) {

    int nk = LL;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + size_quan_per_element * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + size_quan_per_element * n + quan_len_);
      unsigned MaxM = *neighbors;

      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + size_quan_per_element * neighbors[m], _MM_HINT_T0);

      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + size_quan_per_element * id, X, book_);
		int r;
        if (L2 >= LL && dist >= retset[LL - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < LL){
			r = InsertIntoPool(retset.data(), L2, nn);
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), LL, nn);}
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }

  for (size_t i = 0; i < LL; i++) {	
      if(retset[i].flag2 == false){
		  init_point[i] = retset[i].id;
         }
	  else{
		  init_point[i] = trans[retset[i].id];
         }
  } 
  visited_list_pool_->releaseVisitedList(vl); 
}

void SearchWithquanGraph3(float** book_, int X, elem* points, elem* init_point, 
    unsigned int* trans, size_t ef, char* fflag, int level, VisitedListPool *visited_list_pool_) {
  unsigned LL = ef;
  
     VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV;  

  float dist, sum;
   std::vector<Neighbor> retset(LL + 1);
  
	size_t quan_len_ = quan_size_[level];
    size_t size_quan_per_element = size_quan_per_element_[level];
    char* quan_graph_ = quan_level0_memory_[level];	
		
	for(int i = 0; i < LL; i++){
		if(points[i].flag == true){
		    sum = querydistfunc_( (unsigned char*) quan_graph_ + size_quan_per_element * points[i].id, X, book_);
            retset[i] = Neighbor(points[i].id, sum, true, true);
			visited_array[retset[i].id] = visited_array_tag;
		}
	    else{
			retset[i] = Neighbor(points[i].id, points[i].dist, false, false);
		}
	}
	
	std::sort(retset.begin(), retset.begin() + LL);
	int L2= LL;
	
    int k = 0;
    while (k < (int)LL) {
    int nk = LL;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + size_quan_per_element * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + size_quan_per_element * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + size_quan_per_element * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + size_quan_per_element * id, X, book_);
		int r;
        if (L2 >= LL && dist >= retset[LL - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < LL){
			r = InsertIntoPool(retset.data(), L2, nn);
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), LL, nn);}
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
   for (size_t i = 0; i < LL; i++) {
    
    if(retset[i].flag2 == false) {
	    init_point[i].id = retset[i].id;
		init_point[i].dist = retset[i].distance;
	    init_point[i].flag = false;
        continue;		
	}
	
	if(fflag[ retset[i].id ] == 1){
		init_point[i].flag = true;
    }
	
    else{ 
	    init_point[i].flag = false;
    }
	
	init_point[i].id = trans[retset[i].id];
	init_point[i].dist = retset[i].distance;
   } 
   visited_list_pool_->releaseVisitedList(vl);
}


void SearchWithOptGraph(const float *query, size_t K, size_t ef, 
        unsigned *indices, unsigned *init_point, VisitedListPool *visited_list_pool_) {

  unsigned LL = ef;

  std::vector<Neighbor> retset(LL + 1);
  std::vector<unsigned> init_ids(LL);

     VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
  

  for (unsigned i = 0; i < LL; i++) {
	  float* data = (float *) (data_level0_memory_ + init_point[i] * size_data_per_element_ + offsetNorm_);
	float norm = *data;
    data++;
	float dist = norm - 2 * fstdistfunc_( (const void*) query,  (const void*) data, dist_func_param_);
    retset[i] = Neighbor(init_point[i], dist, true, true);
    visited_array[init_point[i]] = visited_array_tag;
  }
  
  std::sort(retset.begin(), retset.begin() + LL);
  int k = 0;
  while (k < (int)LL) {
    int nk = LL;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

	  _mm_prefetch(data_level0_memory_ + size_data_per_element_ * n, _MM_HINT_T0);   	  
	  unsigned *neighbors = (unsigned *)(data_level0_memory_ + size_data_per_element_ * n);
	  
      unsigned MaxM = *neighbors;
      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	   _mm_prefetch(data_level0_memory_ + neighbors[m] * size_data_per_element_ + offsetNorm_, _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float *data = (float *)(data_level0_memory_ + id * size_data_per_element_ + offsetNorm_);
        float norm = *data;
        data++;
        float dist = norm - 2 * fstdistfunc_( (const void*)query, (const void*)data, dist_func_param_);

        if (dist >= retset[LL - 1].distance) continue;
        Neighbor nn(id, dist, true, true);
        int r = InsertIntoPool(retset.data(), LL, nn);

        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
  for (size_t i = 0; i < K; i++) {
    indices[i] = getExternalLabel(retset[i].id, -1);

  }
  visited_list_pool_->releaseVisitedList(vl);
}				

    };

}
