#include "efanna2e/index_nsg.h"
#include <omp.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
#include <queue>

#include "efanna2e/exceptions.h"
#include "efanna2e/parameters.h"

namespace efanna2e {
#define _CONTROL_NUM 100
  // IndexNSG::IndexNSG(){}
  
IndexNSG::IndexNSG(const size_t dimension, const size_t n, Metric m,
                   Index *initializer, const bool flag, float*** book, const int l, float*** quanti, const int level)
  : Index(dimension, n, m, flag, book, l, quanti, level), initializer_{initializer} {}

IndexNSG::~IndexNSG() {}

struct CompareByFirst {
    constexpr bool operator()(std::pair<float, unsigned int> const &a,
        std::pair<float, unsigned int> const &b) const noexcept {
            return a.first < b.first;
        }
};

void IndexNSG::Save(const char *filename) {
  std::ofstream out(filename, std::ios::binary | std::ios::out);
  assert(final_graph_.size() == nd_);

  out.write((char *)&width, sizeof(unsigned));
  out.write((char *)&ep_, sizeof(unsigned));
  for (unsigned i = 0; i < nd_; i++) {
    unsigned GK = (unsigned)final_graph_[i].size();
    out.write((char *)&GK, sizeof(unsigned));
    out.write((char *)final_graph_[i].data(), GK * sizeof(unsigned));
	std::vector<unsigned>().swap(final_graph_[i]); // clear space
  }
  CompactGraph().swap(final_graph_);
  //std::vector<std::vector<unsigned > > ().swap(final_graph_);//clear space
  out.close();
}

void IndexNSG::Load(const char *filename) {
  std::ifstream in(filename, std::ios::binary);
  in.read((char *)&width, sizeof(unsigned));
  in.read((char *)&ep_, sizeof(unsigned));
  // width=100;
  unsigned cc = 0;
  while (!in.eof()) {
    unsigned k;
    in.read((char *)&k, sizeof(unsigned));
    if (in.eof()) break;
    cc += k;
    std::vector<unsigned> tmp(k);
    in.read((char *)tmp.data(), k * sizeof(unsigned));
    final_graph_.push_back(tmp);
  }
  cc /= nd_;
  // std::cout<<cc<<std::endl;
}
void IndexNSG::Load_nn_graph(const char *filename) {
  std::ifstream in(filename, std::ios::binary);
  unsigned k;
  in.read((char *)&k, sizeof(unsigned));
  k = KK;    //k is fixed to 100
  in.seekg(0, std::ios::end);
  std::ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  size_t num = (unsigned)(fsize / (k + 1) / 4);
  in.seekg(0, std::ios::beg);

  final_graph_.resize(num);
  final_graph_.reserve(num);
  unsigned kk = (k + 3) / 4 * 4;
  for (size_t i = 0; i < num; i++) {
    in.seekg(4, std::ios::cur);
    final_graph_[i].resize(k);
    final_graph_[i].reserve(kk);
    in.read((char *)final_graph_[i].data(), k * sizeof(unsigned));
  }
  in.close();
}

void IndexNSG::get_neighbors(const float *query, const Parameters &parameter,
                             std::vector<Neighbor> &retset,
                             std::vector<Neighbor> &fullset) {
  unsigned L = parameter.Get<unsigned>("L");
  int sdim_ = dimension_/length_;
  retset.resize(L + 1);
  std::vector<unsigned> init_ids(L);
  // initializer_->Search(query, nullptr, L, parameter, init_ids.data());

  boost::dynamic_bitset<> flags{nd_, 0};
  L = 0;
  for (unsigned i = 0; i < init_ids.size() && i < final_graph_[ep_].size(); i++) {
    init_ids[i] = final_graph_[ep_][i];
    flags[init_ids[i]] = true;
    L++;
  }
  while (L < init_ids.size()) {
    unsigned id = rand() % nd_;
    if (flags[id]) continue;
    init_ids[L] = id;
    L++;
    flags[id] = true;
  }

  L = 0;
  float* arr = new float[dimension_];
  unsigned char index;
  float dist;
  
  for (unsigned i = 0; i < init_ids.size(); i++) {
    unsigned id = init_ids[i];
    if (id >= nd_) continue;

    if(is_quan == false){
		dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                    (unsigned)dimension_);
	}
	else{
	    for(int j = 0; j < length_; j++){
			index = id* length_ + j;
	        for(int l = 0; l < sdim_; l++){
		        arr[j * sdim_ + l] = quantizer[j][ data2_[index] ][l];
	        }
        }
		dist = distance_->compare(arr, query,
                                    (unsigned)dimension_);		
	}
	
    retset[i] = Neighbor(id, dist, true, true);

    L++;
  }

  std::sort(retset.begin(), retset.begin() + L);
  int k = 0;
  while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      for (unsigned m = 0; m < final_graph_[n].size(); ++m) {
        unsigned id = final_graph_[n][m];
        if (flags[id]) continue;
        flags[id] = 1;

        if(is_quan == false){
		    dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                    (unsigned)dimension_);
	    }
	    else{
	        for(int j = 0; j < length_; j++){
			    index = id* length_ + j;
	            for(int l = 0; l < sdim_; l++){
		            arr[j * sdim_ + l] = quantizer[j][ data2_[index] ][l];
	            }
            }
		    dist = distance_->compare(arr, query,
                                    (unsigned)dimension_);		
	    }
	
        Neighbor nn(id, dist, true, true);
        fullset.push_back(nn);
        if (dist >= retset[L - 1].distance) continue;
        int r = InsertIntoPool(retset.data(), L, nn);

        if (L + 1 < retset.size()) ++L;
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
}

void IndexNSG::get_neighbors(const float *query, const Parameters &parameter,
                             boost::dynamic_bitset<> &flags,
                             std::vector<Neighbor> &retset,
                             std::vector<Neighbor> &fullset) {
  unsigned L = parameter.Get<unsigned>("L");
  int sdim_ = dimension_/length_;
  retset.resize(L + 1);
  std::vector<unsigned> init_ids(L);
  // initializer_->Search(query, nullptr, L, parameter, init_ids.data());

  L = 0;
  float* arr = new float[dimension_];
  unsigned char index;
  float dist;
  
  for (unsigned i = 0; i < init_ids.size() && i < final_graph_[ep_].size(); i++) {
    init_ids[i] = final_graph_[ep_][i];
    flags[init_ids[i]] = true;
    L++;
  }
  while (L < init_ids.size()) {
    unsigned id = rand() % nd_;
    if (flags[id]) continue;
    init_ids[L] = id;
    L++;
    flags[id] = true;
  }

  L = 0;
  for (unsigned i = 0; i < init_ids.size(); i++) {
    unsigned id = init_ids[i];
    if (id >= nd_) continue;

	/*
    float dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                    (unsigned)dimension_);
	*/
    if(is_quan == false){
		dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                    (unsigned)dimension_);
	}
	else{
	    for(int j = 0; j < length_; j++){
			index = id* length_ + j;
	        for(int l = 0; l < sdim_; l++){
		        arr[j * sdim_ + l] = quantizer[j][ data2_[index] ][l];
	        }
        }
		dist = distance_->compare(arr, query,
                                    (unsigned)dimension_);		
	}	
									
    retset[i] = Neighbor(id, dist, true, true);
    fullset.push_back(retset[i]);

    L++;
  }

  std::sort(retset.begin(), retset.begin() + L);
  int k = 0;
  while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      for (unsigned m = 0; m < final_graph_[n].size(); ++m) {
        unsigned id = final_graph_[n][m];
        if (flags[id]) continue;
        flags[id] = 1;
        /*
        float dist = distance_->compare(query, data_ + dimension_ * (size_t)id,
                                        (unsigned)dimension_);
		*/
	    if(is_quan == false){
		dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                    (unsigned)dimension_);
	    }
	    else{
	        for(int j = 0; j < length_; j++){
			    index = id* length_ + j;
	            for(int l = 0; l < sdim_; l++){
		            arr[j * sdim_ + l] = quantizer[j][ data2_[index] ][l];
	            }
            }
		    dist = distance_->compare(arr, query,
                                    (unsigned)dimension_);		
	    }	
        Neighbor nn(id, dist, true, true);
        fullset.push_back(nn);
        if (dist >= retset[L - 1].distance) continue;
        int r = InsertIntoPool(retset.data(), L, nn);

        if (L + 1 < retset.size()) ++L;
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
}

void IndexNSG::init_graph(const Parameters &parameters) {
  float *center = new float[dimension_];
  for (unsigned j = 0; j < dimension_; j++) center[j] = 0;
  for (unsigned i = 0; i < nd_; i++) {
    for (unsigned j = 0; j < dimension_; j++) {
      center[j] += data_[i * dimension_ + j];
    }
  }
  for (unsigned j = 0; j < dimension_; j++) {
    center[j] /= nd_;
  }
  std::vector<Neighbor> tmp, pool;
  ep_ = rand() % nd_;  // random initialize navigating point
  get_neighbors(center, parameters, tmp, pool);
  ep_ = tmp[0].id;
}

void IndexNSG::init_graph2(const Parameters &parameters) {
  float *center = new float[dimension_];
  int sdim_ = dimension_ / length_; 
  
  float* arr = new float[dimension_];
  for (unsigned j = 0; j < dimension_; j++) center[j] = 0;
 
  for (unsigned i = 0; i < nd_; i++) {
	for(int j = 0; j < length_; j++){
	    for(int l = 0; l < sdim_; l++){
		    arr[j * sdim_ + l] = quantizer[j][ data2_[i * length_ + j] ][l];
	    }
    }	  
    for (unsigned j = 0; j < dimension_; j++) {
      center[j] += arr[j];
    }
  }
  for (unsigned j = 0; j < dimension_; j++) {
    center[j] /= nd_;
  }
  std::vector<Neighbor> tmp, pool;
  ep_ = rand() % nd_;  // random initialize navigating point
  get_neighbors(center, parameters, tmp, pool);
  ep_ = tmp[0].id;
}

void IndexNSG::sync_prune(unsigned q, std::vector<Neighbor> &pool,
                          const Parameters &parameter,
                          boost::dynamic_bitset<> &flags,
                          SimpleNeighbor *cut_graph_) {
  unsigned range = parameter.Get<unsigned>("R");
  unsigned maxc = parameter.Get<unsigned>("C");
  width = range;
  unsigned start = 0;

  for (unsigned nn = 0; nn < final_graph_[q].size(); nn++) {
    unsigned id = final_graph_[q][nn];
    if (flags[id]) continue;
	float dist;
	if(is_quan == false)
     dist =
        distance_->compare(data_ + dimension_ * (size_t)q,
                           data_ + dimension_ * (size_t)id, (unsigned)dimension_);
	else{
		dist = 0;
		for(int i = 0; i < length_; i++){
			dist +=  quan_book[i][*(data2_+length_* (size_t)q +i)][*(data2_+length_* (size_t)id +i)];
		}
	}					   
						   
    pool.push_back(Neighbor(id, dist, true, true));
  }

  std::sort(pool.begin(), pool.end());
  std::vector<Neighbor> result;
  if (pool[start].id == q) start++;
  result.push_back(pool[start]);

  while (result.size() < range && (++start) < pool.size() && start < maxc) {
    auto &p = pool[start];
    bool occlude = false;
    for (unsigned t = 0; t < result.size(); t++) {
      if (p.id == result[t].id) {
        occlude = true;
        break;
      }
      float djk; 
	  if(is_quan == false)
	    djk = distance_->compare(data_ + dimension_ * (size_t)result[t].id,
                                     data_ + dimension_ * (size_t)p.id,
                                     (unsigned)dimension_);
	  else{
		djk = 0;
	    for(int i = 0; i < length_; i++){
			djk +=  quan_book[i][*(data2_+length_*(size_t)result[t].id +i)][*(data2_+length_*(size_t)p.id +i)];
		}
	  }							 
									 
      if (djk < p.distance /* dik */) {
        occlude = true;
        break;
      }
    }
    if (!occlude) result.push_back(p);
  }

  SimpleNeighbor *des_pool = cut_graph_ + (size_t)q * (size_t)range;
  for (size_t t = 0; t < result.size(); t++) {
    des_pool[t].id = result[t].id;
    des_pool[t].distance = result[t].distance;
  }
  if (result.size() < range) {
    des_pool[result.size()].distance = -1;
  }
}

void IndexNSG::InterInsert(unsigned n, unsigned range,
                           std::vector<std::mutex> &locks,
                           SimpleNeighbor *cut_graph_) {
  SimpleNeighbor *src_pool = cut_graph_ + (size_t)n * (size_t)range;
  for (size_t i = 0; i < range; i++) {
    if (src_pool[i].distance == -1) break;

    SimpleNeighbor sn(n, src_pool[i].distance);
    size_t des = src_pool[i].id;
    SimpleNeighbor *des_pool = cut_graph_ + des * (size_t)range;

    std::vector<SimpleNeighbor> temp_pool;
    int dup = 0;
    {
      LockGuard guard(locks[des]);
      for (size_t j = 0; j < range; j++) {
        if (des_pool[j].distance == -1) break;
        if (n == des_pool[j].id) {
          dup = 1;
          break;
        }
        temp_pool.push_back(des_pool[j]);
      }
    }
    if (dup) continue;

    temp_pool.push_back(sn);
    if (temp_pool.size() > range) {
      std::vector<SimpleNeighbor> result;
      unsigned start = 0;
      std::sort(temp_pool.begin(), temp_pool.end());
      result.push_back(temp_pool[start]);
      while (result.size() < range && (++start) < temp_pool.size()) {
        auto &p = temp_pool[start];
        bool occlude = false;
        for (unsigned t = 0; t < result.size(); t++) {
          if (p.id == result[t].id) {
            occlude = true;
            break;
          }
		  /*
          float djk = distance_->compare(data_ + dimension_ * (size_t)result[t].id,
                                         data_ + dimension_ * (size_t)p.id,
                                         (unsigned)dimension_);
										 
		  */
          float djk; 
	      if(is_quan == false)
	        djk = distance_->compare(data_ + dimension_ * (size_t)result[t].id,
                                     data_ + dimension_ * (size_t)p.id,
                                     (unsigned)dimension_);
	      else{
		    djk = 0;
	        for(int i = 0; i < length_; i++){
			  djk +=  quan_book[i][*(data2_+length_*(size_t)result[t].id +i)][*(data2_+length_*(size_t)p.id +i)];
		    }
	      }   			  
          if (djk < p.distance /* dik */) {
            occlude = true;
            break;
          }
        }
        if (!occlude) result.push_back(p);
      }
      {
        LockGuard guard(locks[des]);
        for (unsigned t = 0; t < result.size(); t++) {
          des_pool[t] = result[t];
        }
      }
    } else {
      LockGuard guard(locks[des]);
      for (unsigned t = 0; t < range; t++) {
        if (des_pool[t].distance == -1) {
          des_pool[t] = sn;
          if (t + 1 < range) des_pool[t + 1].distance = -1;
          break;
        }
      }
    }
  }
}

void IndexNSG::Link(const Parameters &parameters, SimpleNeighbor *cut_graph_) {
  std::cout << " graph link" << std::endl;
  unsigned range = parameters.Get<unsigned>("R");
  std::vector<std::mutex> locks(nd_);
  int sdim_ = dimension_/length_;
  int index;
  float* arr = new float[dimension_];
//#pragma omp parallel
  {
    std::vector<Neighbor> pool, tmp;
    boost::dynamic_bitset<> flags{nd_, 0};
//#pragma omp for schedule(dynamic, 100)
    for (unsigned n = 0; n < nd_; ++n) {
      pool.clear();
      tmp.clear();
      flags.reset();
	  
      if(is_quan == false){
            get_neighbors(data_ + dimension_ * n, parameters, flags, tmp, pool);
      }
	  else{
		  for(int j = 0; j < length_; j++){
	        index = n * length_ + j;
	        for(int l = 0; l < sdim_; l++){
		       arr[j * sdim_ + l] = quantizer[j][ data2_[index] ][l];
	        }
          }
		  get_neighbors(arr, parameters, flags, tmp, pool);
	  }
	  
      sync_prune(n, pool, parameters, flags, cut_graph_);
    }
  }
//#pragma omp for schedule(dynamic, 100)
  for (unsigned n = 0; n < nd_; ++n) {
    InterInsert(n, range, locks, cut_graph_);
  }
}

void IndexNSG::Build(size_t n, float *data, const Parameters &parameters) {
  std::string nn_graph_path = parameters.Get<std::string>("nn_graph_path");
  unsigned range = parameters.Get<unsigned>("R");
  Load_nn_graph(nn_graph_path.c_str());
  data_ = data;
  init_graph(parameters);    
  SimpleNeighbor *cut_graph_ = new SimpleNeighbor[nd_ * (size_t)range];
  Link(parameters, cut_graph_);
  final_graph_.resize(nd_);

  for (size_t i = 0; i < nd_; i++) {
    SimpleNeighbor *pool = cut_graph_ + i * (size_t)range;
    unsigned pool_size = 0;
    for (unsigned j = 0; j < range; j++) {
      if (pool[j].distance == -1) break;
      pool_size = j;
    }
    pool_size++;
    final_graph_[i].resize(pool_size);
    for (unsigned j = 0; j < pool_size; j++) {
      final_graph_[i][j] = pool[j].id;
    }
  }

  tree_grow(parameters);

  unsigned max = 0, min = 1e6, avg = 0;
  for (size_t i = 0; i < nd_; i++) {
    auto size = final_graph_[i].size();
    max = max < size ? size : max;
    min = min > size ? size : min;
    avg += size;
  }
  avg /= 1.0 * nd_;
  printf("Degree Statistics: Max = %d, Min = %d, Avg = %d\n", max, min, avg);

  has_built = true;
  delete[] cut_graph_;
}

void IndexNSG::Build2(size_t n, unsigned char *data, const Parameters &parameters) {
  std::string nn_graph_path = parameters.Get<std::string>("nn_graph_path");
  unsigned range = parameters.Get<unsigned>("R");
  Load_nn_graph(nn_graph_path.c_str());   //unchanged; read k neighbor ids
  data2_ = data;
  init_graph2(parameters);  // modify get_neighbors
  SimpleNeighbor *cut_graph_ = new SimpleNeighbor[nd_ * (size_t)range];
  Link(parameters, cut_graph_);
  final_graph_.resize(nd_);   //nd_ has been set to count
  for (size_t i = 0; i < nd_; i++) {
    SimpleNeighbor *pool = cut_graph_ + i * (size_t)range;
    unsigned pool_size = 0;
    for (unsigned j = 0; j < range; j++) {
      if (pool[j].distance == -1) break;
      pool_size = j;
    }
    pool_size++;
    final_graph_[i].resize(pool_size);
    for (unsigned j = 0; j < pool_size; j++) {
      final_graph_[i][j] = pool[j].id;
    }
  }
  tree_grow(parameters);
  
  unsigned max = 0, min = 1e6, avg = 0;
  for (size_t i = 0; i < nd_; i++) {
    auto size = final_graph_[i].size();
    max = max < size ? size : max;
    min = min > size ? size : min;
    avg += size;
  }
  avg /= 1.0 * nd_;
  printf("Degree Statistics: Max = %d, Min = %d, Avg = %d\n", max, min, avg);

  has_built = true;
  delete[] cut_graph_;
}

void IndexNSG::Search(const float *query, float *x, size_t K,
                      const Parameters &parameters, unsigned *indices) {
  const unsigned L = parameters.Get<unsigned>("L_search");
  data_ = x;
  std::vector<Neighbor> retset(L + 1);
  std::vector<unsigned> init_ids(L);
  boost::dynamic_bitset<> flags{nd_, 0};
  // std::mt19937 rng(rand());
  // GenRandom(rng, init_ids.data(), L, (unsigned) nd_);

  unsigned tmp_l = 0;
  for (; tmp_l < L && tmp_l < final_graph_[ep_].size(); tmp_l++) {
    init_ids[tmp_l] = final_graph_[ep_][tmp_l];
    flags[init_ids[tmp_l]] = true;
  }

  while (tmp_l < L) {
    unsigned id = rand() % nd_;
    if (flags[id]) continue;
    flags[id] = true;
    init_ids[tmp_l] = id;
    tmp_l++;
  }

  for (unsigned i = 0; i < init_ids.size(); i++) {
    unsigned id = init_ids[i];
    float dist =
        distance_->compare(data_ + dimension_ * id, query, (unsigned)dimension_);
    retset[i] = Neighbor(id, dist, true, true);
    // flags[id] = true;
  }

  std::sort(retset.begin(), retset.begin() + L);
  int k = 0;
  while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      for (unsigned m = 0; m < final_graph_[n].size(); ++m) {
        unsigned id = final_graph_[n][m];
        if (flags[id]) continue;
        flags[id] = 1;
        float dist =
            distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);
        if (dist >= retset[L - 1].distance) continue;
        Neighbor nn(id, dist, true, true);
        int r = InsertIntoPool(retset.data(), L, nn);

        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
  for (size_t i = 0; i < K; i++) {
    indices[i] = retset[i].id;
  }
}

float querydistfunc_(unsigned char* id, int num, float** book){
	float dist = 0;
	for(int i = 0; i < num; i++){
	    dist += book[i][ id[i] ];
	}
	return dist;
}


void IndexNSG::SearchWithsingleGraph(float** book_, int X, unsigned* points, 
    unsigned* enter_point, unsigned int* trans, const Parameters &parameters, VisitedListPool *visited_list_pool_) {
  unsigned L = parameters.Get<unsigned>("L_search");
  
   VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
  
  
  float sum, dist;
  unsigned int currObj, min_id, fan_new;
  if(L < fan){
	fan_new = L;
  }
  else{
	fan_new = fan;
  }
  std::vector<Neighbor> retset(L + 1);
  
  size_t quan_len_ = quan_len[0];
  size_t quan_size_ = quan_size[0];
  char* quan_graph_ = quan_graph[0];
	
    	
    for(int i = 0; i< fan_new; i++){
       sum = querydistfunc_( (unsigned char*) quan_graph_ + quan_size_ * points[i], X, book_);
       retset[i] = Neighbor(points[i], sum, true, true);
	   visited_array[points[i]] = visited_array_tag;
    }
	
    std::sort(retset.begin(), retset.begin() + fan_new);

	int L2 = fan_new;
    int k = 0;

    while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

       _mm_prefetch(quan_graph_ + quan_size_ * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + quan_size_ * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;

      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + quan_size_ * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
        if (visited_array[id] == visited_array_tag) continue;
        visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + quan_size_ * id, X, book_);
		int r;
        if (L2 >= L && dist >= retset[L - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < L && dist >= retset[L2-1].distance){
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

   for (size_t i = 0; i < L; i++) {		
	    enter_point[i] = trans[retset[i].id];
   } 
   visited_list_pool_->releaseVisitedList(vl);
}




void IndexNSG::SearchWithquanGraph(float** book_, int X, unsigned* points, elem* init_point, 
        unsigned int* trans, const Parameters &parameters, char* fflag, int level, VisitedListPool *visited_list_pool_) {
  unsigned L = parameters.Get<unsigned>("L_search");
   
    VisitedList *vl = visited_list_pool_->getFreeVisitedList();
    vl_type *visited_array = vl->mass;
    vl_type visited_array_tag = vl->curV;  
  

  float sum, dist;
  unsigned int currObj, min_id, fan_new;
  if(L < fan){
	fan_new = L;
  }
  else{
	fan_new = fan;
  }
  std::vector<Neighbor> retset(L + 1);

  size_t quan_len_ = quan_len[level];
  size_t quan_size_ = quan_size[level];
  char* quan_graph_ = quan_graph[level];
	
    	
    for(int i = 0; i< fan_new; i++){
       sum = querydistfunc_( (unsigned char*) quan_graph_ + quan_size_ * points[i], X, book_);
       retset[i] = Neighbor(points[i], sum, true, true);
	   visited_array[points[i]] = visited_array_tag;
    }
	
    std::sort(retset.begin(), retset.begin() + fan_new);

	int L2 = fan_new;
    int k = 0;

    while (k < (int)L) {

    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + quan_size_ * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + quan_size_ * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;

      for (unsigned m = 0; m < MaxM; ++m)
	   _mm_prefetch(quan_graph_ + quan_size_ * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + quan_size_ * id, X, book_);
		int r;
        if (L2 >= L && dist >= retset[L - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < L && dist >= retset[L2-1].distance){
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

   for (size_t i = 0; i < L; i++) {	
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

void IndexNSG::SearchWithquanGraph2(float** book_, int X, elem* points, 
    unsigned* init_point, unsigned int* trans, const Parameters &parameters, VisitedListPool *visited_list_pool_) {
  unsigned L = parameters.Get<unsigned>("L_search");

      VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
   
  float dist, sum;
   std::vector<Neighbor> retset(L + 1);

        size_t quan_len_ = quan_len[0];
        size_t quan_size_ = quan_size[0];
        char* quan_graph_ = quan_graph[0];
	
	for(int i = 0; i < L; i++){
		if(points[i].flag == true){
		    sum = querydistfunc_( (unsigned char*) quan_graph_ + quan_size_ * points[i].id, X, book_);
            retset[i] = Neighbor(points[i].id, sum, true, true);
			visited_array[retset[i].id] = visited_array_tag;
		}
	    else{
			retset[i] = Neighbor(points[i].id, points[i].dist, false, false);
		}
	}
	
	std::sort(retset.begin(), retset.begin() + L);
	int L2= L;
	
    int k = 0;
    while (k < (int)L) {

    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + quan_size_ * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + quan_size_ * n + quan_len_);
      unsigned MaxM = *neighbors;

      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + quan_size_ * neighbors[m], _MM_HINT_T0);

      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + quan_size_ * id, X, book_);
		int r;
        if (L2 >= L && dist >= retset[L - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < L){
			r = InsertIntoPool(retset.data(), L2, nn);
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), L, nn);}
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }

  for (size_t i = 0; i < L; i++) {	
      if(retset[i].flag2 == false){
		  init_point[i] = retset[i].id;
         }
	  else{
		  init_point[i] = trans[retset[i].id];
         }
  } 
  visited_list_pool_->releaseVisitedList(vl); 
}

void IndexNSG::SearchWithquanGraph3(float** book_, int X, elem* points, elem* init_point, 
    unsigned int* trans, const Parameters &parameters, char* fflag, int level, VisitedListPool *visited_list_pool_) {
  unsigned L = parameters.Get<unsigned>("L_search");
  
     VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV;  

  float dist, sum;
   std::vector<Neighbor> retset(L + 1);
 
        size_t quan_len_ = quan_len[level];
        size_t quan_size_ = quan_size[level];
        char* quan_graph_ = quan_graph[level];
	
	
	for(int i = 0; i < L; i++){
		if(points[i].flag == true){
		    sum = querydistfunc_( (unsigned char*) quan_graph_ + quan_size_ * points[i].id, X, book_);
            retset[i] = Neighbor(points[i].id, sum, true, true);
			visited_array[retset[i].id] = visited_array_tag;
		}
	    else{
			retset[i] = Neighbor(points[i].id, points[i].dist, false, false);
		}
	}
	
	std::sort(retset.begin(), retset.begin() + L);
	int L2= L;
	
    int k = 0;
    while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(quan_graph_ + quan_size_ * n + quan_len_, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(quan_graph_ + quan_size_ * n + quan_len_);
      unsigned MaxM = *neighbors;
      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	     _mm_prefetch(quan_graph_ + quan_size_ * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
		float dist = querydistfunc_((unsigned char*) quan_graph_ + quan_size_ * id, X, book_);
		int r;
        if (L2 >= L && dist >= retset[L - 1].distance) continue;

		Neighbor nn(id, dist, true, true);
		if(L2 < L){
			r = InsertIntoPool(retset.data(), L2, nn);
			L2++;
		}
		else{ r = InsertIntoPool(retset.data(), L, nn);}
		
        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
   for (size_t i = 0; i < L; i++) {
    
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


void IndexNSG::SearchWithOptGraph(const float *query, size_t K, const Parameters &parameters, 
        unsigned *indices, unsigned *init_point, VisitedListPool *visited_list_pool_) {

  unsigned L = parameters.Get<unsigned>("L_search");
  DistanceFastL2 *dist_fast = (DistanceFastL2 *)distance_;

  std::vector<Neighbor> retset(L + 1);
  std::vector<unsigned> init_ids(L);

     VisitedList *vl = visited_list_pool_->getFreeVisitedList();
   vl_type *visited_array = vl->mass;
   vl_type visited_array_tag = vl->curV; 
  

  for (unsigned i = 0; i < L; i++) {
	float* data = (float *)(opt_graph_ + node_size * init_point[i]);
	float norm = *data;
    data++;
    float dist = dist_fast->compare(query, data, norm, (unsigned)dimension_);
    retset[i] = Neighbor(init_point[i], dist, true, true);
    visited_array[init_point[i]] = visited_array_tag;
  }
  
  std::sort(retset.begin(), retset.begin() + L);
  int k = 0;
  while (k < (int)L) {
    int nk = L;

    if (retset[k].flag) {
      retset[k].flag = false;
      unsigned n = retset[k].id;

      _mm_prefetch(opt_graph_ + node_size * n + data_len, _MM_HINT_T0);
      unsigned *neighbors = (unsigned *)(opt_graph_ + node_size * n + data_len);
      unsigned MaxM = *neighbors;
      neighbors++;
      for (unsigned m = 0; m < MaxM; ++m)
	    _mm_prefetch(opt_graph_ + node_size * neighbors[m], _MM_HINT_T0);
      for (unsigned m = 0; m < MaxM; ++m) {
        unsigned id = neighbors[m];
		
		if (visited_array[id] == visited_array_tag) continue;
            visited_array[id] = visited_array_tag;
		
        float *data = (float *)(opt_graph_ + node_size * id);
        float norm = *data;
        data++;
        float dist = dist_fast->compare(query, data, norm, (unsigned)dimension_);
        if (dist >= retset[L - 1].distance) continue;
        Neighbor nn(id, dist, true, true);
        int r = InsertIntoPool(retset.data(), L, nn);

        if (r < nk) nk = r;
      }
    }
    if (nk <= k)
      k = nk;
    else
      ++k;
  }
  for (size_t i = 0; i < K; i++) {
    indices[i] = retset[i].id;
  }  
  visited_list_pool_->releaseVisitedList(vl);
}

void IndexNSG::OptimizeGraph(float *data) {  // use after build or load

  unsigned cur_width = width;
  if(cur_width > max_width_)
	  cur_width = max_width_;

  if(nd_ > 30000000 && cur_width > 50)
	  cur_width = 50;
  
  data_ = data;
  data_len = (dimension_ + 1) * sizeof(float);
  neighbor_len = (cur_width + 1) * sizeof(unsigned);
  node_size = data_len + neighbor_len;
  opt_graph_ = (char *)malloc(node_size * nd_);
  DistanceFastL2 *dist_fast = (DistanceFastL2 *)distance_;
  for (unsigned i = 0; i < nd_; i++) {
    char *cur_node_offset = opt_graph_ + i * node_size;
    float* tmp3 = data_ + i * dimension_;
    float cur_norm = dist_fast->norm(data_ + i * dimension_, dimension_);

    std::memcpy(cur_node_offset, &cur_norm, sizeof(float));
    std::memcpy(cur_node_offset + sizeof(float), data_ + i * dimension_,
                data_len - sizeof(float));
           
    cur_node_offset += data_len;
    unsigned k = final_graph_[i].size();
	
	if(k > cur_width) k = cur_width;
    std::memcpy(cur_node_offset, &k, sizeof(unsigned));
    std::memcpy(cur_node_offset + sizeof(unsigned), final_graph_[i].data(),
                k * sizeof(unsigned));
    std::vector<unsigned>().swap(final_graph_[i]);
  }
  CompactGraph().swap(final_graph_);
}

void IndexNSG::OptimizeQuan(unsigned char *quan_, int level, int num, int max_level) {  // use after build or load
  int length = pow(2, max_level - level + OFF ); 
  //data_ = data;
  quan_len[level] = length;
  
  unsigned cur_width = width;
  if(cur_width > max_width_)
	  cur_width = max_width_;
	  
  if(nd_ > 30000000 && cur_width > 50)
	  cur_width = 50;  
  
  neighbor_len = (cur_width + 1) * sizeof(unsigned);
  quan_size[level] = quan_len[level] + neighbor_len;
  quan_graph[level] = (char *)malloc(quan_size[level] * num);
  //Distance *dist_fast = (DistanceFastL2 *)distance_;
  for (unsigned i = 0; i < num; i++) {
    char *cur_node_offset = quan_graph[level] + i * quan_size[level];
    std::memcpy(cur_node_offset, quan_ + i * length,
                quan_len[level]);

    cur_node_offset += quan_len[level];
    unsigned k = final_graph_[i].size();
	
	if(k > cur_width) k = cur_width;
    std::memcpy(cur_node_offset, &k, sizeof(unsigned));
    std::memcpy(cur_node_offset + sizeof(unsigned), final_graph_[i].data(),
                k * sizeof(unsigned));
    std::vector<unsigned>().swap(final_graph_[i]);
  }
  CompactGraph().swap(final_graph_);
}

void IndexNSG::DFS(boost::dynamic_bitset<> &flag, unsigned root, unsigned &cnt) {
  unsigned tmp = root;
  std::stack<unsigned> s;
  s.push(root);
  if (!flag[root]) cnt++;
  flag[root] = true;
  while (!s.empty()) {
    unsigned next = nd_ + 1;
    for (unsigned i = 0; i < final_graph_[tmp].size(); i++) {
      if (flag[final_graph_[tmp][i]] == false) {
        next = final_graph_[tmp][i];
        break;
      }
    }
    // std::cout << next <<":"<<cnt <<":"<<tmp <<":"<<s.size()<< '\n';
    if (next == (nd_ + 1)) {
      s.pop();
      if (s.empty()) break;
      tmp = s.top();
      continue;
    }
    tmp = next;
    flag[tmp] = true;
    s.push(tmp);
    cnt++;
  }
}

void IndexNSG::findroot(boost::dynamic_bitset<> &flag, unsigned &root,
                        const Parameters &parameter) {
  int sdim_ = dimension_/length_;
  int index;
  float* arr = new float[dimension_];
  
  unsigned id = nd_;
  for (unsigned i = 0; i < nd_; i++) {
    if (flag[i] == false) {
      id = i;
      break;
    }
  }

  if (id == nd_) return;  // No Unlinked Node

  std::vector<Neighbor> tmp, pool;

  if(is_quan == false)
     get_neighbors(data_ + dimension_ * id, parameter, tmp, pool);
  else{
    for(int j = 0; j <length_; j++){
      index = id * length_ + j;
      for(int l = 0; l < sdim_; l++){
	arr[j * sdim_ + l] = quantizer[j][data2_[index]][l];
      }
    }
    get_neighbors(arr,parameter,tmp,pool);
  }

  
  std::sort(pool.begin(), pool.end());

  unsigned found = 0;
  for (unsigned i = 0; i < pool.size(); i++) {
    if (flag[pool[i].id]) {
      // std::cout << pool[i].id << '\n';
      root = pool[i].id;
      found = 1;
      break;
    }
  }
  if (found == 0) {
    while (true) {
      unsigned rid = rand() % nd_;
      if (flag[rid]) {
        root = rid;
        break;
      }
    }
  }
  final_graph_[root].push_back(id);
}
void IndexNSG::tree_grow(const Parameters &parameter) {
  printf("DFS...");
  unsigned root = ep_;
  boost::dynamic_bitset<> flags{nd_, 0};
  unsigned unlinked_cnt = 0;
  while (unlinked_cnt < nd_) {
    DFS(flags, root, unlinked_cnt);
    //std::cout << unlinked_cnt << '\n';
    if (unlinked_cnt >= nd_) break;
    findroot(flags, root, parameter);
    // std::cout << "new root"<<":"<<root << '\n';
  }
  for (size_t i = 0; i < nd_; ++i) {
    if (final_graph_[i].size() > width) {
      width = final_graph_[i].size();
    }
  }
}

void IndexNSG::Init(size_t dimension, size_t n, Metric m, Index *initializer, 
        const bool flag, float*** book, int l, float*** quanti){
		distance_ = new DistanceL2();	
        dimension_ = dimension; 
        nd_ = n;
        has_built = false; 
        is_quan = flag; 
        quan_book = book;
        length_ = l; 
        quantizer = quanti;			
}

}
