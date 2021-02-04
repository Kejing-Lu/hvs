#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "hnswlib/hnswlib.h"
#include <ctime>
#include <cstdlib>

#include <unordered_set>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace hnswlib;
using namespace cv;

class StopW {
    std::chrono::steady_clock::time_point time_begin;
public:
    StopW() {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro() {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset() {
        time_begin = std::chrono::steady_clock::now();
    }

};


#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))

#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

struct elem_{
	unsigned char* index;
	int id;
	float dist;
};

struct k_elem{
	int id;
	float dist;
};

void kmeans_(
    float** train,
    float** result,
    int n,
    int d
){
        unsigned seed;
	float dist, min_dist;
	int min_id;
	int index;
	int round= 20;
	float** temp_quan = new float*[L];  //temp array for 
	int* tol_quan = new int[L];  //the number of objects on each center 
	for(int i = 0; i < L; i++)
		temp_quan[i] = new float[d];
	
	bool* flag= new bool[n];  // whether to choose it as a center
    for(int i = 0; i < n; i++)
		flag[i] = false;

    double rand_;
	for(int i = 0; i < L; i++){  //generate L intial centers randomly
		 rand_ = double(i) / L;
		 index = int( (n-1) * rand_);
		if(index < 0 || index >= n){
			printf("random_generator error\n");
			exit(0);
		}
				
		if(flag[index] == false){
			for(int j = 0; j < d; j++)
			    result[i][j] = train[index][j];
			flag[index] = true;
		}
		else{
            i--;           //re-generate center
	    }	
	}
    float temp_sum;
	
	for(int l = 0; l < round; l++){
	    temp_sum = 0;
		for(int i = 0; i < L; i++)
			for(int j = 0; j < d; j++)
		        temp_quan[i][j] = 0;
			
		for(int i = 0; i < L; i++)
            tol_quan[i] = 0;			
		
		for(int i = 0; i < n; i++){
			for(int ii = 0; ii < L; ii++){
			    dist = 0;
			    for(int j =0; j < d; j++)
				    dist += (train[i][j] - result[ii][j]) * (train[i][j] - result[ii][j]);
				
				if(ii == 0){
				    min_dist = dist;
				    min_id = 0;
                    continue;					
				}
				
			    if (dist < min_dist){
				    min_dist = dist;
				    min_id = ii;
			    }
		    }
			for(int j = 0; j < d; j++){
			    temp_quan[min_id][j] += train[i][j];
			}
			temp_sum += min_dist;
			tol_quan[min_id]++;
		}

		for(int i = 0; i < L; i++){
			if(tol_quan[i] == 0) continue;
			else{
				for(int j = 0; j < d; j++)
					result[i][j] = temp_quan[i][j] / tol_quan[i];
			}
		}
		
	}

}

void kmeans0_(
    float** train,
    float** result,
    int n,
    int d,
	int cen1
){
        unsigned seed;
	float dist, min_dist;
	int min_id;
	int index;
	int round= 20;
	float** temp_quan = new float*[cen1];  //temp array for 
	int* tol_quan = new int[cen1];  //the number of objects on each center 
	for(int i = 0; i < cen1; i++)
		temp_quan[i] = new float[d];
	
	bool* flag= new bool[n];  // whether to choose it as a center
    for(int i = 0; i < n; i++)
		flag[i] = false;

    double rand_;
	for(int i = 0; i < cen1; i++){  //generate L intial centers randomly
		 rand_ = double(i) / cen1;
		 index = int( (n-1) * rand_);
		if(index < 0 || index >= n){
			printf("random_generator error\n");
			exit(0);
		}
				
		if(flag[index] == false){
			for(int j = 0; j < d; j++)
			    result[i][j] = train[index][j];
			flag[index] = true;
		}
		else{
            i--;           //re-generate center
	    }	
	}
    float temp_sum;

	for(int l = 0; l < round; l++){
	    temp_sum = 0;
		for(int i = 0; i < cen1; i++)
			for(int j = 0; j < d; j++)
		        temp_quan[i][j] = 0;
			
		for(int i = 0; i < cen1; i++)
            tol_quan[i] = 0;			
		
		for(int i = 0; i < n; i++){
			for(int ii = 0; ii < cen1; ii++){
			    dist = 0;
			    for(int j =0; j < d; j++)
				    dist += (train[i][j] - result[ii][j]) * (train[i][j] - result[ii][j]);
				
				if(ii == 0){
				    min_dist = dist;
				    min_id = 0;
                    continue;					
				}
				
			    if (dist < min_dist){
				    min_dist = dist;
				    min_id = ii;
			    }
		    }
			for(int j = 0; j < d; j++){
			    temp_quan[min_id][j] += train[i][j];
			}
			temp_sum += min_dist;
			tol_quan[min_id]++;
		}
		for(int i = 0; i < cen1; i++){
			if(tol_quan[i] == 0) continue;
			else{
				for(int j = 0; j < d; j++)
					result[i][j] = temp_quan[i][j] / tol_quan[i];
			}
		}
		
	}
}

int compare_int(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}	

int compfloats(						// compare two real values
	float v1,							// 1st value (type float)
	float v2)							// 2nd value (type float)
{
	const float FLOATZERO = 1e-6F;
	if (v1 - v2 < -FLOATZERO) return -1;   //?floatzero位置
	else if (v1 - v2 > FLOATZERO) return 1;
	else return 0;
}

int QsortComp(				// compare function for qsort
	const void* e1,						// 1st element
	const void* e2)						// 2nd element
{
	int ret = 0;
	k_elem* value1 = (k_elem *) e1;
	k_elem* value2 = (k_elem *) e2;
	if (value1->dist < value2->dist) {
		ret = -1;
	} else if (value1->dist > value2->dist) {
		ret = 1;
	} else {
		if (value1->id < value2->id) ret = -1;
		else if (value1->id > value2->id) ret = 1;
	}
	return ret;
}
	
void sub_kmeans_(
    float** train, //L * L * d
    float** result,  //L * d
    unsigned char* merge,   //length 2L
    int d,
	float** train1,
	float** train2,
	float** train0) //size_n * d
{

    unsigned seed;
    int n = L*L;
    int index;
	int round= 5;
	int first_num = 10;
	double rand_;
	
	int* weight = new int[L*L];

	for(int i = 0; i < L*L; i++){
        weight[i] = 0;
	}
	
	int min_id1 = 0;
	int min_id2 = 0;
	float min_sum1, min_sum2;	
	float sum1 = 0;
	float sum2 = 0;	
	int half_d = d / 2;
	float** array1 = new float* [L];
	float** array2 = new float* [L];
	for(int i = 0; i < L; i++){
		array1[i] = new float[half_d];
		array2[i] = new float[half_d];
	}	
	for(int i = 0; i < L; i++){
		for(int j = 0; j < half_d; j++){
		    array1[i][j] = train[i*L+i][j];
		    array2[i][j] = train[i*L+i][half_d + j];
		}
	}	
	for(int i = 0; i < size_n; i++){
	    for(int j = 0; j < L; j++){
	          sum1 = 0;
	              sum2 = 0;
			for(int s = 0; s < half_d; s++){
			    sum1 += (train1[i][s] - array1[j][s]) * (train1[i][s] - array1[j][s]);
			    sum2 += (train2[i][s] - array2[j][s]) * (train2[i][s] - array2[j][s]);
			}
			if(j == 0) {min_id1 = 0; min_sum1 = sum1; min_id2 = 0; min_sum2 = sum2;}
			else{
				if(sum1 < min_sum1) {min_id1 = j; min_sum1 = sum1;}
                if(sum2 < min_sum2) {min_id2 = j; min_sum2 = sum2;}
			}
		}
		weight[min_id1 * L + min_id2]++;
	}
	
	float** temp_quan = new float*[L];
	int** count = new int*[L];
	int* tol_quan = new int[L];
	int* w_quan = new int[L];
	
	for(int i = 0; i < L; i++){
		temp_quan[i] = new float[d];
	    count[i] = new int[L*L];
	}
	
	bool * flag= new bool[n];
    for(int i = 0; i < n; i++)
		flag[i] = false;
	
	k_elem sort_array_[L*L];
	int min_id, id1, id2; float min_dist, dist;
	int pointer[L];

	float tmp_dist;
	
	int* ord= new int[size_n];
	
	kmeans0_(train0, result, size_n, d, L);
	
	for(int i = 0; i < size_n; i++){
		for(int j = 0; j < L; j++){
			dist = 0;
			for(int ii = 0; ii < d; ii++)
				dist += (train0[i][ii] - result[j][ii]) * (train0[i][ii] - result[j][ii]);

            if(j == 0){
				min_id = 0;
				min_dist = dist;
			}
            else{
				if(dist < min_dist){
					min_id = j;
					min_dist = dist;
				}		
			}			
		}
		ord[i] = min_id;
	}
	//-------------------moving points------------------------
    for(int ii = 0; ii < L; ii++){
		for(int i = 0; i < n; i++){
			dist = 0;
			for(int j = 0; j < d; j++){
				dist += (train[i][j] - result[ii][j]) * (train[i][j] - result[ii][j]);
			}
			sort_array_[i].id = i;
		    sort_array_[i].dist = dist;
		}
	    qsort(sort_array_, n, sizeof(k_elem), QsortComp);
			
		for(int i = 0; i < first_num; i++){
		    id1 = sort_array_[i].id;
			dist = 0;
			for(int j = 0; j < size_n; j++){
					
			    if(ord[j] != ii) continue;
						
				for(int s = 0; s < d; s++)
					dist += (train[id1][s] - train0[j][s]) * (train[id1][s] - train0[j][s]);
			}
			if(i == 0){
				min_id = id1;
				min_dist = dist;
				continue;
			}
			else{
				if(dist < min_dist){
					min_id = id1;
					min_dist = dist;
				}
			}
				
		}
		tmp_dist += min_dist;
		pointer[ii] = min_id;
		for(int i = 0; i < d; i++)
			result[ii][i] = train[min_id][i];		
	}	

	//---------------------------------------------------------	
	qsort(pointer, L, sizeof(int), compare_int);
	//----------------important sampling-------
	
	float* mean_ = new float[d];
	float* prob = new float[n];
    float* prob2 = new float[n];	
	
	for(int i = 0; i < d; i++) 
		mean_[i] = 0;
	
	for(int i = 0; i < n; i++)
		prob[i] = 0;
	
	int a =0;
	for(int i = 0; i < n; i++){
		if(weight[i] == 0) continue;
		for(int j = 0; j < d; j++){
			mean_[j] +=  train[i][j];
		}
		a++;
	}
	for(int i = 0; i < d; i++)
		mean_[i] = mean_[i] / a;
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < d; j++){
		    prob[i] += weight[i] * ( mean_[j] - train[i][j]) * (mean_[j] - train[i][j]);
        }		
	}	
	float sum_prob = 0;
	
	for(int i = 0; i < n; i++)
	    sum_prob += prob[i];
	
	for(int i = 0; i < n; i++)
		prob[i] = 1 / n / 2.0f +  prob[i] / sum_prob / 2.0f;
	
	prob2[0] = prob[0];
	for(int i = 1; i < n; i++){
		prob2[i] = prob2[i-1] + prob[i];
	}
	
	for(int i = 1; i < L; i++){
		for(int j = 0; j <= i-1; j++){
			if(pointer[i] == pointer[j]){
	            rand_ = rand() / double(RAND_MAX);
			    for(int l = 0; l < n; l++){
				    if(rand_ <= prob2[l]){
					    pointer[i] = l;
			            i--;
						break;
					}
				}
                break;				
			}	
	    }
	}
	qsort(pointer, L, sizeof(int), compare_int);
	//------------------------------------------	
	for(int i = 0; i < L; i++){
		merge[2*i] = pointer[i] / L;
		merge[2*i+1] = pointer[i] % L;
	}
}	
	
void sub_kmeans0_(
    int ii, 
    float*** quantizer,   
    unsigned char* merge,
	float** train0,
	int min_dim,
	int max_dim,
	float** quantizer0)
{
    int index;
	float sum = 0;
	float min_sum;
    int min_id;	
	bool flag;
	
	float** result = new float* [cen];
	for(int i = 0; i < cen; i++)
		result[i] = new float[max_dim];
	
	float* result0 = new float[min_dim];
	
	unsigned char** obj = new unsigned char* [cen];
	for(int i = 0; i < cen; i++)
		obj[i] = new unsigned char[nnum];
	
	kmeans0_(train0, result, size_n, max_dim, cen);
		
    for(int i = 0; i < cen; i++){
		for(int j = 0; j < nnum; j++){			
			for(int l = 0; l < min_dim; l++){
			    result0[l] = result[i][j* min_dim +l];	
			}
            index = ii * nnum + j;
            for(int l = 0; l < L; l++){
				sum = 0;
				for(int jj = 0; jj < min_dim; jj++){
					sum += (quantizer[index][l][jj] - result0[jj]) * (quantizer[index][l][jj] - result0[jj]);
				}
				if(l == 0) {min_id = 0; min_sum = sum;}
				else{
					if(sum < min_sum){
						min_id = l;
						min_sum = sum;
					}
				}
			}
		    obj[i][j] = min_id;
						
			for(int l = 0; l < min_dim; l++)
                quantizer0[i][ j * min_dim + l]	= quantizer[index][min_id][l];		
		}
	}
	
	for(int i = 0; i < cen; i++){
		for(int j = i+1; j < cen; j++ ){
			flag = true;
			for(int l = 0; l < nnum; l++){
				if(obj[i][l] != obj[j][l]){
					flag = false;
					break;
				}
			}
			if(flag == true){
			//	printf("same obj\n");
			}
		}
	}
	
	for(int i = 0; i < cen; i++){
		for(int j = 0; j < nnum; j++){
			merge[i * nnum + j] = obj[i][j];
		}
	}
		
}		
	
void sift_test1B(
    char* path_data,
    int vecsize_,
	int vecdim_,
	int level_,
	float delta_) {

	int efConstruction = 200;
	int M = 16;
    size_t vecsize = vecsize_;
    size_t vecdim = vecdim_;
	float delta = delta_;
	int max_level_ = level_;
	int report_every = 50000;

    printf("vecsize = %d\n",vecsize);
    printf("vecdim = %d\n", vecdim);
    printf("level = %d\n", max_level_);
    printf("delta = %.5f\n", delta);
    
	int dim_;
	int max_num = pow(2, max_level_ + OFF); //OFF = 2 from 64 bits =3 from 128 bits
	int remainder = vecdim % max_num;
	int ratio = vecdim / max_num;
	if(remainder == 0){dim_ = vecdim;}
    else{ dim_ = (ratio + 1)* max_num; }	
		
    float *massb = new float[vecdim];
	float *mass = new float[vecdim];
	
    ifstream input(path_data, ios::binary);
	
    int maxk = 100;
	const float MAXREAL   = 3.402823466e+38F;
	srand(time(NULL));
	
    int in = 0;
    L2Space l2space(vecdim);           
	
	int* length = new int[max_level_];
	int* dim = new int[max_level_];
	for (int i = 0; i < max_level_; i++){
		length[i] = pow(2, max_level_ -i + OFF );  //check
		dim[i] = dim_ / length[i];
	}
	
	 //--------------------some arrays----------------------------------------
	unsigned char*** obj_quantizer = new unsigned char**[max_level_];
	for(int i = 0; i < max_level_; i++){
		obj_quantizer[i] = new unsigned char* [vecsize];
	}
	for(int i = 0; i < max_level_; i++){
		for(int j = 0; j < vecsize; j++)
			obj_quantizer[i][j] = new unsigned char[length[i]];
	}
		
	int** init_obj = new int* [max_level_]; 
    for(int i = 0; i < max_level_; i++){	
        init_obj[i] = new int[vecsize];
	}	 
	 
	unsigned char*** merge = new unsigned char** [max_level_];
	for(int i  = 0; i < max_level_; i++){
		merge[i] = new unsigned char* [length[i]];
	}	
	
	for (int i = 0; i < max_level_; i++){
		for(int ii = 0; ii < length[i]; ii++){
			merge[i][ii] = new unsigned char[2*L];
		}			
	}

	
	int max_dim = dim_ / min_book;
	
	float*** quantizer0 = new float** [min_book];  //for checking start book
	for (int i = 0; i <  min_book; i++){
		quantizer0[i] = new float* [cen];
	}
	for (int i = 0; i < min_book; i++){
		for(int j = 0; j < cen; j++){
			quantizer0[i][j] = new float[max_dim];
		}
	}
	

    HierarchicalNSW<float> *appr_alg;

    //cout << "Building HNSW index:\n";
    appr_alg = new HierarchicalNSW<float>(level_, &l2space, vecsize, M, efConstruction);
        	 

	float**** quantizer = new float*** [max_level_];
	float**** temp_quan = new float*** [max_level_];
		
	for (int i = 0; i < max_level_; i++){
        quantizer[i] = new float** [length[i]];
	    temp_quan[i] = new float** [length[i]];
			
		for(int ii = 0; ii < length[i]; ii++){
			quantizer[i][ii] = new float* [L];
			temp_quan[i][ii] = new float* [L*L];
		}
		for(int ii = 0; ii < length[i]; ii++){
			for(int j = 0; j < L; j++)
				quantizer[i][ii][j] = new float[dim[i]];
			for(int j = 0; j < L*L; j++)
				temp_quan[i][ii][j] = new float[dim[i]];
		}
	}

	int max_book = length[0];
	int min_dim = dim[0];
	
	unsigned int Tol = cen * cen * cen * cen;
	unsigned int** start_book = new unsigned int* [Tol]; 
    for(int i = 0; i < Tol; i++){	
        start_book[i] = new unsigned int[fan];
	}
		
	float*** train = new float** [max_book];
	for(int i = 0; i < max_book; i++)
		train[i] = new float* [size_n];
	for(int i = 0; i < max_book; i++)
		for(int j = 0; j < size_n; j++)
			train[i][j] = new float[min_dim];
		
	float** train_org = new float* [size_n]; //alt
    for(int i = 0; i < size_n; i++)
        train_org[i] = new float[dim_];		
	
    int interval = vecsize / size_n;
    int ind = -interval;
    if(interval < 1) {interval = 1;}

	input.seekg(0, std::ios::beg);	
	for(int i = 0; i < size_n; i++){
		ind += interval;
		if(i > 0) {input.seekg( (4 + 4 * vecdim)* (interval-1), std::ios::cur);}  
		input.read((char*)(&in), 4);
		input.read((char*)(massb), vecdim * 4);
		              
		for(int ii = 0; ii < dim_; ii++){
			if(ii < vecdim)
			{train_org[i][ii] = massb[ii];}
			else {train_org[i][ii] = 0;}
		}
                
	}

	input.seekg(0, std::ios::beg);	
	float** data = new float* [vecsize];   
	for(int i = 0; i < vecsize; i++)
		data[i] = new float[dim_];
	
	for(int i = 0; i < vecsize; i++){
		input.seekg(4, std::ios::cur);
		input.read((char*)(massb), vecdim * 4);
        for(int j= 0; j < dim_; j++){	
            if(j < vecdim) {data[i][j] = massb[j];}
            else{data[i][j] = 0;}
		}    
	}
	
	StopW stopw_full = StopW();
	//----------build base layer------------------------------------
	float* density = new float[vecsize];
	for(int i = 0; i < vecsize; i++) density[i] = 0;
	
	char tmp_path[1024];
	
	std::string data_path;
	data_path = "data_graph.fvecs";
	std::string quan_path[max_level_]; 
	
    for(int i = 0; i < max_level_; i++){
        sprintf(tmp_path, "quan_graph%d.fvecs", i);
		quan_path[i] = tmp_path;
			memset(tmp_path,'\0',sizeof(tmp_path));
	}
    unsigned int* neighbor = new unsigned int[KK-1];  //100NN by default
	unsigned int* neighbor2 = new unsigned int[KK];
	
	appr_alg -> setEf(200);
	
	input.seekg(0, std::ios::beg);

    input.read((char *) &in, 4);
    input.read((char *) massb, 4 * vecdim);

    for (int j = 0; j < vecdim; j++) {
        mass[j] = massb[j] * (1.0f);
    }
    appr_alg->addPoint((void *) (massb), (size_t) 0, NULL, -1, true);
					      
    int j1 = 0;
	
    for (int i = 1; i < vecsize; i++) {
	    if(i % report_every == 0)
		    printf("build HNSW index...%.2f %% complete\n", i * 100 / float(vecsize) );
            
		int j2=0;

        input.read((char *) &in, 4);
        input.read((char *) massb, 4 * vecdim);
        for (int j = 0; j < vecdim; j++) {
            mass[j] = massb[j];
        }
        j1++;
        j2=j1;
				
		appr_alg->addPoint((void *) (mass), (size_t) j2, NULL, -1, false);
    }
    
	//-----------------------build knn graph-----------------------------------------

	std::string tmp_string0;
    tmp_string0 = data_path;
		
	ofstream output_knn(tmp_string0, ios::binary);
		
	int count2_;
	    
	for(int l = 0; l < vecsize; l++){
		if(l > 0 && l % report_every == 0)
		    printf("build kNN graph...%.2f %% complete\n", l * 100 / float(vecsize) );
		
		appr_alg->connect((void *) (data[l]), nullptr, neighbor2, -1, true);				
		count2_ = 0;
		output_knn.write((char*) &l, 4);
		for(int i = 0; i < KK-1; i++){  //remove the same id
			if(l == neighbor2[i]) {count2_++;}
				neighbor[i]  = neighbor2[count2_];
				count2_++;		
			}
			output_knn.write((char*) neighbor, 4 * (KK-1));
			//----------------distity estimation---------------------				
			for(int i = 0; i < dim_; i++){
			    density[l] += (data[l][i] - data[neighbor[0]][i]) * (data[l][i] - data[neighbor[0]][i]);					
			}
			if(density[l] <= 0) { continue;}
			density[l] = -1 * vecdim * log(density[l]) / 2 / log(2);
            //---------------------------------------------------------			
	}
	output_knn.close();

	appr_alg-> deleteLinklist(-1);
	
	//--------------matrix rotation-----------------------------------------
	CvMat* M_X = cvCreateMat(dim_, size_n, CV_32FC1);
	CvMat* M_X2 = cvCreateMat(dim_, size_n, CV_32FC1);
	
	CvMat* M_Y = cvCreateMat(dim_, size_n, CV_32FC1);
	CvMat* M_YT = cvCreateMat(size_n, dim_, CV_32FC1);
	CvMat* M_R = cvCreateMat(dim_, dim_, CV_32FC1);
	CvMat* M_RC = cvCreateMat(dim_, dim_, CV_32FC1);
	CvMat* M_RX = cvCreateMat(dim_, size_n, CV_32FC1);
	CvMat* M_RT = cvCreateMat(dim_, dim_, CV_32FC1);
	
	CvMat* ABt   = cvCreateMat(dim_, dim_, CV_32FC1);
	CvMat* ABt_D = cvCreateMat(dim_, dim_, CV_32FC1);
	CvMat* ABt_U = cvCreateMat(dim_, dim_, CV_32FC1);
	CvMat* ABt_VT = cvCreateMat(dim_, dim_, CV_32FC1);
	
	for (int i = 0; i < M_X->rows; i++) {
        for (int j = 0; j < M_X->cols; j++) {
            cvmSet (M_X, i, j, train_org[j][i]);
        }
    }

    printf("PCA...\n");
	
	CvMat* pMean = cvCreateMat(1, dim_, CV_32FC1);
	CvMat* pEigVals = cvCreateMat(1, dim_, CV_32FC1);
	CvMat* pEigVecs = cvCreateMat(dim_, dim_, CV_32FC1);
	cvCalcPCA(M_X, pMean, pEigVals, pEigVecs, CV_PCA_DATA_AS_COL);
	CvMat* PCA_R = cvCreateMat(dim_, dim_, CV_32FC1);
	
	
	int* ord = new int[max_book];
	int* ord2 = new int[dim_];
	k_elem* prod = new k_elem[max_book];
	
	float ssum;

	for(int i = 0; i < dim_; i++){
		if(i < max_book){
			prod[i].dist = cvmGet (pEigVals, 0, i);
			prod[i].id = i;
			ord[i] = 1;
			ord2[i] = i * min_dim; 
		}
		
		if(i >= max_book){
			
			ssum = cvmGet (pEigVals, 0, i);
			qsort(prod, max_book,sizeof(k_elem), QsortComp);
					
			for(int j = 0; j < max_book; j++){
				if( ord[ prod[j].id ]  < min_dim){
					ord2[i] = prod[j].id * min_dim +  ord[ prod[j].id ];
                    ord[ prod[j].id ]++;
                    prod[j].dist *= ssum;
                    break;					
				}
			}
		}	
	}
	
	float* pca_arr = new float[dim_];
	
	for(int i = 0; i < dim_; i++){
	    for(int j = 0; j < dim_; j++){
	      pca_arr[j] = cvmGet(pEigVecs, i, j);
		}
		ssum = 0;
		for(int j = 0; j < dim_; j++){
			ssum += pca_arr[j] * pca_arr[j];
		}
		for(int j = 0; j < dim_; j++){
		    cvmSet(PCA_R, ord2[i], j, pca_arr[j]/sqrt(ssum));
		}
	}
	
	cvMatMul( PCA_R, M_X, M_X2 );
    for(int i = 0; i < size_n; i++){
		for(int j = 0; j < dim_; j++){
			train_org[i][j] = cvmGet (M_X2, j, i);
		}
	}

	for(int i = 0; i < size_n; i++){
		for(int ii = 0; ii < max_book; ii++){ 
			for(int jj = 0; jj < min_dim; jj++){
				train[ii][i][jj] = train_org[i][ii * min_dim +jj];
			}
		}
	}	
		
	for(int i = 0; i < max_book; i++){ 
		kmeans_(train[i],quantizer[0][i], size_n, min_dim);
	}
    
	
	float*** vec = new float** [max_book];
	for(int i = 0; i < max_book; i++){
		vec[i] = new float* [L];
	}
  
    for(int i = 0; i < max_book; i++){
		for(int l = 0; l < L; l++){ 
	        vec[i][l] = new float[min_dim];     
		}
	}	

	for(int i = 0; i < max_book; i++){
		for(int j = 0; j < L; j++){
			for(int l = 0; l < min_dim; l++){
				vec[i][j][l] = quantizer[0][i][j][l];
			}
		}
	}
	
	
	int* pvec = new int[size_n];
	int ROUND1 = 10;
	int ROUND2 = 20;
	int min_vec;
	float temp, min_temp;
	int* tol_count = new int[L];
	
	printf("rotate matrix(OPQ)\n");
	for(int k1 = 0; k1 < ROUND1; k1++){
	    for(int i = 0; i < max_book; i++){
		    for(int k2 = 0; k2 < ROUND2; k2++){
		        for(int j = 0; j < size_n; j++){  //assignment of train vectors			
			       // max_temp = -1;
			        for(int l = 0; l < L; l++){
				        temp = 0;
                        for(int x = 0; x < min_dim; x++){
			                temp += (train_org[j][i* min_dim +x] - vec[i][l][x]) * (train_org[j][i* min_dim +x] - vec[i][l][x]);
                        }
						
						if( l == 0) {min_temp = temp; min_vec = l;}
                        else if(temp < min_temp) {min_temp = temp; min_vec = l;}				
		            }
                     pvec[j] = min_vec;			
		        }
		
		        for(int j = 0; j < size_n; j++){
			        for(int x = 0; x < min_dim; x++){
			            vec[i][ pvec[j] ][x] = 0; 
			        }	
		        }

			    for(int j = 0; j < L; j++) tol_count[j] = 0;
                for(int j = 0; j < size_n; j++){  // compute new vectors
			        for(int x = 0; x < min_dim; x++){
			            vec[i][ pvec[j] ][x] += train_org[j][i * min_dim + x];	
			        }	
				    tol_count[ pvec[j] ]++;
		        }
			
                for(int j = 0; j < L; j++){
					if(tol_count[j] == 0) continue; 
				    for(int l = 0; l < min_dim; l++){
					    vec[i][j][l] /= tol_count[j];
				    }
			    }
		    }
		
		    for(int j = 0; j < size_n; j++){
                for(int x = 0; x < min_dim; x++){
				    cvmSet (M_Y, i* min_dim + x, j, vec[i][pvec[j]][x]);
			    }
		    }	
	    }

        if(k1 == ROUND1 - 1) break;
	
        cvTranspose(M_Y, M_YT); // transpose

        if(k1 == 0)
        cvMatMul(M_X2, M_YT, ABt);
        else{ cvMatMul(M_RX, M_YT, ABt);}

	    cvSVD(ABt, ABt_D, ABt_U, ABt_VT, CV_SVD_V_T); //SVD
	    cvMatMul( ABt_U, ABt_VT, M_R );
	    cvTranspose(M_R, M_RT);
    
	    if(k1 == 0){			
		    for (int i = 0; i < M_RC->rows; i++) {
                for (int j = 0; j < M_RC->cols; j++) {
				    if(i == j)
                        cvmSet (M_RC, i, j, 1);
			        else{
				        cvmSet (M_RC, i, j, 0);	
				    }
                }
            }
		    cvMatMul( PCA_R, M_RC, M_RC );	
	    }
      
        cvMatMul( M_RT, M_RC, M_RC );

	    cvMatMul( M_RC, M_X, M_RX );
           
        for(int i = 0; i < size_n; i++){
		    for(int j = 0; j < dim_; j++){
			    train_org[i][j] = cvmGet (M_RX, j, i);
		    }
	    }
	}
	
	for(int i = 0; i < max_book; i++){
		for(int j = 0; j < L; j++){
			for(int l = 0; l < min_dim; l++){
				quantizer[0][i][j][l] = vec[i][j][l];
			}
		}
	}
	
	float** train1 = new float* [size_n];
	float** train2 = new float* [size_n];
	float** train0 = new float* [size_n];
	
    float* temp_arr;
	
	float sum3 = 0;
	float min_sum3 = -1;
	int min_id3;

    printf("permutate sub-codebooks\n");	
	for (int i = 1; i < max_level_; i++){   //if max_level_ = 1,skip
		
		float*** quantizer2 = new float** [length[i-1]];
		
		for (int j = 0; j < length[i-1]; j++){
            quantizer2[j] = new float* [L];
	    
		    for(int l = 0; l < L; l++)
			    quantizer2[j][l] = new float [dim[i-1]];
	    }
		
		int** index1 = new int* [length[i-1]];
        for(int j = 0; j < length[i-1]; j++)
            index1[j] = new int[size_n];
				
		int** count1 = new int* [length[i-1]];
		
        for(int j = 0; j < length[i-1]; j++)
            count1[j] = new int[L];
		
		for(int j = 0; j < length[i-1]; j++)
			for(int l = 0; l < L; l++)
				count1[j][l] = 0;
		

        int** con_count = new int*[L];
        for(int j = 0; j < L; j++)
			con_count[j] = new int[L];
	    
        for(int j = 0; j < L; j++)
			for(int l = 0; l < L; l++)
			con_count[j][l] = 0;
		
		temp_arr = new float[dim[i-1]];
		bool* flag3 = new bool[length[i-1]];
		for(int j = 0; j < length[i-1]; j++) flag3[j] = false;
		
		
		for(int j = 0; j < length[i-1]; j++){
			for(int l = 0; l < size_n; l++){
				for(int jj = 0; jj < dim[i-1]; jj++)
					temp_arr[jj] = train_org[l][ j * dim[i-1] +jj ];
			    
				for(int ii = 0; ii < L; ii++){
					sum3 = 0;
					for(int jj = 0; jj < dim[i-1]; jj++)
						sum3 += (temp_arr[jj] - quantizer[i-1][j][ii][jj] ) * (temp_arr[jj] - quantizer[i-1][j][ii][jj]);
					
					if(ii == 0) {min_sum3 = sum3; min_id3 = ii;}
					else if(sum3 < min_sum3) {min_sum3 = sum3; min_id3 = ii;}	
				}
				index1[j][l] = min_id3;
				count1[j][min_id3]++;
			}
		}
		
		k_elem* arr2 =new k_elem[ length[i-1] * length[i-1]];
		for(int ii = 0; ii < length[i-1]*length[i-1]; ii++)
			arr2[ii].id = -1; 
		
		for(int ii = 0; ii < length[i-1]; ii++){
			for(int jj = ii+1; jj < length[i-1]; jj++){

                for(int kk1 = 0; kk1 < L; kk1++ )
                    for(int kk2 = 0; kk2 < L; kk2++ )
                        con_count[kk1][kk2] = 0;

				for(int l = 0; l < size_n; l++)
					con_count[ index1[ii][l] ][ index1[jj][l] ]++;
                                                     
				sum3 = 0;
				for(int i1 = 0; i1 < L; i1++){
					if( count1[ii][i1] == 0) continue;
					
					for(int j1 = 0; j1 < L; j1++){
						if( con_count[i1][j1]==0 || count1[jj][j1] == 0) continue;
						sum3 += con_count[i1][j1] * ( log( double(size_n)* con_count[i1][j1]/ count1[ii][i1]/  count1[jj][j1])  ) / size_n;					    	
					}
				}
            
                arr2[ii* length[i-1] + jj].dist = sum3;
                arr2[ii* length[i-1] + jj].id = ii* length[i-1] + jj;
    
			}			
		}
        
		qsort(arr2, length[i-1] * length[i-1], sizeof(k_elem), QsortComp);
		
		int s1 = 0;
		for (int j = 0; j < M_R->rows; j++) {
            for (int l = 0; l < M_R->cols; l++) {
				cvmSet (M_R, j, l, 0);				
            }
        }		
	
		for(int ii = length[i-1]*length[i-1]-1; ii >= 0; ii--){
			if (arr2[ii].id == -1) continue;
			int e1 = arr2[ii].id / length[i-1];
			int e2 = arr2[ii].id % length[i-1];

			if(flag3[e1] == false && flag3[e2] == false){
                                
				flag3[e1] = true;
				flag3[e2] = true;
				
				int j1 = s1 * dim[i-1];
				for(int i1 = e1 * dim[i-1]; i1 < (e1+1) * dim[i-1]; i1++ ){
					cvmSet (M_R, i1, j1, 1);
				    j1++;	
				}	
				
				for(int l = 0; l < L; l++)
				    for(int jj =0; jj < dim[i-1]; jj++)
		                quantizer2[s1][l][jj] = quantizer[i-1][e1][l][jj];
				
				int j2 = (s1+1) * dim[i-1];
				for(int i2 = e2 * dim[i-1]; i2 < (e2+1) * dim[i-1]; i2++ ){
					cvmSet (M_R, i2, j2, 1);
				    j2++;	
				}

			    for(int l = 0; l < L; l++)
				    for(int jj =0; jj < dim[i-1]; jj++)
		                quantizer2[s1+1][l][jj] = quantizer[i-1][e2][l][jj];				
                
				s1 += 2;
			}	
		}

		cvTranspose(M_R, M_RT) ;
		cvMatMul( M_RT, M_RC, M_RC );
		cvMatMul( M_RC, M_X, M_RX );
		
		for(int i1 = 0; i1 < size_n; i1++){
		    for(int j1 = 0; j1 < dim_; j1++){
			    train_org[i1][j1] = cvmGet (M_RX, j1, i1);
		    }
     	}
		
		if(i == max_level_ -1) break;
        
		for(int j = 0; j < length[i-1]; j++)
			for(int l = 0; l < L; l++)
				for(int jj =0; jj < dim[i-1]; jj++)
		            quantizer[i-1][j][l][jj] = quantizer2[j][l][jj];
		
        //------------------------------------------------------		
		for(int l = 0; l < size_n; l++){
		    train1[l] = new float[dim[i-1]];
			train2[l] = new float[dim[i-1]];
			train0[l] = new float[dim[i]];
		}
		for(int ii = 0; ii < length[i]; ii++){
			
			for(int j1 = 0; j1 < L; j1++){
       			for(int j2 = 0; j2 < L; j2++){
				    for(int l = 0; l < dim[i]; l++){
						if(l < dim[i-1]) {temp_quan[i][ii][j1*L+j2][l] = quantizer[i-1][2*ii][j1][l];} 
						else{temp_quan[i][ii][j1*L+j2][l] = quantizer[i-1][2*ii+1][j2][l-dim[i-1]];}
					}
				}    
			}
			
			for(int l = 0; l < size_n; l++){
			    for(int jj = 0; jj < dim[i-1]; jj++){
				    train1[l][jj] = train_org[l][ 2 * ii * dim[i-1] + jj];
				    train2[l][jj] = train_org[l][ (2 * ii +1) *dim[i-1] +jj ];
			    }
				for(int jj = 0; jj < dim[i]; jj++){
					train0[l][jj] = train_org[l][ii * dim[i] + jj];
				}
			}	
			sub_kmeans_(temp_quan[i][ii],quantizer[i][ii],merge[i][ii], dim[i], train1, train2, train0);
		}
		
		for(int ii = 0; ii < length[i]; ii++){
			for(int j = 0; j < L; j++){
				for(int l = 0; l < dim[i]; l++){
					if(l < dim[i-1])
					    {quantizer[i][ii][j][l] = quantizer[i-1][2*ii][ merge[i][ii][2*j] ][l];}
                    else{quantizer[i][ii][j][l] = quantizer[i-1][2*ii+1][ merge[i][ii][2*j+1] ][l-dim[i-1]];}						
				}
			}
		}		
		
		for(int l = 0; l < size_n; l++){
		    delete[] train1[l];
			delete[] train2[l];
			train1[l] = NULL;
			train2[l] = NULL;
		}
		delete[] arr2; arr2 = NULL;
		for(int j = 0; j < length[i-1]; j++) {delete[] index1[j]; delete[] count1[j];}
        delete[] index1; delete[] count1;
        delete[] flag3;		
	}
		
    printf("merging step\n");
	for(int i = 0; i < size_n; i++){
		for(int ii = 0; ii < max_book; ii++){ 
			for(int jj = 0; jj < min_dim; jj++){
				train[ii][i][jj] = train_org[i][ii * min_dim +jj];
			}
		}
	}	
	
	for(int i = 0; i < max_book; i++){ 
		kmeans_(train[i],quantizer[0][i], size_n, min_dim);
	}
		
	for (int i = 1; i < max_level_; i++){   //if max_level_ = 1,skip
        //------------------------------------------------------		
		for(int l = 0; l < size_n; l++){
		    train1[l] = new float[dim[i-1]];
			train2[l] = new float[dim[i-1]];
			train0[l] = new float[dim[i]];
		}
		for(int ii = 0; ii < length[i]; ii++){
			
			for(int j1 = 0; j1 < L; j1++){
       			for(int j2 = 0; j2 < L; j2++){
				    for(int l = 0; l < dim[i]; l++){
						if(l < dim[i-1]) {temp_quan[i][ii][j1*L+j2][l] = quantizer[i-1][2*ii][j1][l];} 
						else{temp_quan[i][ii][j1*L+j2][l] = quantizer[i-1][2*ii+1][j2][l-dim[i-1]];}
					}
				}    
			}
			
			for(int l = 0; l < size_n; l++){
			    for(int jj = 0; jj < dim[i-1]; jj++){
				    train1[l][jj] = train_org[l][ 2 * ii * dim[i-1] + jj];
				    train2[l][jj] = train_org[l][ (2 * ii +1) *dim[i-1] +jj ];
			    }
				for(int jj = 0; jj < dim[i]; jj++){
					train0[l][jj] = train_org[l][ii * dim[i] + jj];
				}
			}	
			sub_kmeans_(temp_quan[i][ii],quantizer[i][ii],merge[i][ii], dim[i], train1, train2, train0);
		}
		
		for(int ii = 0; ii < length[i]; ii++){
			for(int j = 0; j < L; j++){
				for(int l = 0; l < dim[i]; l++){
					if(l < dim[i-1])
					    {quantizer[i][ii][j][l] = quantizer[i-1][2*ii][ merge[i][ii][2*j] ][l];}
                    else{quantizer[i][ii][j][l] = quantizer[i-1][2*ii+1][ merge[i][ii][2*j+1] ][l-dim[i-1]];}						
				}
			}
		}		
		
		for(int l = 0; l < size_n; l++){
		    delete[] train1[l];
			delete[] train2[l];
			train1[l] = NULL;
			train2[l] = NULL;
		}
	}
	
	unsigned char** merge0 = new unsigned char* [min_book];
	for(int i  = 0; i < min_book; i++){
		merge0[i] = new unsigned char [cen * nnum];  // merge nnum quantizers
	}
	
	float*** train_new = new float** [min_book];
	
	for(int i = 0; i < min_book; i++)
		train_new[i] = new float* [size_n];
	
	for(int i = 0; i < min_book; i++){
     	for(int j = 0; j < size_n; j++){
			train_new[i][j] = new float[max_dim];
	    }
	}
	
	for(int ii = 0; ii < min_book; ii++){				
		for(int j = 0; j < size_n; j++){
			for(int jj = 0; jj < max_dim; jj++){
				train_new[ii][j][jj] = train_org[j][ ii * max_dim + jj];
			}
	    }	
		sub_kmeans0_(ii, quantizer[max_level_-1], merge0[ii], train_new[ii], dim[max_level_-1], max_dim, quantizer0[ii]);
	}


    //----------rotate dataset--------------------------------------

    float* data2 = new float [dim_];
	float* R = new float [dim_ * dim_];
	
	for(int i = 0; i < dim_; i++){
		for(int j = 0; j < dim_; j++){ 
            R[i * dim_ + j] = cvmGet (M_RC, i, j);
		}
	}
	
    for(int i = 0; i < vecsize; i++){
		for(int j = 0; j < dim_; j++){
			data2[j] = 0;
			for(int l = 0; l < dim_; l++){
			    data2[j] += R[j * dim_ + l] * data[i][l];	
			}
		}
		
		for(int l = 0; l < dim_; l++){
			data[i][l] = data2[l];
		}
	}

	//-----------------determine which candidates are eligible------------
		
    std::vector< std::vector<elem_> > object_selection;
		
	float sum = 0;
	float min_sum = 0;
	int min_id = -1;
	int pos_;
    bool flag1 = false;
	bool flag2 = false;
	bool flag;
		
	std::vector<elem_> temp2;	
	for(int i = 0; i < max_level_*L*L; i++){
		object_selection.push_back(temp2);
	}
			
	bool** flag_obj = new bool* [max_level_];
    for(int i = 0; i < max_level_; i++){	
        flag_obj[i] = new bool[vecsize];
	}
	for(int i = 0; i < max_level_; i++){
        for	(int j = 0; j < vecsize; j++){		
            flag_obj[i][j] = false;
		}
	}
	
	elem_ element;
	
	//-----------------determine D_t---------------------------------------------------------
	
	bool** fflag = new bool*[max_level_];
	for(int i = 0; i < max_level_; i++)
		fflag[i] = new bool[vecsize];
	
	for(int j = 0; j < vecsize; j++){
		fflag[max_level_-1][j] = true;
	}
	
	for(int i = 0; i < max_level_-1; i++){
		for(int j = 0; j < vecsize; j++){
			fflag[i][j] = false;
		}
	}
	
	int cut_size = vecsize;
	
	for(int i = max_level_-1; i >= 0; i--){
		
		k_elem* g_t = new k_elem[vecsize];		
	    float* temp_dist = new float [vecsize];
		for(int jj = 0; jj < vecsize; jj++){
			temp_dist[jj] = 0;
		    for(int ii = 0; ii < length[i]; ii++){
			    for(int j = 0; j < L; j++){
				    sum = 0;
				    for(int l = 0; l < dim[i]; l++){
                        sum += (data[jj][ii* dim[i] + l] - quantizer[i][ii][j][l]) * (data[jj][ii* dim[i] + l] - quantizer[i][ii][j][l]);					
				    }
				    if( j == 0) {min_id = 0; min_sum = sum;}
				    else{ if(sum < min_sum) {min_sum = sum; min_id = j;} };
			    }
			    obj_quantizer[i][jj][ii] = min_id;
			    temp_dist[jj] += min_sum;
		    }
			g_t[jj].dist = log(temp_dist[jj]) / 2 + density[jj];
			g_t[jj].id = jj;
		}

		qsort(g_t, vecsize, sizeof(k_elem), QsortComp);
                  	
		int u = 0;

        if(i < max_level_-1){
		    for(int jj = vecsize - 1; jj >=0; jj--){
			    if( fflag[i+1][ g_t[jj].id ] == true){
				    fflag[i][ g_t[jj].id ] = true;
				    u++;
				    if(u >= cut_size) break;
			    }	
		    }	
        }
		cut_size = delta * cut_size;

		for(int jj = 0; jj < vecsize; jj++){
			if(fflag[i][jj] == false) continue; //this points does not exists
			
			element.id = jj;
		    element.index = new unsigned char[length[i]];   //release space
			for(int ii = 0; ii < length[i]; ii++){
			    element.index[ii] =  obj_quantizer[i][jj][ii];
		    }
            element.dist = temp_dist[jj];
		
		    pos_ = i* L * L +  element.index[0] * L +  element.index[1];
		    if(object_selection[pos_].empty()){
			    flag_obj[i][jj] = true;
			    object_selection[pos_].push_back(element);
		    }
		    else{
			    for(int ii = 0; ii < object_selection[pos_].size(); ii++){
				    flag = false;
			        for(int j = 2; j < length[i]; j++){
				        if (object_selection[pos_][ii].index[j] != element.index[j])
				            {flag = true; break;}
			        }
			        if(flag == false){
			            if( element.dist < object_selection[pos_][ii].dist){  //replace and update
					        flag_obj[i][ object_selection[pos_][ii].id ] = false;
				            object_selection[pos_][ii].id = element.id;
      				        object_selection[pos_][ii].dist = element.dist;
                            flag_obj[i][jj] = true;							
			            }
                        break;						
		           }
	            }
			    if(flag == true) {
				    object_selection[pos_].push_back(element);
				    flag_obj[i][jj] = true;
			    }
		    }		
		}
	} 
	
	for(int i = 0; i < vecsize; i++)
		delete[] data[i];
	
	delete[] data;
	
    for(int i = 0; i < size_n; i++)
        delete[] train_org[i];
    delete[] train_org;
	
	   	
	for(int i = max_level_-1; i > 0; i--){ // not deal with level 0

        elem_ element2;
	    element2.index = new unsigned char[length[i-1]];
					
        for(int j = 0; j < vecsize; j++){		
            if(flag_obj[i][j] == true){
				if(fflag[i-1][j] == true){  //the representative exists in the base layer
				    element2.id = j;
				    for(int jj = 0; jj < length[i-1]; jj++)
				        element2.index[jj] = obj_quantizer[i-1][j][jj];
				
				    pos_ = (i-1) * L * L +  element2.index[0] * L +  element2.index[1];
					
				    for(int ii = 0; ii < object_selection[pos_].size(); ii++){
					    flag = false;
					    for(int jj = 2; jj < length[i-1]; jj++){
						    if (object_selection[pos_][ii].index[jj] != element2.index[jj])
						        {flag = true; break;}
					    }
                        if(flag == false){  // find the connected quantizer
                            init_obj[i][j] = object_selection[pos_][ii].id;	 //	
                             break;						
				        }							
			        }
				}
                else{init_obj[i][j] = j;}	// connect to representative in the base layer			
			}			
		}
    }

	//---------------------------------------------
        int count_obj[max_level_];
        for(int i = 0; i < max_level_; i++){
           count_obj[i]=0;
        }

        for(int i = 0; i < max_level_; i++){
            for(int j = 0; j < vecsize; j++){
                if(flag_obj[i][j] == true){
                   count_obj[i]++;
                }
            }
        }



	unsigned char ID[4];
			
    float**** quan_book = new float*** [max_level_];
			
	int* trans = new int[vecsize];
	int* trans2 = new int[vecsize];
	int count_ = 0;
	
	bool** fflag2 = new bool* [max_level_];
    char** fflag3 = new char* [max_level_]; 

    for(int ii = 0; ii < max_level_; ii++){
	    fflag2[ii] = new bool[vecsize]; 
     	fflag3[ii] = new char[vecsize]; 
	}
	
	for(int ii = 0; ii < max_level_; ii++){
		for(int jj = 0; jj < vecsize; jj++){
	        fflag2[ii][jj] = false; 
     	    fflag3[ii][jj] = 0; 
		}
	}	

	for(int ii = 0; ii < max_level_; ii++){              
 
		count_ = 0;
		quan_book[ii] = new float** [length[ii]];
		for(int j = 0; j < length[ii]; j++){
			quan_book[ii][j]= new float* [L];
		}
			
		for(int j = 0; j < length[ii]; j++){
			for(int l = 0; l < L; l++){
				quan_book[ii][j][l]= new float[L];
			}	
		}
		
		for(int j = 0; j < length[ii]; j++){
			for(int l = 0; l < L; l++){
				for(int jj = 0; jj < L; jj++){
					sum = 0;
					for(int s = 0; s < dim[ii]; s++)
						sum += (quantizer[ii][j][l][s] - quantizer[ii][j][jj][s]) * (quantizer[ii][j][l][s] - quantizer[ii][j][jj][s]);
					quan_book[ii][j][l][jj] = sum; 
				}	
			}	
		}
		
		count_ = 0;
		for(int l = 0; l < vecsize; l++){
		    if(flag_obj[ii][l] == true){
			    trans[l] = count_;
                if(ii == 0){
				    init_obj[0][count_] = l; //from internal id to ture id
			    }	
				count_++;
			}		
		}

		if(ii < max_level_ - 1){
			for(int i = 0; i < vecsize; i++){ //from external to internal
			    if(flag_obj[ii+1][i] == true && fflag[ii][i] == true) {
				    init_obj[ii+1][i] = trans[ init_obj[ii+1][i] ];	
				    fflag2[ii+1][i] = true;	
				}	
			}	
		}

		if(ii > 0){ //from internal to internal
			count_ = 0;
			for(int i = 0; i < vecsize; i++){
				if(flag_obj[ii][i] == true){
					if(fflag2[ii][i]  == true){
						fflag3[ii][count_] = 1;
					}
					
			      	init_obj[ii][count_] = init_obj[ii][i];	
					count_++;
				}
			}	
		}	
	}
		
    for(int ii = 0; ii <= max_level_-1; ii++){  
		count_ = 0;
			
		for(int l = 0; l < vecsize; l++){
			if(flag_obj[ii][l] == true){
			    appr_alg->addPoint((void *) (obj_quantizer[ii][l]), (size_t) l, quan_book[ii], ii, true);
				count_ = l;
				break;
			}	
		}	
		  
        int j1 = 0;
	
		for(int l = count_+1; l < vecsize; l++){
			if(l % report_every == 0)
		        printf("Level %d, build HNSW index...%.2f %% complete\n", ii, l * 100 / float(vecsize) );
				
			if(flag_obj[ii][l] == true){
			    appr_alg->addPoint((void *) (obj_quantizer[ii][l]), (size_t) l, quan_book[ii], ii, false);
		    }	
		}
	    
		std::string tmp_string;
		if(ii > -1) tmp_string = quan_path[ii];
		else {tmp_string = data_path;}
		
		ofstream output_knn(tmp_string, ios::binary);
		
		
		count_ = 0;
		int count2_;
	    
		for(int l = 0; l < vecsize; l++){
			if(l > 0 && l % report_every == 0)
		        printf("Level %d, build Knn graph...%.2f %% complete\n", ii, l * 100 / float(vecsize) );
			
		    if(flag_obj[ii][l] == true){
				appr_alg->connect((void *) (obj_quantizer[ii][l]), quan_book[ii], neighbor2, ii, true);  //should return internal IDs			
				count2_ = 0;
				output_knn.write((char*) &count_, 4);
				for(int i = 0; i < KK-1; i++){  //remove the same id
					if(count_ == neighbor2[i]) {count2_++;}
					neighbor[i]  = neighbor2[count2_];
					count2_++;		
				}
				output_knn.write((char*) neighbor, 4 * (KK-1));	
				count_++;
			}			
		}
		output_knn.close();
		if( ii < max_level_-1 )
		    appr_alg-> deleteLinklist(ii);
	}
	
	//-----------------------build top layer----------------------------------------
	unsigned char* tmp = new unsigned char[min_book * nnum]; 
	int tmp_id;
	int eff = 200;
	appr_alg -> setEf(eff);
	
	printf("build top layer\n");
	for(int i = 0; i < cen; i++){
		for(int j = 0; j < nnum; j++)
			tmp[j] = merge0[0][i * nnum + j];
		
		for(int j = 0; j < cen; j++){
			
			for(int l = 0; l < nnum; l++)
				tmp[nnum + l] = merge0[1][j * nnum + l];
			
			for(int l = 0; l < cen; l++){
				
				for(int s = 0; s < nnum; s++)
					tmp[2* nnum + s] = merge0[2][l * nnum +s];
					
				for(int s = 0; s < cen; s++){

				    for(int ii = 0; ii < nnum; ii++)
					    tmp[3* nnum + ii] = merge0[3][s * nnum +ii];
				
					tmp_id = i * cen * cen * cen + j* cen * cen + l * cen + s; 
					appr_alg->connect((void *) (tmp), quan_book[max_level_-1], start_book[tmp_id], max_level_-1, false);
				}	
			}
		}
	}

    appr_alg-> deleteLinklist(max_level_-1);
	
    input.close();
    cout << "Build time:" << 1e-6 * stopw_full.getElapsedTimeMicro() << "  seconds\n";

	//------------------------------store info---------------------------------
	ofstream outputQ("quantizer.gt", ios::binary);
	outputQ.write((char *) length, 4 * max_level_);
	outputQ.write((char *) dim, 4 * max_level_);
	outputQ.write((char *) count_obj, 4 * max_level_);
		
	for(int i = 0; i < max_level_; i++){
		for(int j = 0; j < vecsize; j++){
			if(flag_obj[i][j] == true)
			    outputQ.write( (char *) obj_quantizer[i][j], length[i]);
		}		
	}	
			
	for(int i = 0; i < max_level_; i++){
		for(int j = 0; j < length[i]; j++){
		    for(int l = 0; l < L; l++){
		        outputQ.write( (char *) quan_book[i][j][l], 4 * L);					
		    }
		}	
	}

	for(int i = 0; i < max_level_; i++){
		for(int j = 0; j < length[i]; j++){
		    for(int l = 0; l < L; l++){
		        outputQ.write( (char *) quantizer[i][j][l], 4 * dim[i]);					
		    }
		}	
	}
    outputQ.close();
		
    ofstream outputQ2("searching.gt", ios::binary);	
    outputQ2.write((char *) R, 4 * dim_ * dim_);
		
	for(int i = 0; i < max_level_; i++){    
		outputQ2.write((char *) fflag3[i], vecsize);
	}
		
	for(int i = 0; i < max_level_; i++){
		for(int j = 0; j < length[i]; j++){
		    outputQ2.write((char *) merge[i][j], 2 * L);
		}
	}
	
	for(int i = 0; i < min_book; i++){
		outputQ2.write((char *) merge0[i], nnum * cen);
	}	
	
	for(int i = 0; i < length[0]; i++){
		for(int j = 0; j < L; j++){
		    outputQ2.write((char*) quantizer[0][i][j], 4 * dim[0]);				
		}
    }
	
	for(int i = 0; i < Tol; i++){
		outputQ2.write((char *) start_book[i], 4 * fan);
	}
	
    for(int i = 0; i < max_level_; i++){	
        outputQ2.write((char *) init_obj[i], 4 * count_obj[i]);
	}
    outputQ2.close();		    	
    
    return;

}
