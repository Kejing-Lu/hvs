#include <efanna2e/index_nsg.h>
#include <efanna2e/util.h>
#include <efanna2e/distance.h>
#include <chrono>
#include <string>

using namespace efanna2e;

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


void load_data(char* filename, float*& data, unsigned& num,
               unsigned& dim) {  // load data with sift10K pattern
  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    std::cout << "open file error1" << std::endl;
    exit(-1);
  }
  in.read((char*)&dim, 4);
  // std::cout<<"data dimension: "<<dim<<std::endl;
  in.seekg(0, std::ios::end);
  std::ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = (unsigned)(fsize / (dim + 1) / 4);
  data = new float[(size_t)num * (size_t)dim];

  in.seekg(0, std::ios::beg);
  for (size_t i = 0; i < num; i++) {
    in.seekg(4, std::ios::cur);
    in.read((char*)(data + i * dim), dim * 4);
  }
  in.close();
}




void restore_index2(float *query_data, float* array0, float*** book, 
        float** start_book, unsigned char*** merge, unsigned char** merge0, int* length, int* sdim, int tol_dim, float* temp_arr,
		    DistanceL2New* dist0, DistanceInnerProduct* dist1, int max_level){

     int start_;
     int ll;	
	       
    for(int ii = 1; ii < max_level; ii++){
		for(int i = 0; i < length[ii]; i++){
		    for(int j = 0; j < LL; j++){
			    book[ii][i][j] = book[ii-1][2*i][merge[ii][i][2*j]] + book[ii-1][2*i+1][merge[ii][i][2*j+1]];
		    }	   
	    }	
	}
	
	for(int i = 0; i < min_book; i++){
		for(int j = 0; j < cen; j++){
			ll = i * nnum;
            start_ = j * nnum;				
			for(int l = 0; l < nnum; l++){
			    start_book[i][j] += book[max_level-1][ll][ merge0[i][start_] ];
				ll++;
				start_++;
			}
		}	   

         }
}

void restore_index(float *query_data, float* array0, float* book, 
        float** start_book, unsigned char*** merge, unsigned char** merge0, int* length, int* sdim, int tol_dim, float* temp_arr,
 	DistanceL2New* dist0, DistanceInnerProduct* dist1){
			
    float* temp_arr2 = new float[sdim[0]];

    int base_dim = sdim[0];
		
	float* ind2 = array0 ;
			
	for(int j = 0; j < base_dim; j++){
		temp_arr2[j] = dist1 -> compare(ind2, query_data, tol_dim);
        ind2 += tol_dim;			
	}   

	for(int j = 0; j < LL; j++){
		_mm_prefetch(ind2, _MM_HINT_T0);

		book[j] = dist0 -> compare(ind2, temp_arr2, base_dim);

	    ind2 += base_dim;		
	}
	          			
}

void save_result(const char* filename,
                 std::vector<std::vector<unsigned> >& results) {
  std::ofstream out(filename, std::ios::binary | std::ios::out);

  for (unsigned i = 0; i < results.size(); i++) {
    unsigned GK = (unsigned)results[i].size();
    out.write((char*)&GK, sizeof(unsigned));
    out.write((char*)results[i].data(), GK * sizeof(unsigned));
  }
  out.close();
}
int main(int argc, char** argv) {
  if (argc != 7) {
    std::cout << argv[0]
              << " data_file query_file nsg_path level search_K"
              << std::endl;
    exit(-1);
  }
  float* data_load = NULL;
  unsigned points_num, dim;
  load_data(argv[1], data_load, points_num, dim);
  float* query_load = NULL;
  unsigned query_num, query_dim;
  load_data(argv[2], query_load, query_num, query_dim);  //query add dimension;
  assert(dim == query_dim);

 // unsigned L = (unsigned)atoi(argv[4]);
  unsigned L = 100;
  unsigned K = (unsigned)atoi(argv[5]);
  int max_level = atoi(argv[4]);
  
  std::cout << "Loading GT:\n";
  std::ifstream inputGT(argv[6], std::ios::binary);
  int maxk = 100;
  unsigned int* massQA = new unsigned int[query_num * maxk]; //alt
  for (int i = 0; i < query_num; i++) {
      int t;
      inputGT.read((char *) &t, 4);
      inputGT.read((char *) (massQA + 100 * i), 4 * maxk);
  }
  inputGT.close();
  
  if (L < K) {
    std::cout << "search_L cannot be smaller than search_K!" << std::endl;
    exit(-1);
  }

  // data_load = efanna2e::data_align(data_load, points_num, dim);//one must
  // align the data before build query_load = efanna2e::data_align(query_load,
  // query_num, query_dim);
  efanna2e::IndexNSG index(dim, points_num, efanna2e::FAST_L2, nullptr, false, nullptr, -1, nullptr, max_level);
  index.Load(argv[3]);   

  index.OptimizeGraph(data_load);   
  
  //------------quantizer--------------
  std::string quan_index[max_level];
  char tmp_path[1024];
  
  for(int i = 0; i < max_level; i++){
	  sprintf(tmp_path, "quan_index%d.fvecs", i);
	  quan_index[i] = tmp_path;
	  memset(tmp_path,'\0',sizeof(tmp_path));
  }
  
 	std::ifstream inQ("quantizer.gt", std::ios::binary);
	
	int* length = new int[max_level];
	int* sdim = new int[max_level];
	int* count = new int[max_level];

	inQ.read( (char *) length, 4 * max_level);
	inQ.read( (char *) sdim, 4 * max_level);
	inQ.read( (char *) count, 4 * max_level);
				
	unsigned char** quan_load = new unsigned char*[max_level];
	for(int i = 0; i < max_level; i++){
		quan_load[i] = new unsigned char [count[i] * length[i]];
	}
  	for(int i = 0; i < max_level; i++){
	    inQ.read( (char *) quan_load[i], count[i] * length[i]);		
	}
    inQ.close();
 
    index.quan_graph = (char **)malloc(sizeof(char*) * max_level);
    for(int i = 0; i < max_level; i++){
		index.Load(quan_index[i].c_str()); 
        index.OptimizeQuan(quan_load[i], i, count[i], max_level);	  
    }
  //-------------read array-------------------
  float*** quantizer_Q = new float**[ length[0] ];
  for(int i = 0; i < length[0]; i++)
	quantizer_Q[i] = new float* [LL];
	
  for(int i = 0; i < length[0]; i++){
	for(int j = 0; j < LL; j++){
		quantizer_Q[i][j] = new float[ sdim[0] ];				
	}
  }
  
  unsigned char*** merge = new unsigned char** [max_level];
  for(int i  = 0; i < max_level; i++){
	merge[i] = new unsigned char* [length[i]];
  }	
	
  for (int i = 0; i < max_level; i++){
	for(int ii = 0; ii < length[i]; ii++){
		merge[i][ii] = new unsigned char[2*LL];
	}			
  }
  
  unsigned char** merge0 = new unsigned char* [min_book];
  for(int i  = 0; i < min_book; i++){
	merge0[i] = new unsigned char [cen * nnum];  // merge nnum quantizers
  }
  
  int Tol = cen * cen * cen * cen;
  int** connect = new int* [Tol]; 
  for(int i = 0; i < Tol; i++){	
    connect[i] = new int[fan];
  }

  unsigned int** trans = new unsigned int* [max_level]; 
  for(int i = 0; i < max_level; i++){	
    trans[i] = new unsigned int[count[i]];
  } 
  
  int tol_dim = length[0] * sdim[0];
  float* R = new float [tol_dim * tol_dim];
  
    char** fflag = new char* [points_num]; 

    for(int ii = 0; ii < max_level; ii++){ 
     	fflag[ii] = new char[points_num]; 
	}
  
  //--------------------------------searching parameters------------------------------------
  std::ifstream inQ2("searching.gt", std::ios::binary);
  
  inQ2.read((char *) R, 4 * tol_dim * tol_dim);

  for(int i = 0; i < max_level; i++){    
	 inQ2.read((char *) fflag[i], points_num);
  }  
  
  for(int i = 0; i < max_level; i++){
    for(int j = 0; j < length[i]; j++){
	   inQ2.read((char *) merge[i][j], 2 * LL);
	}
  }
  for(int i = 0; i < min_book; i++){
    inQ2.read((char *) merge0[i], nnum * cen);
  }	
  for(int i = 0; i < length[0]; i++){
    for(int j = 0; j < LL; j++){
	  inQ2.read((char*) quantizer_Q[i][j], 4 * sdim[0]);				
    }
  }		
  for(int i = 0; i < Tol; i++){
    inQ2.read((char *) connect[i], 4 * fan);
  }
  for(int i = 0; i < max_level; i++){
    inQ2.read((char *) trans[i], 4 * count[i]);
  } 
  inQ2.close();  
  
  //-------------initialize book-----------------------------
  float**** book = new float*** [query_num];
  for(int ii = 0; ii < query_num; ii++)
	  book[ii] = new float**[max_level];
		
  for(int ii = 0; ii < query_num; ii++)
	  for(int jj = 0; jj < max_level; jj++)
		  book[ii][jj] = new float* [ length[jj] ];
			
  for(int ii = 0; ii < query_num; ii++){
    for(int jj = 0; jj < max_level; jj++){
	  for(int l = 0; l < length[jj]; l++){
		book[ii][jj][l] = new float [LL];
		for(int ll = 0; ll < LL; ll++)
		  book[ii][jj][l][ll] = 0;
	  } 
    }
  }  

  float*** start_book = new float**[query_num];
  for(int ii = 0; ii < query_num; ii++)
    start_book[ii] = new float* [min_book];
  
  for(int ii = 0; ii < query_num; ii++)
		for(int jj = 0; jj < min_book; jj++)
		    start_book[ii][jj] = new float[cen];  
  for(int ii = 0; ii < query_num; ii++){
	 for(int jj = 0; jj < min_book; jj++){
	    for(int l = 0; l < cen; l++){
            start_book[ii][jj][l] = 0;
		} 
     }
  }	

  efanna2e::Parameters paras;
  paras.Set<unsigned>("L_search", L);
  paras.Set<unsigned>("P_search", L);

  std::vector<std::vector<unsigned> > res(query_num);
  for (unsigned i = 0; i < query_num; i++) res[i].resize(K);

  float* query_load2 = new float[query_num * tol_dim];
  float* tmp;
  float* tmp_query_load2;
  float* tmp_query_load;
  
  for(int i = 0; i < query_num; i++){
    	  tmp_query_load = query_load + i * dim;
     	  tmp_query_load2 = query_load2 + i * tol_dim;
      for(int j = 0; j < tol_dim; j++){
		  if(j < dim){
			  tmp_query_load2[j] = tmp_query_load[j];
		  }
		  else{
			  tmp_query_load2[j] = 0;
		  }
	  }
  }
 

float* array0 = new float[tol_dim * tol_dim + LL * tol_dim]; 
float* ind2; 

  for(int i = 0; i < length[0]; i++){
	ind2 = array0 + i * (sdim[0] * tol_dim + LL * sdim[0]); 
	  
	for(int j = 0; j < sdim[0] * tol_dim; j++)
            ind2[j] = R[i * sdim[0] * tol_dim + j];

    ind2 += sdim[0] * tol_dim;		
	  
	for(int j = 0; j < LL; j++){ 
                ind2[j * sdim[0] ] = 0;
                for(int l = 0; l < sdim[0]; l++)
                    ind2[j * sdim[0]] += quantizer_Q[i][j][l] * quantizer_Q[i][j][l];

		for(int l = 0; l < sdim[0]; l++){
                       // if(l < ssdim)
			    ind2[j * sdim[0] + l] = quantizer_Q[i][j][l];
                       // else ind2[j * sdim[0] + l] = 0;
		}				
	}
  }
  
  	DistanceL2New* dist0 = new DistanceL2New();
	DistanceInnerProduct* dist1 = new DistanceInnerProduct();

  float* temp_arr = new float[sdim[0]];

  VisitedListPool *visited_list_pool_ = new VisitedListPool(1, points_num);
  
StopW stopw2 = StopW();

 int const1 = sdim[0] * tol_dim + LL * sdim[0];

  #pragma omp parallel for num_threads(6)
  for (unsigned i = 0; i < query_num * length[0]; i++) {
          int  j = i / length[0];
          int  l = i % length[0];
                      

	  restore_index(query_load2 + j * tol_dim, array0 + l* const1, book[j][0][l], start_book[j], merge, merge0, length, sdim, tol_dim, temp_arr, dist0, dist1);
  } 

#pragma omp parrallel for num_thread(4)
  for (unsigned i = 0; i < query_num; i++) { 
    restore_index2(query_load2 + i * tol_dim, array0, book[i], start_book[i], merge, merge0, length, sdim, tol_dim, temp_arr, dist0, dist1, max_level);
  }   

  float time2 = stopw2.getElapsedTimeMicro()/query_num;
  std::cout<<"build distance book: "<< time2 << " us\n";

  unsigned int* points = new unsigned int[fan];

  float min_sum;
  int min_id;
  int id[4];

  std::vector<unsigned int> efs;
    
    if( K < 30){
        for (int i = K; i < 30; i++) {
            efs.push_back(i);
        }
		for (int i = 30; i < 100; i += 10) {
            efs.push_back(i);
        }
	}
	else{
		for (int i = K; i < 100; i += 10) {
            efs.push_back(i);
        }		 
	} 
	
    for (int i = 100; i < 1000; i += 50) {
        efs.push_back(i);
    }
  
    for (unsigned int ef : efs) {

        paras.Set<unsigned>("L_search", ef);
        paras.Set<unsigned>("P_search", ef); 
  

        unsigned int *enter_obj = new unsigned[ef];
        elem** init_obj = new elem* [max_level];
	    for(int i = 0; i < max_level; i++)
            init_obj[i] = new elem[ef];		 
	

        StopW stopw = StopW();
        for (unsigned ii = 0; ii < query_num; ii++) {

	    float** tmp_book = start_book[ii];
        for( int i = 0; i < 4; i++){
	        for(int j = 0; j < cen; j++){
		        if(j == 0) {min_id = 0; min_sum = tmp_book[i][j];}
		        else{
			        if(tmp_book[i][j] < min_sum) { min_id = j; min_sum = tmp_book[i][j];}
		        }
	        }
            id[i] = min_id;				
        }
        int temp = cen * cen * cen * id[0] + cen * cen * id[1] + cen * id[2] + id[3];

        for(int i = 0; i < fan; i++){
	        points[i] = connect[temp][i];
	    }
	
		if(max_level == 1){
			index.SearchWithsingleGraph(book[ii][0], length[0], points, enter_obj, trans[0], paras, visited_list_pool_);
		}
		else{
		    index.SearchWithquanGraph(book[ii][max_level-1], length[max_level-1], points, init_obj[max_level - 2], 
			trans[max_level-1], paras, fflag[max_level-1], max_level-1, visited_list_pool_);
	    
	        for(int i = max_level -2; i >= 1; i--)
                index.SearchWithquanGraph3(book[ii][i], length[i], init_obj[i], init_obj[i-1], trans[i], 
			paras, fflag[i], i, visited_list_pool_);
	
            index.SearchWithquanGraph2(book[ii][0], length[0], init_obj[0], enter_obj, trans[0], paras, visited_list_pool_);
        }		
        index.SearchWithOptGraph(query_load + ii * dim, K, paras, res[ii].data(), enter_obj, visited_list_pool_);
        }
		
        float time_us_per_query = stopw.getElapsedTimeMicro() / query_num;
        std::cout <<"ef="<< ef << "\t" << time_us_per_query << " us\n";

        unsigned int** answers = new unsigned int* [query_num];
        for(int i = 0; i < query_num; i++){
	        answers[i] = new unsigned int[K];
        }
  
  
        for (int i = 0; i < query_num; i++) {
            for (int j = 0; j < K; j++) {
                answers[i][j] = massQA[maxk * i + j];
            }
        }
  
        int correct = 0;
        int total = K* query_num;
        for (int i = 0; i < query_num; i++) {
            if(res[i].size() != K){
			    printf("error\n");
			    exit(-1);	
		    }
		    for(int j = 0; j < K; j++){
			    for(int l = 0; l < K; l++){
				    if(answers[i][j] == res[i][l]){
					    correct++;
					    break;
				    }
			    }
		    }	
		
        }
        printf("recall = %.5f\n", 1.0f * correct / total);
    }
    return 0;
}
