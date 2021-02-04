#include <efanna2e/index_nsg.h>
#include <efanna2e/util.h>
#include <iostream>

void load_data(char* filename, float*& data, unsigned& num,
               unsigned& dim) {  // load data with sift10K pattern
  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    std::cout << "open file error" << std::endl;
    exit(-1);
  }
  in.read((char*)&dim, 4);
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

int main(int argc, char** argv) {	
  if (argc != 8) {
    std::cout << argv[0] << " data_file nn_graph_path L R C save_index level"  //data_graph.fvecs
              << std::endl;
    exit(-1);
  }
  float* data_load = NULL;
  unsigned points_num, dim;
  int max_level = atoi(argv[7]);
  
  load_data(argv[1], data_load, points_num, dim); 

//--------------quantizer---------------------------------------	
  std::ifstream inQ("quantizer.gt", std::ios::binary);
	
	int* length = new int[max_level];
	int* sdim = new int[max_level];
	int* count = new int[max_level];

	inQ.read( (char *) length, 4 * max_level);
	inQ.read( (char *) sdim, 4 * max_level);
	inQ.read( (char *) count, 4 * max_level);
				
	unsigned char** obj_quantizer = new unsigned char*[max_level];
	for(int i = 0; i < max_level; i++){
		obj_quantizer[i] = new unsigned char [count[i] * length[i]];
	}

	float**** quantizer = new float*** [max_level];
	for (int i = 0; i < max_level; i++){
        quantizer[i] = new float** [length[i]];
	   			
		for(int ii = 0; ii < length[i]; ii++){
			quantizer[i][ii] = new float* [LL];
		}
		for(int ii = 0; ii < length[i]; ii++){
			for(int j = 0; j < LL; j++)
				quantizer[i][ii][j] = new float[sdim[i]];
		}
	}
	
	float**** quan_book = new float*** [max_level];
	for(int ii = 0; ii < max_level; ii++){
	  //	count_ = 0;
		quan_book[ii] = new float** [length[ii]];
		for(int j = 0; j < length[ii]; j++){
			quan_book[ii][j]= new float* [LL];
		}
			
		for(int j = 0; j < length[ii]; j++){
			for(int l = 0; l < LL; l++){
				quan_book[ii][j][l]= new float[LL];
			}	
		}	
	}
	//--------------------read array----------------------
	for(int i = 0; i < max_level; i++){
	    inQ.read( (char *) obj_quantizer[i], count[i] * length[i]);		
	}	
		
		
	for(int i = 0; i < max_level; i++){
		for(int j = 0; j < length[i]; j++){
		    for(int l = 0; l < LL; l++){
		        inQ.read( (char *) quan_book[i][j][l], 4 * LL);					
		    }
		}	
	}

	for(int i = 0; i < max_level; i++){
		for(int j = 0; j < length[i]; j++){
		    for(int l = 0; l < LL; l++){
		        inQ.read( (char *) quantizer[i][j][l], 4 * sdim[i]);					
		    }
		}	
	}
    inQ.close();
//-------------------------------------------------------------- 

  std::string nn_graph_path(argv[2]);
  unsigned L = (unsigned)atoi(argv[3]);
  unsigned R = (unsigned)atoi(argv[4]);
  unsigned C = (unsigned)atoi(argv[5]);
  //------------------------------------
  std::string quan_path[max_level];
  std::string quan_index[max_level];
  char tmp_path[1024];
  
  
  for(int i = 0; i < max_level; i++){
	  sprintf(tmp_path, "quan_graph%d.fvecs", i);
	  quan_path[i] = tmp_path;
	  memset(tmp_path,'\0',sizeof(tmp_path));
	  
	  sprintf(tmp_path, "quan_index%d.fvecs", i);
	  quan_index[i] = tmp_path;
	  memset(tmp_path,'\0',sizeof(tmp_path));
	  //sprintf(quan_index[i], "quan_index%d.fvecs", i);
  }
  
  efanna2e::IndexNSG index(dim, points_num, efanna2e::L2, nullptr, false, nullptr, -1, nullptr, 1);
  
  auto s = std::chrono::high_resolution_clock::now();
  efanna2e::Parameters paras;
  paras.Set<unsigned>("L", L);
  paras.Set<unsigned>("R", R);
  paras.Set<unsigned>("C", C);

  for(int i = 0; i < max_level; i++){       //for quantizers

	  
    efanna2e::IndexNSG index2(sdim[i]*length[i], count[i], efanna2e::L2, nullptr, true, quan_book[i], length[i], quantizer[i], max_level);
	  paras.Set<std::string>("nn_graph_path", quan_path[i]);  //for reading knn-graph
	  index2.Build2(count[i], obj_quantizer[i], paras);
	  index2.Save(quan_index[i].c_str());
	  delete[] obj_quantizer[i];
	  obj_quantizer[i] = NULL;
  }
  
  paras.Set<std::string>("nn_graph_path", nn_graph_path);
  index.Build(points_num, data_load, paras); //for original data

  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e - s;

  std::cout << "indexing time: " << diff.count() << "\n";
  index.Save(argv[6]);
  
  return 0;
}
