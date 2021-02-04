//
// Copyright (c) 2017 ZJULearning. All rights reserved.
//
// This source code is licensed under the MIT license.
//

#ifndef EFANNA2E_INDEX_H
#define EFANNA2E_INDEX_H

#include <cstddef>
#include <string>
#include <vector>
#include <fstream>
#include "distance.h"
#include "parameters.h"

//#define max_level 1

#define max_width_ 150
#define LL 256
#define OFF 3 //from 64bits
#define min_book 4
#define cen 32
#define nnum 4 //from 128 to 32
#define KK 100
#define fan 10
namespace efanna2e {

class Index {
 public:
  explicit Index(const size_t dimension, const size_t n, Metric metric, const bool flag, float*** book, const int l, float*** quanti, const int level);


  virtual ~Index();

  // virtual void Build0(size_t n, float *data, const Parameters &parameters) = 0;

  // virtual void Search0(
  //     const float *query,
  //     float *x,
  //      size_t k,
	//     const Parameters &parameters,
  //      unsigned *indices) = 0;

  virtual void Save(const char *filename) = 0;

  virtual void Load(const char *filename) = 0;

  inline bool HasBuilt() const { return has_built; }

  inline size_t GetDimension() const { return dimension_; };

  inline size_t GetSizeOfDataset() const { return nd_; }

  inline const float *GetDataset() const { return data_; }

  /*
  virtual void Init(size_t dimension, size_t n, Metric m, Index *initializer, 
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
  */
   protected:
  size_t dimension_;
  float *data_;
  size_t nd_;
  bool has_built;
  Distance* distance_;
  unsigned char* data2_;
  bool is_quan;
  float*** quan_book;
  float*** quantizer;
  //unsigned char** obj_quantizer;
  int length_;
  int max_level;
  size_t* quan_len;
  size_t* quan_size;
  
};

}

#endif //EFANNA2E_INDEX_H
