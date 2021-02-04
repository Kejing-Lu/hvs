//
// Copyright (c) 2017 ZJULearning. All rights reserved.
//
// This source code is licensed under the MIT license.
//
#include <efanna2e/index.h>
namespace efanna2e {
  Index::Index(const size_t dimension, const size_t n, Metric metric, const bool flag, float*** book, const int l, float*** quanti, const int max_level_)
    : dimension_ (dimension), nd_(n), has_built(false), is_quan(flag), quan_book(book), length_(l), quantizer(quanti), max_level(max_level_) {
  /*
    switch (metric) {
      case L2:distance_ = new DistanceL2();
        break;
      default:distance_ = new DistanceL2();
        break;
    }
  */
    quan_len = new size_t[max_level];
    quan_size = new size_t[max_level];
  distance_ = new DistanceL2();
}
Index::~Index() {}
}
