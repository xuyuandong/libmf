#ifndef _MATRIX_FACTORIZATION_VERTEX_H
#define _MATRIX_FACTORIZATION_VERTEX_H

#include <string>
#include <vector>
#include "util.h"

namespace mf {

static const int _MAX_ = 100000;

class Vertex {
  public:
    Vertex(const std::string& id, int D = 100) {
      id_ = id;
      weights_.resize(D);
      for (int i = 0; i < D; ++i) {
        weights_[i] = (double(rand()%(2*_MAX_)) / _MAX_ - 1) / sqrt(D);
      }
    }
    inline bool operator==(const Vertex& v) {
      return this->id_ == v.id_;
    }

  public:
    std::string id_;
    std::vector<double> weights_;
};

};  // end namespace

#endif  // _MATRIX_FACTORIZATION_VERTEX_H
