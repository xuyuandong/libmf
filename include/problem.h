#ifndef _MATRIX_FACTORIZATION_PROBLEM_H_
#define _MATRIX_FACTORIZATION_PROBLEM_H_

#include <map>
#include <string>
#include <iostream>
#include "vertex.h"

namespace mf {

struct ViewClick {
  int view_; 
  int click_;
};

class Problem {
  public:
    Problem();
    ~Problem();

    bool SetMatrix(std::istream& fin);
    bool Optimize(int iter_num = 100);

    void PrintInfo();
    void PrintResult();
    
    int GetDims() {
      return latent_dims_ * (ad_vec_.size() + ps_vec_.size());
    }

  public:
    // for data
    std::vector<Vertex> ad_vec_;
    std::vector<Vertex> ps_vec_;
    std::vector<std::map<int, ViewClick> > graph_;
    int graph_size_;
    // for model
    int latent_dims_;
    double regular_alpha_;
    double regular_beta_;
    // for optimization
    bool opt_flag_;
    double loss_;
    // raw data
    std::map<std::string, int> ad_map_;
    std::map<std::string, int> ps_map_;
};


}  // end namespace

#endif  // _MATRIX_FACTORIZATION_PROBLEM_H_
