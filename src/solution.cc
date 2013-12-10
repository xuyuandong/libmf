#include <cassert>
#include <iostream>
#include "base/logging.h"
#include "solution.h"

using namespace std;
using namespace Ipopt;

namespace mf {

Solution::Solution() {}

Solution::~Solution() {}

// returns the size of the problem
bool Solution::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, IndexStyleEnum& index_style) {
  n = problem_->GetDims();
  m = 0;  //no constraint
  nnz_jac_g = 0;
  nnz_h_lag = 0;  //to use L-BFGS, set hessian to zero
  index_style = TNLP::C_STYLE;
  return true;
}

// returns the variable bounds
bool Solution::get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u) {
  assert(n == problem_->GetDims());
  assert(m == 0);
  for (Index i = 0; i < n; i++) {
    x_l[i] = -2;
    x_u[i] = 2;
  }
  return true;
}

// returns the initial poIndex for the problem
bool Solution::get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda, Number* lambda) {
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  int D = problem_->latent_dims_;
  for(size_t i = 0; i < problem_->ad_vec_.size(); i++){
    for(int j = 0; j < D; j++){
      x[i * D + j] = problem_->ad_vec_[i].weights_[j];
    }
  }

  int offset = problem_->ad_vec_.size() * D;
  for(size_t i = 0; i < problem_->ps_vec_.size(); i++){
    for(int j = 0; j < D; j++){
      x[offset + i * D + j] = problem_->ps_vec_[i].weights_[j];
    }
  }

  return true;
}

// returns the value of the objective function
bool Solution::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
  assert(n == problem_->GetDims());
    
  int D = problem_->latent_dims_;
  int offset = D * problem_->ad_vec_.size();

  double loss = 0;
  for (size_t i = 0; i < problem_->graph_.size(); i++) {
    const map<int, ViewClick>& edge = problem_->graph_[i];
    for (map<int, ViewClick>::const_iterator ia = edge.begin(); 
            ia != edge.end(); ++ia) {
      double click = ia->second.click_;
      double noclick = ia->second.view_ - click;

      int start_offset_ad = i * D;
      int start_offset_ps = ia->first * D + offset;

      double dot = 0;  // product of latent vectors
      for(int j = 0; j < D; j++) {
          dot += x[start_offset_ad + j] * x[start_offset_ps + j];
      }
      if (dot > 100 || dot < -100) {
        VLOG(5) << "debug: dot=" << dot;
      }
      double ctr=sigmod(-dot);
      if (ctr > 1000000) {
        VLOG(5) << "debug: ctr=" << ctr;
      }
      loss += (click * log(ctr) - noclick * log(1 - 1/ctr));
    }   
  } 
  loss = loss / problem_->graph_size_;

  //regulization part
  double loss_alpha = 0, loss_beta = 0;
  for(size_t i = 0; i < problem_->ad_vec_.size(); i++) {
    for(int j = 0; j < D; j++) {
      loss_beta += x[i*D+j] * x[i*D+j];
    }
  }
  for(size_t i = 0; i < problem_->ps_vec_.size(); i++) {
    for(int j = 0; j < D; j++) {
      loss_alpha += x[offset+i*D+j] * x[offset+i*D+j];
    }
  }
  loss += 0.5 * (loss_alpha * problem_->regular_alpha_ + loss_beta * problem_->regular_beta_);
    
  obj_value = loss;
  problem_->loss_ = loss;
  VLOG(4) << "get objective function value " << loss;
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool Solution::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
  VLOG(4) << "get objective function gradient ...";
  assert(n == problem_->GetDims());
  
  int D = problem_->latent_dims_;
  int offset = problem_->ad_vec_.size() * D;
  
  // reset state
  for (int i = 0; i < n; i++) {
    grad_f[i] = 0;
  }

  // gradient from regulized part
  if (problem_->opt_flag_) {
    for (size_t i = 0; i < problem_->ad_vec_.size(); i++) {
      for (int j = 0; j < D; j++){
        grad_f[i*D+j] += x[i*D+j] * problem_->regular_beta_;
      }
    }
  } else {
    for (size_t i = 0; i < problem_->ps_vec_.size(); i++) {
      for (int j = 0; j < D; j++) {
        grad_f[offset+i*D+j] += x[offset+i*D+j] * problem_->regular_alpha_;
      }
    }
  }

  // now do the main part
  for(size_t i = 0; i < problem_->graph_.size(); i++) {
    const map<int, ViewClick>& edge = problem_->graph_[i];
    for (map<int, ViewClick>::const_iterator ia = edge.begin(); 
            ia != edge.end(); ++ia) {
      double click = ia->second.click_;
      double view = ia->second.view_;

      int start_offset_ad = i * D;
      int start_offset_ps = ia->first * D + offset;

      double dot = 0;
      for (int j = 0; j < D; j++) {
        dot += x[start_offset_ad+j] * x[start_offset_ps+j];
      }

      double ctr = 1 / sigmod(-dot);
      if (problem_->opt_flag_) {
        for(int j = 0; j < D; j++){
          grad_f[start_offset_ad+j] += (view * ctr - click) * x[start_offset_ps+j] / problem_->graph_size_;
        }
      } else {
        for(int j = 0; j < D; j++){
          grad_f[start_offset_ps+j] += (view * ctr - click) * x[start_offset_ad+j] / problem_->graph_size_;
        }
      }
    } // end for map iter
  } // end for graph vec
 
  return true;
}


// return the value of the constraIndexs: g(x)
bool Solution::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
  return true;
}

// return the structure or values of the jacobian
bool Solution::eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow, Index *jCol,
        Number* values) {
  return true;
}

void Solution::finalize_solution(SolverReturn status,
        Index n, const Number* x, const Number* z_L, const Number* z_U,
        Index m, const Number* g, const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq) {
  VLOG(3) << "get final result values ...";
  int D = problem_->latent_dims_;
  int offset = problem_->ad_vec_.size() * D;
  for (size_t i = 0; i < problem_->ad_vec_.size(); i++) {
    for (int j = 0; j < D; j++) {
      problem_->ad_vec_[i].weights_[j] = x[i*D+j];
    }   
  }
  for (size_t i = 0; i < problem_->ps_vec_.size(); i++) {
    for (int j = 0; j < D; j++) {
      problem_->ps_vec_[i].weights_[j] = x[offset+i*D+j];
    }
  }
}




}  /// end namespace
