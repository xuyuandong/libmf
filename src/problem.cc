#include <iostream>
#include "IpIpoptApplication.hpp"
#include "base/flags.h"
#include "base/logging.h"
#include "base/string_util.h"
#include "problem.h"
#include "solution.h"
#include "util.h"

using namespace std;
using namespace Ipopt;

DEFINE_string(ipopt, "ipopt.opt", "config file for ipopt optimization options");
DEFINE_bool(filter, false, "filter out click = 0 samples");

namespace mf {

Problem::Problem() {
  opt_flag_ = false;
  latent_dims_ = 40;
  regular_alpha_ = .5;
  regular_beta_ = .5;
}

Problem::~Problem() {
}

bool Problem::SetMatrix(istream& fin) {
  string line;
  graph_size_ = 0;
  while(getline(fin, line)) {
    vector<string> tokens;
    SplitString(line, '\t', &tokens);
    if (tokens.size() != 3) {
      VLOG(1) << "unexpected token size " << line;
      continue;
    }
    string name = tokens[0];
    int click = StringToInt(tokens[1]);
    int view = StringToInt(tokens[2]);
    if (FLAGS_filter && click == 0) {
      VLOG(2) << "filter out click=0 sample";
      continue;
    }
    
    tokens.clear();
    SplitString(name, '#', &tokens);
    if (tokens.size() != 2) {
      VLOG(1) << "unexpected token size " << name;
      continue;
    }
    string ad = tokens[0];
    string ps = tokens[1];

    if (ad_map_.find(ad) == ad_map_.end()) {
      ad_map_[ad] = ad_vec_.size();
      Vertex v(ad, latent_dims_);
      ad_vec_.push_back(v);
    } 
    if (ps_map_.find(ps) == ps_map_.end()) {
      ps_map_[ps] = ps_vec_.size();
      Vertex v(ps, latent_dims_);
      ps_vec_.push_back(v);
    }

    int ad_id = ad_map_[ad];
    int ps_id = ps_map_[ps];
    while(graph_.size() <= (size_t)ad_id) {
      map<int, ViewClick> edge;
      graph_.push_back(edge);
    } 
    map<int, ViewClick>& edge = graph_[ad_id];
    edge[ps_id].view_ = view;
    edge[ps_id].click_ = click;

    graph_size_ += 1;
  }
  return true;
}

bool Problem::Optimize(int iter_num) {
  int t = 0;
  double last_loss = 1e+10;
  vector<double> last_vec;
  last_vec.resize(latent_dims_ * (ad_vec_.size() + ps_vec_.size()));
  
  for (int n = 0; n < iter_num; ++n) {
    // backup
    VLOG(3) << "backup latent variables ...";
    t = 0;
    for (size_t i = 0; i < ad_vec_.size(); ++i) {
      for (int j = 0; j < latent_dims_; ++j) {
        last_vec[t++] = ad_vec_[i].weights_[j];
      }
    }
    for (size_t i = 0; i < ps_vec_.size(); ++i) {
      for (int j = 0; j < latent_dims_; ++j) {
        last_vec[t++] = ps_vec_[i].weights_[j];
      }
    }
    
    // optimization
    VLOG(3) << "optimize by ipopt ...";
    {
      SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
      ApplicationReturnStatus status = app->Initialize(FLAGS_ipopt);
      app->Options()->SetNumericValue("tol", 1e-2);
      app->Options()->SetStringValue("mu_strategy", "adaptive");
      app->Options()->SetStringValue("hessian_approximation", "limited-memory"); 
      if (status != Solve_Succeeded) {
        VLOG(1) << "initialize ipopt error";
        return (int) status;
      }
      
      SmartPtr<Solution> mynlp = new Solution();
      mynlp->SetProblem(this);
      status = app->OptimizeTNLP(mynlp);
      if (status == Solve_Succeeded) {
        cerr << "Converge at opt_flag=" << opt_flag_ << endl;
      } else {
        cerr << "!!!Warning: Unconverge at opt_flag=" << opt_flag_ << endl;
      }
    }
    
    // switch flag
    VLOG(3) << "switch optimization flag ...";
    opt_flag_ = !opt_flag_;
    
    // check step
    VLOG(3) << "check optimization step ...";
    double max_step = 0.0;
    t = 0;
    for (size_t i = 0; i < ad_vec_.size(); ++i) {
      for (int j = 0; j < latent_dims_; ++j) {
        double step = fabs(last_vec[t++] - ad_vec_[i].weights_[j]);
        max_step = max_step < step ? step : max_step;
      }
    }
    for (size_t i = 0; i < ps_vec_.size(); ++i) {
      for (int j = 0; j < latent_dims_; ++j) {
        double step = fabs(last_vec[t++] - ps_vec_[i].weights_[j]);
        max_step = max_step < step ? step : max_step;
      }
    }

    // output info
    cerr << "Iteration=" << n 
        << " , Loss=" << loss_ 
        << " , step=" << max_step << endl << endl;

    // step too small
    double threshold = 1 / (latent_dims_ * graph_size_);
    if (max_step < threshold) {
      cerr << "Step goes too small: " << max_step <<", terminate" << endl;
      break;
    }

    // loss too small
    if (fabs(loss_ - last_loss) < 10) {
      cerr << "Loss converge <" 
          << loss_ << " - " << last_loss 
          << ">, terminate" << endl;
      break;
    }
    last_loss = loss_;
  } 

  return true;
}

void Problem::PrintInfo() {
  VLOG(2) << "Matrix Factorization Info:";
  cerr << "ad size: " << ad_vec_.size() << endl;
  cerr << "ps size: " << ps_vec_.size() << endl;
  cerr << "graph size: " << graph_size_ << endl;
  cerr << "latent dimension size: " << latent_dims_ << endl;
  cerr << "regular alpha: " << regular_alpha_ << endl;
  cerr << "regular beta: " << regular_beta_ << endl;
}

void Problem::PrintResult() {
  VLOG(2) << "Matrix Factorization Result:";
  map<string, int>::iterator it;
  for (it = ad_map_.begin(); it != ad_map_.end(); ++it) {
    cout << it->first;
    const Vertex& v = ad_vec_[it->second];
    for (int j = 0; j < latent_dims_; ++j) {
      cout << " " << v.weights_[j];
    }
    cout << endl;
  }
  for (it = ps_map_.begin(); it != ps_map_.end(); ++it) {
    cout << it->first;
    const Vertex& v = ps_vec_[it->second];
    for (int j = 0; j < latent_dims_; ++j) {
      cout << " " << v.weights_[j];
    }
    cout << endl;
  }
}


}  // end namespace 
