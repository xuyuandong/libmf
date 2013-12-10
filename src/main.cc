#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include "base/flags.h"
#include "base/logging.h"
#include "problem.h"

DEFINE_string(matrix, "", "matrix for factorization");
DEFINE_int32(iteration, 20, "optimization iteration number");

int main(int argc, char* argv[]) {
  base::ParseCommandLineFlags(&argc, &argv, true);    
  srand( (unsigned)time( NULL ) );
 
  mf::Problem ptr;
  
  if (FLAGS_matrix.empty()) {
    ptr.SetMatrix(std::cin);
  } else {
    std::ifstream fin(FLAGS_matrix.c_str());
    ptr.SetMatrix(fin);
    fin.close();
  }

  ptr.PrintInfo();
  ptr.Optimize(FLAGS_iteration);
  ptr.PrintResult();  
  
  return 0;
}
