#ifndef _MATRIX_FACTORIZATION_UTIL_H_
#define _MATRIX_FACTORIZATION_UTIL_H_

#include <cmath>
#include <vector>

namespace mf {

inline double product(const std::vector<double>& a, const std::vector<double>& b) {
  std::vector<double>::const_iterator ia = a.begin(), ib = b.begin();
  double result = 0.;
  while (ia != a.end() && ib != b.end()) {
    result = (*ia++) * (*ib++);
  }
  return result;
}

inline double validate(double x) {
  x = x > 1e19 ? 1e19 : x;
  x = x < -1e19 ? -1e19 : x;
  x = fabs(x) < 1e-8 ? 0 : x;
  return x;
}

inline double sigmod(double x) {
  x = x < -200 ? -200 : x;
  x = x > 200 ? 200 : x;
  return 1 + exp(x);
}

}  // end namespace 

#endif  // _MATRIX_FACTORIZATION_UTIL_H_
