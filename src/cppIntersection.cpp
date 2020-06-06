#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> cppIntersection(const IntegerVector &A, const IntegerVector &B) {
  std::vector<int> output;
  std::set_intersection(A.begin(), A.end(), B.begin(), B.end(),
                        std::back_inserter(output));
  return output;
}
