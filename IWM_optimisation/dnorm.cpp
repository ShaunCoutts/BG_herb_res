#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector my_dnorm(NumericVector x, NumericVector means, double sd) {
  int n = x.size();
  NumericVector res(n);
  double x0;
  double c = 1/sqrt(2*PI);
  for (int i=0; i<n; i++) {
    x0 = x[i] - means[i];
    res[i] = exp(-x0*x0/(2*sd*sd))*c/sd;
  }
  return res;
}
