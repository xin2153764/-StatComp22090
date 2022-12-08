#include <Rcpp.h>
using namespace Rcpp;

//' @name gibbs
//' @title gibbs sample by Rcpp
//' @description Generate the gibbs sample by Rcpp,using the method of gibbs sampling.
//' @param N the sample size
//' @return a Numeric Matrix,sp stands for random sample of multivariate distribution.
//' @examples
//' \dontrun{
//' res <- gibbsC(5000)
//' }
//' @export
//[[Rcpp::export]]
NumericMatrix gibbsC(int N) {
  int burn = 1000;
  double mu1 = 0;
  double mu2 = 0;
  double sigma1 = 1;
  double sigma2 = 1;
  double rho = 0.9;
  NumericMatrix mat(N, 2);
  double x = mu1, y = mu2;
  mat(1,0) = x;
  mat(1,1) = y;
  double s1 = sqrt(1-pow(rho,2))*sigma1;
  double s2 = sqrt(1-pow(rho,2))*sigma2;
  for(int i = 1; i < N; i++) {
    for(int j = 0; j < burn; j++) {
      y = mat(j-1, 2);
      double m1 = 0 + rho * (y - 0) * sigma1/sigma2;
      x = rnorm(1, m1, s1)[0];
      double m2 = 0 + rho * (x - 0) * sigma2/sigma1;
      y = rnorm(1, m2, s2)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}  