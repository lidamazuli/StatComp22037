#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples (int)
//' @param burn burn-in length, the length of deserted samples from the first
//' @param mu1 the mean of first marginal value 
//' @param mu2 the mean of second marginal value 
//' @param sigma1 the variance of first marginal value 
//' @param sigma2 the variance of second marginal value 
//' @param rho correlation
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rnC <-gibbsC(5000,2,0,0,1,1,.9)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin, double mu1, double mu2,double sigma1,double sigma2,double rho) {
  NumericMatrix mat(N, 2);
  double m1,m2;
  double x = mu1, y = mu2;
  mat(1,0) = x;
  mat(1,1) = y;
  double s1 = sqrt(1-pow(rho,2))*sigma1;
  double s2 = sqrt(1-pow(rho,2))*sigma2;
  for(int i = 1; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      y = mat(j-1, 2);
      m1 = 0 + rho * (y - 0) * sigma1/sigma2;
      x = rnorm(1, m1, s1)[0];
      m2 = 0 + rho * (x - 0) * sigma2/sigma1;
      y = rnorm(1, m2, s2)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
