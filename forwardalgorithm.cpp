/* Forward algorithm for an N-state 
 * hidden Markov model 
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double foralg(int n, int N, arma::mat delta, arma::mat gamma, arma::mat emprobs) {


  
  
  return(loglike);
  
}
  
