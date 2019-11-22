/* Forward/backward algorithm using the log_sum_exp function to ensure
 * numerical stability
 * 
 * N = # of states in HMM
 * 
 * n = length of time series
 * 
 * log_gamma = transition probability matrix for a homogeneous state process
 *             log transformed 
 *             :: NxN matrix
 * 
 * log_tr_gamma = transition probability matrix for a homogeneous state process
 *                transposed and log transformed 
 *                :: NxN matrix
 * 
 * log_allprobs = logarithm of the state-dependent densities evaluated 
 *                at each observation 
 *                :: n x N matrix
 *                
 * log_delta  = log initial distribution 
 *              :: vector of length N
 *              
 */


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_sum_exp(arma::vec foo1, arma::vec foo2){

  double msf;
  arma::vec sumfoos = foo1 + foo2;
  double exps;
  double lse;
  int n = foo1.size();

  msf = max(sumfoos);
  exps = 0.0; 
  for(int i=0; i< n; i++){
    exps+=exp(sumfoos[i] -msf);
  }

  lse = msf + log(exps);
  
  return(lse);
}

// [[Rcpp::export]]
arma::mat forbackalg(int n, int N, arma::mat log_foo, arma::vec log_tr_gamma, arma::vec log_gamma, arma::mat log_allprobs) {

  
  arma::cube log_tr_gamma_array(log_tr_gamma.begin(), N, N, n-1, false);
  arma::cube log_gamma_array(log_gamma.begin(), N, N, n-1, false);
  
  double lscale;
  arma::vec zero(N); /* defaults to a vector of zeros */
  arma::mat logalpha(n, N); /* a numeric matrix of zeros */
  arma::mat logbeta(n, N); /* a numeric matrix of zeros */
  arma::mat stateprobs(n, N); /* a numeric matrix of zeros */
  
  /*NumericVector ltga(N);
  NumericVector ltg(N);*/
  
  /* forward algorithm */
  
  for(int j=0; j<N; j++){ 
    logalpha(0,j)= log_foo[j] + log_allprobs(0,j);
  }
  
  for (int i=1; i < n; i++){
    for(int j=0; j<N; j++){ 
      logalpha(i,j) = log_sum_exp(logalpha.row(i-1), log_tr_gamma_array.slice(j)) + log_allprobs(i,j);
    }
  }
  
  
  /* likelihood */
  
  lscale = log_sum_exp(logalpha.row(n-1), zero);
  
  /*backward algorithm */
  
   for(int j=0; j<N; j++){ 
     logbeta(n-1,j) = 0;
  }
  
  for (int i = n-2 ; i > -1; i--){
    for(int j=0; j<N; j++){ 
      logbeta(i,j) = log_sum_exp(logbeta.row(i+1), log_gamma_array.slice(j)) + log_allprobs(i+1,j);
    }
    
  }
  
  /* evaluating state probabilities */
  
  for(int i=0; i < n; i++){
    for(int j=0; j< N; j++){
      stateprobs(i,j) = exp(logalpha(i,j) + logbeta(i,j) - lscale);
    }
  }
  
  
  return(stateprobs);
  
}


