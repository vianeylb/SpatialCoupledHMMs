/* Forward algorithm using the log_sum_exp function to ensure
 * numerical stability
 * 
 * log_tr_gamma = transition probability matrix for a non-homogeneous process
 *                transposed and log transformed at each time 't'
 * 
 * log_allprobs = logarithm of the density evaluated for at observation, under
 *                each state process
 *                
 * log_delta  = log initial distribution 
 * 
 */


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_sum_exp(NumericVector foo1, NumericVector foo2){

  double msf;
  NumericVector sumfoos = foo1 + foo2;
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
double foralg_alt(int n, int N, NumericVector log_foo, NumericMatrix log_tr_gamma, NumericMatrix log_allprobs) {

  double lscale;
  NumericVector foo1(N);
  NumericVector foo2(N);
  NumericVector zero(N);

  for(int j=0; j<N; j++){ 
    foo1[j] = log_foo[j] + log_allprobs(0,j);
  }

  for (int i=1; i < n; i++){
    for(int j=0; j<N; j++){ 
      foo2[j] = log_sum_exp(foo1, log_tr_gamma.row(j)) + log_allprobs(i,j);
    }
    
    for(int j=0; j < N; j++){
      foo1[j] = foo2[j];
    }
  }
    

  lscale = log_sum_exp(foo1, zero);
  
  return(lscale);
  
}


