// Basic HMM in Stan

data {
  int<lower=1> N; // number of states
  int<lower=1> T; // length of data set
  real y[T]; // observations
  
  //priors for tpm
  vector[N] alpha[N];
  
}

parameters {
  simplex[N] tpm[N]; // N x N tpm
  vector<lower=0>[N] sigma;
 //ordered[N] mu; // state-dependent parameters
  simplex[N] init;
  real rho[N];
}  

model {

  //log of the transposed tpm
  vector[N] log_tpm_tr[N];
  
  //forward variables
  vector[N] lp;
  vector[N] lp_p1;
  
  // prior for mu - non-exchangeable preferred
   //mu ~ student_t(3, 0, 1);
   // mu[4] ~ normal(15, 1);
  //mu[6] ~ normal(4, 0.3);
  
  //prior for sigma - non-exchangeable preferred
  //sigma ~ student_t(3, 0, 1);
  sigma ~ normal(0.1, 0.3);
  
  rho ~ normal(0,0.3);
  
  for(n in 1:N){
    tpm[n] ~ dirichlet(alpha[n]);
  }

  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);

  
  // forward algorithm implementation
  
  //alpha_1
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[1] | 0, sigma[n]);

  //alpha_t
  for (t in 2:T) { // looping over observations
    for (n in 1:N){ // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp);
      lp_p1[n] = lp_p1[n] +  normal_lpdf(y[t] |  rho[n]*y[t-1], sigma[n]); 
    }
    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}


