// Basic HMM in Stan

data {
  int<lower=1> N; // number of states
  int<lower=1> T; // length of data set
  real y[T]; // observations
}

parameters {
  simplex[N] tpm[N]; // N x N tpm
 // vector[N] mu; // state-dependent parameters
  vector<lower=0>[N] sigma;
  
  simplex[N] init;
  vector[4] rho[N];
}  

model {

  //log of the transposed tpm
  vector[N] log_tpm_tr[N];
  
  //forward variables
  vector[N] lp;
  vector[N] lp_p1;
  
  // prior for mu - non-exchangeable preferred
  //mu ~ normal(1 ,1);
  //mu[2] ~ normal(1.5, 0.3);
  //mu[3] ~ normal(2, 0.3);
  //mu[3] ~ normal(3, 0.3);
  //mu[5] ~ normal(3, 0.3);
  //mu[5] ~ normal(3, 0.3);
  //mu[6] ~ normal(4, 0.3);
  
  //prior for sigma - non-exchangeable preferred
  //sigma ~ student_t(3, 0, 1);
  sigma ~ normal(0, 1);
  
  for(n in 1:N)  rho[n] ~ normal(0,0.3);

  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);

  
  // forward algorithm implementation
  
  //alpha_1
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[1] | 0, sigma[n]);

  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[2] | rho[n,1]*(y[1]-0), sigma[n]);

  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[3] | rho[n,1]*(y[2]-0) + 
                                              rho[n,2]*(y[1]-0), sigma[n]);
  
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[4] | rho[n,1]*(y[3]-0) + 
                                              rho[n,2]*(y[2]-0) + 
                                              rho[n,3]*(y[1]-0), sigma[n]);

  //alpha_t
  for (t in 5:T) { // looping over observations
    for (n in 1:N){ // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp);
      
       lp_p1[n] = lp_p1[n] +  normal_lpdf(y[t] | 0 + rho[n,1]*(y[t-1]-0) + 
                                                           rho[n,2]*(y[t-2]-0)+ 
                                                           rho[n,3]*(y[t-3]-0)+ 
                                                           rho[n,4]*(y[t-4]-0), sigma[n]); 
    }
    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}


