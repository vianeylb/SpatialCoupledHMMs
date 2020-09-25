// Basic AR(1) HMM in Stan

data {

  int<lower=1> N; // number of states
  int<lower=1> T; // length of data set
  real y[T]; // observations

} 

parameters {
  simplex[N] tpm[N]; // N x N tpm
  vector<lower=0>[N] sigma;
  vector[N] mu; // state-dependent parameters
  simplex[N] init;
  real rho[N];
}  

model {

  //log of the transposed tpm
  vector[N] log_tpm_tr[N];
  
  //forward variables
  vector[N] lp;
  vector[N] lp_p1;
  
  // prior for mu 
   mu ~ normal(5, 3);
   
  //non-exchangeable means lead to multiple posterior modes
  // a prior that 'shrinks' the larger and smaller means toward
  // a center point (like an avg of the data/process) works well

  //prior for sigma 
  //sigma ~ student_t(3, 0, 1);
  sigma ~ normal(0.1, 0.3);
  
  // placing most of the probability mass on within state 
  // stationarity -- doesn't necessarily lead to stationary HMM
  rho ~ normal(0,0.3);
  
  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);

  
  // forward algorithm implementation
  
  //alpha_1
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[1] | mu[n], sigma[n]);

  //alpha_t
  for (t in 2:T) { // looping over observations
    for (n in 1:N){ // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp);
      lp_p1[n] = lp_p1[n] +  normal_lpdf(y[t] | mu[n] + rho[n]*(y[t-1]-mu[n]), sigma[n]); 
    }
    lp = lp_p1;
  }

  target += log_sum_exp(lp);

}

generated quantities{
  
  //TRYING TO GENERATE FORECASTS FOR THE NEXT 24 OBSERVATIONS FROM 
  // A TIME-HOMOGENEOUS HMM
  
  vector[24] y_pred;
  //forward variables
  vector[N] lp;
  vector[N] lp_p1;
  //log of the transposed tpm
  vector[N] log_tpm_tr[N];
  vector[N] stateProbT;
  
  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);
  
  //repeating the forward algorithm -- this is probably a dumb way to do this
  
  //alpha_1
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + normal_lpdf(y[1] | mu[n], sigma[n]);

  //alpha_t
  for (t in 2:T) { // looping over observations
    for (n in 1:N){ // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp);
      lp_p1[n] = lp_p1[n] +  normal_lpdf(y[t] | mu[n] + rho[n]*(y[t-1]-mu[n]), sigma[n]); 
    }
    lp = lp_p1;
  }
  
  for(n in 1:N){
    stateProbT[n] = lp[n]/sum(lp);
  }
  
  
  
  
}




