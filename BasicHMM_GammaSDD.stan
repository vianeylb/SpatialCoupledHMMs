// Basic HMM in Stan

data {
  int<lower=1> N; // number of states
  int<lower=1> T; // length of data set
  real y[T]; // observations
}

parameters {
  simplex[N] tpm[N]; // N x N tpm
  positive_ordered[N] mu; // state-dependent parameters
  vector<lower=0>[N] sigma;
  
  simplex[N] init;
}  

model {

  //log of the transposed tpm
  vector[N] log_tpm_tr[N];
  
  //forward variables
  vector[N] lp;
  vector[N] lp_p1;

  //transformed parameters
  vector[N] shape = mu .* mu ./ (sigma .* sigma);
  vector[N] rate =  mu ./ (sigma .* sigma);

  
  // prior for mu - non-exchangeable preferred
  mu[1] ~ normal(2, 1);
  mu[2] ~ normal(10, 1);
  mu[3] ~ normal(20, 1);
  //mu[4] ~ normal(15, 2);
  
  //prior for sigma - non-exchangeable preferred
  sigma ~ student_t(3, 0, 1);

  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);

  
  // forward algorithm implementation
  
  //alpha_1
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + gamma_lpdf(y[1] | shape[n], rate[n]);

  //alpha_t
  for (t in 2:T) { // looping over observations
    for (n in 1:N){ // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp);
      
      if(y[t] < 900)
       lp_p1[n] = lp_p1[n] +  normal_lpdf(y[t] | shape[n], rate[n]); 
    }
    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}


