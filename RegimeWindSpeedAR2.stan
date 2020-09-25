functions{
  
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb;
    real p_ub;
    real u;
    real y;
    
    p_lb = normal_cdf(lb, mu, sigma);
    p_ub = normal_cdf(ub, mu, sigma);
    if(p_ub == 0) p_ub = p_ub + 0.00001;
    if(p_lb == 1) p_lb = p_lb - 0.00001;
  
    u = uniform_rng(p_lb, p_ub);
    y = mu + sigma * inv_Phi(u);
    
    return y;
  }
}

data {
  int<lower=1> TT;
  vector[TT] obs;
  int regime[TT];
  int<lower=1> Nregime;

  real regimelimits[Nregime,2];
  
  int counts[Nregime, Nregime];
  
}


parameters {
  //ordered[Nregime] mu;
  vector<lower=0>[Nregime] sigma;
  vector[2] rho[Nregime];  
  
  simplex[Nregime] tpm[Nregime];
  
}


model {
  
  sigma ~ normal(0, 1);
  
  for(j in 1:Nregime)
   rho[j] ~ normal(0, 0.3);
  
  //obs[1] ~ normal(, sigma[regime[1]])T[regimelimits[regime[1], 1], regimelimits[regime[1], 2]];
  
  for(t in 3:TT){
    if(obs[t] < 1000){
      if(obs[t-1] < 1000){
        
          if(regime[t] == 1){
            obs[t] ~ normal(rho[regime[t],1]*obs[t-1] + rho[regime[t],2]*obs[t-2], 
                            sigma[regime[t]])T[,regimelimits[regime[t],2]];            
          } else {
            obs[t] ~ normal(rho[regime[t],1]*obs[t-1] + rho[regime[t],2]*obs[t-2], 
                            sigma[regime[t]])T[regimelimits[regime[t],1],regimelimits[regime[t],2]];
          }
      }
    } 
  }
  
  for(n in 1:Nregime){
    counts[n] ~ multinomial_lpmf(tpm[n]);
  }
  
}


generated quantities{

  // forecasts for the next hour
  vector[12] obs_forecast;
  int regime_forecast[12];

  regime_forecast[1] = categorical_rng(tpm[regime[TT]]);
  //obs_forecast[1] = normal_rng(rho[regime[TT]]*obs[TT], sigma[regime[TT]])T(regimelimits[TT,1], regimelimits[TT,2]);
  obs_forecast[1] = normal_lub_rng(rho[regime_forecast[1],1]*obs[TT] + rho[regime_forecast[1],2]*obs[TT-1], 
                                   sigma[regime_forecast[1]],
                                   regimelimits[regime_forecast[1],1],
                                   regimelimits[regime_forecast[1],2]);

  
  regime_forecast[2] = categorical_rng(tpm[regime_forecast[1]]);
  //obs_forecast[1] = normal_rng(rho[regime[TT]]*obs[TT], sigma[regime[TT]])T(regimelimits[TT,1], regimelimits[TT,2]);
  obs_forecast[2] = normal_lub_rng(rho[regime_forecast[2],1]*obs_forecast[1] + 
                                   rho[regime_forecast[2],2]*obs[TT], 
                                   sigma[regime_forecast[2]],
                                   regimelimits[regime_forecast[2],1],
                                   regimelimits[regime_forecast[2],2]);

  
  for(j in 3:12){

    regime_forecast[j] = categorical_rng(tpm[regime_forecast[j-1]]);

    obs_forecast[j] = normal_lub_rng(rho[regime_forecast[j],1]*obs_forecast[j-1] + 
                                     rho[regime_forecast[j],2]*obs_forecast[j-2], 
                                     sigma[regime_forecast[j]],
                                     regimelimits[regime_forecast[j],1],
                                     regimelimits[regime_forecast[j],2]);
  }

}



