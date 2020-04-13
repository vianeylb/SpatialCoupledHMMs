data {
  int<lower=1> TT;
  vector[TT] obs;
  int regime[TT];
  int<lower=1> Nregime;

  real regimelimits[Nregime,2];
}


parameters {
  //ordered[Nregime] mu;
  vector<lower=0>[Nregime] sigma;
  vector[Nregime] rho;  
}


model {
  
  //mu ~ normal(0, 1);
  sigma ~ student_t(3, 0, .5);
  rho ~ normal(0, 0.3);
  
  //obs[1] ~ normal(, sigma[regime[1]])T[regimelimits[regime[1], 1], regimelimits[regime[1], 2]];
  
  for(t in 2:TT){
    if(obs[t] < 5){
      if(obs[t-1] < 5){
        
          if(regime[t] == 1){
            obs[t] ~ normal(rho[regime[t]]*obs[t-1], sigma[regime[t]])T[,regimelimits[regime[t],2]];            
          } else {
            obs[t] ~ normal(rho[regime[t]]*obs[t-1], 
                            sigma[regime[t]])T[regimelimits[regime[t],1],regimelimits[regime[t],2]];
          }
      }
    } 
  }
  
  
  
}
