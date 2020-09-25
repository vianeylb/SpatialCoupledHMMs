// Spatial Coupled HMMs - 
//   Trying out the composite likelihood appraoch


//refer back to coupledHMM_complike_2x2.stan -- this is a test version
data {
  
  //observations
  int<lower=1> TT;
  vector[TT] y1;
  vector[TT] y2;
  
  //states
  int<lower=1> J; 
  
}


parameters {
  
  //state process
  simplex[J] tpm11[J];
  simplex[J] tpm12[J];
  simplex[J] tpm21[J];
  simplex[J] tpm22[J];
  
  vector<lower=0, upper=1>[J] weights;

  //observation process
  vector<lower=0>[J] sigma1;
  ordered[J] mean1;
  vector<lower=0>[J] sigma2;
  ordered[J] mean2;
}


model {

  vector[J] log_tr_tpm11[J];
  vector[J] log_tr_tpm12[J];
  vector[J] log_tr_tpm21[J];
  vector[J] log_tr_tpm22[J];
  
  vector[J] lp1;
  vector[J] lp1_temp;
  vector[J] lp2_temp;
  vector[J] lp2;

  mean1 ~ normal(0, 2);
  mean2 ~ normal(0, 2); 
  sigma1 ~ normal(1, .1);
  sigma2 ~ normal(1, .1);
  
  //taking the logarithm and transposing the tpm
  for(j in 1:J){
    for(k in 1:J){
      log_tr_tpm11[j,k] = log(tpm11[k,j]);
      log_tr_tpm12[j,k] = log(tpm12[k,j]);
      log_tr_tpm21[j,k] = log(tpm21[k,j]);
      log_tr_tpm22[j,k] = log(tpm22[k,j]);
      
      // log_tr_etpm[j,k] = log(ext_tpm[k,j]);
    }
  }

   
  lp1 = rep_vector(0, J);
  lp2 = rep_vector(0, J);
  
  for(j in 1:J){
    lp1[j] = -log(J) + normal_lpdf(y1[1]| mean1[j], sigma1[j]);
    lp2[j] = -log(J) + normal_lpdf(y2[1]| mean2[j], sigma2[j]);
  }
    
  for(t in 2:TT){
    for(j in 1:J){
      
      lp1_temp[j] = log_mix(weights[1], 
                        log_sum_exp(log_tr_tpm11[j] + lp1),
                        log_sum_exp(log_tr_tpm21[j] + lp2)) +
                    normal_lpdf(y1[t]| mean1[j], sigma1[j]);
      
      lp2_temp[j] = log_mix(weights[2], 
                        log_sum_exp(log_tr_tpm22[j] + lp2),
                        log_sum_exp(log_tr_tpm12[j] + lp1)) +
                    normal_lpdf(y2[t]| mean2[j], sigma2[j]);      
  
         
    }
    
    lp1 = lp1_temp;
    lp2 = lp2_temp;

  }


  target += log_sum_exp(lp1);
  target += log_sum_exp(lp2);
 
}

