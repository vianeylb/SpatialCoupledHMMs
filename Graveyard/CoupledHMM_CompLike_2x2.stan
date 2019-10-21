// Spatial Coupled HMMs - 
//   Trying out the composite likelihood appraoch


//possibly garbage approximate likelihood approach
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
  
  simplex[J] weights[J];

  //observation process
  vector<lower=0>[J] sigma1;
  ordered[J] mean1;
  vector<lower=0>[J] sigma2;
  ordered[J] mean2;
}


model {

  vector[2*J] ext_tpm[2*J];
  vector[2*J] log_tr_etpm[2*J];
  vector[2*J] lp_temp;
  vector[2*J] lp;
  vector[J] lp1;
  vector[J] lp2;
  //vector[J] lp2_temp;
  //vector[J] lp2;

  //constructing the overall matrix with transitions across locations 
  // all possible states
  ext_tpm[1,1] = weights[1,1]*tpm11[1,1];
  ext_tpm[2,1] = weights[1,1]*tpm11[2,1];
  ext_tpm[3,1] = weights[1,2]*tpm21[1,1];
  ext_tpm[4,1] = weights[1,2]*tpm21[2,1];

  ext_tpm[1,2] = weights[1,1]*tpm11[1,2];  
  ext_tpm[2,2] = weights[1,1]*tpm11[2,2];
  ext_tpm[3,2] = weights[1,2]*tpm21[1,2];
  ext_tpm[4,2] = weights[1,2]*tpm21[2,2];
  
  ext_tpm[1,3] = weights[2,1]*tpm12[1,1];
  ext_tpm[2,3] = weights[2,1]*tpm12[2,1];
  ext_tpm[3,3] = weights[2,2]*tpm22[1,1];
  ext_tpm[4,3] = weights[2,2]*tpm22[2,1];

  ext_tpm[1,4] = weights[2,1]*tpm12[1,2];  
  ext_tpm[2,4] = weights[2,1]*tpm12[2,2];
  ext_tpm[3,4] = weights[2,2]*tpm22[1,2];
  ext_tpm[4,4] = weights[2,2]*tpm22[2,2];


  mean1 ~ normal(0, 2);
  mean2 ~ normal(0, 2); 
  sigma1 ~ normal(1, .1);
  sigma2 ~ normal(1, .1);
  
  //taking the logarithm and transposing the tpm
  for(j in 1:(2*J)){
    for(k in 1:(2*J)){
      log_tr_etpm[j,k] = log(ext_tpm[k,j]);
    }
  }

   
  lp = rep_vector(0, 2*J);
  for(j in 1:J){
    lp[j] = -log(J) + normal_lpdf(y1[1]| mean1[j], sigma1[j]);
    lp[j+J] = -log(J) + normal_lpdf(y2[1]| mean2[j], sigma2[j]);
  }
    
  for(t in 2:TT){
    for(j in 1:J){
      lp_temp[j] = log_sum_exp(log_tr_etpm[j] + lp) + 
         normal_lpdf(y1[t]| mean1[j], sigma1[j]);
      lp_temp[j+J] = log_sum_exp(log_tr_etpm[j+J] + lp) + 
         normal_lpdf(y2[t]| mean2[j], sigma2[j]);
    }
    
    lp = lp_temp;

  }

  for(j in 1:J){
    lp1[j] = lp[j];
    lp2[j] = lp[j + J];
  }


  target += log_sum_exp(lp1);
  target += log_sum_exp(lp2);
 
}

