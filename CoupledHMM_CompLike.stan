// Spatial Coupled HMMs - 
//   Trying out the composite likelihood appraoch

data {
  
  //observations
  int<lower=1> N;
  vector[N] y1;
  vector[N] y2;
  
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
  vector<lower=0>[J] mean1;
  vector<lower=0>[J] sigma2;
  vector<lower=0>[J] mean2;
}


model {

  matrix[2*J , J] ext_tpm1; 
  matrix[2*J , J] ext_tpm2;
  matrix[J , 2*J] log_tr_etpm1;
  matrix[J , 2*J] log_tr_etpm2;
  vector[2*J] lp1_temp;
  vector[2*J] lp1;
  vector[2*J] lp2_temp;
  vector[2*J] lp2;

  
  ext_tpm1[1,1] = weights[1,1]*tpm11[1,1] + weights[2,1]*tpm21[1,1];
  ext_tpm1[2,1] = weights[1,1]*tpm11[1,1] + weights[2,1]*tpm21[2,1];
  ext_tpm1[3,1] = weights[1,1]*tpm11[2,1] + weights[2,1]*tpm21[1,1];
  ext_tpm1[4,1] = weights[1,1]*tpm11[2,1] + weights[2,1]*tpm21[2,1];

  ext_tpm1[1,2] = weights[1,1]*tpm11[2,2] + weights[2,1]*tpm21[1,2];  
  ext_tpm1[2,2] = weights[1,1]*tpm11[2,2] + weights[2,1]*tpm21[2,2];
  ext_tpm1[3,2] = weights[1,1]*tpm11[1,2] + weights[2,1]*tpm21[1,2];
  ext_tpm1[4,2] = weights[1,1]*tpm11[1,2] + weights[2,1]*tpm21[2,2];
  
  ext_tpm2[1,1] = weights[1,2]*tpm12[1,1] + weights[2,2]*tpm22[1,1];
  ext_tpm2[2,1] = weights[1,2]*tpm12[1,1] + weights[2,2]*tpm22[2,1];
  ext_tpm2[3,1] = weights[1,2]*tpm12[2,1] + weights[2,2]*tpm22[1,1];
  ext_tpm2[4,1] = weights[1,2]*tpm12[2,1] + weights[2,2]*tpm22[2,1];

  ext_tpm2[1,2] = weights[1,2]*tpm12[2,2] + weights[2,2]*tpm22[1,2];  
  ext_tpm2[2,2] = weights[1,2]*tpm12[2,2] + weights[2,2]*tpm22[2,2];
  ext_tpm2[3,2] = weights[1,2]*tpm12[1,2] + weights[2,2]*tpm22[1,2];
  ext_tpm2[4,2] = weights[1,2]*tpm12[1,2] + weights[2,2]*tpm22[2,2];

 
  for(j in 1:J){
    for(k in 1:(2*J)){
      log_tr_etpm1[j,k] = log(ext_tpm1[k,j]);
      log_tr_etpm1[j,k] = log(ext_tpm1[k,j]);
    }
  }


  lp = rep_vector(-log(J), J);
  for(j in 1:J){
    lp[j] = lp[j] + normal_lpdf(y1[1]| mean1[j], sigma1[j]);
    lp[j+J] = lp[j+J] + normal_lpdf(y1[1]| mean2[j], sigma2[j]);
  }
    
  for(n in 2:N){
    for(j in 1:J){
      lp1_temp[j] = log_sum_exp(log_tr_etpm1[j])
    }
  }
 

 
}

