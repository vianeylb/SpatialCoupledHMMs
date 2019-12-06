// paper: Bayesian regularisation in structured additive regression: a unifying perspective on shrinkage, smoothing and predictor selection

data {
int<lower=0> N1; // length of data set
int<lower=0> N2; // length of data set
int<lower=0> N3; // length of data set
int<lower=0> N4; // length of data set
int<lower=1> K; // number of basis
int<lower=1> M; // number of states
vector[K] B1[N1];
vector[K] B2[N2];
vector[K] B3[N3];
vector[K] B4[N4];
}

parameters {
vector[K] as1[M]; // coefficients
vector[K] as2[M]; // coefficients
vector[K] as3[M]; // coefficients
vector[K] as4[M]; // coefficients
simplex[M] theta[M]; // tpm
vector<lower=0>[M] sigma1;
vector<lower=0>[M] sigma2;
vector<lower=0>[M] sigma3;
vector<lower=0>[M] sigma4;
}  

transformed parameters {
  
  simplex[K] as1mat[M]; // coefficients
  simplex[K] as2mat[M]; // coefficients 
  simplex[K] as3mat[M]; // coefficients 
  simplex[K] as4mat[M]; // coefficients

  for(m in 1:M)
    as1mat[m] = softmax(as1[m]);

  for(m in 1:M)
    as2mat[m] = softmax(as2[m]);

  for(m in 1:M)
    as3mat[m] = softmax(as3[m]);

  for(m in 1:M)
    as4mat[m] = softmax(as4[m]);

}

model {

vector[M] log_theta_tr[M];
vector[M] lp;
vector[M] lp_p1;

//sigma ~ beta(1, 10);
sigma1 ~ normal(0, 0.1);
sigma2 ~ normal(0, 0.1);
sigma3 ~ normal(0, 0.1);
sigma4 ~ normal(0, 0.1);

for(i in 1:3){
  for(j in 1:3){
    as1[j,i] ~ student_t(3, 0, 1);
    as2[j,i] ~ student_t(3, 0, 1);
    as3[j,i] ~ student_t(3, 0, 1);
    as4[j,i] ~ student_t(3, 0, 1);
  }
}

for(m in 1:M) {
  for(k in 4:K)
    as1[m,k] ~ normal(3*as1[m, k-1]-3*as1[m,k-2]+as1[m,k-3], sigma1[m]);
}
for(m in 1:M) {
  for(k in 4:K)
    as2[m,k] ~ normal(3*as2[m, k-1]-3*as2[m,k-2]+as2[m,k-3], sigma2[m]);
}
for(m in 1:M) {
  for(k in 4:K)
    as3[m,k] ~ normal(3*as3[m, k-1]-3*as3[m,k-2]+as3[m,k-3], sigma3[m]);
}
for(m in 1:M) {
  for(k in 4:K)
    as4[m,k] ~ normal(3*as4[m, k-1]-3*as4[m,k-2]+as4[m,k-3], sigma4[m]);
}

// entries of the t.p.m.
for (m_from in 1:M)
for (m in 1:M)
log_theta_tr[m, m_from] = log(theta[m_from, m]);

lp = rep_vector(-log(M), M);

for (n in 1:N1) {
for (m in 1:M)
lp_p1[m] = log_sum_exp(log_theta_tr[m] + lp) + log(sum(B1[n] .* as1mat[m]));
             
lp = lp_p1;
}

target += log_sum_exp(lp) ;

for (n in 1:N2) {
for (m in 1:M)
lp_p1[m] = log_sum_exp(log_theta_tr[m] + lp) + log(sum(B2[n] .* as2mat[m]));
             
lp = lp_p1;
}

target += log_sum_exp(lp) ;

for (n in 1:N3) {
for (m in 1:M)
lp_p1[m] = log_sum_exp(log_theta_tr[m] + lp) + log(sum(B3[n] .* as3mat[m]));
             
lp = lp_p1;
}

target += log_sum_exp(lp) ;

for (n in 1:N4) {
for (m in 1:M)
lp_p1[m] = log_sum_exp(log_theta_tr[m] + lp) + log(sum(B4[n] .* as4mat[m]));
             
lp = lp_p1;
}

target += log_sum_exp(lp) ;

}
