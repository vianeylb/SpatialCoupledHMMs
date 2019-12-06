// paper: Bayesian regularisation in structured additive regression: a unifying perspective on shrinkage, smoothing and predictor selection

data {

int<lower=0> Tlen; // length of data set

int<lower=1> K; // number of basis
int<lower=1> M; // number of states
vector[K] Bmat[Tlen];

}

parameters {
vector[K] as[M]; // coefficients

simplex[M] theta[M]; // tpm
vector<lower=0>[M] sigma;

}  

transformed parameters {
  
  simplex[K] asmat[M]; // coefficients


  for(m in 1:M)
    asmat[m] = softmax(as[m]);

}

model {

vector[M] log_theta_tr[M];
vector[M] lp;
vector[M] lp_p;

//sigma ~ beta(1, 10);
sigma ~ normal(0, 0.1);

for(i in 1:2){
  for(j in 1:2){
    as[j,i] ~ student_t(3, 0, 1);
  }
}

for(m in 1:M) {
  for(k in 3:K)
    as[m,k] ~ normal(2*as[m, k-1]-as[m,k-2], sigma[m]);
}


// entries of the t.p.m.
for (m_from in 1:M)
for (m in 1:M)
log_theta_tr[m, m_from] = log(theta[m_from, m]);

lp = rep_vector(-log(M), M);

for (n in 1:Tlen) {
for (m in 1:M)
lp_p[m] = log_sum_exp(log_theta_tr[m] + lp) + log(sum(Bmat[n] .* asmat[m]));
             
lp = lp_p;
}

target += log_sum_exp(lp) ;

}
