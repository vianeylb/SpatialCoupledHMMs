log_sum_exp <- function(x1, x2){
  
  msf = max(x1 + x2)
  exps = sum(exp(x1+x2 - msf))
  
  lse = msf + log(exps)
  
  return(lse)
  
}

ffbs_chmm <- function(N, TT, log_gamma_tr, gamma, log_delta, steps, shape, rate){
  
  stateDraws <- rep(NA, TT)
  logalpha <- matrix(NA, nrow=TT, ncol=N)
  
  ## forward variables
  for(n in 1:N){
    logalpha[1,n] = log_delta[n] 
    if(steps[1]>0)
      logalpha[1,n] = logalpha[1,n] + log(dgamma(x=steps[1], shape = shape[n], rate=rate[n])) 
  }
  
  for(t in 2:TT){
    for(n in 1:N){
      logalpha[t,n] = log_sum_exp(log_gamma_tr[n,], logalpha[t-1,]) 
      if(steps[t]>0)
        logalpha[t,n] = logalpha[t,n] + log(dgamma(x=steps[t], shape = shape[n], rate=rate[n])) 
    }
  }
  
  ## log likelihood
  llk = log_sum_exp(logalpha[TT,], rep(0, N))  
  
  ## state draws from the joint posterior distribution
  stateDraws[TT] = sample(x = 1:N, size=1, prob = exp(logalpha[TT,]) - llk)
  
  for(t0 in (TT-1):1){
    
    t0prob_unnorm <- exp(logalpha[t0,])*gamma[,stateDraws[t0+1]]
    stateDraws[t0] = sample(x=1:N, size=1, prob = t0prob_unnorm/sum(t0prob_unnorm))
    
  }
  
  
  return(stateDraws)
  
}