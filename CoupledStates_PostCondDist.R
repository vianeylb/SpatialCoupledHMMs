##----------------------------------------------------------------------------------------------
## FUNCTIONS TO BE CALLED DURING MCMC STEPS
## 
## inv.logit == inverse multinomial logit for tpm construction
## tpm.construction == construct non-homogeneous tpm using neighboring state information
## ffbs == forward-filtering backward-sampling algorithm to get a draw from joint 
##         posterior distribution of the states
##
##----------------------------------------------------------------------------------------------

## inv.multinomial.logit

inv.logit <- function(x, group, G){ 
  
  y <- numeric(G)
  y[-group] <- x
  
  rowprobs <- exp(y)/(sum(exp(y)))
  
  return(rowprobs)
  
}

# @ G:    number of states - integer
# @ obs:  vector of observation of the s_n location
# @ sdd.param: list of length G of state-dependent parameters - mu, sigma, rho
# @ no.neighbors: number of neighborbs 
# @ state.neighbors: list of length no.neighbors with vectors of neighboring state sequences
# @ tpm.b0s: matrix of intercept terms for the transition probability matrix (t.p.m.) 
# @ tpm.spatialeff - one value for the spatialeff

tpm.construction <- function(G, obs, no.neighbors, state.neighbors, tpm.b0s, tpm.spatialeff){

  tpm <- array(NA, dim=c(G, G, length(obs)-1))

  for(t in 1:(length(obs)-1)){
  
    weights <- numeric(G)
  
    ## double sum over the neighbors
    ## and computing the r.v. I(Z_t = k)
  
    for(w in 1:G){
       for(q in 1:no.neighbors){
         weights[w] <- weights[w] + sum(state.neighbors[[q]][t] == w)
       }
     }   
  
    for(g in 1:G){
      tpm[g,,t] <- inv.logit(x = tpm.b0s[g,-g] + weights[-g], group = g, G = G)
    }
   
  }
  
  return(tpm)
}


log_sum_exp <- function(x1, x2){
  
  msf = max(x1 + x2)
  exps = sum(exp(x1+x2 - msf))
  
  lse = msf + log(exps)
  
  return(lse)
  
}


ffbs <- function(G, TT, log_gamma_tr, gamma, log_delta, obs, sdd.param){
  
  stateDraws <- rep(NA, TT)
  logalpha <- matrix(NA, nrow=TT, ncol=G)
  
  ## forward variables
  for(g in 1:G){
    logalpha[1,g] = log_delta[g] + dnorm(x=obs[1], mean = sdd.param$mu[g], sd=sdd.param$sigma[g], log = T)
  }
  
  for(t in 2:TT){
    for(g in 1:G){
      logalpha[t,g] = log_sum_exp(log_gamma_tr[g,,t-1], logalpha[t-1,]) + 
          dnorm(x=obs[t], mean = sdd.param$mu[g] + 
                      sdd.param$rho[g]*(obs[t-1] - sdd.param$mu[g]), 
                    sd=sdd.param$sigma[g], log = T) 
    }
  }
  
  ## log likelihood
  llk = log_sum_exp(logalpha[TT,], rep(0, G))  
  
  ## state draws from the joint posterior distribution
  stateDraws[TT] = sample(x = 1:G, size=1, prob = exp(logalpha[TT,] - llk))
  
  
  ## need to account for underflow
  for(t0 in (TT-1):1){
    #t0prob_unnorm <- exp(logalpha[t0,])*gamma[,stateDraws[t0+1],t0]
    t0prob_unnorm <- exp(logalpha[t0,]  + log(gamma[,stateDraws[t0+1],t0]))
    stateDraws[t0] = sample(x=1:G, size=1, prob = t0prob_unnorm/sum(t0prob_unnorm))
  }
  
  return(stateDraws)
  
}
