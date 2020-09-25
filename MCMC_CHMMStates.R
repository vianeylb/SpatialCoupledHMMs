#install.packages("spdep")
library(spdep)

## constructing rook neighbros for a 5x5 grid cell
neighbors <- cell2nb(5, 5, type = "rook")
position <- data.frame(index = 1:25, x = rep(1:5, each=5), y=rep(1:5, 5))


## initial state for each location

states0 <- array(NA, c(5,5, 200))

for(k in 1:200){
  for(j in 1:5){
    states0[j,,k] <- sample(x=1:3, size=5, prob=rep(1/3,3), replace = T)
  }
}


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

tpm.construction <- function(G, obs, sdd.param, no.neighbors, state.neighbors, tpm.b0s, tpm.spatialeff){

  tpm <- array(NA, dim=c(G, G, length(obs)-1))

  for(t in 1:(length(obs)-1)){
  
    weights <- numeric(G)
  
    ## double sum over the neighbors
    ## and computing the r.v. I(Z_t = k)
  
    for(w in 1:3){
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
  for(n in 1:G){
    logalpha[1,n] = log_delta[n] + log(dnorm(x=obs[1], mean = sdd.param$mu[n], sd=sdd.param$sigma[n])) 
  }
  
  for(t in 2:TT){
    for(n in 1:G){
      logalpha[t,n] = log_sum_exp(log_gamma_tr[n,,t-1], logalpha[t-1,]) + 
          log(dnorm(x=obs[t], mean = sdd.param$mu[n] + 
                      sdd.param$rho[n]*(obs[t-1] - sdd.param$mu[n]), 
                    sd=sdd.param$sigma[n])) 
    }
  }
  
  ## log likelihood
  llk = log_sum_exp(logalpha[TT,], rep(0, G))  
  
  ## state draws from the joint posterior distribution
  stateDraws[TT] = sample(x = 1:G, size=1, prob = exp(logalpha[TT,]) - llk)
  
  for(t0 in (TT-1):1){
    t0prob_unnorm <- exp(logalpha[t0,])*gamma[,stateDraws[t0+1],t0]
    stateDraws[t0] = sample(x=1:G, size=1, prob = t0prob_unnorm/sum(t0prob_unnorm))
  }
  
  return(stateDraws)
  
}


##----------------------------------------------------------------------------------------------
##
## Sampling Steps
##
##----------------------------------------------------------------------------------------------


# @ S: Number of samples
# @ cstate.draws: Array of coupled HMM state draws




S <- 200
G <- 3
mu <- c(0, 5, 10)
rho <- c(-0.5, 0.5, 0.9)
sigma <- c(0.1, 0.5, 1)

tpmb0s <- diag(3)
b0s <- rep(-1, 6)
tpmb0s[!tpmb0s] <- b0s
spatialparam <- -1

sdd.param <- list(mu, rho, sigma)
names(sdd.param) <- c("mu", "rho", "sigma")

# for(s in 1:S){
# 
#   for(k in 1:25){
# 
#     state.neighbors <- list()
#     for(n in 1:length(neighbors[[k]])){
#       state.neighbors[[n]] <- states[position[neighbors[[k]][n],2], position[neighbors[[k]][n],3],]
#     }
# 
#     tpm <- tpm.construction(G=3,
#                             obs = obs[position[k,2],position[k,3],],
#                             sdd.param = sdd.param,
#                             no.neighbors = length(neighbors[[k]]),
#                             state.neighbors = state.neighbors,
#                             tpm.b0s = tpm.b0s, tpm.spatialeff = -1)
# 
# 
#    log_tpm_tr <- array(NA, c(G, G, 199))
#    for(j in 1:199){
#      log_tpm_tr[,,j] <- log(t(tpm[,,j]))
#    }
# 
#    cstate.draws[position[k,2],position[k,3],,s] <- ffbs(G = 3, TT = 200, log_gamma_tr = log_tpm_tr,
#                                                         gamma = tpm, log_delta = log(rep(1/3, 3)),
#                                                         obs = obs[position[k,2],position[k,3],], sdd.param = sdd.param)
# 
#   }
# }

cstate.draws <- array(NA, c(5, 5, 200, S))
cstate.draws[,,,1] <- states0

for(s in 2:S){
  
  index.reord <- sample(x = 1:25, size = 25, prob = rep(1/25, 25))
  cstate.draws[,,,s] <- cstate.draws[,,,s-1]

  for(k in index.reord){
    
    state.neighbors <- list()
    for(n in 1:length(neighbors[[k]])){
      state.neighbors[[n]] <- cstate.draws[position[neighbors[[k]][n],2], position[neighbors[[k]][n],3],,s]
    }
    
    tpm <- tpm.construction(G=3, 
                            obs = obs[position[k,2],position[k,3],], 
                            no.neighbors = length(neighbors[[k]]), 
                            state.neighbors = state.neighbors, 
                            tpm.b0s = tpmb0s, tpm.spatialeff = -1)
    
    
    log_tpm_tr <- array(NA, c(G, G, 199))
    for(j in 1:199){
      log_tpm_tr[,,j] <- log(t(tpm[,,j]))
    } 
    
    cstate.draws[position[k,2],position[k,3],,s] <- ffbs(G = 3, TT = 200, log_gamma_tr = log_tpm_tr, 
                                                         gamma = tpm, log_delta = log(rep(1/3, 3)), 
                                                         obs = obs[position[k,2],position[k,3],], sdd.param = sdd.param)
    
  }
}




