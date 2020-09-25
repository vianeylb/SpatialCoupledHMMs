## Running the MCMC 20 times for 1000 iterations each

## -----------------------------------------------------------------------------
## Simulate data + read in functions for MCMC

library(spdep)
library(mvtnorm)

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

## simulating data
source("Simulate_StateGrid_Observations.R")

## sourcing cond. post. distributions for the states
source("CoupledStates_PostCondDist.R")

## sourcing cond. post. distributions for mu, sigma, rho
source("SDD_PostCondDist.R")

## sourcing cond. post. distributions for beta, psi
source("TPM_PostCondDist.R")

## -----------------------------------------------------------------------------
## Setup data for MCMC (start at 180 -- end at 199 obs across the grid)


# @ S: Number of samples
# @ cstate.draws: Array of coupled HMM state draws

S <- 1000 ## number of posterior draws
G <- 3 ## number of states
TT <- 180:199 ## length of the data set
N <- 25

mu <- c(0, 50, 100)
rho <- c(-0.5, 0.5, 0.9)
sigma <- c(0.1, 0.5, 1)

tpmb0s <- diag(3)
b0s <- rep(-3, 6)
tpmb0s[!tpmb0s] <- b0s
spatialparam <- 1

sdd.param <- list(mu, rho, sigma)
names(sdd.param) <- c("mu", "rho", "sigma")


## Objects to save the posterior draws
mu.draws <- list()
sigma.draws <- list()
rho.draws <- list()
beta.draws <- list()
acc.beta <- list()
psi.draws <- list()
acc.psi <- list()
cstate.draws <- list()


for(e in 1:length(TT)){
  
  mu.draws[[e]] <- matrix(NA, ncol=G, nrow=S)
  sigma.draws[[e]] <- matrix(NA, ncol=G, nrow=S)
  rho.draws[[e]] <- matrix(NA, ncol=G, nrow=S)
  beta.draws[[e]] <- array(NA, dim=c(S, G, G-1))
  acc.beta[[e]] <- matrix(NA, ncol=G, nrow=S)
  psi.draws[[e]] <- numeric(S)
  acc.psi[[e]] <- numeric(S)

  ## Initialize MCMC
  cstate.draws[[e]] <- array(NA, c(5, 5, TT[e], S))
  cstate.draws[[e]][,,,1] <- states[,,1:TT[e]]
  mu.draws[[e]][1,] <- mu
  sigma.draws[[e]][1,] <- sigma
  rho.draws[[e]][1,] <- rho
  beta.draws[[e]][1,,] <- b0s
  psi.draws[[e]][1] <- spatialparam

}


## -----------------------------------------------------------------------------
## MCMC Loop (start at 180 -- end at 199 obs across the grid)

for(e in 1:length(TT)){
  for(s in 2:S){
    
    print(c(e,s))
    
    sdd.param <- list(mu.draws[[e]][s-1,], rho.draws[[e]][s-1,], sigma.draws[[e]][s-1,])
    names(sdd.param) <- c("mu", "rho", "sigma")
    
    ## Sampling from states
    
    index.reord <- sample(x = 1:N, size = N, prob = rep(1/N, N), replace = F)
    cstate.draws[[e]][,,,s] <- cstate.draws[[e]][,,,s-1]
    
    for(k in index.reord){
      
      state.neighbors <- list()
      for(n in 1:length(neighbors[[k]])){
        state.neighbors[[n]] <- cstate.draws[[e]][position[neighbors[[k]][n],2], position[neighbors[[k]][n],3],,s]
      }
      
      tpmb0s <- diag(G)
      for(g in 1:G){
        tpmb0s[g,-g] <- beta.draws[[e]][s-1,g,] 
      }
      tpm <- tpm.construction(G=3, 
                              obs = obs[position[k,2],position[k,3],1:TT[e]], 
                              no.neighbors = length(neighbors[[k]]), 
                              state.neighbors = state.neighbors, 
                              tpm.b0s = tpmb0s, tpm.spatialeff = psi.draws[[e]][s-1])
      
      ## need to fix beta.draws
      
      log_tpm_tr <- array(NA, c(G, G, TT[e]-1))
      for(j in 1:(TT[e]-1)){
        log_tpm_tr[,,j] <- log(t(tpm[,,j]))
      } 
      
      cstate.draws[[e]][position[k,2],position[k,3],,s] <- ffbs(G = 3, TT = TT[e], log_gamma_tr = log_tpm_tr, 
                                                           gamma = tpm, log_delta = log(rep(1/3, 3)), 
                                                           obs = obs[position[k,2],position[k,3],1:TT[e]], 
                                                           sdd.param = sdd.param)
      
    }
    
    ## Splitting the observations into G groups
    obs.mat <- list()
    for(g in 1:G){
      
      obs.mat[[g]] <- cbind(obs[1, 1, which(cstate.draws[[e]][1,1,2:TT[e],s]==g)+1], 
                            obs[1, 1, which(cstate.draws[[e]][1,1,2:TT[e],s]==g)])
      
      for(n in 2:N){
        obs.mat[[g]] <- rbind(obs.mat[[g]], 
                              cbind(obs[position[n,2], position[n,3], 
                                        which(cstate.draws[[e]][position[n,2], position[n,3],2:TT[e],s]==g)+1], 
                                    obs[position[n,2], position[n,3], 
                                        which(cstate.draws[[e]][position[n,2], position[n,3],2:TT[e],s]==g)]))
      }
    }
    
    ## Setting up the current values of the parameters before sampling
    mu.current <- mu.draws[[e]][s-1,]
    sigma.current <- sigma.draws[[e]][s-1,]
    rho.current <- rho.draws[[e]][s-1,]
    
    
    ## Sampling from state-dependent distribution
    
    for(g in 1:G){
      mu.draws[[e]][s,g] <- mu.condpost.draw(obs.mat[[g]], mu.current[g], 
                                             prior.mean=c(0, 50, 100)[g], prior.sd=1, 
                                             sigma.current[g], rho.current[g])   
    }
    
    for(g in 1:G){
      sigma.draws[[e]][s,g] <- sqrt(sigmasq.condpost.draw(obs.mat[[g]], mu.current[g], 
                                                          rho.current[g], prior.df = 1, prior.scale = 1))
      rho.draws[[e]][s,g] <- rho.condpost.draw(obs.mat[[g]], mu.draws[[e]][s,g], 
                                               sigma.draws[[e]][s,g], 
                                               prior.mean.rho = 0, prior.sd.rho = 0.3)
    }
    
    
    ## Sampling from tpm
    covmat.candidate = matrix(c(0.01, 0, 0, 0.01), nrow=G-1, byrow=T)
    
    for(g in 1:G){
      beta.draws.acc <- beta.condpost.draw(N = N, TT = TT[e], G=G, state.draws = cstate.draws[[e]][,,,s], 
                                           position = position, neighbors = neighbors, 
                                           covmat.candidate = covmat.candidate, 
                                           beta.current = beta.draws[[e]][s-1, g,], 
                                           spatialparam.current = psi.draws[[e]][s-1], 
                                           beta.row = g, prior.mean = c(-3, -3), 
                                           prior.covmat = matrix(c(0.1, 0, 0, 0.1), byrow=T, ncol=2))
      beta.draws[[e]][s,g,] <- beta.draws.acc[[1]]
      acc.beta[[e]][s,g] <- beta.draws.acc[[2]] 
    }  
    
    psi.draws.acc <- psi.condpost.draw(N = N, TT = TT[e], G=G, state.draws = cstate.draws[[e]][,,,s],
                                       position = position, neighbors = neighbors,
                                       cand.sd = 0.1, beta.current = beta.draws[[e]][s,,],
                                       spatialparam.current = psi.draws[[e]][s-1], prior.mean = 1, prior.sd = 0.2)
    psi.draws[[e]][s] <- psi.draws.acc[[1]]
    acc.psi[[e]][s] <- psi.draws.acc[[2]]
    
  }
}
