# -----------------------------------------------------------------------
# Metropolis Step
# -----------------------------------------------------------------------

# library(mvtnorm)
# library(dplyr)

# N -- number of neighbors
# TT -- length of the time series
# G -- number of states (same across each location)
# state.neighbors -- a list of each of the state sequences of the neighbors for each location

tpmfun <- function(TT, G, beta.current, spatialparam.current, state.draws, position, neighbors, beta.row){
  
  working.tpm <- matrix(diag(G)[beta.row,], nrow=G, ncol=TT-1)
  natural.tpm <- matrix(NA, nrow=G, ncol=TT-1)
  
  state.neighbors <- list()
  for(j in 1:length(neighbors)){
    state.neighbors[[j]] <- state.draws[position[neighbors[j],2], position[neighbors[j],3],]
  }
  
  state.neighbors <- dplyr::bind_cols(state.neighbors)
  
  ## counting the transitions into state g
  weights <- matrix(NA, ncol=TT-1, nrow=G)
  
  for(g in (1:G)[-beta.row]){
    #weights[g,] = sapply(1:(TT-1), function(w) sum(state.neighbors[w,] == g))
    weights[g,] = rowSums(state.neighbors[1:(TT-1),] == g)
  }
  
  working.tpm[-beta.row, ] <- exp(matrix(beta.current, ncol=TT-1, nrow=G-1) + spatialparam.current*weights[-beta.row,])
  
  for(t in 1:(TT-1)){
    natural.tpm[,t] <- working.tpm[,t]/sum(working.tpm[,t])
  }

  return(natural.tpm)
}


beta.clp <- function(TT, state.draws, tpm, beta.row){

  clp <- 0
  
  for(t in 1:(TT-1)){
    if(state.draws[t] == beta.row){
      clp <- clp + log(tpm[state.draws[t+1],t])
    }
  }
  
  #clp <- sum(log(tpm[,state.draws[2:TT]][state.draws[1:(TT-1)] == beta.row]))

  return(clp)
  
}


beta.condpost.draw <- function(N, TT, G, state.draws, position, neighbors, covmat.candidate, beta.current, spatialparam.current, beta.row, prior.mean, prior.covmat){
  
  #---------------------------------------------------------------------------------------------------
  # COMPUTING LOG POSTERIORS FOR CURRENT AND CANDIDATE
  
  lpcurrent <- 0
  lpproposal <- 0
  
  # log posterior for the current value
  for(n in 1:N){
    tpm.current <- tpmfun(TT, G, beta.current, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
    lpcurrent <- lpcurrent + beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.current, beta.row)
  }
  
  
  lpcurrent <- lpcurrent + dmvnorm(beta.current, prior.mean, prior.covmat, log=T)
  #print(lpcurrent)
  
  # generating a new candidate value
  beta.candidate <- rmvnorm(1, mean=beta.current, sigma = covmat.candidate)
  # print(beta.candidate)
  # log posterior for the candidate value
  for(n in 1:N){
    tpm.candidate <- tpmfun(TT, G, beta.candidate, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
    lpproposal <- lpproposal + beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.candidate, beta.row)
  }
  
    
  lpproposal <- lpproposal + dmvnorm(beta.candidate, prior.mean, prior.covmat, log = T)
  #print(lpproposal)
  
  #---------------------------------------------------------------------------------------------------
  # METROPOLIS-HASTINGS STEP
  
  acc <- 0
  
  # Compute acceptance ratio: (symmetric proposal)

  R <- exp(lpproposal - lpcurrent)  
  #print(R)
  U <- runif(1)                          
  if(U<R){       
    beta   <- beta.candidate
    #lpcurrent  <- lpproposal
    acc <- 1
  } else {
    beta <- beta.current
  }
  
  return(list(beta, acc))
    
}


psi.condpost.draw <- function(N, TT, G, state.draws, position, neighbors, cand.sd, beta.current, spatialparam.current, prior.mean, prior.sd){
  
  
  #---------------------------------------------------------------------------------------------------
  # COMPUTING LOG POSTERIORS FOR CURRENT AND CANDIDATE
  
  lpcurrent <- 0
  lpproposal <- 0
  
  # log posterior for the current value
  for(n in 1:N){
    for(g in 1:G){
      tpm.current <- tpmfun(TT, G, beta.current[g,], spatialparam.current, state.draws, position, neighbors[[n]], beta.row=g)
      lpcurrent <- lpcurrent + beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.current, beta.row=g) 
    }
  }
  
  lpcurrent <- lpcurrent + dnorm(spatialparam.current, mean = prior.mean, sd = prior.sd)
  
  # generating a new candidate value
  spatialparam.candidate <- rnorm(n= 1, mean=spatialparam.current, sd = cand.sd)
  # log posterior for the candidate value
  for(n in 1:N){
    for(g in 1:G){
      tpm.candidate <- tpmfun(TT, G, beta.current[g,], spatialparam.candidate, state.draws, position, neighbors[[n]], beta.row=g)
      lpproposal <- lpproposal + beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.candidate, beta.row=g)
    }
  }
  
  lpproposal <- lpproposal + dnorm(spatialparam.candidate, mean=prior.mean, sd=prior.sd)
  
  #---------------------------------------------------------------------------------------------------
  # METROPOLIS-HASTINGS STEP
  
  acc <- 0
  
  # Compute acceptance ratio: (symmetric proposal)
  R <- exp(lpproposal - lpcurrent)  
  #print(R)
  U <- runif(1)                          
  if(U<R){       
    spatialparam   <- spatialparam.candidate
    lpcurrent  <- lpproposal
    acc <- 1
  } else {
    spatialparam <- spatialparam.current
  }
  
  return(list(spatialparam, acc))
  
}



## Extras: 

# lpcurrent <- foreach(n = 1:N, .combine='sum') %do% {
#   tpm.current <- tpmfun(TT, G, beta.current, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
#   beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.current, beta.row)
# }
# 

# lpcurrent <- sum(future_map_dbl(1:N, function(n){    
#   tpm.current <- tpmfun(TT, G, beta.current, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
#   beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.current, beta.row)
# }))

# lpproposal <- foreach(n = 1:N, .combine='sum') %do% {
# tpm.candidate <- tpmfun(TT, G, beta.candidate, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
# beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.candidate, beta.row)
# }

# lpproposal <- sum(future_map_dbl(1:N, function(n){  
#   tpm.candidate <- tpmfun(TT, G, beta.candidate, spatialparam.current, state.draws, position, neighbors[[n]], beta.row)
#   beta.clp(TT, state.draws[position[n,2], position[n,3],], tpm.candidate, beta.row)
# }))

# for(t in 2:(TT-1)){
#   if(state.draws[t-1] == beta.row){
#     clp <- clp + tpm[state.draws[t],t] 
#   }
# }