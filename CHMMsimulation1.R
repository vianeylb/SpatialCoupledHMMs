## Simulation to produce two time series

# @ G - number of states
# @ len - length of each time series
# @ initprob - initial distribution 

G <- 2
len <- 200
initprob <- c(0.5, 0.5)

# @ effect - \beta_0^{ij} terms for each process
# @ spatialeff - effect of spatial dependence

effect <- list()
effect[[1]] <- rep(-2, 2)
effect[[2]] <- rep(-2.5, 2)
spatialeff <- c(-0.5, -0.5)


tpm <- function(effect, spatialeff, nghborstate){
  
  indicator1 <- as.numeric(nghborstate == 1)
  indicator2 <- as.numeric(nghborstate == 2)
  
  tpm12 <- exp(effect[1] + spatialeff[1]*indicator1)/(1+exp(effect[1] + spatialeff[1]*indicator1))
  tpm21 <- exp(effect[2] + spatialeff[2]*indicator2)/(1+exp(effect[2] + spatialeff[2]*indicator2))
  
  tpm.mat <- matrix(c(1-tpm12, tpm12, tpm21, 1-tpm21), nrow=2, byrow=T)
  
}


states <- list()

for(g in 1:G){
  states[[g]] <- numeric(len)
}

## need to initialize one of the state sequences in order to start the simulation process
states[[2]] <- sample(x=1:2, size=len, replace = T, prob=c(0.5, 0.5))


for(g in 1:G){
  
  if(g == 1){
    nghbor <- 2
  } 
  
  if(g == 2){
    nghbor <- 1
  }
  
  states[[g]][1] <- sample(x=1:2, size=1, prob=initprob)
  for(t in 2:len){
    
    temptpm <- tpm(effect = effect[[g]], spatialeff = spatialeff, nghborstate = states[[nghbor]][t-1])
    
    states[[g]][t] <- sample(x=1:2, size=1, prob=temptpm[states[[g]][t-1],])
  }
}


# state-dependent distributions
mu <- list()
sigma <- list()

mu[[1]] <- c(-10, -5)
mu[[2]] <- c(-9, -4)

sigma[[1]] <- c(1, 1)
sigma[[2]] <- c(1, 1)


obs <- list()
obs[[1]] <- rnorm(len, mean=mu[[1]][states[[1]]], sd=1)
obs[[2]] <- rnorm(len, mean=mu[[2]][states[[2]]], sd=1)



