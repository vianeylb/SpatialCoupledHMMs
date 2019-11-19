##------------------------------------------
## Simulating data for a coupled HMM 
##  using two time series
##------------------------------------------

TT <- 1000 #number of observations

init1 <- rep(0.5, 2)
init2 <- rep(0.5, 2)

tpm11 <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2, byrow=T)
tpm12 <- matrix(c(.9, 0.1, 0.1, 0.9), nrow=2, byrow=T)
tpm22 <- matrix(c(0.95, 0.05, 0.01, 0.99), nrow=2, byrow=T)
tpm21 <- matrix(c(.9, 0.1, 0.1, 0.9), nrow=2, byrow=T)
## columns must sum to one for nghbor_weights
nghbor_weights <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow=2, byrow=T)

states1 <- numeric(TT)
states2 <- numeric(TT)

states1[1] <- sample(x = 1:2, size = 1, prob = init1)
states2[1] <- sample(x = 1:2, size = 1, prob = init2) 

for(tt in 2:TT){
  states1[tt] = sample(x=1:2, size=1, 
                       prob=nghbor_weights[1,1]*tpm11[states1[tt-1],] + 
                            nghbor_weights[2,1]*tpm21[states2[tt-1],])
  states2[tt] = sample(x=1:2, size=1, 
                       prob= nghbor_weights[2,2]*tpm22[states2[tt-1],] + 
                             nghbor_weights[1,2]*tpm12[states1[tt-1],])
}


mu <- c(-3, 1)
sdev <- c(1, 1)

obs1 <- rnorm(n=TT, mean=mu[states1], sd=sdev[states1])
obs2 <- rnorm(n=TT, mean=mu[states2], sd=sdev[states2])


##------------------------------------------
## Fitting the model in Stan
##------------------------------------------


library(rstan)

fit.data <- list(TT=TT, 
                 y1=obs1, 
                 y2=obs2, 
                 J = 2)

fit <- stan(file="CoupledHMM_CompLike_2x2.stan", data = fit.data, chains=1,
            init = list(list(mean1=mu1,mean2=mu2)))











