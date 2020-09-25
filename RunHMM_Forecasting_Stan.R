## 

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##-------------------------------------------------------------------
## FITTING HMM TO INDIVIDUAL TIME SERIES

TT <- seq(200, 390, by=10)

fit11 <- list()

for(j in 1:length(TT)){

  stan.data <- list(T = TT[j], 
                    y = obs[1,1,1:TT[j]], 
                    N = 3)
  
  fit11[[j]] <- stan("AR_HMM.stan", data=stan.data, chains=3)  

}


fit33 <- list()

for(j in 1:length(TT)){
  
  stan.data <- list(T = TT[j], 
                    y = obs[3,3,1:TT[j]], 
                    N = 3)
  
  fit33[[j]] <- stan("AR_HMM.stan", data=stan.data, chains=3)  
  
}


##-------------------------------------------------------------------
## FITTING HMM TO ALL TIME SERIES (SAME STATE ASSUMPTIONS)

obs.mat <- matrix(NA, ncol=25, nrow=400)

TT <- 1000

stan.data <- list(T = TT[1], 
                  y = obs[3,3,1:TT[1]],# - mean(obs[3,3,1:TT[1]]), 
                  N = 3, 
                  alpha = matrix(c(1, 1, 1, 
                                   1, 1, 1, 
                                   1, 1, 1), 
                                 ncol=3, byrow=T))


fit <- rstan::stan("AR_HMM.stan", data=stan.data, chains = 3)

fit <- rstan::stan("AR_HMM2.stan", data=stan.data, chains = 3)

##------------------------------------------------------------------

TT <- 4000

stan.data <- list(T = TT[1], 
                  y = TurbineDataSet$Turbine1_Speed[1:4000] - mean(TurbineDataSet$Turbine1_Speed), 
                  N = 3, 
                  alpha = matrix(c(10, 1, 1, 
                                   1, 10, 1, 
                                   1, 1, 10), 
                                 ncol=3, byrow=T))

fit <- rstan::stan("AR_HMM.stan", data=stan.data, chains = 3)

##------------------------------------------------------------------

stan.data <- list(T = TT[1], 
                  y = D[[1]]$Speed[1:4000] - mean(D[[1]]$Speed), 
                  N = 3, 
                  alpha = matrix(c(10, 1, 1, 
                                   1, 10, 1, 
                                   1, 1, 10), 
                                 ncol=3, byrow=T))

fit.bpa <- rstan::stan("AR_HMM.stan", data=stan.data, chains = 3)
