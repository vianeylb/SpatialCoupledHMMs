# Fitting models to Augspurger

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

augs <- wind.data%>%filter(location=="Augspurger")
augs$`Wind Speed (MPH)`[is.na(augs$`Wind Speed (MPH)`)] <- 1000

stan.data <- list(N=4, 
                  T=dim(augs)[1], 
                  y = augs$`Wind Speed (MPH)`)

fit <- stan(file="BasicHMM.stan", data=stan.data, chains=1, iter = 1000) 
