# read in Augspurger data - Jan 2015

library(lubridate)

smonth_year <- month_year[1,]
names <- "Augspurger"
augs$`Date/Time (UTC)` <- mdy_hms(augs$`Date/Time (UTC)`)

# Fitting models to Augspurger

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

augs <- wind.data%>%filter(location=="Augspurger")
augs$`Wind Speed (MPH)`[is.na(augs$`Wind Speed (MPH)`)] <- 1000
augs$`Wind Speed (MPH)`[augs$`Wind Speed (MPH)` == 0] <- 1000

stan.data <- list(N=2, 
                  T=dim(augs)[1], 
                  y = log(augs$`Wind Speed (MPH)`))

fit <- stan(file="BasicHMM.stan", data=stan.data, chains=2, iter = 1000) 


N <- 4
mu <- rstan::extract(fit, pars=c("mu"))[[1]]
sigma <- rstan::extract(fit, pars=c("sigma"))[[1]]
tpm <- rstan::extract(fit, pars=c("tpm"))[[1]]
delta <- sapply(1:100, function(j) solve(t(diag(N) - tpm[j,,] + 1), rep(1, N)))

dx <- log(seq(0, 40, length.out=100))
hist(log(augs$`Wind Speed (MPH)`), breaks=100, freq = F, xlim=c(-2, 4), ylim=c(0, 1))

for(j in 1:100){
lines(dx, delta[1,j]*dnorm(dx, mu[j,1], sigma[j,1]))
lines(dx, delta[2,j]*dnorm(dx, mu[j,2], sigma[j,2]))
lines(dx, delta[3,j]*dnorm(dx, mu[j,3], sigma[j,3]))
lines(dx, delta[4,j]*dnorm(dx, mu[j,4], sigma[j,4]))
#lines(dx, delta[5]*dnorm(dx, mu[j,5], sigma[j,5]))
lines(dx,  delta[1,j]*dnorm(dx, mu[j,1], sigma[j,1]) + 
        delta[2,j]*dnorm(dx, mu[j,2], sigma[j,2]) + 
        delta[3,j]*dnorm(dx, mu[j,3], sigma[j,3]) + 
        delta[4,j]*dnorm(dx, mu[j,4], sigma[j,4]),#+ 
        #delta[5]*dnorm(dx, mu[j,5], sigma[j,5]), 
        col="red")
}


#--------------------------------------------------------------------------------


stan.data <- list(N=4, 
                  T=dim(augs)[1], 
                  y = log(augs$`Wind Speed (MPH)`))

fit_ar4 <- stan(file="AR_HMM.stan", data=stan.data, chains=2, iter = 1000) 
#fit_drift3 <- stan(file="AR_HMM_Ys.stan", data=stan.data, chains=2, iter = 1000) 


N <- 5
mu <- rstan::extract(fit_ar, pars=c("mu"))[[1]]
sigma <- rstan::extract(fit_ar, pars=c("sigma"))[[1]]
rho <- rstan::extract(fit_ar, pars=c("rho"))[[1]]
tpm <- rstan::extract(fit_ar, pars=c("tpm"))[[1]]
delta <- solve(t(diag(N) - tpm[1,,] + 1), rep(1, N))

dx <- seq(0, 40, length.out=100)
hist(wind.data$`Wind Speed (MPH)`, breaks=50, freq = F)

for(j in 1:100){
  lines(dx, delta[1]*dnorm(dx, mu[j,1], sigma[j,1]/sqrt((1-rho[j]))))
  lines(dx, delta[2]*dnorm(dx, mu[j,2], sigma[j,2]/sqrt((1-rho[j]))))
  lines(dx, delta[3]*dnorm(dx, mu[j,3], sigma[j,3]/(1-rho[j])))
  lines(dx, delta[4]*dnorm(dx, mu[j,4], sigma[j,4]/(1-rho[j])))
  lines(dx,  delta[1]*dnorm(dx, mu[j,1], sigma[j,1]/(1-rho[j])) + 
          delta[2]*dnorm(dx, mu[j,2], sigma[j,2]/(1-rho[j])) + 
          delta[3]*dnorm(dx, mu[j,3], sigma[j,3]/(1-rho[j])) + 
          delta[4]*dnorm(dx, mu[j,4], sigma[j,4]/(1-rho[j])), 
        col="red")
}
