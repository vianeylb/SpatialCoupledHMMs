## Playing around with multiple temporal scales

library(tidyverse)
library(lubridate)

month_year <- month_year[1,]
names <- "Augspurger"
augs <- wind.data%>%filter(location=="Augspurger")
augs$`Date/Time (UTC)` <- mdy_hms(augs$`Date/Time (UTC)`)

augs.avg.hourly <- augs%>%mutate(
  month = month(`Date/Time (UTC)`),
  day = day(`Date/Time (UTC)`),
  hour = hour(`Date/Time (UTC)`)
)%>%group_by(month, day, hour)%>%
  summarize(houravg = mean(`Wind Speed (MPH)`),
            hourvar = var(`Wind Speed (MPH)`), 
            DT = `Date/Time (UTC)`[1])

augs.avg.daily <- augs%>%mutate(
  month = month(`Date/Time (UTC)`),
  day = day(`Date/Time (UTC)`),
  hour = hour(`Date/Time (UTC)`)
)%>%group_by(month, day)%>%
  summarize(dayavg = mean(`Wind Speed (MPH)`),
            dayvar = var(`Wind Speed (MPH)`))


library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

augs.avg.hourly$houravg[is.na(augs.avg.hourly$houravg)] <- 1000
augs.avg.hourly$houravg[augs.avg.hourly$houravg== 0] <- 1000

stan.data <- list(N=5, 
                  T=dim(augs.avg.hourly)[1], 
                  y = log(augs.avg.hourly$houravg))

fit <- stan(file="BasicHMM.stan", data=stan.data, chains=1, iter = 1000) 

N <- 5
mu <- rstan::extract(fit, pars=c("mu"))[[1]]
sigma <- rstan::extract(fit, pars=c("sigma"))[[1]]
tpm <- rstan::extract(fit, pars=c("tpm"))[[1]]

dx <- seq(0, 40, length.out=100)
hist(augs.avg$houravg, breaks=1000, freq = F, xlim=c(0, 50), ylim=c(0, 0.08))


for(j in 1:100){
  delta <- solve(t(diag(N) - tpm[j,,] + 1), rep(1, N))
  
  lines(dx, delta[1]*dnorm(dx, mu[j,1], sigma[j,1]))
  lines(dx, delta[2]*dnorm(dx, mu[j,2], sigma[j,2]))
  lines(dx, delta[3]*dnorm(dx, mu[j,3], sigma[j,3]))
  lines(dx, delta[4]*dnorm(dx, mu[j,4], sigma[j,4]))
#  lines(dx, delta[5]*dnorm(dx, mu[j,5], sigma[j,5]))
  lines(dx,  delta[1]*dnorm(dx, mu[j,1], sigma[j,1]) + 
          delta[2]*dnorm(dx, mu[j,2], sigma[j,2]) + 
          delta[3]*dnorm(dx, mu[j,3], sigma[j,3]) + 
          delta[4]*dnorm(dx, mu[j,4], sigma[j,4]),
#          delta[5]*dnorm(dx, mu[j,5], sigma[j,5]), 
          col="red")
}



stan.data <- list(N=4, 
                  T=dim(augs.avg)[1], 
                  y = augs.avg$houravg)

fit_ar <- stan(file="AR_HMM.stan", data=stan.data, chains=1, iter = 1000) 

N <- 4
mu <- rstan::extract(fit_ar, pars=c("mu"))[[1]]
sigma <- rstan::extract(fit_ar, pars=c("sigma"))[[1]]
rho <- rstan::extract(fit_ar, pars=c("rho"))[[1]]
tpm <- rstan::extract(fit_ar, pars=c("tpm"))[[1]]
delta <- solve(t(diag(N) - tpm[1,,] + 1), rep(1, N))

dx <- seq(0, 40, length.out=100)
hist(augs.avg$houravg, freq = F, 
     xlim=c(0, 50), breaks=800, ylim=c(0, 0.08))

for(j in 90:100){
  lines(dx, delta[1]*dnorm(dx, mu[j,1], sigma[j,1]/sqrt((1-rho[j,1]^2))))
  lines(dx, delta[2]*dnorm(dx, mu[j,2], sigma[j,2]/sqrt((1-rho[j,2]^2))))
  lines(dx, delta[3]*dnorm(dx, mu[j,3], sigma[j,3]/sqrt((1-rho[j,3]^2))))
  lines(dx, delta[4]*dnorm(dx, mu[j,4], sigma[j,4]/sqrt((1-rho[j,4]^2))))
  lines(dx, delta[5]*dnorm(dx, mu[j,5], sigma[j,5]/sqrt((1-rho[j,5]^2))))
  lines(dx,  delta[1]*dnorm(dx, mu[j,1], sigma[j,1]/sqrt((1-rho[j,1]^2))) + 
          delta[2]*dnorm(dx, mu[j,2], sigma[j,2]/sqrt((1-rho[j,2]^2))) + 
          delta[3]*dnorm(dx, mu[j,3], sigma[j,3]/sqrt((1-rho[j,3]^2))) + 
          delta[4]*dnorm(dx, mu[j,4], sigma[j,4]/sqrt((1-rho[j,4]^2))) + 
          delta[5]*dnorm(dx, mu[j,5], sigma[j,5]/sqrt((1-rho[j,5]^2))),
        col="red")
}




