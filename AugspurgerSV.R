## Playing around with multiple temporal scales

library(tidyverse)
library(lubridate)

month_year <- month_year[1:3,]
names <- "Augspurger"

augs <- wind.data%>%filter(location=="Augspurger")
augs$`Date/Time (UTC)` <- mdy_hms(augs$`Date/Time (UTC)`)

augs <- augs%>%mutate(wsdiff = lead(`Wind Speed (MPH)`) - `Wind Speed (MPH)`)


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

augs$wsdiff[is.na(augs$wsdiff)] <- 100

stan.data <- list(N=2, 
                  T=dim(augs)[1], 
                  y = augs$wsdiff)

fit2 <- stan(file="BasicHMM_SV.stan", data=stan.data, chains=2, iter = 1000) 


#------------------------------------------------------------------------------------

tpm <- rstan::extract(fit, pars=c("tpm"))[[1]]
delta <- solve(t(diag(2)- tpm[1,,] + 1), rep(1, 2))

hist(augs$wsdiff, breaks=500, xlim=c(-6, 6), freq = FALSE)
xs <- seq(-6, 6, length.out=100)
lines(xs, delta[1]*dnorm(xs, 0, 0.72), col="blue")
lines(xs, delta[2]*dnorm(xs, 0, 1.88), col="red")
lines(xs, delta[1]*dnorm(xs, 0, 0.72) + 
        delta[2]*dnorm(xs, 0, 1.88), col="purple", lwd=2)

#------------------------------------------------------------------------------------

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

augs$`Wind Speed (MPH)`[is.na(augs$`Wind Speed (MPH)`)] <- 1000
augs$`Wind Speed (MPH)`[augs$`Wind Speed (MPH)` == 0] <- 1000

stan.data <- list(N=2, 
                  T=dim(augs)[1], 
                  y = augs$`Wind Speed (MPH)`)

fit <- stan(file="MSDrift_HMM.stan", data=stan.data, chains=2, iter = 1000) 


