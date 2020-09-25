library(lubridate)


#names <- c("NaselleRidge", "Megler")

wind.data$`Date/Time (UTC)` <- mdy_hms(wind.data$`Date/Time (UTC)`)

# Fitting models to Augspurger

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

augs <- wind.data%>%filter(location=="HoodRiver")
augs$`Wind Speed (MPH)`[is.na(augs$`Wind Speed (MPH)`)] <- 1000
augs$`Wind Speed (MPH)`[augs$`Wind Speed (MPH)` == 0] <- 1000


 
stan.data <- list(N=3, 
                  T=dim(augs)[1], 
                  y = log(augs$`Wind Speed (MPH)`))


stan.data <- list(N=4, 
                  T=3000, 
                  y = log(augs$`Wind Speed (MPH)`)[1:3000])


fit_ws4 <- stan(file="AR_HMM.stan", data=stan.data, chains=3, iter = 1000) 

## AR(2)

stan.data <- list(N=4, 
                  T=3000, 
                  y = log(augs$`Wind Speed (MPH)`)[1:3000])


#fitar2_ws2 <- stan(file="ARP_HMM.stan", data=stan.data, chains=3, iter = 1000) 

fithr4_ws4 <- stan(file = "AR4_HMM.stan", data=stan.data, chains=3, iter=2000)



###---------------------------------------------------------------------------------
# FITTING HOURLY DATA

augs <- wind.data%>%filter(location=="Augspurger")

augs.avg.hourly <- augs%>%mutate(
  month = month(`Date/Time (UTC)`),
  day = day(`Date/Time (UTC)`),
  hour = hour(`Date/Time (UTC)`)
)%>%group_by(month, day, hour)%>%
  summarize(houravg = mean(`Wind Speed (MPH)`, na.rm = T),
            hourvar = var(`Wind Speed (MPH)`, na.rm = T), 
            DT = `Date/Time (UTC)`[1])

#augs.avg.hourly$houravg[is.na(augs.avg.hourly$houravg)] <- 1000
augs.avg.hourly$houravg[augs.avg.hourly$houravg == 0] <- 1000


stan.avg.data <- list(N=2, 
                  T=dim(augs.avg.hourly)[1], 
                  y = log(augs.avg.hourly$houravg))


fitar4avg_ws2 <- stan(file="AR4_HMM.stan", data=stan.avg.data, chains=1, iter = 2000) 

