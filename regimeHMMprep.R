## Regime-specific HMM type of model

#--------------------------------------------------------------------------
# R PACKAGES
#--------------------------------------------------------------------------

library(rvest)
library(tidyverse)
library(readr)

#--------------------------------------------------------------------------
# USING RVEST TO EXTRACT HTML CODE AND NAMES FOR CSV FILES
#--------------------------------------------------------------------------

bpahtml <- read_html("https://transmission.bpa.gov/business/operations/wind/MetData/default.aspx")
stationTF <- unlist(html_nodes(bpahtml, "a")%>%html_attrs())%>%str_detect(pattern="MetMonthly.aspx")
all_names <- ((unlist(html_nodes(bpahtml, "a")%>%
                        html_attrs())[stationTF]%>%
                 str_split(pattern="=", simplify = T))[,2]%>%
                str_split(pattern="&", simplify=T))[,1]
## don't want the last four stations as they're no longer in service: 
names <-all_names[!(all_names %in% c("BrowningDepotLvl1", "BrowningDepotLvl2", "CapeBlanco", "Chinook"))]

## start: feb 2010
## end: april 2018
month_year <- data.frame(year = c(rep(2010, 11),
                                  rep(2011:2017, each=12), 
                                  rep(2018, 4)), 
                         month = c(2:12, rep(1:12, times=length(2011:2017)), 
                                   1:4))


month_year <- data.frame(year=rep(2015, each=12), 
                         month=c(1:12))

## For only one location, during a certain period of the year use: 

#names <- c("NaselleRidge", "Megler") #-- where j indicates the indices desired
names <- c("Augspurger", "HoodRiver") 
month_year <- month_year[1:6,] #-- where k indicates the indices desired


#--------------------------------------------------------------------------
# LOOPING OVER STATION NAMES AND MONTH/YEARS OF DATA
#--------------------------------------------------------------------------

datalist <- list() 
for(i in 1:(length(names))){
  for(j in 1:dim(month_year)[1]){
    
    temp_path <- paste("https://transmission.bpa.gov/business/operations/wind/MetData/Monthly/", 
                       names[i], "/", names[i], "_", month_year$year[j], "_", 
                       formatC(month_year$month[j], width=2, flag="0"), ".csv", sep="")
    
    #datalist[[(i-1)*dim(month_year)[1] + j]] <- read_csv(temp_path, skip=6)%>%
    #    mutate(location=names[i])
    tryCatch({ 
      datalist[[(i-1)*dim(month_year)[1] + j]] <- 
        read_csv(temp_path, skip=6, col_types = cols(
          `Date/Time (UTC)` = col_character(),
          `Barometric Pressure (INHG)` = col_double(),
          `Relative Humidity (PT)` = col_double(),
          `Temperature (F)` = col_double(),
          `Wind Direction (DEG)` = col_double(),
          `Wind Speed (MPH)` = col_double(),
          `Peak Time (HMS)` = col_character(),
          `Peak Speed (MPH)` = col_double(),
          `Peak Direction (DEG)` = col_double()
        ))%>%
        mutate(location=names[i]) 
    }, error = function (e) { return(NULL) }) 
    
  }
}

wind.data <- bind_rows(datalist)

#-------------------------------------------------------------------------------------------
## DATA PREP FOR MODEL FITTING

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#--------------------------------------------------------------------------
# AUGSPURGER
#--------------------------------------------------------------------------


augs <- wind.data%>%filter(location == "Augspurger")
# augs1000 <- augs[1:1000,]
# augs$`Wind Speed (MPH)`[which(augs$`Wind Speed (MPH)` == 0)] <- 1000
augs$`Wind Speed (MPH)`[is.na(augs$`Wind Speed (MPH)`)] <- 8

regime <- as.numeric(cut(augs$`Wind Speed (MPH)`,
                             breaks= c(-Inf, quantile(augs$`Wind Speed (MPH)`, probs=c(1/4, 3/4)), Inf)))

# regime <- as.numeric(cut(augs$`Wind Speed (MPH)`, 
#                          breaks= c(-Inf, 5, 15, Inf)))


q25 <- quantile(augs$`Wind Speed (MPH)`, probs=c(1/4))
q75 <- quantile(augs$`Wind Speed (MPH)`, probs=c(3/4))
rawcts <- table(paste(regime, lead(regime), sep=""))
counts <- matrix(rawcts[-4], ncol=3, nrow=3, byrow=T)

# counts <- matrix(table(paste(regime, lead(regime), sep=""))[1:9], 
#                  nrow=3, byrow = T)

# stan.data <- list(TT = dim(augs)[1],
#                       obs = log(augs$`Wind Speed (MPH)`),
#                       regime = regime,
#                       Nregime= 3,
#                       counts = counts,
#                       regimelimits = matrix(c(-0.001, q25+0.1,
#                                               q25-0.1, q75+0.1,
#                                               q75-0.1, 5), nrow=3, ncol=2, byrow=T))

stan.data <- list(TT = dim(augs[1:4000,])[1],
                  obs = augs$`Wind Speed (MPH)`[1:4000],
                  regime = regime[1:4000],
                  Nregime= 3,
                  counts = counts,
                  regimelimits = matrix(c(-0.001, q25,
                                          q25, q75,
                                          q75, 45), nrow=3, ncol=2, byrow=T))

fit <- stan(file="RegimeWindSpeed.stan", data=stan.data, chains = 3)

fitar2 <- stan(file="RegimeWindSpeedAR2.stan", data=stan.data, chains = 3)




