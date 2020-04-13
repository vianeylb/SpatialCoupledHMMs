## Reading in data from BPA - across various locations, from feb 2010 to april 2018


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
month_year <- month_year[1,] #-- where k indicates the indices desired


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
