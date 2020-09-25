library(tidyverse)
library(lubridate)
library(fda)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(splines)


month_year <- month_year[1:3,]
names <- "Augspurger"

augs <- wind.data%>%filter(location=="Augspurger")
augs$`Date/Time (UTC)` <- mdy_hms(augs$`Date/Time (UTC)`)

augs <- augs%>%mutate(wsdiff = lead(`Wind Speed (MPH)`) - `Wind Speed (MPH)`)
misindex <- augs$`Wind Speed (MPH)` == 0 | is.na(augs$`Wind Speed (MPH)`)
augs$`Wind Speed (MPH)`[misindex] = runif(sum(misindex), min=0, max=20)

## -------------------------------------------------------------------------------------

n <- dim(augs)[1]
dat <- log(augs$`Wind Speed (MPH)`)

# calculation of the B-spline design matrix 
ord<-4
degree<-ord-1
K<-10
nb<-2*K+1
dr <- range(dat)
knots <- seq(dr[1]-0.001, dr[2]+0.001, length.out=25)
#knots<- c(.0699, sort( sample(unique(round(x, 2)), nrknots+2*degree-2)), 3)
B0<-spline.des(knots,seq(knots[1],knots[length(knots)],length=100000),degree+1,outer.ok=T)$design 
ws<-rep(NA,nb)
for (k in 1:nb){
  ws[k]<-(diff(c(knots[1],knots[length(knots)]))/100000*sum(B0[,k]))^(-1) # this computes the integrals of the B-spline basis functions (which are then standardized below)
} 
B<-t(t(spline.des(knots,dat,degree+1,outer.ok=T)$design)*ws) 

## -------------------------------------------------------------------------------------

stan.data<-list(K=nb,Tlen=n,M=2,Bmat=B)

fit.np<-stan(file="HMM_Psplines_SecondOrderDiffPrior.stan",data=stan.data,chains = 1)
