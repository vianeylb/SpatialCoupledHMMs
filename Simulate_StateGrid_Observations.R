set.seed(17)

########################################################################
##  Simulating 2 weeks of initial 'training' data and 100 days of 'test'
##  Training data will be a rolling window of 2 weeks
##  Assuming 24 obs per day, 1 obs per hour (generally an average)
##
########################################################################

#install.packages("spdep")
library(spdep)

TT <- 14*24 + 100*24

M <- 5

## constructing rook neighbros for a 5x5 grid cell
neighbors <- cell2nb(M, M, type = "rook")

position <- data.frame(index = 1:(M^2), x = rep(1:M, each=M), y=rep(1:M, M))

## simulating from a coupled HMM with state-dependent distributions
## y|y_t-1 = N(mu_g + rho_g(y_{t-1} - mu_g), sigma_g)

mu <- c(0, 5, 10)
#mu <- rep(0, 3)
rho <- c(0.5, -0.5, 0.9)
sigma <- c(0.1, 1, 2)

tpmb0s <- diag(3)
b0s <- rep(-2, 6)
tpmb0s[!tpmb0s] <- b0s
spatialparam <- -1
rawtpm <- tpmb0s

#-----------------------------------------------------------------

inv.logit <- function(x, group, G){ 
  
  y <- numeric(G)
  y[-group] <- x
  
  rowprobs <- exp(y)/(sum(exp(y)))
  
  return(rowprobs)
  
}

#-------------------------------------------------------------------------------------------

# number of states
G <- 3

states <- array(NA, dim=c(M, M, TT))
init <- rep(1/3, 3)

tpm.array <- array(NA, c(3, 3, TT-1, M^2))


for(j in 1:TT){
  
  weights <- rep(0, 3)
  
  for(k in 1:(M^2)){
    
    if(j == 1){
      
      states[position[k,2], position[k,3],j] <- sample(x = 1:3, size=1, prob=rep(1/3, 3))
    
    } else {
      
      weights <- rep(0, 3)
      
      for(w in 1:3){
        for(q in 1:length(neighbors[[k]])){
          weights[w] <- weights[w] + sum(states[position[neighbors[[k]][q],2],
                                                position[neighbors[[k]][q],3],j-1] == w)
        }
      }     
      
      tpmbs <- rawtpm[states[position[k,2], position[k,3], j-1], ] + spatialparam*weights
      
      for(g in 1:G){    
        tpm.array[g,,j-1, k] <- inv.logit(x = rawtpm[g,-g] + spatialparam*weights[-g], group = g, G=3)
      }
      
     states[position[k,2], position[k,3], j] <- 
        sample(x=1:3, size=1, 
               prob = inv.logit(x = tpmbs[-states[position[k,2], position[k,3], j-1]],
                                group = states[position[k,2], position[k,3], j-1], 
                                G=3))
      
    }
  }
}

#-------------------------------------------------------------------------------------------

## observations

obs <- array(NA, dim=c(M, M, TT))

for(j in 1:TT){
  for(k in 1:(M^2)){
    if(j == 1){
      obs[position[k,2], position[k,3],j] <- rnorm(1, 
                                                   mean=mu[states[position[k,2], position[k,3], j]], 
                                                   sd=sigma[states[position[k,2], position[k,3], j]])
    } else {
      
      tempmean <- mu[states[position[k,2], position[k,3], j]] + 
        rho[states[position[k,2], position[k,3], j]] * 
        (obs[position[k,2], position[k,3], j-1] - mu[states[position[k,2], position[k,3], j]])
      
      obs[position[k,2], position[k,3],j] <- rnorm(1, 
                                                   mean=tempmean, 
                                                   sd=sigma[states[position[k,2], position[k,3], j]])
    }
  }                                                                  
}

# #-------------------------------------------------------------------------------------------
# #  Visulizing simulated states
# #
# 
# position.states <- data.frame(time = 1, x=position$x, y=position$y, state = c(t(states[,,1])))
# 
# for(j in 2:200){
#   position.states <- rbind(position.states,  data.frame(time = j, x=position$x, y=position$y, state = c(t(states[,,j]))))
# }
# 
# library(ggplot2)
# library(tidyverse)
# library(gganimate)
# 
# 
# state.animation <- ggplot(data=position.states%>%filter(time %in% c(1:20)), aes(x, y)) + 
#   geom_point(aes(color=state), size=15) + 
#   scale_colour_viridis_c() + theme_classic() + 
#   labs(title = 'Time: {frame_time}') + 
#   transition_time(time)
# 
# position.obs <- data.frame(time = 1, x=position$x, y=position$y, obs = c(t(obs[,,1])))
# 
# for(j in 2:200){
#   position.obs <- rbind(position.obs, data.frame(time = j, x=position$x, y=position$y, obs= c(t(obs[,,j]))))
# }
# 
# ggplot(data=position.obs%>%filter(time %in% c(1:100)), aes(x, y)) + geom_point(aes(color=obs), size=3) + 
#   scale_color_viridis_c() + facet_wrap(~time)
# 
# cor(c(obs[3,3,]),c(obs[4,3,])) 
# 
# ggplot(data=position.obs%>%filter(time %in% c(1:100)), aes(x, y)) + geom_point(aes(color=obs), size=3) + 
#   scale_color_viridis_c() + facet_wrap(~time)
# 
# animation <- ggplot(position.obs, aes(x, y)) + 
#   geom_point(aes(color=obs), size=15) + 
#   scale_color_viridis_c() + theme_classic() + 
#   labs(title = 'Time: {frame_time}') + 
#   transition_time(time)
# 
# animate(animation, nframes = 200, fps = 10)
# 
# position.states <- data.frame(time = 1, x=position$x, y=position$y, states = c(t(states[,,1])))
# 
# for(j in 2:200){
#   position.states <- rbind(position.states, data.frame(time = j, x=position$x, y=position$y, states= c(t(states[,,j]))))
# }
# 
# animation <- ggplot(position.states, aes(x, y)) + 
#   geom_point(aes(color=factor(states)), size=15) + 
#   scale_color_viridis_d() + theme_classic() + 
#   labs(title = 'Time: {frame_time}') + 
#   transition_time(time)
# 
# animate(animation, nframes = 200, fps = 5)
#         
