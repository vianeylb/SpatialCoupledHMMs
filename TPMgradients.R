# @ N - number of locations/time series
# @ states_primary -- a list of the state sequence of the state sequences at location `n`
# @ states_nghbors -- the characteristics of the neighbors of location `n`
# @ spatialeff -- spatial effect in the tpm
# @ tpmeffects -- \beta_0^(ij) terms - baselines 


beta.cllike <- function(Tlen, N, states_primary, states_nghbors, spatialeff, tpmeffects){

  
  ## initial distributions for each of the time series - fixed at (0.5, 0.5)
  llike <- N*log(0.5) 
  
  for(n in 1:N){
    for(t in 2:Tlen){
      ## summing the numerator
    
      if(states_primary[[n]][t] != states_primary[[n]][t-1] & 
        states_primary[[n]][t] == states_nghbors[[n]][t-1]){
      
      
        ## beta_0^(ij)[1] + \spatialbeta
        ## order of tpmeffects = \beta_0^(12), \beta_0^(21)
        llike <- llike + tpmeffects[[n]][states_primary[[n]][t-1]] + spatialeff
      }
    
      if(states_primary[[n]][t] != states_primary[[n]][t-1] & 
       states_primary[[n]][t] != states_nghbors[[n]][t-1]){
      
      
        ## beta_0^(ij)[1] 
        llike <- llike + tpmeffects[[n]][states_primary[[n]][t-1]] 
      }
    
      ## summing the denominator
  
      if(states_primary[[n]][t] == states_nghbors[[n]][t-1]){
  
        ## log(1 + exp(\beta_0^(ij)[1] + \spatialbeta))
        llike <- llike - 
          log(1+ exp(tpmeffects[[n]][states_primary[[n]][t-1]] + spatialeff))
      }
    
      if(states_primary[[n]][t] != states_nghbors[[n]][t-1]){
      
        ## log(1 + exp(\beta_0^(ij)[1] + \spatialbeta))
        llike <- llike - 
          log(1+ exp(tpmeffects[[n]][states_primary[[n]][t-1]]))
      }
    }
  }
  
  lliks <- llike + dnorm(x=tpmeffects[[1]][1:2], mean = -2, sd = 0.5, log = T) +
    dnorm(x=tpmeffects[[2]][1:2], mean = -2, sd = 0.5, log = T) + 
    dnorm(x=spatialeff, mean=0, sd=0.5, log=T)
  
  return(llike)
    
}



beta.gradients <- function(Tlen, states_primary, states_nghbors, spatialeff, tpmeffects){
  
  b12_1 <- sum((lead(states_primary[[1]]) - states_primary [[1]]) == 1, na.rm=T) - 
    sum(states_primary[[1]][2:Tlen]==1 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[1]][1] + spatialeff)/( 1 + exp(tpmeffects[[1]][1] + spatialeff) ) - 
    (Tlen - 1 - sum(states_primary[[1]][2:Tlen]==1 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0))*
    exp(tpmeffects[[1]][1])/(1+exp(tpmeffects[[1]][1])) + 1/(2*0.5)

  
  b21_1 <- sum((lead(states_primary[[1]]) - states_primary[[1]]) == -1, na.rm=T) - 
    sum(states_primary[[1]][2:Tlen]==2 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[1]][2] + spatialeff)/( 1 + exp(tpmeffects[[1]][2] + spatialeff)) - 
    (Tlen - 1 - sum(states_primary[[1]][2:Tlen]==2 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0))*
    exp(tpmeffects[[1]][2])/(1+exp(tpmeffects[[1]][2])) + 1/(2*0.5)
  
  b12_2 <- sum(lead(states_primary[[2]]) - states_primary[[2]] == 1, na.rm=T) - 
    sum(states_primary[[2]][2:Tlen]==1 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0)*
    exp(tpmeffects[[2]][1] + spatialeff)/( 1 + exp(tpmeffects[[2]][1] + spatialeff) ) - 
    (Tlen - 1 - sum(states_primary[[2]][2:Tlen]==1 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0))*
    exp(tpmeffects[[2]][1])/(1+exp(tpmeffects[[2]][1])) + 1/(2*0.5)
  
  b21_2 <- sum(lead(states_primary[[2]]) - states_primary[[2]] == -1, na.rm=T) - 
    sum(states_primary[[2]][2:Tlen]==2 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[2]][2] + spatialeff)/( 1 + exp(tpmeffects[[2]][2] + spatialeff) ) - 
    (Tlen - 1 - sum(states_primary[[2]][2:Tlen]==2 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0))*
    exp(tpmeffects[[2]][2])/(1+exp(tpmeffects[[2]][2])) + 1/(2*0.5)
  
  spatialeff <- sum((lead(states_primary[[1]]) - states_primary [[1]]) == 1 & lead(states_primary[[1]]) == states_nghbors[[1]], na.rm=T) + 
    sum((lead(states_primary[[1]]) - states_primary [[1]]) == -1 & lead(states_primary[[1]]) == states_nghbors[[1]], na.rm=T) + 
    sum((lead(states_primary[[2]]) - states_primary [[2]]) == 1 & lead(states_primary[[2]]) == states_nghbors[[2]], na.rm=T) + 
    sum((lead(states_primary[[2]]) - states_primary [[2]]) == -1 & lead(states_primary[[2]]) == states_nghbors[[2]], na.rm=T) - 
    sum(states_primary[[1]][2:Tlen]==1 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[1]][1] + spatialeff)/( 1 + exp(tpmeffects[[1]][1] + spatialeff) ) - 
    sum(states_primary[[1]][2:Tlen]==2 & (states_primary[[1]][2:Tlen] - states_nghbors[[1]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[1]][2] + spatialeff)/( 1 + exp(tpmeffects[[1]][2] + spatialeff)) - 
    sum(states_primary[[2]][2:Tlen]==1 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0)*
    exp(tpmeffects[[2]][1] + spatialeff)/( 1 + exp(tpmeffects[[2]][1] + spatialeff) ) - 
    sum(states_primary[[2]][2:Tlen]==2 & (states_primary[[2]][2:Tlen] - states_nghbors[[2]][1:(Tlen-1)])==0, na.rm=T)*
    exp(tpmeffects[[2]][2] + spatialeff)/( 1 + exp(tpmeffects[[2]][2] + spatialeff) ) + 1/(2*0.5)
    

 gradients <- c(b12_1, b21_1, b12_2, b21_2, spatialeff)
 
 return(gradients)
    
}

