states

# time  = 3
# mcmc = 4
cstate.draws[1,1,1,]

correctpostdraws <- array(NA, c(5, 5, 1000))

for(j in 1:200){
  for(k in 1:5){
    for(l in 1:5){
      correctpostdraws[k,l,j] <- 
        sum(states[k,l,j] == cstate.draws[k,l,j,1001:2000])/1000
    }
    
  }
}