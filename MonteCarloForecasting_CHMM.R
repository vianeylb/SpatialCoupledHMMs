## statesTTpred -- a 5x5 array of predicted states at time TT
## obsTT -- a 5x5 array of observed values at time TT


forecasting.chmm <- function(timesteps, obsTT, statesTTpred, G, sdd.param, beta.draw, psi.draw, position, neighbors){
  
  K <- dim(position)[1]
  
  obs.forecast <- array(NA, c(5, 5,timesteps))
  states.forecast <- array(NA, c(5, 5, timesteps))

  weights <- matrix(0, nrow=K, ncol=G)
  
  for(k in 1:K){
    
    state.neighbors <- list()
    for(n in 1:length(neighbors[[k]])){
      state.neighbors[[n]] <- statesTTpred[position[neighbors[[k]][n],2], position[neighbors[[k]][n],3]]
    }
    
    for(w in 1:G){
      for(q in 1:length(state.neighbors)){
        weights[k,w] <- weights[k,w] + sum(state.neighbors[[q]][t] == w)
      }
    }
    
    tpm.row <- diag(G)[statesTTpred[position[k,2], position[k,3]]]
    tpm.row[-statesTTpred[position[k,2], position[k,3]]] <-   
      exp(beta.draw[statesTTpred[position[k,2], position[k,3]]] + weights[k,])/
      (1 + exp(beta.draw[statesTTpred[position[k,2], position[k,3]]] + weights[k,]))
    
  }
  
  
  
  
}
