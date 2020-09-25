#------------------------------------------------------------------------
# Gradient for transformed autoregressive component in the
# spatial coupled HMM framework with AR(1) state-dependent distributions
#------------------------------------------------------------------------

tpm_intercept_gradient <- function(TransCount, rowInd, colInd,  beta.int, spatial.eff, state.neighbors){
  
  val.num <- exp(beta.int[j] + spatial.eff*state.neighbors[,j])
    
  val.denom <- 1 + exp(beta.int[,1] + spatial.eff*state.neighbors[,1])

  for(i in 2:(G-1)){
    val.denom <- val.denom + exp(beta.int[,i] + spatial.eff*state.neighbors[,i])
  }

  gradient <- Nobs + val.num/val.denom
  
  return(gradient)
  
}

