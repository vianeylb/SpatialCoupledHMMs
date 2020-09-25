#------------------------------------------------------------------------
# Gradient for transformed autoregressive component in the
# spatial coupled HMM framework with AR(1) state-dependent distributions
#------------------------------------------------------------------------

arparam_gradient <- function(y, lagy,  eta, mu, rho){
  
  gradient <- e^(-2*eta)*(y - (mu + rho*(lagy - mu)))*(lagy - mu)
  
  return(sum(gradient))
  
}
