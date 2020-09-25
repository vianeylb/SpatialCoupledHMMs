#-------------------------------------------------------------
# Gradient for transformed means in the spatial coupled HMM
#  framework with AR(1) state-dependent distributions
#-------------------------------------------------------------

mean_gradient <- function(y, lagy,  eta, mu, rho){
  
  gradient <- exp(-2*eta)*(y - mu*(1+rho) + rho*lagy)*(1 + rho)
  
  return(sum(gradient))
  
}