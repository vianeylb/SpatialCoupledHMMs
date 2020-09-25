#-------------------------------------------------------------
# Gradient for transformed std dev in the spatial coupled HMM
#  framework with AR(1) state-dependent distributions
#-------------------------------------------------------------

stddev_gradient <- function(y, lagy, eta, mu, rho){
  
  gradient <- -eta + exp(-2*eta)*(y - (mu + rho*(lagy - mu)))^2
  
  return(sum(gradient))
  
}