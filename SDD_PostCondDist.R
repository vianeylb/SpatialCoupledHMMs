library(geoR) #for the inverse chi squred distribution

# -----------------------------------------------------------------------
# Gibbs Step
# -----------------------------------------------------------------------

### ASSUMING mu ~ N(prior.mean, prior.sd) (explicitly setting k0 = 1 -- BDA notation)

mu.condpost.draw <- function(obs, mu.current, prior.mean, prior.sd, sigma.current, rho.current) {

  #obs[,1] = y_t
  #obs[,2] = y_{t-1}  
  
  n <- dim(obs)[1]
  
  update.sd <- 1/sqrt(n*(1-rho.current)^2/sigma.current^2 + 1/prior.sd^2)
  update.mean <- ((1-rho.current)/(sigma.current^2) * (sum(obs[,1]) - rho.current*sum(obs[,2])) +  (1/prior.sd^2) *prior.mean)*update.sd^2
  
  
  mu.new <- rnorm(1, mean = update.mean, sd = update.sd)  
    
  return(mu.new)
}

sigmasq.condpost.draw <- function(obs, mu.current, rho.current, prior.df, prior.scale) {

  #obs[,1] = y_t
  #obs[,2] = y_{t-1}  
  
  n <- dim(obs)[1]
  
  update.df <- (n + prior.df)
  update.scale1 <- prior.df*prior.scale 
  update.scale2 <- sum((obs[,1] - (mu.current + rho.current*(obs[,2]-mu.current)))^2) 
    
  update.scale <- (update.scale1 + update.scale2)/(update.df)
    
  sigma.new <- rinvchisq(1, df=update.df, scale=update.scale)
  
  return(sigma.new)  
  
}

rho.condpost.draw <- function(obs, mu.current, sigma.current, prior.mean.rho, prior.sd.rho){
 
  #obs[,1] = y_t
  #obs[,2] = y_{t-1}
  
  update.sd <- 1/sqrt(sum((obs[,2] - mu.current)^2)/sigma.current^2 + 1/prior.sd.rho^2)
  update.mean <- (1/sigma.current^2 * sum((obs[,2] - mu.current)*(obs[,1] - mu.current)) + prior.mean.rho/prior.sd.rho^2)*update.sd^2
  
  rho.new <- rnorm(1, mean=update.mean, sd= update.sd)
   
}

