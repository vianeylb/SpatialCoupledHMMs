library(geoR) #for the inverse chi squred distribution

# -----------------------------------------------------------------------
# Gibbs Step
# -----------------------------------------------------------------------

### ASSUMING mu ~ N(prior.mean, prior.sd) (explicitly setting k0 = 1 -- BDA notation)

mu.condpost.draw <- function(obs, mu.current, prior.mean, prior.sd, sigma.current) {

  n <- length(obs) 
  
  update.mean <- (prior.mean + sum(y))/(1 + n)
  update.sd <- sigma.current/(1 + n)
  
  mu.new <- rnorm(1, mean = update.mean, sd = update.sd)  
    
  return(mu.new)
}



sigma.condpost.draw <- function(obs, sigma.current, prior.mean, prior.sd, prior.df, prior.scale) {

  n <- length(obs)
  
  update.df <- prior.df + n
  update.scale <- prior.scale*prior.sd^2 + (n-1)var(obs) + (n/(1+n))*(mean(obs) - prior.mean)^2
  
  sigma.new <- rinvchisq(1, df=update.df, scale=update.scale)
  
  return(sigma.new)  
  
}


