# -----------------------------------------------------------------------
# Metropolis Step within Gibbs
# -----------------------------------------------------------------------

Bayes.mlogistic <- function(tpm.row.obs, 
                            X, 
                            curr.beta, 
                            can.sd=0.1,
                            prior.mean.int, 
                            prior.sd.int,
                            prior.mean.slope,
                            prior.sd.slope,
                            row, 
                            G){
  
  #-------------------------------------------------------------------
  # Set-up
  #-------------------------------------------------------------------
  #Initial values:
  ## curr.beta -- G-1x2 matrix
  
  ## acceptance 
  acc   <- numeric(G)
  #-------------------------------------------------------------------

  
  #-------------------------------------------------------------------
  # MH Step
  #------------------------------------------------------------------- 
  ## conditional posterior
  curlp <- cond_log_post(tpm.row.obs,
                         X,
                         betarow =beta, 
                         row, 
                         prior.mean.int, prior.sd.int, 
                         prior.mean.slope, prior.sd.slope,
                         G) 
  
  
  beta <- currbeta
  
  #Update betas using MH sampling:
  for(j in 1:G){
    
    # Draw candidate:
    if(j != G){
      beta[j,1] <- rnorm(n = G-1, mean = currbeta[j,1], sd = can.sd)
    }
    if(j == G){
      beta[,2] <- rnorm(n = 1, mean = currrbeta[1,2], sd=can.sd)
    }
    
    canbeta    <- beta
    canlp      <- cond_log_post(tpm.row.obs,
                                X,
                                betarow = canbeta, 
                                row, 
                                prior.mean.int, prior.sd.int, 
                                prior.mean.slope, prior.sd.slope,
                                G)
    
    # Compute acceptance ratio:
    R <- exp(canlp-curlp)  
    U <- runif(1)                          
    if(U<R){       
      beta   <- canbeta
      curlp  <- canlp
      acc[j] <- acc[j]+1
    }
  }

  
  # Return the posterior samples of beta and
  # the Metropolis acceptance rates
  list(beta=beta,acc=acc)
  
}


