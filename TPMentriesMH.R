# -----------------------------------------------------------------------
# Metropolis Step within Gibbs -- sob
# -----------------------------------------------------------------------

Bayes.mlogistic <- function(tpm.row, 
                            X, 
                            n.samples=10000, 
                            can.sd=0.1, 
                            ncovs, 
                            prior.mean, 
                            prior.sd,
                            row){
  
  #Initial values:
  beta <- matrix(rnorm((ncovs+1)*(G-1)), nrow=G-1, ncol=(ncovs + 1))
  
  # Keep track of the samples   
  keep.beta     <- array(0,dim=c(G-1, (ncovs+1), n.samples))
  keep.beta[,,1] <- beta
  
  acc   <- att <- rep(0,2)
  
  curlp <- cond_log_post(tpm.row,X,betarow, row, prior.mean, prior.sd) # log posterior at current beta
  
  for(i in 2:n.samples){
    
    #Update beta using MH sampling:
    for(j in 1:2){
      
      att[j] <- att[j] + 1
      
      # Draw candidate:
      canbeta    <- beta
      canbeta[j] <- rnorm(1,beta[j],can.sd)
      canlp      <- log_post(Y,X,canbeta)
      
      # Compute acceptance ratio:
      
      R <- exp(canlp-curlp)  
      U <- runif(1)                          
      if(U<R){       
        beta   <- canbeta
        curlp  <- canlp
        acc[j] <- acc[j]+1
      }
    }
    keep.beta[i,]<-beta
    
  }
  # Return the posterior samples of beta and
  # the Metropolis acceptance rates
  list(beta=keep.beta,acc.rate=acc/att)
  
}


