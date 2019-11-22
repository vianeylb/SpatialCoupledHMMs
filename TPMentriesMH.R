# -----------------------------------------------------------------------
# Metropolis Step for TPM Entries
# -----------------------------------------------------------------------

Bayes.mlogistic <- function(tpm.row.obs, 
                            X, 
                            n.samples=10000, 
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
  # Want G-1 beta intercept terms and one spatial effect term per row
  beta.int <- rnorm(G-1)
  beta.slope <- rnorm(1)
  #G-1X2 matrix
  beta <- matrix(c(beta.int, 
                   rep(beta.slope, G-1)), 
                 nrow=G-1, ncol = 2, 
                 byrow=F)
  
  # Keep track of the samples   
  keep.beta     <- array(0,dim=c(G-1, 2, n.samples))
  keep.beta[,,1] <- beta
  
  
  ## acceptance rates  
  acc   <- att <- rep(0,G)
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
  
  for(i in 2:n.samples){
  
    #Update betas using MH sampling:
    for(j in 1:G){
      
      att[j] <- att[j] + 1
      
      # Draw candidate:
      if(j != G){
        beta[j,1] <- rnorm(n = G-1, mean = keep.beta[j,1,i-1], sd = can.sd)
      }
      if(j == G){
        beta[,2] <- rnorm(n = 1, mean = keep.beta[1,2, i-1], sd=can.sd)
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
    keep.beta[,,i]<-beta
    
  }
  # Return the posterior samples of beta and
  # the Metropolis acceptance rates
  list(beta=keep.beta,acc.rate=acc/att)
  
}


