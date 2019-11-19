###--------------------------------------------------------------
## Conditional distributions for the parameters in the tpm
###--------------------------------------------------------------

# -----------------------------------------------------------------------
# multinomial logistic function
#
# @ `probs`: a vector of probabilities
# @ `group`: an integer in {1,...,G}
# -----------------------------------------------------------------------

mult.logit <- function(probs, group){
  
  m.logit <-  log(probs[-group]/probs[group])
  
  return(m.logit)
  
}

# -----------------------------------------------------------------------
# inverse multinomial logistic function
#
# @ `x`: a vector of X%*%beta values
# @ `group`: an integer in {1,...,G}
# @ `G`: an integer denoting the number of groups
# -----------------------------------------------------------------------

inv.logit <- function(x, group, G){ 
  
  y <- numeric(G)
  y[-group] <- x
  
  rowprobs <- exp(y)/(sum(exp(y)))
  
  return(rowprobs)
  
}

# -----------------------------------------------------------------------
# Computing the conditional log posterior
# 
# @ `Y_row`: A matrix of 1s and 0s denoting the 
#            state process and transitions
# @ `X`: A matrix of covariates
# @ `beta`: A matrix of parameters (estimating)
# @ `prior.mean`: prior means for the beta terms
# @ `prior.sd`: prior standard deviations for the beta terms
# -----------------------------------------------------------------------

cond_log_post<-function(tpm.row, X, betarow, row, prior.mean, prior.sd){
  
  ## beta is a (ncovs + 1)*(G-1) matrix
  ## X is a [T, ncovs +1] matrix
  
  rowprobs  <- sapply(1:dim(X)[1], function(j) inv.logit(X[j,]%*%betarow, group=row, G=G))
  
  like   <- sum(dmultinom(tpm.row, size = 1,prob = rowprobs,log=TRUE))
  
  prior  <- sum(dnorm(betarow ,prior.mean, prior.sd,log=TRUE))
  
  return(like+prior)
  
}

# -----------------------------------------------------------------------
# Computing the conditional log posterior for each row of the TPM
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


