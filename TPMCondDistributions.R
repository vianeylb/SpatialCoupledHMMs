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

