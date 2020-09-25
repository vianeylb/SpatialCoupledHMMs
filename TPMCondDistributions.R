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

cond_log_post<-function(tpm.row.obs, X, betarow, row, 
                        prior.mean.int, prior.sd.int, 
                        prior.mean.slope, prior.sd.slope, G){
  
  ## beta is a (G-1)x2 matrix
  ## X is a [T, 2] matrix
  
  Xbeta <- X%*%t(betarow)
  
  rowprobs  <- t(sapply(1:dim(X)[1], function(j) inv.logit(Xbeta[j,], group=row, G=G)))
  
  like   <- sum(sapply(1:dim(X)[1], function(j) dmultinom(tpm.row.obs[j,], size = 1,prob = rowprobs[j,],log=TRUE)))
  
  prior  <- sum(dnorm(betarow[,1] ,prior.mean.int, prior.sd.int,log=TRUE)) + 
    dnorm(betarow[1,2], prior.mean.slope, prior.sd.slope, log=TRUE)
  
  return(like+prior)
  
}

