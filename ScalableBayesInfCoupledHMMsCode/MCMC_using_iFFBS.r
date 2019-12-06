# 
# MCMC for Markov Model using iFFBS algorithm to update the hidden states
# 
#

require("compiler")
enableJIT(3)
require(doParallel)
require(foreach) # parallel computing
registerDoParallel(2) # here we can use up to number of pens depending on how many cores our machine has

##### Calling C from R #####

dyn.load("EcoliFFBScondmod.so")

##############################

# Input data
source("SimulatedDataMarkovModel.r") 

source("iFFBSalgorithm.r")


#########################################################
#########################################################


### Find number of true positives

Positives <- function(R,F,Y)
{
    R <- as.vector(R)
    F <- as.vector(F)
    Y <- as.vector(Y)
    remove_NA <- cbind(R,F,Y)
    remove_NA <- subset(remove_NA, R>-1 & F>-1 & Y >-1) # removed NA
    R <- remove_NA[,1]
    F <- remove_NA[,2]
    Y <- remove_NA[,3]
    data_n11 <- subset(remove_NA, R==1 & F ==1 & Y==1)
    Y1 <- dim(data_n11)[1]
    data_n10 <- subset(remove_NA, R==1 & F ==0 & Y==1)
    Y2 <- dim(data_n10)[1]
    data_n01 <- subset(remove_NA, R==0 & F ==1 & Y==1)
    Y3 <- dim(data_n01)[1]
    data_n00 <- subset(remove_NA, R==0 & F ==0 & Y==1)
    Y4 <- dim(data_n00)[1]
    return(cbind(Y1,Y2,Y3,Y4))
}


## Hamiltonian

grad <- function(a, b, m0, c, s,  p, Pens)
{
    
    likea <- 0
    likeb <- 0
    likem0 <- 0
    
    a <- exp(a)
    b <- exp(b)
    
    
    if ( m0 > 0)
    {
        for (pen in 1:p)
        {
            z <- Pens[,,pen]
            
            for (j in 1:(s-1)) {
                avoid<-length(which(z[,j]==0 & z[,j+1]==0))
                infec<-length(which(z[,j]==0 & z[,j+1]==1))
                inf<-length(which(z[,j]==1))
                ratio <- exp(-a-inf*b)/(1-exp(-a-inf*b))
                likea<-likea -avoid*a + a*infec*ratio
                likeb<-likeb -avoid*inf*b + infec*ratio*inf*b
                n10 <-length(which(z[,j]==1 & z[,j+1]==0))
                n11 <-length(which(z[,j]==1 & z[,j+1]==1))
                likem0 <-likem0  + n11/m0 - (n11+n10)/(m0+1)
            }
        }
        likea <- likea + a1 - b1 * a
        likeb <- likeb + a2 - b2 * b
        likem0 <- likem0 + (a3 - 1)/m0 - b3
        
        return(c(likea, likeb,  likem0))
    } else {return(c(log(likea), log(likeb), log(likem0)))}
}

log_post <- function(a, b, m0, c, s, p, Pens)
{
    a <- exp(a)
    b <- exp(b)
    
    like <- 0
    if (m0 > 0)
    {
        for (pen in 1:p)
        {
            z <- Pens[,,pen]
            
            for (j in 1:(s-1)) {
                avoid<-length(which(z[,j]==0 & z[,j+1]==0))
                infec<-length(which(z[,j]==0 & z[,j+1]==1))
                inf<-length(which(z[,j]==1))
                like<-like+avoid*(-a-inf*b)+infec*log(1-exp(-a-inf*b))
                n10 <-length(which(z[,j]==1 & z[,j+1]==0))
                n11 <-length(which(z[,j]==1 & z[,j+1]==1))
                like <- like + n11*log(m0) -(n11+n10)*log(m0+1)
            }
        }
        a_prior <- dgamma(a, shape = a1, rate = b1, log = TRUE) +log(a)
        b_prior <-  dgamma(b, shape = a2, rate = b2, log = TRUE) +log(b)
        m0_prior <- dgamma(m0, shape = a3, rate = b3, log = TRUE)
        logprior <- a_prior + b_prior + m0_prior
        return(like + logprior)
    } else {return(log(like))}
}



HMC <- function(cura, curb,  curm0, c, s, pen, Pens, epsilon, L)
{
    
    q <- c(cura, curb, curm0)
    p <- rnorm(length(q), 0,1)
    curp <- p
    
    p <- p + epsilon *grad(q[1], q[2], q[3], c, s,  pen, Pens)/2
    
    L <- ceiling(runif(1)*L)-1
    
    for (i in 1:L)
    {
        q <- q + epsilon*p
        if (i!=L) {p <- p + epsilon *grad(q[1], q[2], q[3], c, s,  pen, Pens) }
    }
    
    p <- p + epsilon *grad(q[1], q[2], q[3], c, s,  pen, Pens)/2
    p <- -p
    
    ProposedH  <-  -log_post(q[1], q[2], q[3], c, s, pen, Pens) + 0.5*sum(p^2)
    CurrentH <-  -log_post(cura, curb, curm0, c, s, pen, Pens) + 0.5*sum(curp ^2)
    ap <- -ProposedH + CurrentH
    alpha <- min(1,exp(ap))
    u <- runif(1)

    if (u < alpha)
    { return(c(q, ProposedH))
    } else {return(c(cura, curb, curm0, CurrentH))}
}

# Update hidden states

updateXFFBSmod <- function(x, Pens, tNA, thetaR, thetaF, ObservedR, ObservedF, maxs)
{
    
    for (jj in 1:c)
    {
        rr <- FFBS(exp(x[1]), exp(x[2]), 1/(x[3]+1), x[4], jj, tNA[jj], Pens, thetaR, thetaF, ObservedR[jj,], ObservedF[jj,], maxs)
        Pens[jj, ] <- rr[[1]]
    }
    return(Pens)
}


#################### MCMC algorithm #####################
#########################################################


MCMCMarkovModelusingiFFBS <- function(N, burn, init1, init2, init3, init4,c, s, init5, init6, ObservedR, ObservedF, Pensinit, tNA, maxs)
{
  # N: number of iterations
  # burn: burn in period
  # init*: initial values
  # s: sampling period
  # c: number of cows in a pen
  # ObservedR, ObservedF: Observed RAMS and faecal data, respectively
  # Pensinit: initial infection status matrix
  # tNA: matrix with the day that a cow withdrawn from the study
  # maxs: maximum day of the study
  
  # we use the transformation log(a), log(b)
    
  out <- rep(NA, 6)
  sum_X <- rep(NA,p)
  
  x <- cbind(init1, init2, init3, init4)
  thetaR <- init5
  thetaF <- init6
  Pens <- Pensinit
  epsilon <- 0.03
  L <- 30


  for (i in 1:N)
  { 
    
    

    ### Sample augmented data
    
    
    rr1 <- foreach(kk=1:p) %dopar% {
        updateXFFBSmod(x, Pens[,,kk], tNA[,kk], thetaR, thetaF, ObservedR[,,kk], ObservedF[,,kk], maxs)
    }
    
    
    for (kk in 1:p)
    {
        Pens[ , ,kk] <- rr1[[kk]]
    }
    
    
    # Hamiltonian for updating a, b, m0
    
    x[1:3] <- HMC(x[1], x[2], x[3], c, s, p, Pens, epsilon, L)[1:3]
    
    # Gibbs Sampling for parameter mu
    
    Ininf <- sum(Pens[ ,1,])
    x[4] <- rbeta(n = 1, shape1 = Ininf + a4, shape2 = (c*p) - Ininf + b4, ncp = 0)
    
    
    # Gibbs Sampling for thetaR, thetaF
    
    Posit <- Positives(ObservedR,ObservedF,Pens)
    Y1 <- Posit[1]
    Y2 <- Posit[2]
    Y3 <- Posit[3]
    Y4 <- Posit[4]
    thetaR <- rbeta(n = 1, shape1 = Y1 + Y2 + a5, shape2 = Y3 + Y4 + b5, ncp = 0)
    thetaF <- rbeta(n = 1, shape1 = Y1 + Y3 + a6, shape2 = Y2 + Y4 + b6, ncp = 0)
    
    out <- rbind(out, c(exp(x[1:2]), x[3:4], thetaR, thetaF) )
    
    sum_X<- rbind(sum_X, apply(Pens, 3, function(x) sum(x, na.rm=T)) )
    
    # adjust epsilon
    
    if (((i%%100) == 0) & (i<=burn))  {
        acc <- length(which(diff(out[(i-99):i,1])!=0))/99
        y <- 1 + 1000*(acc-0.7)*(acc-0.7)*(acc-0.7)
        if (y < 0.9) {epsilon <- 0.9*epsilon
        } else if (y > 1.1) { epsilon <- 1.1*epsilon}
    }
    
  }
  list(out, sum_X, epsilon)
}


# prior parameters

a1 <- 1 # for parameter a
b1 <- 1 # for parameter a
a2 <- 1 # for parameter b
b2 <- 1 # for parameter b
a3 <- 0.01 # for parameter m0
b3 <- 0.01 # for parameter m0
a4 <- 1 # for parameter mu
b4 <- 1 # for parameter mu
a5 <- 1 # for parameter thetaR
b5 <- 1 # for parameter thetaR
a6 <- 1 # for parameter thetaF
b6 <- 1 # for parameter thetaF

# initial values

init1 <- init2 <- init4 <- 0.02
init3 <- 11
init5 <- init6 <- 0.5

maxs <- 99


Pens_starold <- array(0, c(c,s,p))
for (i in 1:p)
{
    Pens_starold[ , ,i] <- simul(c, s, init1, init2, init3, init4)
}

Pens_star <- array(NA, c(c,s,p))
for (i in 1:p)
{
    for (j in 1:c)
    {
        Pens_star[j,1:tNA[j,i],i] <- Pens_starold[j,1:tNA[j,i],i]
    }
}


Pens_new <- Pens_star

for (ii in 1:p)
{
    for (jj in 1:c)
    {
        rr <- FFBS(init1, init2, 1/(init3+1), init4, jj, tNA[jj,ii], Pens_new[,,ii], init5, init6, ObservedR[jj,,ii], ObservedF[jj,,ii], maxs)
        Pens_new[jj, ,ii] <- rr[[1]]
    }
    
}

Pensinit <- Pens_new


N <- 11000 # number of iterations
burn <- 1000 # burn in period


Output1 <- MCMCMarkovModelusingiFFBS(N, burn,log(init1), log(init2), init3, init4,c, s, init5, init6, ObservedR, ObservedF, Pensinit, tNA, maxs)
Output <- Output1[[1]][-1, ]
sum_X <- Output1[[2]][-1, ]
epsilon <- Output1[[3]]



