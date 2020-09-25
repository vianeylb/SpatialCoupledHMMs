########## MHiFFBS algorithm  ##########
########################################


FFBS <- function(a, b, pr, mu, cow, s, X, thetaR, thetaF, ObservedR, ObservedF)
{
    # cow: index of cow to be updated
    # s: the day that a cow from the pen withdrawn
    # a: external infection rate
    # b: with-in pen infection rate
    # m0: mean infectious period - 1
    # pr: 1/(1+m0)
    # mu: initial probability of infection
    # ObservedR, ObservedF: Observed test results for RAMS and faecal respectively
    # thetaR, thetaF : sensitivity of the RAMS and faecal test respectively
    # X: infection status matrix

CP1 <- matrix(0,s,2)
CP2 <- matrix(0,s,2)
Binary <- c(0,1)
I <- apply(X[-cow, ], 2, function(x) sum(x, na.rm = T))
likelihood <- rep(NA,s)

# for t=1
CP1[1, ] <- mu^Binary * (1-mu)^(1-Binary)
prob_star <- c(1,1)
prob_star <- prob_star* 0^(((ObservedR[1] == 1) | (ObservedF[1] == 1)) * (Binary == 0))
prob_star <- prob_star* (thetaR * thetaF)^(((ObservedR[1]== 1) & (ObservedF[1] == 1)) * (Binary == 1))
prob_star <- prob_star* (thetaR * (1-thetaF))^(((ObservedR[1] == 1) & (ObservedF[1] == 0)) * (Binary == 1))
prob_star <- prob_star* (thetaF * (1-thetaR))^(((ObservedR[1] == 0) & (ObservedF[1] == 1)) * (Binary == 1))
prob_star <- prob_star* ((1-thetaF) * (1-thetaR))^(((ObservedR[1] == 0) & (ObservedF[1] == 0)) * (Binary == 1))
prob_star <- prob_star *  CP1[1,]
CP2[1,] <- prob_star / sum(prob_star)

for (t in 2:s)
{
    p00 <- exp(-a - b*I[t-1])
    p10 <- pr
    CP1[t, 1] <- p00*(CP2[t-1, 1]) + p10*CP2[t-1,2]
    CP1[t, 2] <- (1-p00)*(CP2[t-1, 1]) + (1-p10)*CP2[t-1,2]
    
    if ( is.na(ObservedR[t])==TRUE)
    {
        CP2[t, ] <- CP1[t, ]
    } else {
        prob_star <- c(1,1)
        prob_star <- prob_star* 0^(((ObservedR[t] == 1) | (ObservedF[t] == 1)) * (Binary == 0))
        prob_star <- prob_star* (thetaR * thetaF)^(((ObservedR[t]== 1) & (ObservedF[t] == 1)) * (Binary == 1))
        prob_star <- prob_star* (thetaR * (1-thetaF))^(((ObservedR[t] == 1) & (ObservedF[t] == 0))* (Binary == 1))
        prob_star <- prob_star* (thetaF * (1-thetaR))^(((ObservedR[t] == 0) & (ObservedF[t] == 1))* (Binary == 1))
        prob_star <- prob_star* ((1-thetaF) * (1-thetaR))^(((ObservedR[t] == 0) & (ObservedF[t] == 0)) *(Binary == 1))
        prob_star <- prob_star *  CP1[t,]
        CP2[t,] <- prob_star / sum(prob_star)
    }
}

# Backward Sampling

pbs <- CP2[s, 2]
X[cow,s] <- rbinom(1, 1, prob=pbs)
likelihood[s] <- ifelse(X[cow,s] == 1, pbs, 1-pbs)

for (t in (s-1):1)
{
    p11 <- 1-pr
    p10 <- pr
    
    if (X[cow,t+1] == 1)
    { pbs <- (p11*CP2[t,2])/(CP1[t+1, 2])
        X[cow, t] <- rbinom(1, 1, prob=pbs)
        likelihood[t] <- ifelse(X[cow,t] == 1, pbs, 1-pbs)
    } else {
        pbs <- (p10*CP2[t, 2])/(CP1[t+1,1])
        X[cow, t] <- rbinom(1, 1, prob=pbs)
        likelihood[t] <- ifelse(X[cow,t] == 1, pbs, 1-pbs)
    }
}
list(X[cow, ], sum(log(likelihood)))
}


