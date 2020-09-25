simulid = 1 # for example

seedid <- simulid*6 # Simulate same data across methods

################ Simulations of the Data ################
####################  Markov Model  #####################

simul <- function(c, s, a, b, m0, mu)
{
    # c: number of cows in a pen
    # s: number of days
    # a: external infection rate
    # b: with-in pen infection rate
    # m0: mean infectious period - 1
    # mu: initial probability of infection
    
    
    X <- matrix(0, c, s)
    I <- rep(0, s) # number of infected animals at time t
    
    pr <- 1/(m0 + 1)
    p <- mu
    for (i in 1:c)
    {
        U <- runif(1)
        if (U < p ) { X[i,1] <- 1}
    }
    
    for (t in 2:s)
    {
        for (i in 1:c)
        {
            I[t-1] <- sum(X[ ,t-1])
            if (X[i,t-1] == 0)
            { U <- runif(1)
                p <- 1-exp(-a - b*I[t-1])
                if (U < p ) { X[i,t] <- 1}
            } else {
                U <- runif(1)
                p <- 1 - pr
                if (U < p ) { X[i,t] <- 1}
            }    
        }
    }
    return(X)
}



# Simulated values

a0 <- 0.009        # external infection rate
b0 <- 0.01         # with-in pen infection rate
m00 <- 8           # mean infectious period - 1
mu0 <- 0.1         # initial probability of infection
theta0R <- 0.8     # sensitivity of RAMS test
theta0F <- 0.5     # sensitivity of faecal test

c <- cows <- 8      # number of cows in a pen
s <- days <- 99     # number of days
p <- pens <- 20     # number of pens


set.seed(seedid) # Simulate same data across methods

# Construct hidden infection states

Pens_star <- array(0, c(c,s,p))
for (i in 1:p)
{
    Pens_star[ , ,i] <- simul(c, s, a0, b0, m00, mu0)
}


# Construct observed data

ObservedR <- array(0,c(c,s,p))
ObservedF <- array(0,c(c,s,p))

for (j in 1:p)
{
    for (i in 1:c)
    {
        for (t in 1:s)
        {
            if (Pens_star[i,t,j] == 0)
            {
                ObservedR[i,t,j] <- 0
                ObservedF[i,t,j] <- 0
            } else {
                U <- runif(1)
                if (U < theta0R) { ObservedR[i,t,j] <- 1}
                U <- runif(1)
                if (U < theta0F) { ObservedF[i,t,j] <- 1}
            }
        }
    }
}


# Actual sampling frame employed in the real dataset

sample_new <- rep(1, 29)
for (i in 2:29)
{
    if ( (i %% 2) == 0 )
    {
        sample_new[i] <- sample_new[i-1] + 3
    } else { sample_new[i] <-sample_new[i-1] + 4 }
}
sample_new <- sample_new[-18]
sample_new[13] <- sample_new[13] + 1
sample_new[24] <- sample_new[24] + 1

ObservedR[,-sample_new,]  <- NA
ObservedF[,-sample_new,]  <- NA

cNA = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 8, 0, 0, 3, 7, 0, 5, 0, 0, 0, 8) # cow index withdrawn from the study before the completion
tNA_vector = c(100, 100, 100, 100, 57, 100, 100, 100, 100, 57, 100, 100, 71, 50, 100, 81, 100, 100, 100, 99)

for (i in 1:p)
{
    if (tNA_vector[i] < 100)
    {
        ObservedR[cNA[i], tNA_vector[i] : s,i]  <- NA
        ObservedF[cNA[i], tNA_vector[i] : s,i]  <- NA
    }
    
}

tNA = matrix( c(99,99,99,99,56,99,99,99,99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 70, 99, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 99, 99, 99, 80, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 99, 99, 99, 99, 49, 99, 99, 99, 99, 99, 99,
99,99,99,99,99,99,99,99,99, 56, 99, 99, 99, 99, 99, 99, 99, 99, 99, 98), nrow=8, byrow=T)

