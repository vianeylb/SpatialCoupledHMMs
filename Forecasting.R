## need alpha_T, then just use the equations from the model to forecast (should be simple enough)
## alpha_T%*%Gamma^h%*%P(x) (for x from some range to another)
## P(x) should just be the density functions

Tlen <- 3000
N <- 3
y <- log(augs$`Wind Speed (MPH)`)[1:3000]


logap <- function(y, rho, sigma, Tlen, N){

  log_all_probs <- matrix(NA, nrow=Tlen, ncol=N)
  
  for(n in 1:N)
    log_all_probs[1,n] =  dnorm(y[1], mean = 0, sd = sigma[n], log = TRUE);
  
  for(j in 1:N) 
    log_all_probs[2,n] = dnorm(y[2], mean =  rho[n,1]*(y[1]-0), sd =  sigma[n], log=TRUE);
  
  for(j in 1:N) 
    log_all_probs[3,n] = dnorm(y[3], mean = rho[n,1]*(y[2]-0) + 
                                       rho[n,2]*(y[1]-0), sd = sigma[n], log = TRUE);
  
  for(j in 1:N) 
    log_all_probs[4,n] = dnorm(y[4], mean = rho[n,1]*(y[3]-0) + 
                                       rho[n,2]*(y[2]-0) + 
                                       rho[n,3]*(y[1]-0), sd = sigma[n], log = TRUE);
  

  for (t in 5:Tlen) { 
    for (n in 1:N){ 
      if(y[t] < 5){
        if(y[t-1] < 5){
          if(y[t-2] < 5){
            if(y[t-3] < 5){
              if(y[t-4] < 5)
                log_all_probs[t,n] = dnorm(y[t], mean = 0 + rho[n,1]*(y[t-1]-0) + 
                                                     rho[n,2]*(y[t-2]-0)+ 
                                                     rho[n,3]*(y[t-3]-0)+ 
                                                     rho[n,4]*(y[t-4]-0), sd = sigma[n], log = TRUE); 
            }
          }
        }
      }
    }
  }
  
  log_all_probs[which(is.na(log_all_probs))] <- 1
  
  return(log_all_probs)
}



rhohr3 <- rstan::extract(fithr4_ws3, pars=c("rho"))[[1]]
sigmahr3 <- rstan::extract(fithr4_ws3, pars=c("sigma"))[[1]]
foo <- rstan::extract(fithr4_ws3, pars=c("init"))[[1]]
gamma <- rstan::extract(fithr4_ws3, pars=c("tpm"))[[1]]

log_all_probs <- logap(y = y, rho = rhohr3[1,,], sigma = sigmahr3[1,], Tlen = 3000, N = 3)
alphahr3 <- foralg(n = 3000, N = 3, log_foo = log(foo[1,]), log_tr_gamma = log(t(gamma[1,,])), log_allprobs = log_all_probs)



hmm_forecast <- function(alpha, gamma, sigma, rho, ys){
  
  alpha%*%gamma
  
}


rhos1o <- data.frame(rhos1[,c(2, 5, 8, 11, 3, 6, 9, 12, 1, 4, 7, 10)], location="Augspurger")
rhos2o <- data.frame(rhos2[,c(2, 5, 8, 11,  1, 4, 7, 10, 3, 6, 9, 12)], location="HoodRiver")
colnames(rhos1o)[1:12] <- paste("Rho ", rep(1:3, each=4), rep(1:4, 3), sep="")
colnames(rhos2o)[1:12] <- paste("Rho ", rep(1:3, each=4), rep(1:4, 3), sep="")

rhoscomb <- rbind(rhos1o, rhos2o)



231 au

16:27



213 hr



