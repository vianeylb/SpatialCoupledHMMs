
leapfrog_step <- function(gradient, step_size, position, momentum, d, Tlen, states_primary, states_nghbors) {
  
  gradient <- beta.gradients(Tlen = Tlen, states_primary = states_primary, states_nghbors = states_nghbors, spatialeff = position[5], tpmeffects = list(position[1:2], position[3:4]))
  
  momentum1 <- momentum + gradient * 0.5 * step_size
  position1 <- position + step_size * momentum1
  
  gradient1 <- beta.gradients(Tlen = Tlen, states_primary = states_primary, states_nghbors = states_nghbors, spatialeff = position1[5], tpmeffects = list(position1[1:2], position1[3:4]))
  
  momentum2 <- momentum1 + gradient1 * 0.5 * step_size
  
  matrix(c(position1, momentum2), ncol = d*2)
}

leapfrogs <- function(gradient, step_size, l, position, momentum, d, Tlen, states_primary, states_nghbors) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum, d, Tlen, states_primary, states_nghbors)
    position <- pos_mom[seq_len(d)]
    momentum <- pos_mom[-seq_len(d)]
  }
  pos_mom
}

log_acceptance <- function(propPosition,
                           propMomentum,
                           position,
                           momentum,
                           log_posterior, Tlen, N, states_primary, states_nghbors) {
  
  
  beta.cllike(Tlen, N, states_primary, states_nghbors, spatialeff=propPosition[5], tpmeffects=list(propPosition[1:2], propPosition[3:4])) + sum(dnorm(propMomentum, log = T)) - 
    beta.cllike(Tlen, N, states_primary, states_nghbors, spatialeff=propPosition[5], tpmeffects=list(position[1:2], position[3:4])) - sum(dnorm(momentum, log = T))
  
  #log_posterior(propPosition) + sum(dnorm(propMomentum, log = T)) - 
  #  log_posterior(position) - sum(dnorm(momentum, log = T))
}

hmc_step <- function(log_posterior, gradient, step_size, l, position, Tlen, N, states_primary, states_nghbors) {
  d <- length(position)
  momentum <- rnorm(d)
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum, d, Tlen, states_primary, states_nghbors)
  propPosition <- pos_mom[seq_len(d)]
  propMomentum <- pos_mom[-seq_len(d)]
  a <- log_acceptance(propPosition, propMomentum, position, momentum, log_posterior, Tlen, N, states_primary, states_nghbors)
  #print(a)
  if (log(runif(1)) < a) {
    propPosition
  } else {
    position
  }
}



hmc <- function(log_posterior, gradient, step_size, l, initP, m, Tlen, N, states_primary, states_nghbors){
  out <- matrix(NA_real_, nrow = m, ncol = length(initP))
  out[1, ] <- initP
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i-1,], Tlen = Tlen, N=N, states_primary = states_primary, states_nghbors= states_nghbors)
    
    ## new steps
    log_posterior <- beta.cllike(Tlen = Tlen, N = N, states_primary = states_primary, states_nghbors = states_nghbors, spatialeff = out[i,5], tpmeffects = list(out[i,1:2], out[i,3:4])) 
    gradient <- beta.gradients(Tlen = Tlen, states_primary = states_primary, states_nghbors = states_nghbors, spatialeff = out[i,5], tpmeffects = list(out[i,1:2], out[i,3:4]))
    
  }
  out
}




