log_posterior <- beta.cllike(Tlen = 200, N = 2, states_primary = states, states_nghbors = list(states[[2]], states[[1]]), spatialeff = spatialeff[1], tpmeffects = effect) 

gradient <- beta.gradients(Tlen = 200, states_primary = states, states_nghbors = list(states[[2]], states[[1]]), spatialeff = spatialeff[1], tpmeffects = effect)


testhmc <- hmc(log_posterior = log_posterior, gradient = gradient, step_size=0.005, l=10, initP = c(-1.5, -2.5, -1, -3, 0), m = 10000, Tlen = 200, N=2, states_primary = states, states_nghbors = list(states[[2]], states[[1]]))

# effect[[1]] <- rep(-2, 2)
# effect[[2]] <- rep(-2.5, 2)
# spatialeff <- c(-0.5, -0.5)

#----------------------------------------------------------------------------
# visualizing results


