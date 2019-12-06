README file for SupplementaryMaterial.zip

Supplement.pdf -- Additional text and figures from the paper:
  A) Toy example: three interacting individuals
  B) MCMC details for the Markov model
  C) Additional results for the Markov model
  D) MCMC details for the semi-Markov model
  E) Additional results for the semi-Markov model
  F) Application on real data

Outline of all files for the epidemic example in Section 4.2 of: Scalable Bayesian inference for coupled hidden Markov and semi-Markov models

Simulated data set 
SimulatedData.r -- Section 4.2 for Markov epidemic model 

Subroutines: algorithms 
iFFBSalgorithm.r -- iFFBS algorithm
MHiFFBSalgorithm.r -- MHiFFBS algorithm

Subroutines Calling Cpp from R
EcoliFFBScondmod.cpp -- Function for calculating the transition probabilities required by iFFBS method
EcoliFFBScondmod.so -- Compiled file using the command: R CMD SHLIB EcoliFFBScondmod.cpp 

Implementation 
MCMC_using_iFFBS.r -- MCMC code to estimate parameters of Markov epidemic model using the proposed iFFBS algorithm for updating the hidden states
MCMC_using_MHiFFBS.r -- MCMC code to estimate parameters of Markov epidemic model using the proposed MHiFFBS algorithm for updating the hidden states

