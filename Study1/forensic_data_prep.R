#
rm(list=ls())
setwd("C:/matlab_files/fiance/forensic_beads_pub_repo/Forensic-beads-paper-1/Study1")

library(tidyverse)

# load data
d <- read_csv("forensic_beads_study1_data.csv")

# select 1 participant
d_i <- d[1:11 + 11,]
data.frame(d_i)

# prepare data for stan
d_stan <- list(
  prob_rating = d_i$HumanProbability / 100,
  claim = d_i$Claim,
  q = 0.6
)

# load stan 
library(rstan)
options(mc.cores = parallel::detectCores())

# run MCMC sampling
stan_fit <- stan(file="beads_single.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))

# visualize posterior
pairs(stan_fit, pars=c("prior_p", "prior_z", "beta","var_beta"))

# trace plot for chain convergence
traceplot(stan_fit)

# more parameters and credible intervals
print(stan_fit, pars=c("prior_p", "beta","var_beta"))
