#############behavior_LO, start#################
behavior_LO <- function(parameter_means,d_stan){
  
  prior_p = parameter_means["prior_p"]
  lambda = parameter_means["lambda"]
  p0 = parameter_means["p0"]
  q = parameter_means["q"]
  
  claim = d_stan$claim
  
  prob <- rep(NaN,length(claim))
  
  #q = .6
  
  
  for(i in 1:11){
    
    if(i==1){
      # subjective probability is equal to prior in first response
      prob[i] = prior_p;
      
      # I am unsure whether this probability should be "distorted"...
      lo_prob = lambda*log(prob[i]/(1-prob[i])) + (1-lambda)*p0;
      prob[i] = exp(lo_prob)/(1+exp(lo_prob));
      
    }else{
      
      # keep track of claims
      n_d = i-1;
      n_g = sum(claim[2:i]);
      
      # compute optimal probability
      prob[i] = prior_p/(prior_p + (q/(1-q)^(n_d-n_g) * (1-prior_p)));
      
      # LO
      lo_prob = lambda*log(prob[i]/(1-prob[i]))+ (1-lambda)*p0;
      prob[i] = exp(lo_prob)/(1+exp(lo_prob));
    }  #if first claim or not

  }   #for loop through claims
  
  return(prob)

}   #behavior_LO function definition
#############behavior_LO, start#################






#############behavior_split, start#################
behavior_split <- function(parameter_means,d_stan){
  
  
  prior_p = parameter_means["prior_p"]
  q = parameter_means["q"]
  intercept = parameter_means["intercept"]
  slope = parameter_means["slope"]
  
  claim = d_stan$claim
  
  optimal_prob <- rep(NaN,length(claim))
  
  # likelihood
  for(i in 1:11){
    
    if(i==1){
      optimal_prob[1] = prior_p;
    }else{
      
      # keep track of claims
      n_d = i-1;
      n_g = sum(claim[2:i]);
      
      # compute optimal probability
      optimal_prob[i] = prior_p/(prior_p + (q/(1-q)^(n_d-n_g) * (1-prior_p)));
      
      # scale optimal probs
      optimal_prob[i] = intercept + slope*optimal_prob[[i]];
      
    }   #first claim or not?
  }   #loop through claims
  
  return(optimal_prob)
  
} #behavior_split, end
#############behavior_split, end#################






#############behavior_delta, start#################
behavior_delta <- function(parameter_means,d_stan){
  
  prior_p = parameter_means["prior_p"]
  alpha = parameter_means["alpha"]
  invtemp_beta = parameter_means["invtemp_beta"]
  
  claim = d_stan$claim
  
  prob <- rep(NaN,length(claim))

  for(i in 1:11){
    
    if(i==1){
      qhat = prior_p;
      prob[1] = qhat;
    }else{
      
      qhat = qhat + alpha*(claim[i] - qhat);
      prob[i] = exp(invtemp_beta*qhat)/(exp(invtemp_beta*qhat) + exp(invtemp_beta*(1-qhat)));
      
    }   #first claim or not?
    
  }#loop through claims

  return(prob)
    
} #behavior_delta, end
#############behavior_delta, start#################








############################################
#MAIN SCRIPT BODY
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
#options(mc.cores = parallel::detectCores())

# run MCMC sampling
#stan_fit <- stan(file="beads_single_nick.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))
#stan_fit <- stan(file="beads_delta_nick.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))
#stan_fit <- stan(file="beads_split_nick.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))
#stan_fit <- stan(file="beads_single_LO.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))
stan_fit <- stan(file="beads_q_LO.stan", data=d_stan, iter=8000, control=list(adapt_delta=0.9))

# visualize posterior
#pairs(stan_fit, pars=c("prior_p", "prior_z", "alpha","invtemp_beta","var_beta"))
#pairs(stan_fit, pars=c("prior_p", "prior_z", "q","beta","var_beta"))
#pairs(stan_fit, pars=c("prior_p", "lambda","sigma","p0"))
pairs(stan_fit)

# trace plot for chain convergence
traceplot(stan_fit)

# more parameters and credible intervals
#print(stan_fit, pa0rs=c("prior_p", "alpha","q","beta","var_beta"))
print(stan_fit)

#GET PARAMETER MEANS#####

stan_summary <- summary(stan_fit)    # Get the summary of the stan_fit object
parameter_names <- rownames(stan_summary$summary)   # Get the parameter names

parameter_means <- numeric(length(parameter_names))    # Initialize a vector to store the means

# Loop through the parameters and assign each mean to the vector
for (i in seq_along(parameter_names)) {
  parameter_means[i] <- stan_summary$summary[i, "mean"]
}

names(parameter_means) <- parameter_names # Optionally, name the vector elements for clarity

#print(parameter_means) # View the result

#probabilities <- behavior_delta(parameter_means,d_stan)
#probabilities <- behavior_split(parameter_means,d_stan)
probabilities <- behavior_LO(parameter_means,d_stan)

plot(probabilities, type="o", col="blue", pch=16, xlab="Sequence position", ylab="Probability", ylim=c(0, 1))   #Create the line plot for probabilities_LO
points(d_stan$prob_rating, type="o", col="red", pch=16)    ## Add a second line for d_stan$prob_rating with circular markers
legend("topleft", legend=c("Model", "Human"), col=c("blue", "red"), pch=16, lty=1)

