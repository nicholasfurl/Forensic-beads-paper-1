data {
  vector[11] prob_rating;
  //int<lower=0, upper=3> claim[11];
  array[11] int<lower=0, upper=3> claim;
  //real<lower=0, upper=1> q; 
}

parameters {
  real prior_z;
  //vector[2] beta;
  real<lower=0, upper=1> alpha;
  real<lower=0> invtemp_beta;
  real<lower=0> var_beta;
}

transformed parameters {
 real prior_p;
 prior_p = Phi(prior_z);
}

model {
  //real optimal_prob[11];
  array[11] real optimal_prob;
  //real n_d;
  //real n_g;
  real qhat;
  real nu;
  real a;
  real b;
  
  // priors
  alpha ~ normal(.5, 1);
  invtemp_beta ~ normal(1, 1);
  prior_z ~ normal(0, 2);
  // sigma ~ normal(0,1) T[0,];
  var_beta ~ normal(0,1) T[0,];
  
  // likelihood
  
  for(i in 1:11){
    
    if(i==1){
      qhat = prior_p;
      optimal_prob[1] = qhat;
    }else{
    
      // keep track of claims
      //n_d = i-1;
      //n_g = sum(claim[2:i]);
    
      // compute optimal probability
      //optimal_prob[i] = prior_p/(prior_p + (q/(1-q)^(n_d-n_g) * (1-prior_p)));
      qhat = qhat + alpha*(claim[i] - qhat);
    
      // scale optimal probs
      //optimal_prob[i] = beta[1] + beta[2]*optimal_prob[i];
      optimal_prob[i] = exp(invtemp_beta*qhat)/(exp(invtemp_beta*qhat) + exp(invtemp_beta*(1-qhat)));
    
    }
    
    // compute beta_parameters
    nu = (optimal_prob[i]*(1-optimal_prob[i]))/var_beta;
    a = optimal_prob[i]*nu;
    b = (1-optimal_prob[i])*nu;
    
    // likelihood
    prob_rating[i] ~ beta(a, b);

  }
  // prob_rating ~ normal(beta[1] + beta[2]*optimal_prob, sigma);

}
