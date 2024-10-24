data {
  vector[11] prob_rating;
  //int<lower=0, upper=3> claim[11];
  array[11] int<lower=0, upper=3> claim;
//  real<lower=0, upper=1> q; 
}

parameters {
  real prior_z;
  real<lower=0, upper=1> q; 
  real intercept;
  real slope;
  //vector[2] beta;
  // real<lower=0> sigma;
  real<lower=0> var_beta;
}

transformed parameters {
 real prior_p;
 prior_p = Phi(prior_z);
}

model {
  //real optimal_prob[11];
  array[11] real optimal_prob;
  real n_d;
  real n_g;
  real nu;
  real a;
  real b;
  
  // priors
  intercept ~ normal(0, 1);
  slope ~ normal(1, 1);
  prior_z ~ normal(0, 2);
  q ~ normal(0.6, 1);
  // sigma ~ normal(0,1) T[0,];
  var_beta ~ normal(0,1) T[0,];
  
  // likelihood
  for(i in 1:11){
    
    if(i==1){
      optimal_prob[1] = prior_p;
    }else{
    
      // keep track of claims
      n_d = i-1;
      n_g = sum(claim[2:i]);
    
      // compute optimal probability
      optimal_prob[i] = prior_p/(prior_p + (q/(1-q)^(n_d-n_g) * (1-prior_p)));
    
      // scale optimal probs
      optimal_prob[i] = intercept + slope*optimal_prob[i];
    
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
