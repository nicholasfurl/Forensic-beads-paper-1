data {
  vector[11] prob_rating;
  //int<lower=0, upper=3> claim[11];
  array[11] int<lower=0, upper=3> claim;
  real<lower=0, upper=1> q; 
}

parameters {
  real prior_z;
  real p0_z;
  real log_lambda;
  real log_sigma;
}

transformed parameters {
 real prior_p;
 real p0;
 real lambda;
 real sigma;
 prior_p = Phi(prior_z);
 lambda = exp(log_lambda);
 sigma = exp(log_sigma);
 p0 = Phi(p0_z);
}

model {
  //real prob[11];
  array[11] real prob;
  real lo_prob;
  real n_d;
  real n_g;
  real nu;
  real a;
  real b;
  
  // priors
  log_lambda ~ normal(0,0.2);
  prior_z ~ normal(0, 2);
  p0_z ~ normal(0, 2);
  log_sigma ~ normal(-3,1);
  
  // likelihood
  for(i in 1:11){
    
    if(i==1){
      // subjective probability is equal to prior in first response
      prob[i] = prior_p;
      
      // I am unsure whether this probability should be "distorted"...
      lo_prob = lambda*log(prob[i]/(1-prob[i])) + (1-lambda)*p0;
      prob[i] = exp(lo_prob)/(1+exp(lo_prob));
        
    }else{
    
      // keep track of claims
      n_d = i-1;
      n_g = sum(claim[2:i]);
    
      // compute optimal probability
      prob[i] = prior_p/(prior_p + (q/(1-q)^(n_d-n_g) * (1-prior_p)));
      
      // LO
      lo_prob = lambda*log(prob[i]/(1-prob[i]))+ (1-lambda)*p0;
      prob[i] = exp(lo_prob)/(1+exp(lo_prob));
    }
    
    // compute beta_parameters
    nu = (prob[i]*(1-prob[i]))/sigma^2;
    a = prob[i]*nu;
    b = (1-prob[i])*nu;
    
    // likelihood
    prob_rating[i] ~ beta(a, b);

  }
}
