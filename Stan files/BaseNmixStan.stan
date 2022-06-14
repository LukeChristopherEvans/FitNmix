
data{
 int n;
 int nsite;
 int nyear;

 int site[n];
 int year[n];
 int maxvisits;
 int count[n,maxvisits]; 
 int visits[n]; // visits
 int<lower=0> k; // max at site
 int<lower=0> minpop[n];
 }

 parameters{
  // parameters on the ecology
  real beta1;
  vector[nsite] beta0; // intercepts on the ecology

  // parameters on the observation
   real delta0; // intercepts on the observation 

 }

//
model {
  to_vector(beta0) ~ normal(2,1);
  beta1 ~ normal(0,0.05);
  delta0 ~ normal(0,1.6);

  for (i in 1:n) {
    real loglamda = beta0[site[i]] + year[i] * beta1; // I later move this to transformed parameters - linear model on the lamba
    real logitp = delta0; // shared p fixed through time
    // new vector for marginalsing over N - here j
    vector[k - minpop[i] + 1] lp; // this is a vector of length minpop to k but starting at 1
    for(j in minpop[i]:k) {
      lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp);
    }

   target += log_sum_exp(lp);
  }

 }

generated quantities {
 vector[n] N;
 real p;
 vector[n] log_lik;

// repeat the loop here to track the likelihood for waic statistics 
 for(i in 1:n) {
 real loglamda = beta0[site[i]] + year[i] * beta1;
 real logitp = delta0;
 vector[k - minpop[i] + 1] lp; 
 for(j in minpop[i]:k) {
   lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp);
 }

  log_lik[i] = log_sum_exp(lp);
  N[i] = poisson_log_rng(loglamda);

  }

   p = inv_logit(delta0);

 }
