
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

 }

 transformed data {
 int veclength[n] ;// vector of lengths of lp

 for (i in 1:n) {
  int veclengthhere = 0;
  for(l in 1:visits[i]) {
  veclengthhere +=  k - count[i,l] + 1;
  }
  veclength[i] = veclengthhere;
 }
 }

 parameters{
  // parameters on the ecology
  real beta1;
  vector[nsite] beta0; // intercepts on the ecology

  // parameters on the observation
   real delta0; // intercepts on the observation only to begin
 }


//
model {
  to_vector(beta0) ~ normal(2,1);
  beta1 ~ normal(0,0.05);
  delta0 ~ normal(0,1.6);

  for (i in 1:n) {
  // new vector for marginalsing over k
    real loglamda = beta0[site[i]] + year[i] * beta1;
    real logitp = delta0;

    vector[veclength[i]] lp; 
    int increment = 1;

    for(m in 1:visits[i]) {
     for(j in count[i,m]:k) { // evaluate by each visit
      lp[increment] = poisson_log_lpmf( j |loglamda) + binomial_logit_lpmf(count[i,m] | j, logitp);
      increment +=  1;
      }
     }
   target += log_sum_exp(lp);
  }

 }

generated quantities {
 vector[n] N;
 real p;
 vector[n] log_lik;

 for(i in 1:n) {
   real loglamda = beta0[site[i]] + year[i] * beta1;
   real logitp = delta0;

   vector[veclength[i]] lp; 
   int increment = 1;

   for(m in 1:visits[i]) {
    for(j in count[i,m]:k) { // evaluate by each visit the minimum population 
     lp[increment] = poisson_log_lpmf( j |loglamda) + binomial_logit_lpmf(count[i,m] | j, logitp);//vectorize binomal
     increment +=  1;
     }
    }

  log_lik[i] = log_sum_exp(lp);
  N[i] = poisson_log_rng(loglamda);

  }

   p = inv_logit(delta0);

 }
