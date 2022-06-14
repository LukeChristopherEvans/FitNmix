
data{
 int n;
 int nsite;
 int nyear;
 int nsplines;
 int maxvisits;
 int nDetectors;

 int site[n];
 int year[n];
 int count[n,maxvisits];
 int Detector[n,maxvisits];
 int SpotSections[n,maxvisits];

 int visits[n]; // visits
 int<lower=0> k[n]; // max at site
 int<lower=0> minpop[n];

 matrix[nyear,nsplines] splineMat; // spline matrix
 }

 parameters{
  // parameters on the ecology
  vector[nsplines] beta1; // regression coefficient for each spline
  vector[nsite] beta0; // intercepts on the ecology

  // parameters on the observation
   real delta0; 
   vector[nDetectors] delta1;
   real delta2;

 }

// move to transpars
 transformed parameters {
 vector[n] loglamda;
 matrix[n,maxvisits] logitp;
 vector[nyear] yearweights = splineMat * beta1; // matrix multiplication

 for (i in 1:n) {
  loglamda[i] = beta0[site[i]] + sum(yearweights[1:year[i]]);
      for(c in 1:visits[i]) {
    logitp[i,c] = delta0 + delta2 * SpotSections[i,c] + delta1[Detector[i,c]] ;
      }
 }

 }

//
model {
  beta0 ~ normal(2,1);
  beta1 ~ normal(0,0.05);
  delta0 ~ normal(0,1.6);
  delta2  ~ normal(0,1.6);
  to_vector(delta1)  ~ normal(0,1.6); 

  for (i in 1:n) {
  
    vector[k[i] - minpop[i] + 1] lp; // this is of length minpop to k but starting at 1
    for(j in minpop[i]:k[i]) {
      lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp[i,1:visits[i]]);//vectorize binomal
    }

   target += log_sum_exp(lp);
  }

 }

generated quantities {
 vector[n] N;
 vector[n] log_lik;

 for(i in 1:n) {
 vector[k[i] - minpop[i] + 1] lp;
 for(j in minpop[i]:k[i]) {
   lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp[i,1:visits[i]]);
 }

  log_lik[i] = log_sum_exp(lp);
  N[i] = poisson_log_rng(loglamda[i]);

  }


 }
