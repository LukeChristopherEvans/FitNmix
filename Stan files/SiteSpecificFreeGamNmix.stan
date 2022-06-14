
data{
 int n;
 int nsite;
 int nyear;
 int nsplines;
 int maxvisits;

 int site[n];
 int year[n];
 int count[n,maxvisits]; 
 int visits[n]; // visits
 int<lower=0> k; // max at site
 int<lower=0> minpop[n];

 matrix[nyear,nsplines] splineMat; // spline matrix
 }

 parameters{
  // parameters on the ecology
  matrix[nsplines,nsite] zbeta1;
  vector[nsite] zbeta0; // intercepts on the ecology

  // parameters on the observation
   real zdelta0; // intercepts on the observation only to begin
 }

// move to transpars
 transformed parameters {
 vector[n] loglamda;
 real logitp;
 vector[nsite] beta0;
 matrix[nsplines,nsite] beta1;

 beta0 = 2 + zbeta0 * 1;

 for(s in 1:nsplines){
  beta1[s,] =  zbeta1[s,] * 0.05 ;
 }

 for (i in 1:n) {
  vector[nyear] yearweights = splineMat * beta1[,site[i]];
  loglamda[i] = beta0[site[i]] + sum(yearweights[1:year[i]]);
 }
 logitp = zdelta0 * 1.6;
 }

//
model {
  // priors
  to_vector(zbeta0) ~ std_normal();
  zdelta0 ~ std_normal();
  to_vector(zbeta1) ~ std_normal();


 for(i in 1:n) {
   vector[k - minpop[i] + 1] lp; // this is of length minpop to k but starting at 1
    
    for(j in minpop[i]:k) {
        lp[j - minpop[i] + 1]= poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]]| j, logitp);
    }
   
    target += log_sum_exp(lp);
   }


 }

generated quantities {
 vector[n] N;
 real p;
 vector[n] log_lik;

 for(i in 1:n) {
   vector[k - minpop[i] + 1] lp; 
    
    for(j in minpop[i]:k) {
        lp[j - minpop[i] + 1]= poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]]| j, logitp);
    }
   
   log_lik[i] = log_sum_exp(lp);
   N[i] = poisson_log_rng(loglamda[i]);

  }

   p = inv_logit(logitp);

 }
