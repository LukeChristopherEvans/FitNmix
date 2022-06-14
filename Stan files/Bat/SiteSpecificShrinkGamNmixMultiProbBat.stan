
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
//

 parameters{
  // parameters on the ecology
  matrix[nsplines,nsite] zbeta1;
  vector[nsite] zbeta0; // intercepts on the ecology

  // parameters on the observation
  real zdelta0; 
  vector [nDetectors] zdelta1;
  real zdelta2;

  vector[nsplines] hyperbeta; // hyper parameters
  vector<lower=0>[nsplines] sigmabeta;
 }


 transformed parameters {

 vector[n] loglamda;
 matrix[n,maxvisits] logitp;
 vector[nsite] beta0;
 matrix[nsplines,nsite] beta1;
 real delta2;
 real delta0;
 vector[nDetectors] delta1;

// re centre
 delta2 =  zdelta2 * 1.6;
 delta0 =  zdelta0 * 1.6;

 for(d in 1:nDetectors){
 delta1[d] =  zdelta1[d] * 1.6;
 }

 beta0 = 2 + zbeta0 * 1;

 for(s in 1:nsplines){
  beta1[s,] = hyperbeta[s] + zbeta1[s,] *sigmabeta[s] ;
 }


 for (i in 1:n) {
  vector[nyear] yearweights = splineMat * beta1[,site[i]];
  loglamda[i] = beta0[site[i]] + sum(yearweights[1:year[i]]);
  for(c in 1:visits[i]) {
       logitp[i,c] = delta0 + delta2 * SpotSections[i,c] + delta1[Detector[i,c]] ;
     }
 }



 }


model {
   // priors
  zdelta0 ~ std_normal();
  to_vector(zdelta1) ~ std_normal();
  zdelta2 ~ std_normal();
  zbeta0 ~ std_normal();
  to_vector(zbeta1) ~ std_normal();
  to_vector(hyperbeta)~normal(0,0.05);
  to_vector(sigmabeta)~exponential(5);


 for (i in 1:n) {
     // new vector for marginalsing over k
   vector[k[i] - minpop[i] + 1] lp;

    for(j in minpop[i]:k[i]) { // evaluate by each possible lambda
     lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp[i,1:visits[i]]);
    }

   target += log_sum_exp(lp);
   }


 }

 generated quantities {
  vector[n] N;
  vector[n] log_lik;

  for (i in 1:n) {
   vector[k[i] - minpop[i] + 1] lp;

    for(j in minpop[i]:k[i]) { 
     lp[j - minpop[i] + 1] = poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]] | j, logitp[i,1:visits[i]]);
    }

    log_lik[i] = log_sum_exp(lp);
    N[i] = poisson_log_rng(loglamda[i]);
   }

  }
