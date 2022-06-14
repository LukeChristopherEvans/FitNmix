functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}

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
 matrix[nsite,nsite] DMat; // distance matrix
 }


 parameters{
  // parameters on the ecology
  matrix[nsplines,nsite] zbeta1;
  vector[nsite] zbeta0; // intercepts on the ecology

  // parameters on the observation
   real zdelta0; // intercepts on the observation only to begin

   // gaussian params
    real<lower=0> etasq;
    real<lower=0> rhosq;

 }


 transformed parameters {
 vector[n] loglamda;
 real logitp;
 vector[nsite] beta0;
 matrix[nsplines,nsite] beta1;

 matrix[nsite, nsite] SigmaDIST;
 matrix[nsite, nsite] LSigmaDIST;

// gaussian componet
 SigmaDIST = cov_GPL2(DMat, etasq, rhosq, 0.1);
 LSigmaDIST= cholesky_decompose(SigmaDIST);

 beta0 = 2 + zbeta0 * 1;
 beta1 =   zbeta1 * LSigmaDIST   ;//matrix multiplication
 
 for (i in 1:n) {
  vector[nyear] yearweights = splineMat * beta1[,site[i]];
  loglamda[i] = beta0[site[i]] + sum(yearweights[1:year[i]]);
 }
 
 logitp = zdelta0 * 1.6;
 }


model {
  to_vector(zbeta0) ~ normal(0,1);
  zdelta0 ~ normal(0,1);

  // Gaussian params
  rhosq ~ exponential(1);
  etasq ~ exponential(1);

 // put priors on the columns
  to_vector(zbeta1) ~ normal(0,1);
 

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

 for (i in 1:n) {
   vector[k - minpop[i] + 1] lp; 
    
    for(j in minpop[i]:k) {
        lp[j - minpop[i] + 1]= poisson_log_lpmf( j |loglamda[i]) + binomial_logit_lpmf(count[i,1:visits[i]]| j, logitp);
    }
   
    log_lik[i] = log_sum_exp(lp);
    N[i] = poisson_log_rng(loglamda[i]);
   }
 
   p = inv_logit(logitp);

  }




