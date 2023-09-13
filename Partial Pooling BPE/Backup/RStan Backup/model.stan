functions {
  // DDE is declared here
  real[] dv_dt
}

data {
  int m; // number of time series to fit
  int n; // length of the times series
  int time0; // initial timepoints of the time series
  real initial_values[5, 10] // 2D array of the initial values for the 10 trajectfries
  int time[m,n]; // placeholder for times of observation
  real vl[m,n]; // placeholder for the observed living tumour
  real dv[m,n]; // placeholder for the observed dead tumour
}

transformed data {
  real x_r[0]
  int x_i[0]; 
}

parameters {
  real mu1_k6;
  real mu2_k6;
  real<lower=0> sigma1_k6
  real<lower=0> sigma2_k6;
  real<lower=0, upper=1> alpha_k6
  real k6[m];

  real mu1_d1;
  real mu2_d1;
  real<lower=0> sigma1_d1
  real<lower=0> sigma2_d1;
  real<lower=0, upper=1> alpha_d1
  real d1[m];

  real mu1_s2;
  real mu2_s2;
  real<lower=0> sigma1_s2
  real<lower=0> sigma2_s2;
  real<lower=0, upper=1> alpha_s2
  real s2[m];
}

transformed parameters {
  // Where the differential equations are solved

}

model {
 // Hyperpriors (ϕ vector)
 mu1_k6 ~ normal(,);
 mu2_k6 ~ normal(,);
 sigma1_k6 ~ normal(,);
 sigma2_k6 ~ normal(,);
 
 mku1_d1 ~ normal(,);
 mu2_d1 ~ normal(,);
 sigma1_d1 ~ normal(,);
 sigma2_d1 ~ normal(,);
 
 mu1_s2 ~ normal(,);
 mu2_s2 ~ normal(,);
 sigma1_s2 ~ normal(,);
 sigma2_s2 ~ normal(,);

  // Priors (θ vector)
  k6 ~ alpha_k6 * normal(mu1_k6, sigma1_k6) + (1 - alpha_k6) * normal(mu2_k6, sigma_2_k6)
  d1 ~ alpha_d1 * normal(mu1_d1, sigma1_d1) + (1 - alpha_d1) * normal(mu2_d1, sigma_2_d1)
  s2 ~ alpha_s2 * normal(mu1_s2, sigma1_s2) + (1 - alpha_s2) * normal(mu2_s2, sigma_2_s2)

  // Likelihood

}


