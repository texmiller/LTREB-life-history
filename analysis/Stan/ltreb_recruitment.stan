data {
  int<lower=0> N;
  int<lower=0> y[N];
  int<lower=0> f1[N];
  int<lower=0> f2[N];
  int<lower=0> n_spp;
  int<lower=0,upper=n_spp> species[N];
  int<lower=0,upper=1> endo[N];
  int<lower=0> n_plots;
  int<lower=0,upper=n_plots> plot[N];
  int<lower=0> n_years;
  int<lower=0,upper=n_years> year[N];
}

parameters {
  vector[n_spp] alpha1;
  vector[n_spp] beta1;
  vector[n_spp] alpha2;
  vector[n_spp] beta2;
  vector[n_plots] rho;
  matrix[n_spp,n_years] tau;
  real<lower=0> sigma_year;
  real<lower=0> sigma_plot;
}

transformed parameters {
  real r1[N];
  real r2[N];
  real<lower=0> lambda[N];
  for(i in 1:N){
    r1[i]=exp(alpha1[species[i]] + beta1[species[i]]*endo[i] + rho[plot[i]] + tau[species[i],year[i]]);
    r2[i]=exp(alpha2[species[i]] + beta2[species[i]]*endo[i] + rho[plot[i]] + tau[species[i],year[i]]);
    lambda[i]=r1[i]*f1[i] + r2[i]*f2[i] + 0.01; //there are some instances of f1=0 and f2=0; lambda must be positive
  }
}

model {
  to_vector(tau) ~ normal(0,sigma_year);
  rho ~ normal(0,sigma_plot);
  y ~ poisson(lambda);
}

generated quantities {
  int<lower = 0> y_sim[N];
  y_sim = poisson_rng(lambda);
}

