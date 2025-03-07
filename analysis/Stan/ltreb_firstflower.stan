data {
  int<lower=0> N;
  int<lower=0> y[N];
  int<lower=0> n_spp;
  int<lower=0,upper=n_spp> species[N];
  int<lower=0,upper=1> endo[N];
  int<lower=0> n_plots;
  int<lower=0,upper=n_plots> plot[N];
  int<lower=0> n_years;
  int<lower=0,upper=n_years> year[N];
}

parameters {
  vector[n_spp] alpha;
  vector[n_spp] beta;
  vector[n_plots] rho;
  matrix[n_spp,n_years] tau;
  real<lower=0> sigma_year;
  real<lower=0> sigma_plot;
}

transformed parameters {
real<lower=0> mu[N];
  for(i in 1:N){
    mu[i]=exp(alpha[species[i]] + beta[species[i]]*endo[i] + rho[plot[i]] + tau[species[i],year[i]]);
  }
}

model {
  to_vector(tau) ~ normal(0,sigma_year);
  rho ~ normal(0,sigma_plot);
  y ~ poisson(mu);
}

generated quantities {
  int<lower = 0> y_sim[N];
  y_sim = poisson_rng(mu);
}

