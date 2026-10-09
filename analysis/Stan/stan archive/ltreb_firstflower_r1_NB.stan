data {
  int<lower=0> N;
  int<lower=1> y[N];
  int<lower=0> n_spp;
  int<lower=0,upper=n_spp> species[N];
  int<lower=0,upper=1> endo[N];
  int<lower=0> n_plots;
  int<lower=0,upper=n_plots> plot[N];
  int<lower=0> n_cohorts;
  int<lower=0,upper=n_cohorts> cohort[N];
}

parameters {
  vector[n_spp] alpha;
  vector[n_spp] beta;
  vector[n_plots] z_p;
  matrix[n_spp,n_cohorts] z_c;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_cohort;
  real<lower=0> phi;
}

transformed parameters{
  vector[N] logmu;
  for(i in 1:N){
    logmu[i]=alpha[species[i]] + beta[species[i]]*endo[i] + 
    z_p[plot[i]]*sigma_plot + z_c[species[i],cohort[i]]*sigma_cohort;
  }
}

model {
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  sigma_cohort ~ normal(0,1);
  sigma_plot ~ normal(0,1);
  to_vector(z_c) ~ normal(0,1);
  z_p ~ normal(0,1);
  phi ~ normal(0,1);
  
  for(i in 1:N){
    real log_prob_zero;
    log_prob_zero = neg_binomial_2_log_lpmf(0 | logmu[i], phi);
    target += neg_binomial_2_log_lpmf(y[i] | logmu[i], phi) - 
    log1m_exp(log_prob_zero);
  }

}
