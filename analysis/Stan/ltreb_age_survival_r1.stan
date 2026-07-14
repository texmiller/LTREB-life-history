data {
  int<lower=0> n_spp;               
  int<lower=0> n_years;
  int<lower=0> n_cohorts;
  int<lower=0> n_plots;
  int<lower=0> n_total;
  //Agrostis perennans
  int<lower=0> n_Ap;               
  int<lower=0, upper=1> y_Ap[n_Ap]; 
  int<lower=0> beta_Ap_dim;
  matrix[n_Ap,beta_Ap_dim] X_Ap;
  int<lower=0,upper=n_years> year_Ap[n_Ap];
  int<lower=0,upper=n_cohorts> cohort_Ap[n_Ap];
  int<lower=0,upper=n_plots> plot_Ap[n_Ap];
  //Elymus villosus
  int<lower=0> n_Er;               
  int<lower=0, upper=1> y_Er[n_Er]; 
  int<lower=0> beta_Er_dim;
  matrix[n_Er,beta_Er_dim] X_Er;
  int<lower=0,upper=n_years> year_Er[n_Er];
  int<lower=0,upper=n_cohorts> cohort_Er[n_Er];
  int<lower=0,upper=n_plots> plot_Er[n_Er];
  //Elymus virginicus
  int<lower=0> n_Ev;               
  int<lower=0, upper=1> y_Ev[n_Ev]; 
  int<lower=0> beta_Ev_dim;
  matrix[n_Ev,beta_Ev_dim] X_Ev;
  int<lower=0,upper=n_years> year_Ev[n_Ev];
  int<lower=0,upper=n_cohorts> cohort_Ev[n_Ev];
  int<lower=0,upper=n_plots> plot_Ev[n_Ev];
  //Festuca subverticillata
  int<lower=0> n_Fs;               
  int<lower=0, upper=1> y_Fs[n_Fs]; 
  int<lower=0> beta_Fs_dim;
  matrix[n_Fs,beta_Fs_dim] X_Fs;
  int<lower=0,upper=n_years> year_Fs[n_Fs];
  int<lower=0,upper=n_cohorts> cohort_Fs[n_Fs];
  int<lower=0,upper=n_plots> plot_Fs[n_Fs];
  //Poa alsodes
  int<lower=0> n_Pa;               
  int<lower=0, upper=1> y_Pa[n_Pa]; 
  int<lower=0> beta_Pa_dim;
  matrix[n_Pa,beta_Pa_dim] X_Pa;
  int<lower=0,upper=n_years> year_Pa[n_Pa];
  int<lower=0,upper=n_cohorts> cohort_Pa[n_Pa];
  int<lower=0,upper=n_plots> plot_Pa[n_Pa];
  //Poa autumnalis
  int<lower=0> n_Pu;               
  int<lower=0, upper=1> y_Pu[n_Pu]; 
  int<lower=0> beta_Pu_dim;
  matrix[n_Pu,beta_Pu_dim] X_Pu;
  int<lower=0,upper=n_years> year_Pu[n_Pu];
  int<lower=0,upper=n_cohorts> cohort_Pu[n_Pu];
  int<lower=0,upper=n_plots> plot_Pu[n_Pu];
  //Poa sylvestris
  int<lower=0> n_Ps;               
  int<lower=0, upper=1> y_Ps[n_Ps]; 
  int<lower=0> beta_Ps_dim;
  matrix[n_Ps,beta_Ps_dim] X_Ps;
  int<lower=0,upper=n_years> year_Ps[n_Ps];
  int<lower=0,upper=n_cohorts> cohort_Ps[n_Ps];
  int<lower=0,upper=n_plots> plot_Ps[n_Ps];
}
parameters {
  //species-specific beta vectors
  vector[beta_Ap_dim] beta_Ap;
  vector[beta_Er_dim] beta_Er;
  vector[beta_Ev_dim] beta_Ev;
  vector[beta_Fs_dim] beta_Fs;
  vector[beta_Pa_dim] beta_Pa;
  vector[beta_Pu_dim] beta_Pu;
  vector[beta_Ps_dim] beta_Ps;
  //random intercepts
  matrix[n_spp,n_years] alpha_y;
  matrix[n_spp,n_cohorts] alpha_c;
  vector[n_plots] alpha_p;
  //random effect variances shared across species
  real<lower=0> sigma_year;
  real<lower=0> sigma_cohort;
  real<lower=0> sigma_plot;
}
model {
  beta_Ap ~ normal(0,1);
  beta_Er ~ normal(0,1);
  beta_Ev ~ normal(0,1);
  beta_Fs ~ normal(0,1);
  beta_Pa ~ normal(0,1);
  beta_Pu ~ normal(0,1);
  beta_Ps ~ normal(0,1);
  to_vector(alpha_y) ~ normal(0,sigma_year);
  to_vector(alpha_c) ~ normal(0,sigma_cohort);
  alpha_p ~ normal(0,sigma_plot);
  sigma_year ~ normal(0,1);
  sigma_cohort ~ normal(0,1);
  sigma_plot ~ normal(0,1);
  y_Ap ~ bernoulli_logit_glm(X_Ap,to_vector(alpha_y[1,year_Ap])+to_vector(alpha_c[1,cohort_Ap])+alpha_p[plot_Ap],beta_Ap);
  y_Er ~ bernoulli_logit_glm(X_Er,to_vector(alpha_y[2,year_Er])+to_vector(alpha_c[2,cohort_Er])+alpha_p[plot_Er],beta_Er);
  y_Ev ~ bernoulli_logit_glm(X_Ev,to_vector(alpha_y[3,year_Ev])+to_vector(alpha_c[3,cohort_Ev])+alpha_p[plot_Ev],beta_Ev);
  y_Fs ~ bernoulli_logit_glm(X_Fs,to_vector(alpha_y[4,year_Fs])+to_vector(alpha_c[4,cohort_Fs])+alpha_p[plot_Fs],beta_Fs);
  y_Pa ~ bernoulli_logit_glm(X_Pa,to_vector(alpha_y[5,year_Pa])+to_vector(alpha_c[5,cohort_Pa])+alpha_p[plot_Pa],beta_Pa);
  y_Pu ~ bernoulli_logit_glm(X_Pu,to_vector(alpha_y[6,year_Pu])+to_vector(alpha_c[6,cohort_Pu])+alpha_p[plot_Pu],beta_Pu);
  y_Ps ~ bernoulli_logit_glm(X_Ps,to_vector(alpha_y[7,year_Ps])+to_vector(alpha_c[7,cohort_Ps])+alpha_p[plot_Ps],beta_Ps);
}
generated quantities {
  vector[n_total] log_lik;
  vector[n_Ap] eta_Ap = to_vector(alpha_y[1, year_Ap]) + to_vector(alpha_c[1,cohort_Ap]) + alpha_p[plot_Ap] + X_Ap * beta_Ap;
  vector[n_Er] eta_Er = to_vector(alpha_y[2, year_Er]) + to_vector(alpha_c[2,cohort_Er]) + alpha_p[plot_Er] + X_Er * beta_Er;
  vector[n_Ev] eta_Ev = to_vector(alpha_y[3, year_Ev]) + to_vector(alpha_c[3,cohort_Ev]) + alpha_p[plot_Ev] + X_Ev * beta_Ev;
  vector[n_Fs] eta_Fs = to_vector(alpha_y[4, year_Fs]) + to_vector(alpha_c[4,cohort_Fs]) + alpha_p[plot_Fs] + X_Fs * beta_Fs;
  vector[n_Pa] eta_Pa = to_vector(alpha_y[5, year_Pa]) + to_vector(alpha_c[5,cohort_Pa]) + alpha_p[plot_Pa] + X_Pa * beta_Pa;
  vector[n_Pu] eta_Pu = to_vector(alpha_y[6, year_Pu]) + to_vector(alpha_c[6,cohort_Pu]) + alpha_p[plot_Pu] + X_Pu * beta_Pu;
  vector[n_Ps] eta_Ps = to_vector(alpha_y[7, year_Ps]) + to_vector(alpha_c[7,cohort_Ps]) + alpha_p[plot_Ps] + X_Ps * beta_Ps;
  //calculate log likelihoods in a single vector
  for (n in 1:n_Ap) {
    log_lik[n] = bernoulli_logit_lpmf(y_Ap[n] | eta_Ap[n]);
  }
  for (n in 1:n_Er) {
    log_lik[n_Ap + n] = bernoulli_logit_lpmf(y_Er[n] | eta_Er[n]);
  }
  for (n in 1:n_Ev) {
    log_lik[n_Ap + n_Er + n] = bernoulli_logit_lpmf(y_Ev[n] | eta_Ev[n]);
  }
  for (n in 1:n_Fs) {
    log_lik[n_Ap + n_Er + n_Ev + n] = bernoulli_logit_lpmf(y_Fs[n] | eta_Fs[n]);
  }
  for (n in 1:n_Pa) {
    log_lik[n_Ap + n_Er + n_Ev + n_Fs + n] = bernoulli_logit_lpmf(y_Pa[n] | eta_Pa[n]);
  }
  for (n in 1:n_Pu) {
    log_lik[n_Ap + n_Er + n_Ev + n_Fs + n_Pa + n] = bernoulli_logit_lpmf(y_Pu[n] | eta_Pu[n]);
  }
  for (n in 1:n_Ps) {
    log_lik[n_Ap + n_Er + n_Ev + n_Fs + n_Pa + n_Pu + n] = bernoulli_logit_lpmf(y_Ps[n] | eta_Ps[n]);
  }
}