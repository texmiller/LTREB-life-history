data {
  int<lower=0> n_spp;               
  int<lower=0> n_years;
  int<lower=0> n_plots;
  //Agrostis perennans
  int<lower=0> n_Ap;               
  int<lower=0, upper=1> y_Ap[n_Ap]; 
  int<lower=0> beta_Ap_dim;
  matrix[n_Ap,beta_Ap_dim] X_Ap;
  int<lower=0,upper=n_years> year_Ap[n_Ap];
  int<lower=0,upper=n_plots> plot_Ap[n_Ap];
  //Elymus villosus
  int<lower=0> n_Er;               
  int<lower=0, upper=1> y_Er[n_Er]; 
  int<lower=0> beta_Er_dim;
  matrix[n_Er,beta_Er_dim] X_Er;
  int<lower=0,upper=n_years> year_Er[n_Er];
  int<lower=0,upper=n_plots> plot_Er[n_Er];
  //Elymus virginicus
  int<lower=0> n_Ev;               
  int<lower=0, upper=1> y_Ev[n_Ev]; 
  int<lower=0> beta_Ev_dim;
  matrix[n_Ev,beta_Ev_dim] X_Ev;
  int<lower=0,upper=n_years> year_Ev[n_Ev];
  int<lower=0,upper=n_plots> plot_Ev[n_Ev];
  //Festuca subverticillata
  int<lower=0> n_Fs;               
  int<lower=0, upper=1> y_Fs[n_Fs]; 
  int<lower=0> beta_Fs_dim;
  matrix[n_Fs,beta_Fs_dim] X_Fs;
  int<lower=0,upper=n_years> year_Fs[n_Fs];
  int<lower=0,upper=n_plots> plot_Fs[n_Fs];
  //Poa alsodes
  int<lower=0> n_Pa;               
  int<lower=0, upper=1> y_Pa[n_Pa]; 
  int<lower=0> beta_Pa_dim;
  matrix[n_Pa,beta_Pa_dim] X_Pa;
  int<lower=0,upper=n_years> year_Pa[n_Pa];
  int<lower=0,upper=n_plots> plot_Pa[n_Pa];
  //Poa autumnalis
  int<lower=0> n_Pu;               
  int<lower=0, upper=1> y_Pu[n_Pu]; 
  int<lower=0> beta_Pu_dim;
  matrix[n_Pu,beta_Pu_dim] X_Pu;
  int<lower=0,upper=n_years> year_Pu[n_Pu];
  int<lower=0,upper=n_plots> plot_Pu[n_Pu];
  //Poa sylvestris
  int<lower=0> n_Ps;               
  int<lower=0, upper=1> y_Ps[n_Ps]; 
  int<lower=0> beta_Ps_dim;
  matrix[n_Ps,beta_Ps_dim] X_Ps;
  int<lower=0,upper=n_years> year_Ps[n_Ps];
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
  //parameters shared across species
  matrix[n_spp,n_years] alpha_y;
  vector[n_plots] alpha_p;
  real<lower=0> sigma_year;
  real<lower=0> sigma_plot;
}
model {
  //this prior has a mean survival probability of 25%
  //one SD above and below the mean correspond to 0.24% and 98% survival, respectively
  //so, a very accommodating prior
  beta_Ap ~ normal(-1,5);
  beta_Er ~ normal(-1,5);
  beta_Ev ~ normal(-1,5);
  beta_Fs ~ normal(-1,5);
  beta_Pa ~ normal(-1,5);
  beta_Pu ~ normal(-1,5);
  beta_Ps ~ normal(-1,5);
  to_vector(alpha_y) ~ normal(0,sigma_year);
  alpha_p ~ normal(0,sigma_plot);
  y_Ap ~ bernoulli_logit_glm(X_Ap,to_vector(alpha_y[1,year_Ap])+alpha_p[plot_Ap],beta_Ap);
  y_Er ~ bernoulli_logit_glm(X_Er,to_vector(alpha_y[2,year_Er])+alpha_p[plot_Er],beta_Er);
  y_Ev ~ bernoulli_logit_glm(X_Ev,to_vector(alpha_y[3,year_Ev])+alpha_p[plot_Ev],beta_Ev);
  y_Fs ~ bernoulli_logit_glm(X_Fs,to_vector(alpha_y[4,year_Fs])+alpha_p[plot_Fs],beta_Fs);
  y_Pa ~ bernoulli_logit_glm(X_Pa,to_vector(alpha_y[5,year_Pa])+alpha_p[plot_Pa],beta_Pa);
  y_Pu ~ bernoulli_logit_glm(X_Pu,to_vector(alpha_y[6,year_Pu])+alpha_p[plot_Pu],beta_Pu);
  y_Ps ~ bernoulli_logit_glm(X_Ps,to_vector(alpha_y[7,year_Ps])+alpha_p[plot_Ps],beta_Ps);
}