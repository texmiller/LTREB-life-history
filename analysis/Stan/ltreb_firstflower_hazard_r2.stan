data {
  int<lower=1> n_obs;                   // plant-age observations
  int<lower=0,upper=1> flowered[n_obs];    // 1 = flowered in this age, else 0
  int<lower=1> age[n_obs];
  int<lower=1> n_ages;                     // global max age, used only for array sizing
  int<lower=1> n_spp;
  int<lower=1,upper=n_spp> species[n_obs];
  int<lower=0,upper=1> endo[n_obs];     // 0 = endo-, 1 = endo+
  int<lower=1> n_plots;
  int<lower=1,upper=n_plots> plot[n_obs];
  int<lower=1> n_years;
  int<lower=1,upper=n_years> year[n_obs];
  int<lower=1,upper=n_ages> max_age[n_spp];  // species-specific max observed flowering age
  // --- one row per individual for the PPC ---
  int<lower=1> n_cohorts;
  int<lower=1> n_ind;
  int<lower=1,upper=n_spp> species_ind[n_ind];
  int<lower=0,upper=1> endo_ind[n_ind];
  int<lower=1,upper=n_plots> plot_ind[n_ind];
  int<lower=1,upper=n_years> year_ind[n_ind];
  int<lower=1,upper=n_cohorts> cohort_ind[n_ind];
}
parameters {
  real intercept;
  array[n_spp] matrix[2, n_ages] age_eff;  // age_eff[s][e+1, a]: hazard of flowering, species s, endo status e, age a
  // real<lower=0> sigma_age; //age effects drawn from sigma age, mean zero
  vector[n_plots] z_p;
  matrix[n_spp, n_years] z_y;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_year;
}
transformed parameters {
  vector[n_obs] logit_h;
  for (i in 1:n_obs) {
    logit_h[i] = intercept + age_eff[species[i]][endo[i]+1,age[i]]
               + z_p[plot[i]]*sigma_plot
               + z_y[species[i],year[i]]*sigma_year;
  }
}
model {
  intercept~normal(0,1);
  for (s in 1:n_spp){
  to_vector(age_eff[s])~normal(0,1);  
  }
  // sigma_age~normal(0,1);
  sigma_plot~normal(0,1);
  sigma_year~normal(0,1);
  z_p~normal(0,1);
  to_vector(z_y)~normal(0,1);
  flowered ~ bernoulli_logit(logit_h);
}
generated quantities {
  // --- (1) derived mean age of first flowering, per species x endo group ---
  matrix[n_spp,2] mean_age;
  for (s in 1:n_spp) {
    for (e in 1:2) {
      real p_unflowered_prior=1;  // P(not yet flowered through age a-1); starts at 1 before age 1
      real m=0;                   // accumulates the weighted sum: age * P(flowers at that age)
      for (a in 1:max_age[s]) {          // ragged: stop at this species' observed max age
        real h_a = inv_logit(intercept + age_eff[s][e, a]);     // hazard of flowering at this age, this species/endo group
        real p_flower_at_a = p_unflowered_prior * h_a;  // P(flowers exactly at age a) =
                                                          //   P(unflowered through a-1) * P(flowers at a | unflowered)
        m += a * p_flower_at_a;                     // add this age's contribution to the mean
        p_unflowered_prior *= (1 - h_a);             // update to P(unflowered through age a), for next iteration
      }
      mean_age[s, e] = m;
    }
  }

  // --- (2) posterior predictive check: simulate age-at-flowering for each observed individual ---
  // For each individual, start at age 1 and sequentially flip a weighted coin using that
  // individual's own species/endo/plot/cohort hazard: 
  int<lower=1,upper=n_ages> age_rep[n_ind];
  for (i in 1:n_ind) {
    int s = species_ind[i];
    int e = endo_ind[i] + 1;
    int a = 1; //ages start at 1
    int has_flowered = 0;
    while (has_flowered == 0 && a <= n_ages) {
      int yr = cohort_ind[i] + a - 1; //year must advance with age
      if (yr > n_years) {
      break;  // simulated age has run past the years you have data/effects for
      }
      real h_a = inv_logit(intercept+age_eff[s][e,a]+
                            z_p[plot_ind[i]]*sigma_plot+
                            z_y[s,yr]*sigma_year);
      if (bernoulli_rng(h_a) == 1) {
        has_flowered = 1;
      } else {
        a += 1;
      }
    }
    age_rep[i] = min(a,n_ages);  // capped at n_ages if it never flowers in the simulation
  }

}