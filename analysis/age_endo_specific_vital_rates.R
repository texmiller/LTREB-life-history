## Purpose: fit models for age- and endo-specific survival and fertility
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(scales)
library(xtable)

## functions
invlogit<-function(x){exp(x)/(1+exp(x))}
prop_zero <- function(x) mean(x == 0)

## set working directory
tom<-"C:/Users/tm9/Dropbox/github/LTREB-life-history"
setwd(tom)

## read in QAQC'd data
## as of 4.18.2024 there are still some data issues (see comments) but this is 
## far enough along to proceed
ltreb<-read.csv("./data prep/ltreb_allspp_qaqc.csv") %>% 
  ## drop original plants
  filter(original==0) %>% 
  ## drop LOAR (so few recruits)
  filter(species!="LOAR") %>% 
  ## drop rows with -1 ages (this is first appearance in)
  filter(age>=0) %>% 
  ## convert species to factor
  mutate(species=factor(species))

## what are the age limits we can use for each species
table(ltreb$age,ltreb$species,ltreb$endo_01)
em_tab <- table(ltreb$age[ltreb$endo_01==0],ltreb$species[ltreb$endo_01==0])
ep_tab <- table(ltreb$age[ltreb$endo_01==1],ltreb$species[ltreb$endo_01==1])
tab<-cbind(c(rep(0,times=nrow(em_tab)),rep(1,times=nrow(ep_tab))),
           as.integer(c(rownames(em_tab),rownames(ep_tab))),
           rbind(em_tab,ep_tab))
colnames(tab)[1:2]<-c("Endo","Age")
## Appendix Table
print(xtable(tab,digits=0),include.rownames=FALSE)

## we need a criterion for the max age that we will try to model
## before lumping tail ages as "old"
## for each species and endophyte status, what is the most advanced age
## with at least N observations?
N<-10
ltreb %>% 
  group_by(species,endo_01,age) %>% 
  summarise(count=n()) %>% 
  filter(count>=N) %>% 
  summarise(lump_age=max(age)+1) %>% 
  ## if E+ and E- differ, use the younger age -- *controversial*
  group_by(species) %>% 
  summarise(lump_age=min(lump_age))->age_limits
## lump age is the terminal group that combines ages >=lump age

## create a new age variable that assigns a max age value for all ages over the lower sample size limit
## notice that while this ensures a sample size of N for max_age, then next older group may have < N
ltreb %>% left_join(.,age_limits,by=c("species")) %>% 
  mutate(age_lump=ifelse(age<lump_age,age,lump_age))->ltreb_age_lump

# survival model ----------------------------------------------------------
surv_data<-ltreb_age_lump %>% select(species,endo_01,id,plot,year_t,age,age_lump,surv_t1) %>% drop_na() %>% 
  mutate(species_index = as.numeric(species),
         year_index = year_t-(min(year_t)-1),
         endo_index = endo_01+1,
         ind_index = as.numeric(factor(id)))
## note the species index is alphabetical (1=AGPE, 7=POSY)

## take each species as a data subset -- this helps keep track of the parameters
Ap_surv <- surv_data %>% filter(species=="AGPE") %>% droplevels()
Xs_Ap<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ap_surv)

Er_surv <- surv_data %>% filter(species=="ELRI") %>% droplevels()
Xs_Er<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Er_surv)

Ev_surv <- surv_data %>% filter(species=="ELVI") %>% droplevels()
Xs_Ev<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ev_surv)

Fs_surv <- surv_data %>% filter(species=="FESU") %>% droplevels()
Xs_Fs<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Fs_surv)

Pa_surv <- surv_data %>% filter(species=="POAL") %>% droplevels()
Xs_Pa<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pa_surv)

Pu_surv <- surv_data %>% filter(species=="POAU") %>% droplevels()
Xs_Pu<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pu_surv)

Ps_surv <- surv_data %>% filter(species=="POSY") %>% droplevels()
Xs_Ps<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ps_surv)

stan_dat_surv <- list(n_spp=7,
                 n_years=max(surv_data$year_index),
                 n_plots=max(surv_data$plot),
                 ##AGPE
                 y_Ap=Ap_surv$surv_t1, 
                 n_Ap=length(Ap_surv$surv_t1),
                 beta_Ap_dim=ncol(Xs_Ap),
                 X_Ap=Xs_Ap,
                 year_Ap=Ap_surv$year_index,
                 plot_Ap=Ap_surv$plot,
                 ##ELRI
                 y_Er=Er_surv$surv_t1, 
                 n_Er=length(Er_surv$surv_t1),
                 beta_Er_dim=ncol(Xs_Er),
                 X_Er=Xs_Er,
                 year_Er=Er_surv$year_index,
                 plot_Er=Er_surv$plot,
                 ##ELVI
                 y_Ev=Ev_surv$surv_t1, 
                 n_Ev=length(Ev_surv$surv_t1),
                 beta_Ev_dim=ncol(Xs_Ev),
                 X_Ev=Xs_Ev,
                 year_Ev=Ev_surv$year_index,
                 plot_Ev=Ev_surv$plot,
                 ##FESU
                 y_Fs=Fs_surv$surv_t1, 
                 n_Fs=length(Fs_surv$surv_t1),
                 beta_Fs_dim=ncol(Xs_Fs),
                 X_Fs=Xs_Fs,
                 year_Fs=Fs_surv$year_index,
                 plot_Fs=Fs_surv$plot,
                 ##POAL
                 y_Pa=Pa_surv$surv_t1, 
                 n_Pa=length(Pa_surv$surv_t1),
                 beta_Pa_dim=ncol(Xs_Pa),
                 X_Pa=Xs_Pa,
                 year_Pa=Pa_surv$year_index,
                 plot_Pa=Pa_surv$plot,
                 ##POAU
                 y_Pu=Pu_surv$surv_t1, 
                 n_Pu=length(Pu_surv$surv_t1),
                 beta_Pu_dim=ncol(Xs_Pu),
                 X_Pu=Xs_Pu,
                 year_Pu=Pu_surv$year_index,
                 plot_Pu=Pu_surv$plot,
                 ##POAU
                 y_Ps=Ps_surv$surv_t1, 
                 n_Ps=length(Ps_surv$surv_t1),
                 beta_Ps_dim=ncol(Xs_Ps),
                 X_Ps=Xs_Ps,
                 year_Ps=Ps_surv$year_index,
                 plot_Ps=Ps_surv$plot)

survival_model <- stan_model("analysis/Stan/ltreb_age_survival.stan")
surv_fit<-sampling(survival_model,data = stan_dat_surv,
                       chains=1,
                       #control = list(adapt_delta=0.99,stepsize=0.1),
                       iter=8000,thin=2,
                       pars = c("beta_Ap","beta_Er",
                                "beta_Ev","beta_Fs",
                                "beta_Pa","beta_Pu",
                                "beta_Ps","sigma_year","sigma_plot"), 
                       save_warmup=F)
#write_rds(surv_fit,"analysis/Stan/surv_fit.rds")
surv_fit<-readRDS("analysis/Stan/surv_fit.rds")

## check a few trace plots
bayesplot::mcmc_trace(surv_fit,pars = c("sigma_year","sigma_plot"))
bayesplot::mcmc_trace(surv_fit,pars = c("beta_Ap[1]","beta_Er[1]"))

## wrangle parameter indices to get age- and endo-specific survival
quantile_probs<-c(0.1,0.25,0.5,0.75,0.9)
## Agrostis perennans
Ap_par <- rstan::extract(surv_fit,pars="beta_Ap")
colnames(Ap_par$beta_Ap)<-colnames(Xs_Ap)
age_limits %>% filter(species=="AGPE")## AGPE goes to lump age 6
Ap_em_surv <- invlogit(apply(cbind(Ap_par$beta_Ap[,"(Intercept)"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)2"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)3"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)4"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)5"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)6"]),
                2,quantile,probs=quantile_probs))
Ap_ep_surv <- invlogit(apply(cbind(Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)2"]+Ap_par$beta_Ap[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)3"]+Ap_par$beta_Ap[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)4"]+Ap_par$beta_Ap[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)5"]+Ap_par$beta_Ap[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)6"]+Ap_par$beta_Ap[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Elymus villosus
Er_par <- rstan::extract(surv_fit,pars="beta_Er")
colnames(Er_par$beta_Er)<-colnames(Xs_Er)
age_limits %>% filter(species=="ELRI")## ELRI goes to lump age 4
Er_em_surv <- invlogit(apply(cbind(Er_par$beta_Er[,"(Intercept)"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)2"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)3"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Er_ep_surv <- invlogit(apply(cbind(Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)1"]+Er_par$beta_Er[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)2"]+Er_par$beta_Er[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)3"]+Er_par$beta_Er[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)4"]+Er_par$beta_Er[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


## Elymus virginicus
Ev_par <- rstan::extract(surv_fit,pars="beta_Ev")
colnames(Ev_par$beta_Ev)<-colnames(Xs_Ev)
age_limits %>% filter(species=="ELVI")## ELVI goes to lump age 4
Ev_em_surv <- invlogit(apply(cbind(Ev_par$beta_Ev[,"(Intercept)"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)2"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)3"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Ev_ep_surv <- invlogit(apply(cbind(Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)2"]+Ev_par$beta_Ev[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)3"]+Ev_par$beta_Ev[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)4"]+Ev_par$beta_Ev[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Festuca subverticillata
Fs_par <- rstan::extract(surv_fit,pars="beta_Fs")
colnames(Fs_par$beta_Fs)<-colnames(Xs_Fs)
age_limits %>% filter(species=="FESU")## FESU goes to lump age 6
Fs_em_surv <- invlogit(apply(cbind(Fs_par$beta_Fs[,"(Intercept)"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)2"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)3"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)4"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)5"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)6"]),
                             2,quantile,probs=quantile_probs))
Fs_ep_surv <- invlogit(apply(cbind(Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)2"]+Fs_par$beta_Fs[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)3"]+Fs_par$beta_Fs[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)4"]+Fs_par$beta_Fs[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)5"]+Fs_par$beta_Fs[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)6"]+Fs_par$beta_Fs[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa alsodes
Pa_par <- rstan::extract(surv_fit,pars="beta_Pa")
colnames(Pa_par$beta_Pa)<-colnames(Xs_Pa)
age_limits %>% filter(species=="POAL")## POAL goes to lump age 2
Pa_em_surv <- invlogit(apply(cbind(Pa_par$beta_Pa[,"(Intercept)"],
                                   Pa_par$beta_Pa[,"(Intercept)"]+Pa_par$beta_Pa[,"as.factor(age_lump)1"],
                                   Pa_par$beta_Pa[,"(Intercept)"]+Pa_par$beta_Pa[,"as.factor(age_lump)2"]),
                             2,quantile,probs=quantile_probs))
Pa_ep_surv <- invlogit(apply(cbind(Pa_par$beta_Pa[,"(Intercept)"]+Pa_par$beta_Pa[,"as.factor(endo_01)1"],
                                   Pa_par$beta_Pa[,"(Intercept)"]+Pa_par$beta_Pa[,"as.factor(endo_01)1"]+Pa_par$beta_Pa[,"as.factor(age_lump)1"]+Pa_par$beta_Pa[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pa_par$beta_Pa[,"(Intercept)"]+Pa_par$beta_Pa[,"as.factor(endo_01)1"]+Pa_par$beta_Pa[,"as.factor(age_lump)2"]+Pa_par$beta_Pa[,"as.factor(age_lump)2:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa autumnalis
Pu_par <- rstan::extract(surv_fit,pars="beta_Pu")
colnames(Pu_par$beta_Pu)<-colnames(Xs_Pu)
age_limits %>% filter(species=="POAU")## POAU goes to lump age 5
Pu_em_surv <- invlogit(apply(cbind(Pu_par$beta_Pu[,"(Intercept)"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)2"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)3"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)4"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)5"]),
                             2,quantile,probs=quantile_probs))
Pu_ep_surv <- invlogit(apply(cbind(Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)2"]+Pu_par$beta_Pu[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)3"]+Pu_par$beta_Pu[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)4"]+Pu_par$beta_Pu[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)5"]+Pu_par$beta_Pu[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa sylvestris
Ps_par <- rstan::extract(surv_fit,pars="beta_Ps")
colnames(Ps_par$beta_Ps)<-colnames(Xs_Ps)
age_limits %>% filter(species=="POSY")## POSY goes to lump age 7
Ps_em_surv <- invlogit(apply(cbind(Ps_par$beta_Ps[,"(Intercept)"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)2"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)3"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)4"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)5"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)6"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)7"]),
                             2,quantile,probs=quantile_probs))
Ps_ep_surv <- invlogit(apply(cbind(Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)2"]+Ps_par$beta_Ps[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)3"]+Ps_par$beta_Ps[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)4"]+Ps_par$beta_Ps[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)5"]+Ps_par$beta_Ps[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)6"]+Ps_par$beta_Ps[,"as.factor(age_lump)6:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)7"]+Ps_par$beta_Ps[,"as.factor(age_lump)7:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


## nice figure
pdf("manuscript/figures/age_specific_survival.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(2,4),mar=c(4,4,2,1))
plot(Ap_surv$age_lump,Ap_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,6.5),axes=F)
points(jitter(Ap_surv$age_lump[Ap_surv$endo_01==0])-0.25,
     jitter(Ap_surv$surv_t1[Ap_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ap_surv$age_lump[Ap_surv$endo_01==1])+0.25,
       jitter(Ap_surv$surv_t1[Ap_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Ap_em_surv[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Ap_em_surv[2,1:7],
      (0:6)-.1,Ap_em_surv[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Ap_em_surv[1,1:7],
       (0:6)-.1,Ap_em_surv[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Ap_ep_surv[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Ap_ep_surv[2,1:7],
       (0:6)+.1,Ap_ep_surv[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Ap_ep_surv[1,1:7],
       (0:6)+.1,Ap_ep_surv[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Agrostis perennans",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Er_surv$age_lump,Er_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,4.5),axes=F)
points(jitter(Er_surv$age_lump[Er_surv$endo_01==0])-0.25,
       jitter(Er_surv$surv_t1[Er_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Er_surv$age_lump[Er_surv$endo_01==1])+0.25,
       jitter(Er_surv$surv_t1[Er_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Er_em_surv[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Er_em_surv[2,1:5],
       (0:4)-.1,Er_em_surv[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Er_em_surv[1,1:5],
       (0:4)-.1,Er_em_surv[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Er_ep_surv[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Er_ep_surv[2,1:5],
       (0:4)+.1,Er_ep_surv[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Er_ep_surv[1,1:5],
       (0:4)+.1,Er_ep_surv[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Elymus villosus",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Ev_surv$age_lump,Ev_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,4.5),axes=F)
points(jitter(Ev_surv$age_lump[Ev_surv$endo_01==0])-0.25,
       jitter(Ev_surv$surv_t1[Ev_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ev_surv$age_lump[Ev_surv$endo_01==1])+0.25,
       jitter(Ev_surv$surv_t1[Ev_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Ev_em_surv[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Ev_em_surv[2,1:5],
       (0:4)-.1,Ev_em_surv[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Ev_em_surv[1,1:5],
       (0:4)-.1,Ev_em_surv[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Ev_ep_surv[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Ev_ep_surv[2,1:5],
       (0:4)+.1,Ev_ep_surv[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Ev_ep_surv[1,1:5],
       (0:4)+.1,Ev_ep_surv[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Elymus virginicus",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Fs_surv$age_lump,Fs_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,6.5),axes=F)
points(jitter(Fs_surv$age_lump[Fs_surv$endo_01==0])-0.25,
       jitter(Fs_surv$surv_t1[Fs_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Fs_surv$age_lump[Fs_surv$endo_01==1])+0.25,
       jitter(Fs_surv$surv_t1[Fs_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Fs_em_surv[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Fs_em_surv[2,1:7],
       (0:6)-.1,Fs_em_surv[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Fs_em_surv[1,1:7],
       (0:6)-.1,Fs_em_surv[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Fs_ep_surv[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Fs_ep_surv[2,1:7],
       (0:6)+.1,Fs_ep_surv[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Fs_ep_surv[1,1:7],
       (0:6)+.1,Fs_ep_surv[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Festuca subverticillata",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Pa_surv$age_lump,Pa_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,2.5),axes=F)
points(jitter(Pa_surv$age_lump[Pa_surv$endo_01==0])-0.25,
       jitter(Pa_surv$surv_t1[Pa_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pa_surv$age_lump[Pa_surv$endo_01==1])+0.25,
       jitter(Pa_surv$surv_t1[Pa_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:2)-.1,Pa_em_surv[3,1:3],pch=16,cex=2,col="tomato")
arrows((0:2)-.1,Pa_em_surv[2,1:3],
       (0:2)-.1,Pa_em_surv[4,1:3],length=0,lwd=3,col="tomato")
arrows((0:2)-.1,Pa_em_surv[1,1:3],
       (0:2)-.1,Pa_em_surv[5,1:3],length=0,lwd=1,col="tomato")
points((0:2)+.1,Pa_ep_surv[3,1:3],pch=16,cex=2,col="cornflowerblue")
arrows((0:2)+.1,Pa_ep_surv[2,1:3],
       (0:2)+.1,Pa_ep_surv[4,1:3],length=0,lwd=3,col="cornflowerblue")
arrows((0:2)+.1,Pa_ep_surv[1,1:3],
       (0:2)+.1,Pa_ep_surv[5,1:3],length=0,lwd=1,col="cornflowerblue")
title("Poa alsodes",font.main=3,adj=0)
axis(1,at=0:2,labels=c("0","1","2+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Pu_surv$age_lump,Pu_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,5.5),axes=F)
points(jitter(Pu_surv$age_lump[Pu_surv$endo_01==0])-0.25,
       jitter(Pu_surv$surv_t1[Pu_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pu_surv$age_lump[Pu_surv$endo_01==1])+0.25,
       jitter(Pu_surv$surv_t1[Pu_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Pu_em_surv[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Pu_em_surv[2,1:6],
       (0:5)-.1,Pu_em_surv[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Pu_em_surv[1,1:6],
       (0:5)-.1,Pu_em_surv[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Pu_ep_surv[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Pu_ep_surv[2,1:6],
       (0:5)+.1,Pu_ep_surv[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Pu_ep_surv[1,1:6],
       (0:5)+.1,Pu_ep_surv[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Poa autumnalis",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Ps_surv$age_lump,Ps_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,7.5),axes=F)
points(jitter(Ps_surv$age_lump[Ps_surv$endo_01==0])-0.25,
       jitter(Ps_surv$surv_t1[Ps_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ps_surv$age_lump[Ps_surv$endo_01==1])+0.25,
       jitter(Ps_surv$surv_t1[Ps_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:7)-.1,Ps_em_surv[3,1:8],pch=16,cex=2,col="tomato")
arrows((0:7)-.1,Ps_em_surv[2,1:8],
       (0:7)-.1,Ps_em_surv[4,1:8],length=0,lwd=3,col="tomato")
arrows((0:7)-.1,Ps_em_surv[1,1:8],
       (0:7)-.1,Ps_em_surv[5,1:8],length=0,lwd=1,col="tomato")
points((0:7)+.1,Ps_ep_surv[3,1:8],pch=16,cex=2,col="cornflowerblue")
arrows((0:7)+.1,Ps_ep_surv[2,1:8],
       (0:7)+.1,Ps_ep_surv[4,1:8],length=0,lwd=3,col="cornflowerblue")
arrows((0:7)+.1,Ps_ep_surv[1,1:8],
       (0:7)+.1,Ps_ep_surv[5,1:8],length=0,lwd=1,col="cornflowerblue")
title("Poa sylvestris",font.main=3,adj=0)
axis(1,at=0:7,labels=c("0","1","2","3","4","5","6","7+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(0,0,type="n",axes=F,xlab=" ",ylab=" ")
legend("left",legend=c("E-","E+"),col=c("tomato","cornflowerblue"),pch=16,cex=2)
dev.off()


# fertility model ---------------------------------------------------------
## do zero-year-olds ever flower?
ltreb_age_lump %>% 
  filter(age_lump==0) %>% 
  group_by(species) %>% 
  summarise(mean(flw_count_t,na.rm=T))
## yes, occasionally, and for now I will believe that

fert_data<-ltreb_age_lump %>% select(species,endo_01,id,plot,year_t,age,age_lump,flw_count_t) %>% drop_na() %>% 
  mutate(species_index = as.numeric(species),
         year_index = year_t-(min(year_t)-1),
         endo_index = endo_01+1,
         ind_index = as.numeric(factor(id)))

## take each species as a data subset -- this helps keep track of the parameters
Ap_fert <- fert_data %>% filter(species=="AGPE") %>% droplevels()
Xf_Ap<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ap_fert)

Er_fert <- fert_data %>% filter(species=="ELRI") %>% droplevels()
Xf_Er<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Er_fert)

Ev_fert <- fert_data %>% filter(species=="ELVI") %>% droplevels()
Xf_Ev<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ev_fert)

Fs_fert <- fert_data %>% filter(species=="FESU") %>% droplevels()
Xf_Fs<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Fs_fert)

Pa_fert <- fert_data %>% filter(species=="POAL") %>% droplevels()
Xf_Pa<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pa_fert)

Pu_fert <- fert_data %>% filter(species=="POAU") %>% droplevels()
Xf_Pu<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pu_fert)

Ps_fert <- fert_data %>% filter(species=="POSY") %>% droplevels()
Xf_Ps<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ps_fert)

stan_dat_fert <- list(n_spp=7,
                      n_years=max(fert_data$year_index),
                      n_plots=max(fert_data$plot),
                      ##AGPE
                      y_Ap=Ap_fert$flw_count_t, 
                      n_Ap=length(Ap_fert$flw_count_t),
                      beta_Ap_dim=ncol(Xf_Ap),
                      X_Ap=Xf_Ap,
                      year_Ap=Ap_fert$year_index,
                      plot_Ap=Ap_fert$plot,
                      ##ELRI
                      y_Er=Er_fert$flw_count_t, 
                      n_Er=length(Er_fert$flw_count_t),
                      beta_Er_dim=ncol(Xf_Er),
                      X_Er=Xf_Er,
                      year_Er=Er_fert$year_index,
                      plot_Er=Er_fert$plot,
                      ##ELVI
                      y_Ev=Ev_fert$flw_count_t, 
                      n_Ev=length(Ev_fert$flw_count_t),
                      beta_Ev_dim=ncol(Xf_Ev),
                      X_Ev=Xf_Ev,
                      year_Ev=Ev_fert$year_index,
                      plot_Ev=Ev_fert$plot,
                      ##FESU
                      y_Fs=Fs_fert$flw_count_t, 
                      n_Fs=length(Fs_fert$flw_count_t),
                      beta_Fs_dim=ncol(Xf_Fs),
                      X_Fs=Xf_Fs,
                      year_Fs=Fs_fert$year_index,
                      plot_Fs=Fs_fert$plot,
                      ##POAL
                      y_Pa=Pa_fert$flw_count_t, 
                      n_Pa=length(Pa_fert$flw_count_t),
                      beta_Pa_dim=ncol(Xf_Pa),
                      X_Pa=Xf_Pa,
                      year_Pa=Pa_fert$year_index,
                      plot_Pa=Pa_fert$plot,
                      ##POAU
                      y_Pu=Pu_fert$flw_count_t, 
                      n_Pu=length(Pu_fert$flw_count_t),
                      beta_Pu_dim=ncol(Xf_Pu),
                      X_Pu=Xf_Pu,
                      year_Pu=Pu_fert$year_index,
                      plot_Pu=Pu_fert$plot,
                      ##POAU
                      y_Ps=Ps_fert$flw_count_t, 
                      n_Ps=length(Ps_fert$flw_count_t),
                      beta_Ps_dim=ncol(Xf_Ps),
                      X_Ps=Xf_Ps,
                      year_Ps=Ps_fert$year_index,
                      plot_Ps=Ps_fert$plot)

fertility_model <- stan_model("analysis/Stan/ltreb_age_fertility.stan")
fert_fit<-sampling(fertility_model,data = stan_dat_fert,
                   chains=3,
                   control = list(adapt_delta=0.99,stepsize=0.1),
                   iter=10000,thin=2,
                   pars = c("beta_Ap","beta_Er",
                            "beta_Ev","beta_Fs",
                            "beta_Pa","beta_Pu",
                            "beta_Ps","phi_spp",
                            "sigma_year","sigma_plot",
                            "sim_Ap","sim_Er","sim_Ev","sim_Fs","sim_Pa","sim_Pu","sim_Ps"), 
                   save_warmup=F)
write_rds(fert_fit,"analysis/Stan/fert_fit.rds")
fert_fit<-readRDS("analysis/Stan/fert_fit.rds")

## check a few trace plots
bayesplot::mcmc_trace(fert_fit,pars = c("sigma_year","sigma_plot"))
bayesplot::mcmc_trace(fert_fit,pars = c("beta_Ap[1]","beta_Er[1]"))

## could the negbin generate the zero-heavy data? YES
y_Ap_sim <- rstan::extract(fert_fit,pars="sim_Ap")
ppc_stat(stan_dat_fert$y_Ap, y_Ap_sim$sim_Ap, stat = "prop_zero", binwidth = 0.005)

y_Er_sim <- rstan::extract(fert_fit,pars="sim_Er")
ppc_stat(stan_dat_fert$y_Er, y_Er_sim$sim_Er, stat = "prop_zero", binwidth = 0.005)

y_Ev_sim <- rstan::extract(fert_fit,pars="sim_Ev")
ppc_stat(stan_dat_fert$y_Ev, y_Ev_sim$sim_Ev, stat = "prop_zero", binwidth = 0.005)

y_Fs_sim <- rstan::extract(fert_fit,pars="sim_Fs")
ppc_stat(stan_dat_fert$y_Fs, y_Fs_sim$sim_Fs, stat = "prop_zero", binwidth = 0.005)

y_Pa_sim <- rstan::extract(fert_fit,pars="sim_Pa")
ppc_stat(stan_dat_fert$y_Pa, y_Pa_sim$sim_Pa, stat = "prop_zero", binwidth = 0.005)

y_Pu_sim <- rstan::extract(fert_fit,pars="sim_Pu")
ppc_stat(stan_dat_fert$y_Pu, y_Pu_sim$sim_Pu, stat = "prop_zero", binwidth = 0.005)

y_Ps_sim <- rstan::extract(fert_fit,pars="sim_Ps")
ppc_stat(stan_dat_fert$y_Ps, y_Ps_sim$sim_Ps, stat = "prop_zero", binwidth = 0.005)

## wrangle parameter indices to get age- and endo-specific survival
## Agrostis perennans
Ap_fert_par <- rstan::extract(fert_fit,pars="beta_Ap")
colnames(Ap_fert_par$beta_Ap)<-colnames(Xf_Ap)
age_limits %>% filter(species=="AGPE")## AGPE goes to lump age 6
Ap_em_fert <- exp(apply(cbind(Ap_fert_par$beta_Ap[,"(Intercept)"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)6"]),
                             2,quantile,probs=quantile_probs))
Ap_ep_fert <- exp(apply(cbind(Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)6"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Elymus villosus
Er_par_fert <- rstan::extract(fert_fit,pars="beta_Er")
colnames(Er_par_fert$beta_Er)<-colnames(Xf_Er)
age_limits %>% filter(species=="ELRI")## ELRI goes to lump age 4
Er_em_fert <- invlogit(apply(cbind(Er_par_fert$beta_Er[,"(Intercept)"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)3"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Er_ep_fert <- invlogit(apply(cbind(Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)3"]+Er_par_fert$beta_Er[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)4"]+Er_par_fert$beta_Er[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Elymus virginicus
Ev_par_fert <- rstan::extract(fert_fit,pars="beta_Ev")
colnames(Ev_par_fert$beta_Ev)<-colnames(Xf_Ev)
age_limits %>% filter(species=="ELVI")## ELRI goes to lump age 4
Ev_em_fert <- invlogit(apply(cbind(Ev_par_fert$beta_Ev[,"(Intercept)"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Ev_ep_fert <- invlogit(apply(cbind(Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)4"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Festuca subverticillata
Fs_par_fert <- rstan::extract(fert_fit,pars="beta_Fs")
colnames(Fs_par_fert$beta_Fs)<-colnames(Xf_Fs)
age_limits %>% filter(species=="FESU")## ELRI goes to lump age 4
Fs_em_fert <- invlogit(apply(cbind(Fs_par_fert$beta_Fs[,"(Intercept)"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)6"]),
                             2,quantile,probs=quantile_probs))
Fs_ep_fert <- invlogit(apply(cbind(Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)6"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa alsodes
Pa_par_fert <- rstan::extract(fert_fit,pars="beta_Pa")
colnames(Pa_par_fert$beta_Pa)<-colnames(Xf_Pa)
age_limits %>% filter(species=="POAL")## POAL goes to lump age 2
Pa_em_fert <- invlogit(apply(cbind(Pa_par_fert$beta_Pa[,"(Intercept)"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2"]),
                             2,quantile,probs=quantile_probs))
Pa_ep_fert <- invlogit(apply(cbind(Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa autumnalis
Pu_par_fert <- rstan::extract(fert_fit,pars="beta_Pu")
colnames(Pu_par_fert$beta_Pu)<-colnames(Xf_Pu)
age_limits %>% filter(species=="POAU")## POAU goes to lump age 5
Pu_em_fert <- invlogit(apply(cbind(Pu_par_fert$beta_Pu[,"(Intercept)"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)5"]),
                             2,quantile,probs=quantile_probs))
Pu_ep_fert <- invlogit(apply(cbind(Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)5"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa sylvestris
Ps_par_fert <- rstan::extract(fert_fit,pars="beta_Ps")
colnames(Ps_par_fert$beta_Ps)<-colnames(Xf_Ps)
age_limits %>% filter(species=="POSY")## POSY goes to lump age 7
Ps_em_fert <- invlogit(apply(cbind(Ps_par_fert$beta_Ps[,"(Intercept)"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)7"]),
                             2,quantile,probs=quantile_probs))
Ps_ep_fert <- invlogit(apply(cbind(Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)7"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)7:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


pdf("manuscript/figures/age_specific_fertility.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(2,4),mar=c(4,4,2,1))
ylim_quantile<-0.95
plot(Ap_fert$age_lump,Ap_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,6.5),ylim=c(0,quantile(Ap_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ap_fert$age_lump[Ap_fert$endo_01==0])-0.25,
       jitter(Ap_fert$flw_count_t[Ap_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ap_fert$age_lump[Ap_fert$endo_01==1])+0.25,
       jitter(Ap_fert$flw_count_t[Ap_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Ap_em_fert[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Ap_em_fert[2,1:7],
       (0:6)-.1,Ap_em_fert[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Ap_em_fert[1,1:7],
       (0:6)-.1,Ap_em_fert[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Ap_ep_fert[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Ap_ep_fert[2,1:7],
       (0:6)+.1,Ap_ep_fert[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Ap_ep_fert[1,1:7],
       (0:6)+.1,Ap_ep_fert[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Agrostis perennans",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
axis(2,at=0:round(quantile(Ap_fert$flw_count_t,ylim_quantile)))

plot(Er_fert$age_lump,Er_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,4.5),ylim=c(0,quantile(Er_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Er_fert$age_lump[Er_fert$endo_01==0])-0.25,
       jitter(Er_fert$flw_count_t[Er_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Er_fert$age_lump[Er_fert$endo_01==1])+0.25,
       jitter(Er_fert$flw_count_t[Er_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Er_em_fert[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Er_em_fert[2,1:5],
       (0:4)-.1,Er_em_fert[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Er_em_fert[1,1:5],
       (0:4)-.1,Er_em_fert[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Er_ep_fert[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Er_ep_fert[2,1:5],
       (0:4)+.1,Er_ep_fert[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Er_ep_fert[1,1:5],
       (0:4)+.1,Er_ep_fert[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Elymus villosus",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=0:round(quantile(Er_fert$flw_count_t,ylim_quantile)))

plot(Ev_fert$age_lump,Ev_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,4.5),ylim=c(0,quantile(Ev_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ev_fert$age_lump[Ev_fert$endo_01==0])-0.25,
       jitter(Ev_fert$flw_count_t[Ev_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ev_fert$age_lump[Ev_fert$endo_01==1])+0.25,
       jitter(Ev_fert$flw_count_t[Ev_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Ev_em_fert[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Ev_em_fert[2,1:5],
       (0:4)-.1,Ev_em_fert[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Ev_em_fert[1,1:5],
       (0:4)-.1,Ev_em_fert[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Ev_ep_fert[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Ev_ep_fert[2,1:5],
       (0:4)+.1,Ev_ep_fert[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Ev_ep_fert[1,1:5],
       (0:4)+.1,Ev_ep_fert[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Elymus virginicus",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=0:round(quantile(Ev_fert$flw_count_t,ylim_quantile)))

plot(Fs_fert$age_lump,Fs_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,6.5),ylim=c(0,quantile(Fs_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Fs_fert$age_lump[Fs_fert$endo_01==0])-0.25,
       jitter(Fs_fert$flw_count_t[Fs_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Fs_fert$age_lump[Fs_fert$endo_01==1])+0.25,
       jitter(Fs_fert$flw_count_t[Fs_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Fs_em_fert[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Fs_em_fert[2,1:7],
       (0:6)-.1,Fs_em_fert[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Fs_em_fert[1,1:7],
       (0:6)-.1,Fs_em_fert[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Fs_ep_fert[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Fs_ep_fert[2,1:7],
       (0:6)+.1,Fs_ep_fert[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Fs_ep_fert[1,1:7],
       (0:6)+.1,Fs_ep_fert[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Festuca subverticillata",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
axis(2,at=0:round(quantile(Fs_fert$flw_count_t,ylim_quantile)))

plot(Pa_fert$age_lump,Pa_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,2.5),ylim=c(0,1),axes=F)
points(jitter(Pa_fert$age_lump[Pa_fert$endo_01==0])-0.25,
       jitter(Pa_fert$flw_count_t[Pa_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pa_fert$age_lump[Pa_fert$endo_01==1])+0.25,
       jitter(Pa_fert$flw_count_t[Pa_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:2)-.1,Pa_em_fert[3,1:3],pch=16,cex=2,col="tomato")
arrows((0:2)-.1,Pa_em_fert[2,1:3],
       (0:2)-.1,Pa_em_fert[4,1:3],length=0,lwd=3,col="tomato")
arrows((0:2)-.1,Pa_em_fert[1,1:3],
       (0:2)-.1,Pa_em_fert[5,1:3],length=0,lwd=1,col="tomato")
points((0:2)+.1,Pa_ep_fert[3,1:3],pch=16,cex=2,col="cornflowerblue")
arrows((0:2)+.1,Pa_ep_fert[2,1:3],
       (0:2)+.1,Pa_ep_fert[4,1:3],length=0,lwd=3,col="cornflowerblue")
arrows((0:2)+.1,Pa_ep_fert[1,1:3],
       (0:2)+.1,Pa_ep_fert[5,1:3],length=0,lwd=1,col="cornflowerblue")
title("Poa alsodes",font.main=3,adj=0)
axis(1,at=0:2,labels=c("0","1","2+"))
axis(2,at=0:1)

plot(Pu_fert$age_lump,Pu_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,5.5),ylim=c(0,quantile(Pu_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Pu_fert$age_lump[Pu_fert$endo_01==0])-0.25,
       jitter(Pu_fert$flw_count_t[Pu_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pu_fert$age_lump[Pu_fert$endo_01==1])+0.25,
       jitter(Pu_fert$flw_count_t[Pu_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Pu_em_fert[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Pu_em_fert[2,1:6],
       (0:5)-.1,Pu_em_fert[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Pu_em_fert[1,1:6],
       (0:5)-.1,Pu_em_fert[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Pu_ep_fert[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Pu_ep_fert[2,1:6],
       (0:5)+.1,Pu_ep_fert[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Pu_ep_fert[1,1:6],
       (0:5)+.1,Pu_ep_fert[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Poa autumnalis",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
axis(2,at=0:round(quantile(Pu_fert$flw_count_t,ylim_quantile)))

plot(Ps_fert$age_lump,Ps_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,7.5),ylim=c(0,quantile(Ps_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ps_fert$age_lump[Ps_fert$endo_01==0])-0.25,
       jitter(Ps_fert$flw_count_t[Ps_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ps_fert$age_lump[Ps_fert$endo_01==1])+0.25,
       jitter(Ps_fert$flw_count_t[Ps_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:7)-.1,Ps_em_fert[3,1:8],pch=16,cex=2,col="tomato")
arrows((0:7)-.1,Ps_em_fert[2,1:8],
       (0:7)-.1,Ps_em_fert[4,1:8],length=0,lwd=3,col="tomato")
arrows((0:7)-.1,Ps_em_fert[1,1:8],
       (0:7)-.1,Ps_em_fert[5,1:8],length=0,lwd=1,col="tomato")
points((0:7)+.1,Ps_ep_fert[3,1:8],pch=16,cex=2,col="cornflowerblue")
arrows((0:7)+.1,Ps_ep_fert[2,1:8],
       (0:7)+.1,Ps_ep_fert[4,1:8],length=0,lwd=3,col="cornflowerblue")
arrows((0:7)+.1,Ps_ep_fert[1,1:8],
       (0:7)+.1,Ps_ep_fert[5,1:8],length=0,lwd=1,col="cornflowerblue")
title("Poa sylvestris",font.main=3,adj=0)
axis(1,at=0:7,labels=c("0","1","2","3","4","5","6","7+"))
axis(2,at=0:round(quantile(Ps_fert$flw_count_t,ylim_quantile)))

plot(0,0,type="n",axes=F,xlab=" ",ylab=" ")
legend("left",legend=c("E-","E+"),col=c("tomato","cornflowerblue"),pch=16,cex=2)
dev.off()

# recruitment model -------------------------------------------------------
## estimate the rate of recruitment per inflorescence (which is the unit of fertility)



# basement: data simulation -----------------------------------------------
# simulate a simple bernoulli mixed model and see if I can correctly use stan's bernoulli_logit_glm()

## toy simulation
## three levels of the fixed effect
invlogit<-function(x){exp(x)/(1+exp(x))}
N<-1000
trt<-factor(rep(c("A","B"),each=N/2))
beta <- c(-2,0.5)
x<-model.matrix(~trt-1)
y<-rbinom(N,size=1,prob=invlogit(x%*%beta))

