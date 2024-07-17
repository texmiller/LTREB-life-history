## Purpose: fit models for age- and endo-specific survival and fertility
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(scales)
library(xtable)
library(popbio)

## functions
invlogit<-function(x){exp(x)/(1+exp(x))}
prop_zero <- function(x) mean(x == 0)

## set working directory
tom<-"C:/Users/tm9/Dropbox/github/LTREB-life-history"
setwd(tom)

## read in QAQC'd data
## as of 4.18.2024 there are still some data issues (see comments) but this is 
## far enough along to proceed
ltreb_allplants<-read.csv("./data prep/ltreb_allspp_qaqc.csv") %>% 
  ## drop LOAR (so few recruits)
  filter(species!="LOAR") %>% 
  ## convert species to factor
  mutate(species=factor(species))

##drop original plants
ltreb<-ltreb_allplants %>% 
  ## drop rows with -1 ages (this is first appearance in year t1 -- note this also drops originals! [age is NA])
  filter(age>=0) %>% 
  ## drop original plants
  filter(original==0)
  
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
N<-20
ltreb %>% 
  group_by(species,endo_01,age) %>% 
  summarise(count=n()) %>% 
  filter(count>=N) %>% 
  summarise(lump_age=max(age)+1)%>% 
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
age_limits %>% filter(species=="AGPE")## AGPE goes to lump age 5
Ap_em_surv <- invlogit(apply(cbind(Ap_par$beta_Ap[,"(Intercept)"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)2"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)3"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)4"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(age_lump)5"]),
                2,quantile,probs=quantile_probs))
Ap_ep_surv <- invlogit(apply(cbind(Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)2"]+Ap_par$beta_Ap[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)3"]+Ap_par$beta_Ap[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)4"]+Ap_par$beta_Ap[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                Ap_par$beta_Ap[,"(Intercept)"]+Ap_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_par$beta_Ap[,"as.factor(age_lump)5"]+Ap_par$beta_Ap[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
                2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Elymus villosus
Er_par <- rstan::extract(surv_fit,pars="beta_Er")
colnames(Er_par$beta_Er)<-colnames(Xs_Er)
age_limits %>% filter(species=="ELRI")## ELRI goes to lump age 2
Er_em_surv <- invlogit(apply(cbind(Er_par$beta_Er[,"(Intercept)"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(age_lump)2"]),
                             2,quantile,probs=quantile_probs))
Er_ep_surv <- invlogit(apply(cbind(Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)1"]+Er_par$beta_Er[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Er_par$beta_Er[,"(Intercept)"]+Er_par$beta_Er[,"as.factor(endo_01)1"]+Er_par$beta_Er[,"as.factor(age_lump)2"]+Er_par$beta_Er[,"as.factor(age_lump)2:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


## Elymus virginicus
Ev_par <- rstan::extract(surv_fit,pars="beta_Ev")
colnames(Ev_par$beta_Ev)<-colnames(Xs_Ev)
age_limits %>% filter(species=="ELVI")## ELVI goes to lump age 3
Ev_em_surv <- invlogit(apply(cbind(Ev_par$beta_Ev[,"(Intercept)"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)2"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(age_lump)3"]),
                             2,quantile,probs=quantile_probs))
Ev_ep_surv <- invlogit(apply(cbind(Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)2"]+Ev_par$beta_Ev[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ev_par$beta_Ev[,"(Intercept)"]+Ev_par$beta_Ev[,"as.factor(endo_01)1"]+Ev_par$beta_Ev[,"as.factor(age_lump)3"]+Ev_par$beta_Ev[,"as.factor(age_lump)3:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Festuca subverticillata
Fs_par <- rstan::extract(surv_fit,pars="beta_Fs")
colnames(Fs_par$beta_Fs)<-colnames(Xs_Fs)
age_limits %>% filter(species=="FESU")## FESU goes to lump age 5
Fs_em_surv <- invlogit(apply(cbind(Fs_par$beta_Fs[,"(Intercept)"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)2"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)3"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)4"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(age_lump)5"]),
                             2,quantile,probs=quantile_probs))
Fs_ep_surv <- invlogit(apply(cbind(Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)2"]+Fs_par$beta_Fs[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)3"]+Fs_par$beta_Fs[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)4"]+Fs_par$beta_Fs[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Fs_par$beta_Fs[,"(Intercept)"]+Fs_par$beta_Fs[,"as.factor(endo_01)1"]+Fs_par$beta_Fs[,"as.factor(age_lump)5"]+Fs_par$beta_Fs[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
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
age_limits %>% filter(species=="POAU")## POAU goes to lump age 4
Pu_em_surv <- invlogit(apply(cbind(Pu_par$beta_Pu[,"(Intercept)"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)2"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)3"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Pu_ep_surv <- invlogit(apply(cbind(Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)2"]+Pu_par$beta_Pu[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)3"]+Pu_par$beta_Pu[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Pu_par$beta_Pu[,"(Intercept)"]+Pu_par$beta_Pu[,"as.factor(endo_01)1"]+Pu_par$beta_Pu[,"as.factor(age_lump)4"]+Pu_par$beta_Pu[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa sylvestris
Ps_par <- rstan::extract(surv_fit,pars="beta_Ps")
colnames(Ps_par$beta_Ps)<-colnames(Xs_Ps)
age_limits %>% filter(species=="POSY")## POSY goes to lump age 6
Ps_em_surv <- invlogit(apply(cbind(Ps_par$beta_Ps[,"(Intercept)"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)2"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)3"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)4"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)5"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(age_lump)6"]),
                             2,quantile,probs=quantile_probs))
Ps_ep_surv <- invlogit(apply(cbind(Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)2"]+Ps_par$beta_Ps[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)3"]+Ps_par$beta_Ps[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)4"]+Ps_par$beta_Ps[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)5"]+Ps_par$beta_Ps[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Ps_par$beta_Ps[,"(Intercept)"]+Ps_par$beta_Ps[,"as.factor(endo_01)1"]+Ps_par$beta_Ps[,"as.factor(age_lump)6"]+Ps_par$beta_Ps[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


## nice figure
pdf("manuscript/figures/age_specific_survival.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(2,4),mar=c(4,4,2,1))
plot(Ap_surv$age_lump,Ap_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,5.5),axes=F)
points(jitter(Ap_surv$age_lump[Ap_surv$endo_01==0])-0.25,
     jitter(Ap_surv$surv_t1[Ap_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ap_surv$age_lump[Ap_surv$endo_01==1])+0.25,
       jitter(Ap_surv$surv_t1[Ap_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Ap_em_surv[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Ap_em_surv[2,1:6],
      (0:5)-.1,Ap_em_surv[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Ap_em_surv[1,1:6],
       (0:5)-.1,Ap_em_surv[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Ap_ep_surv[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Ap_ep_surv[2,1:6],
       (0:5)+.1,Ap_ep_surv[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Ap_ep_surv[1,1:6],
       (0:5)+.1,Ap_ep_surv[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Agrostis perennans",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Er_surv$age_lump,Er_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,2.5),axes=F)
points(jitter(Er_surv$age_lump[Er_surv$endo_01==0])-0.25,
       jitter(Er_surv$surv_t1[Er_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Er_surv$age_lump[Er_surv$endo_01==1])+0.25,
       jitter(Er_surv$surv_t1[Er_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:2)-.1,Er_em_surv[3,1:3],pch=16,cex=2,col="tomato")
arrows((0:2)-.1,Er_em_surv[2,1:3],
       (0:2)-.1,Er_em_surv[4,1:3],length=0,lwd=3,col="tomato")
arrows((0:2)-.1,Er_em_surv[1,1:3],
       (0:2)-.1,Er_em_surv[5,1:3],length=0,lwd=1,col="tomato")
points((0:2)+.1,Er_ep_surv[3,1:3],pch=16,cex=2,col="cornflowerblue")
arrows((0:2)+.1,Er_ep_surv[2,1:3],
       (0:2)+.1,Er_ep_surv[4,1:3],length=0,lwd=3,col="cornflowerblue")
arrows((0:2)+.1,Er_ep_surv[1,1:3],
       (0:2)+.1,Er_ep_surv[5,1:3],length=0,lwd=1,col="cornflowerblue")
title("Elymus villosus",font.main=3,adj=0)
axis(1,at=0:2,labels=c("0","1","2+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Ev_surv$age_lump,Ev_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,3.5),axes=F)
points(jitter(Ev_surv$age_lump[Ev_surv$endo_01==0])-0.25,
       jitter(Ev_surv$surv_t1[Ev_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ev_surv$age_lump[Ev_surv$endo_01==1])+0.25,
       jitter(Ev_surv$surv_t1[Ev_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:3)-.1,Ev_em_surv[3,1:4],pch=16,cex=2,col="tomato")
arrows((0:3)-.1,Ev_em_surv[2,1:4],
       (0:3)-.1,Ev_em_surv[4,1:4],length=0,lwd=3,col="tomato")
arrows((0:3)-.1,Ev_em_surv[1,1:4],
       (0:3)-.1,Ev_em_surv[5,1:4],length=0,lwd=1,col="tomato")
points((0:3)+.1,Ev_ep_surv[3,1:4],pch=16,cex=2,col="cornflowerblue")
arrows((0:3)+.1,Ev_ep_surv[2,1:4],
       (0:3)+.1,Ev_ep_surv[4,1:4],length=0,lwd=3,col="cornflowerblue")
arrows((0:3)+.1,Ev_ep_surv[1,1:4],
       (0:3)+.1,Ev_ep_surv[5,1:4],length=0,lwd=1,col="cornflowerblue")
title("Elymus virginicus",font.main=3,adj=0)
axis(1,at=0:3,labels=c("0","1","2","3+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Fs_surv$age_lump,Fs_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,5.5),axes=F)
points(jitter(Fs_surv$age_lump[Fs_surv$endo_01==0])-0.25,
       jitter(Fs_surv$surv_t1[Fs_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Fs_surv$age_lump[Fs_surv$endo_01==1])+0.25,
       jitter(Fs_surv$surv_t1[Fs_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Fs_em_surv[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Fs_em_surv[2,1:6],
       (0:5)-.1,Fs_em_surv[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Fs_em_surv[1,1:6],
       (0:5)-.1,Fs_em_surv[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Fs_ep_surv[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Fs_ep_surv[2,1:6],
       (0:5)+.1,Fs_ep_surv[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Fs_ep_surv[1,1:6],
       (0:5)+.1,Fs_ep_surv[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Festuca subverticillata",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
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
     xlim=c(-0.5,4.5),axes=F)
points(jitter(Pu_surv$age_lump[Pu_surv$endo_01==0])-0.25,
       jitter(Pu_surv$surv_t1[Pu_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pu_surv$age_lump[Pu_surv$endo_01==1])+0.25,
       jitter(Pu_surv$surv_t1[Pu_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Pu_em_surv[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Pu_em_surv[2,1:5],
       (0:4)-.1,Pu_em_surv[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Pu_em_surv[1,1:5],
       (0:4)-.1,Pu_em_surv[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Pu_ep_surv[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Pu_ep_surv[2,1:5],
       (0:4)+.1,Pu_ep_surv[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Pu_ep_surv[1,1:5],
       (0:4)+.1,Pu_ep_surv[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Poa autumnalis",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=c(0,0.25,0.5,0.75,1))

plot(Ps_surv$age_lump,Ps_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,6.5),axes=F)
points(jitter(Ps_surv$age_lump[Ps_surv$endo_01==0])-0.25,
       jitter(Ps_surv$surv_t1[Ps_surv$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ps_surv$age_lump[Ps_surv$endo_01==1])+0.25,
       jitter(Ps_surv$surv_t1[Ps_surv$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Ps_em_surv[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Ps_em_surv[2,1:7],
       (0:6)-.1,Ps_em_surv[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Ps_em_surv[1,1:7],
       (0:6)-.1,Ps_em_surv[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Ps_ep_surv[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Ps_ep_surv[2,1:7],
       (0:6)+.1,Ps_ep_surv[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Ps_ep_surv[1,1:7],
       (0:6)+.1,Ps_ep_surv[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Poa sylvestris",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
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
                   iter=8000,thin=2,
                   pars = c("beta_Ap","beta_Er",
                            "beta_Ev","beta_Fs",
                            "beta_Pa","beta_Pu",
                            "beta_Ps","phi_spp",
                            "sigma_year","sigma_plot",
                            "sim_Ap","sim_Er","sim_Ev","sim_Fs","sim_Pa","sim_Pu","sim_Ps"), 
                   save_warmup=F)
#write_rds(fert_fit,"analysis/Stan/fert_fit.rds")
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
age_limits %>% filter(species=="AGPE")## AGPE goes to lump age 5
Ap_em_fert <- exp(apply(cbind(Ap_fert_par$beta_Ap[,"(Intercept)"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5"]),
                             2,quantile,probs=quantile_probs))
Ap_ep_fert <- exp(apply(cbind(Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ap_fert_par$beta_Ap[,"(Intercept)"]+Ap_fert_par$beta_Ap[,"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5"]+Ap_fert_par$beta_Ap[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
                             2,quantile,probs=quantile_probs))

## Elymus villosus
Er_par_fert <- rstan::extract(fert_fit,pars="beta_Er")
colnames(Er_par_fert$beta_Er)<-colnames(Xf_Er)
age_limits %>% filter(species=="ELRI")## ELRI goes to lump age 2
Er_em_fert <- exp(apply(cbind(Er_par_fert$beta_Er[,"(Intercept)"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2"]),
                             2,quantile,probs=quantile_probs))
Er_ep_fert <- exp(apply(cbind(Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Er_par_fert$beta_Er[,"(Intercept)"]+Er_par_fert$beta_Er[,"as.factor(endo_01)1"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2"]+Er_par_fert$beta_Er[,"as.factor(age_lump)2:as.factor(endo_01)1"]),
                             2,quantile,probs=quantile_probs))

## Elymus virginicus
Ev_par_fert <- rstan::extract(fert_fit,pars="beta_Ev")
colnames(Ev_par_fert$beta_Ev)<-colnames(Xf_Ev)
age_limits %>% filter(species=="ELVI")## ELRI goes to lump age 3
Ev_em_fert <- exp(apply(cbind(Ev_par_fert$beta_Ev[,"(Intercept)"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3"]),
                             2,quantile,probs=quantile_probs))
Ev_ep_fert <- exp(apply(cbind(Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ev_par_fert$beta_Ev[,"(Intercept)"]+Ev_par_fert$beta_Ev[,"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3"]+Ev_par_fert$beta_Ev[,"as.factor(age_lump)3:as.factor(endo_01)1"]),
                             2,quantile,probs=quantile_probs))

## Festuca subverticillata
Fs_par_fert <- rstan::extract(fert_fit,pars="beta_Fs")
colnames(Fs_par_fert$beta_Fs)<-colnames(Xf_Fs)
age_limits %>% filter(species=="FESU")## FESU goes to lump age 5
Fs_em_fert <- exp(apply(cbind(Fs_par_fert$beta_Fs[,"(Intercept)"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5"]),
                             2,quantile,probs=quantile_probs))
Fs_ep_fert <- exp(apply(cbind(Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Fs_par_fert$beta_Fs[,"(Intercept)"]+Fs_par_fert$beta_Fs[,"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5"]+Fs_par_fert$beta_Fs[,"as.factor(age_lump)5:as.factor(endo_01)1"]),
                             2,quantile,probs=quantile_probs))

## Poa alsodes
Pa_par_fert <- rstan::extract(fert_fit,pars="beta_Pa")
colnames(Pa_par_fert$beta_Pa)<-colnames(Xf_Pa)
age_limits %>% filter(species=="POAL")## POAL goes to lump age 2
Pa_em_fert <- exp(apply(cbind(Pa_par_fert$beta_Pa[,"(Intercept)"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2"]),
                             2,quantile,probs=quantile_probs))
Pa_ep_fert <- exp(apply(cbind(Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pa_par_fert$beta_Pa[,"(Intercept)"]+Pa_par_fert$beta_Pa[,"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2"]+Pa_par_fert$beta_Pa[,"as.factor(age_lump)2:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))

## Poa autumnalis
Pu_par_fert <- rstan::extract(fert_fit,pars="beta_Pu")
colnames(Pu_par_fert$beta_Pu)<-colnames(Xf_Pu)
age_limits %>% filter(species=="POAU")## POAU goes to lump age 4
Pu_em_fert <- exp(apply(cbind(Pu_par_fert$beta_Pu[,"(Intercept)"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4"]),
                             2,quantile,probs=quantile_probs))
Pu_ep_fert <- exp(apply(cbind(Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Pu_par_fert$beta_Pu[,"(Intercept)"]+Pu_par_fert$beta_Pu[,"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4"]+Pu_par_fert$beta_Pu[,"as.factor(age_lump)4:as.factor(endo_01)1"]),
                             2,quantile,probs=quantile_probs))

## Poa sylvestris
Ps_par_fert <- rstan::extract(fert_fit,pars="beta_Ps")
colnames(Ps_par_fert$beta_Ps)<-colnames(Xf_Ps)
age_limits %>% filter(species=="POSY")## POSY goes to lump age 6
Ps_em_fert <- exp(apply(cbind(Ps_par_fert$beta_Ps[,"(Intercept)"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6"]),
                             2,quantile,probs=quantile_probs))
Ps_ep_fert <- exp(apply(cbind(Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)1:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)2:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)3:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)4:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)5:as.factor(endo_01)1"],
                                   Ps_par_fert$beta_Ps[,"(Intercept)"]+Ps_par_fert$beta_Ps[,"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6"]+Ps_par_fert$beta_Ps[,"as.factor(age_lump)6:as.factor(endo_01)1"]),
                             2,quantile,probs=c(0.1,0.25,0.5,0.75,0.9)))


pdf("manuscript/figures/age_specific_fertility.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(2,4),mar=c(4,4,2,1))
ylim_quantile<-0.975
plot(Ap_fert$age_lump,Ap_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,5.5),ylim=c(0,quantile(Ap_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ap_fert$age_lump[Ap_fert$endo_01==0])-0.25,
       jitter(Ap_fert$flw_count_t[Ap_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ap_fert$age_lump[Ap_fert$endo_01==1])+0.25,
       jitter(Ap_fert$flw_count_t[Ap_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Ap_em_fert[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Ap_em_fert[2,1:6],
       (0:5)-.1,Ap_em_fert[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Ap_em_fert[1,1:6],
       (0:5)-.1,Ap_em_fert[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Ap_ep_fert[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Ap_ep_fert[2,1:6],
       (0:5)+.1,Ap_ep_fert[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Ap_ep_fert[1,1:6],
       (0:5)+.1,Ap_ep_fert[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Agrostis perennans",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
axis(2,at=0:round(quantile(Ap_fert$flw_count_t,ylim_quantile)))

plot(Er_fert$age_lump,Er_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,2.5),ylim=c(0,quantile(Er_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Er_fert$age_lump[Er_fert$endo_01==0])-0.25,
       jitter(Er_fert$flw_count_t[Er_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Er_fert$age_lump[Er_fert$endo_01==1])+0.25,
       jitter(Er_fert$flw_count_t[Er_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:2)-.1,Er_em_fert[3,1:3],pch=16,cex=2,col="tomato")
arrows((0:2)-.1,Er_em_fert[2,1:3],
       (0:2)-.1,Er_em_fert[4,1:3],length=0,lwd=3,col="tomato")
arrows((0:2)-.1,Er_em_fert[1,1:3],
       (0:2)-.1,Er_em_fert[5,1:3],length=0,lwd=1,col="tomato")
points((0:2)+.1,Er_ep_fert[3,1:3],pch=16,cex=2,col="cornflowerblue")
arrows((0:2)+.1,Er_ep_fert[2,1:3],
       (0:2)+.1,Er_ep_fert[4,1:3],length=0,lwd=3,col="cornflowerblue")
arrows((0:2)+.1,Er_ep_fert[1,1:3],
       (0:2)+.1,Er_ep_fert[5,1:3],length=0,lwd=1,col="cornflowerblue")
title("Elymus villosus",font.main=3,adj=0)
axis(1,at=0:2,labels=c("0","1","2+"))
axis(2,at=0:round(quantile(Er_fert$flw_count_t,ylim_quantile)))

plot(Ev_fert$age_lump,Ev_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,3.5),ylim=c(0,quantile(Ev_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ev_fert$age_lump[Ev_fert$endo_01==0])-0.25,
       jitter(Ev_fert$flw_count_t[Ev_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ev_fert$age_lump[Ev_fert$endo_01==1])+0.25,
       jitter(Ev_fert$flw_count_t[Ev_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:3)-.1,Ev_em_fert[3,1:4],pch=16,cex=2,col="tomato")
arrows((0:3)-.1,Ev_em_fert[2,1:4],
       (0:3)-.1,Ev_em_fert[4,1:4],length=0,lwd=3,col="tomato")
arrows((0:3)-.1,Ev_em_fert[1,1:4],
       (0:3)-.1,Ev_em_fert[5,1:4],length=0,lwd=1,col="tomato")
points((0:3)+.1,Ev_ep_fert[3,1:4],pch=16,cex=2,col="cornflowerblue")
arrows((0:3)+.1,Ev_ep_fert[2,1:4],
       (0:3)+.1,Ev_ep_fert[4,1:4],length=0,lwd=3,col="cornflowerblue")
arrows((0:3)+.1,Ev_ep_fert[1,1:4],
       (0:3)+.1,Ev_ep_fert[5,1:4],length=0,lwd=1,col="cornflowerblue")
title("Elymus virginicus",font.main=3,adj=0)
axis(1,at=0:3,labels=c("0","1","2","3+"))
axis(2,at=0:round(quantile(Ev_fert$flw_count_t,ylim_quantile)))

plot(Fs_fert$age_lump,Fs_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,5.5),ylim=c(0,quantile(Fs_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Fs_fert$age_lump[Fs_fert$endo_01==0])-0.25,
       jitter(Fs_fert$flw_count_t[Fs_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Fs_fert$age_lump[Fs_fert$endo_01==1])+0.25,
       jitter(Fs_fert$flw_count_t[Fs_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:5)-.1,Fs_em_fert[3,1:6],pch=16,cex=2,col="tomato")
arrows((0:5)-.1,Fs_em_fert[2,1:6],
       (0:5)-.1,Fs_em_fert[4,1:6],length=0,lwd=3,col="tomato")
arrows((0:5)-.1,Fs_em_fert[1,1:6],
       (0:5)-.1,Fs_em_fert[5,1:6],length=0,lwd=1,col="tomato")
points((0:5)+.1,Fs_ep_fert[3,1:6],pch=16,cex=2,col="cornflowerblue")
arrows((0:5)+.1,Fs_ep_fert[2,1:6],
       (0:5)+.1,Fs_ep_fert[4,1:6],length=0,lwd=3,col="cornflowerblue")
arrows((0:5)+.1,Fs_ep_fert[1,1:6],
       (0:5)+.1,Fs_ep_fert[5,1:6],length=0,lwd=1,col="cornflowerblue")
title("Festuca subverticillata",font.main=3,adj=0)
axis(1,at=0:5,labels=c("0","1","2","3","4","5+"))
axis(2,at=0:round(quantile(Fs_fert$flw_count_t,ylim_quantile)))

plot(Pa_fert$age_lump,Pa_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,2.5),c(0,quantile(Pa_fert$flw_count_t,ylim_quantile)),axes=F)
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
axis(2,at=0:round(quantile(Pa_fert$flw_count_t,ylim_quantile)))
##note different y lim for Pu
plot(Pu_fert$age_lump,Pu_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,4.5),ylim=c(0,10),axes=F)
points(jitter(Pu_fert$age_lump[Pu_fert$endo_01==0])-0.25,
       jitter(Pu_fert$flw_count_t[Pu_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Pu_fert$age_lump[Pu_fert$endo_01==1])+0.25,
       jitter(Pu_fert$flw_count_t[Pu_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:4)-.1,Pu_em_fert[3,1:5],pch=16,cex=2,col="tomato")
arrows((0:4)-.1,Pu_em_fert[2,1:5],
       (0:4)-.1,Pu_em_fert[4,1:5],length=0,lwd=3,col="tomato")
arrows((0:4)-.1,Pu_em_fert[1,1:5],
       (0:4)-.1,Pu_em_fert[5,1:5],length=0,lwd=1,col="tomato")
points((0:4)+.1,Pu_ep_fert[3,1:5],pch=16,cex=2,col="cornflowerblue")
arrows((0:4)+.1,Pu_ep_fert[2,1:5],
       (0:4)+.1,Pu_ep_fert[4,1:5],length=0,lwd=3,col="cornflowerblue")
arrows((0:4)+.1,Pu_ep_fert[1,1:5],
       (0:4)+.1,Pu_ep_fert[5,1:5],length=0,lwd=1,col="cornflowerblue")
title("Poa autumnalis",font.main=3,adj=0)
axis(1,at=0:4,labels=c("0","1","2","3","4+"))
axis(2,at=0:10)

plot(Ps_fert$age_lump,Ps_fert$flw_count_t,type="n",xlab="Age group",ylab="Fertility (# infs)",
     xlim=c(-0.5,6.5),ylim=c(0,quantile(Ps_fert$flw_count_t,ylim_quantile)),axes=F)
points(jitter(Ps_fert$age_lump[Ps_fert$endo_01==0])-0.25,
       jitter(Ps_fert$flw_count_t[Ps_fert$endo_01==0],factor=0.1),col=alpha("tomato",0.25))
points(jitter(Ps_fert$age_lump[Ps_fert$endo_01==1])+0.25,
       jitter(Ps_fert$flw_count_t[Ps_fert$endo_01==1],factor=0.1),col=alpha("cornflowerblue",0.25))
points((0:6)-.1,Ps_em_fert[3,1:7],pch=16,cex=2,col="tomato")
arrows((0:6)-.1,Ps_em_fert[2,1:7],
       (0:6)-.1,Ps_em_fert[4,1:7],length=0,lwd=3,col="tomato")
arrows((0:6)-.1,Ps_em_fert[1,1:7],
       (0:6)-.1,Ps_em_fert[5,1:7],length=0,lwd=1,col="tomato")
points((0:6)+.1,Ps_ep_fert[3,1:7],pch=16,cex=2,col="cornflowerblue")
arrows((0:6)+.1,Ps_ep_fert[2,1:7],
       (0:6)+.1,Ps_ep_fert[4,1:7],length=0,lwd=3,col="cornflowerblue")
arrows((0:6)+.1,Ps_ep_fert[1,1:7],
       (0:6)+.1,Ps_ep_fert[5,1:7],length=0,lwd=1,col="cornflowerblue")
title("Poa sylvestris",font.main=3,adj=0)
axis(1,at=0:6,labels=c("0","1","2","3","4","5","6+"))
axis(2,at=0:round(quantile(Ps_fert$flw_count_t,ylim_quantile)))

plot(0,0,type="n",axes=F,xlab=" ",ylab=" ")
legend("left",legend=c("E-","E+"),col=c("tomato","cornflowerblue"),pch=16,cex=2)
dev.off()

# recruitment model -------------------------------------------------------
## estimate the rate of recruitment per inflorescence (which is the unit of fertility)

## first summarise total number of infs produced per plot per year -- including originals!
ltreb_allplants %>% 
  group_by(species,plot,endo_01,year_t) %>% 
  summarise(total_infs = sum(flw_count_t,na.rm=T)) %>% 
  mutate(year_tp1 = year_t+1,
         year_tp2 = year_t+2)->inf_counts

## second summarise total number of recruits (0-yo's) per plot per year
ltreb_age_lump %>% 
  filter(age_lump==0) %>% 
  group_by(species,plot,endo_01,year_t) %>% 
  summarise(recruits=n())->recruit_counts

## join two data frames -- this uses two joins to get the 2-year lag
left_join(left_join(recruit_counts,
          inf_counts %>% select(species,plot,endo_01,year_tp1,total_infs),
          by=c("species","plot","endo_01","year_t"="year_tp1")),
          inf_counts %>% select(species,plot,endo_01,year_tp2,total_infs),
          by=c("species","plot","endo_01","year_t"="year_tp2")) %>% 
  rename(recruits_t = recruits,
         total_infs_tm1 = total_infs.x,
         total_infs_tm2 = total_infs.y) -> recruits_infs

plot(recruits_infs$total_infs_tm1,recruits_infs$recruits_t)
plot(recruits_infs$total_infs_tm2,recruits_infs$recruits_t)
plot(recruits_infs$total_infs_tm2+recruits_infs$total_infs_tm1,recruits_infs$recruits_t)

## are there any cases of non-zero recruitment with zero infs in the last two years?
recruits_infs %>% filter(recruits_t>0 & total_infs_tm1==0 & total_infs_tm2==0)
## yeah this happened in 21 plot-years -- so there is potentially a longer-lived seed bank
## but I am going to limit this to a two-year seed bank

## there are some NAs in the inf data...why?
recruits_infs %>% filter(is.na(total_infs_tm1))
inf_counts %>% filter(species=="ELVI",plot==95)
## there are some plots that have no inf data for some years
## not sure why, but proceeding for now to drop NAs
recruits_infs %>% 
  mutate(species_index = as.numeric(species),
         year_index = year_t-(min(year_t)-1)) %>% 
  drop_na() -> recruits_infs

## if I just wanted raw ratios of recruits:infs, what would that look like?
hist(recruits_infs$recruits_t/recruits_infs$total_infs_tm1)
hist(recruits_infs$recruits_t/recruits_infs$total_infs_tm2)

## prep data
recruit_dat<-list(N=nrow(recruits_infs),
                  y=recruits_infs$recruits_t,
                  f1=recruits_infs$total_infs_tm1,
                  f2=recruits_infs$total_infs_tm2,
                  n_spp=length(levels(recruits_infs$species)),
                  n_years=max(recruits_infs$year_index),
                  n_plots=max(recruits_infs$plot),
                  species=recruits_infs$species_index,
                  endo=recruits_infs$endo_01,
                  plot=recruits_infs$plot,
                  year=recruits_infs$year_index
  )

## fit stan model
recruit_model <- stan_model("analysis/Stan/ltreb_recruitment.stan")
recruit_fit<-sampling(recruit_model,data = recruit_dat,
                   chains=3,
                   control = list(adapt_delta=0.99,stepsize=0.1),
                   iter=5000,thin=2,
                   pars = c("alpha1","beta1",
                            "alpha2","beta2",
                            "sigma_year","sigma_plot",
                            "y_sim"), 
                   save_warmup=F)
#write_rds(recruit_fit,"analysis/Stan/recruit_fit.rds")
recruit_fit<-read_rds("analysis/Stan/recruit_fit.rds")

## check a few trace plots
bayesplot::mcmc_trace(recruit_fit,pars = c("sigma_year","sigma_plot"))
bayesplot::mcmc_trace(recruit_fit,pars = c("alpha1[4]","beta1[4]",
                                           "alpha2[4]","beta2[4]"))
recruit_sim<-rstan::extract(recruit_fit,"y_sim")
ppc_dens_overlay(recruit_dat$y,recruit_sim$y_sim)

recruit_params<-rstan::extract(recruit_fit,c("alpha1","beta1",
                                             "alpha2","beta2"))
r1_em<-exp(apply(recruit_params$alpha1,2,quantile,probs=quantile_probs))
r1_ep<-exp(apply(recruit_params$alpha1+recruit_params$beta1,2,quantile,probs=quantile_probs))
r2_em<-exp(apply(recruit_params$alpha2,2,quantile,probs=quantile_probs))
r2_ep<-exp(apply(recruit_params$alpha2+recruit_params$beta2,2,quantile,probs=quantile_probs))

##make a plot of parameter estimates
spp_list<-c("A.p.","E.vil.","E.vir.","F.s.","P.al.","P.au.","P.s.")
pdf("manuscript/figures/recruitment.pdf",height = 6, width = 6,useDingbats = F)
plot(1:7,r1_em[3,],type="n",ylim=c(0,max(c(r1_em,r1_ep,r2_em,r2_ep))),xlim=c(0.8,7.2),
     ylab="Recruits / Inflorescence",xlab="Species",axes=F,cex.lab=1.2)
arrows((1:7)-.2,r1_em[2,],
       (1:7)-.2,r1_em[4,],length=0,lwd=3,col="tomato")
arrows((1:7)-.2,r1_em[1,],
       (1:7)-.2,r1_em[5,],length=0,lwd=1,col="tomato")
points((1:7)-.2,r1_em[3,],pch=16,cex=2,col="tomato")

arrows((1:7)-.1,r2_em[2,],
       (1:7)-.1,r2_em[4,],length=0,lwd=3,col="tomato")
arrows((1:7)-.1,r2_em[1,],
       (1:7)-.1,r2_em[5,],length=0,lwd=1,col="tomato")
points((1:7)-.1,r2_em[3,],pch=21,cex=2,col="tomato",bg="white")

arrows((1:7)+.1,r1_ep[2,],
       (1:7)+.1,r1_ep[4,],length=0,lwd=3,col="cornflowerblue")
arrows((1:7)+.1,r1_ep[1,],
       (1:7)+.1,r1_ep[5,],length=0,lwd=1,col="cornflowerblue")
points((1:7)+.1,r1_ep[3,],pch=16,cex=2,col="cornflowerblue")

arrows((1:7)+.2,r2_ep[2,],
       (1:7)+.2,r2_ep[4,],length=0,lwd=3,col="cornflowerblue")
arrows((1:7)+.2,r2_ep[1,],
       (1:7)+.2,r2_ep[5,],length=0,lwd=1,col="cornflowerblue")
points((1:7)+.2,r2_ep[3,],pch=21,cex=2,col="cornflowerblue",bg="white")
box()
axis(1,at=1:7,labels=spp_list,font=3)
axis(2,at=pretty(c(r1_em,r1_ep,r2_em,r2_ep)))
legend("topright",legend=c("E-","E+",expression(Infs[t-1]),expression(Infs[t-2])),ncol=2,
       pch=c(16,16,16,1),col=c("tomato","cornflowerblue","black","black"),cex=1.2)
dev.off()


# Assemble list of matrices from posterior samples ------------------------

##grab random samples of indices
n_post<-100
## the recruit model has fewer posterior samples
set.seed(12291980)
survfert_i<-sample(1:dim(Ap_par$beta_Ap)[1],n_post,replace = F)
recruit_i<-sample(1:dim(recruit_params$alpha1)[1],n_post,replace = F)

## set up matrix dimensions for each species
Ap_dim<-age_limits$lump_age[age_limits$species=="AGPE"]+2
Er_dim<-age_limits$lump_age[age_limits$species=="ELRI"]+2
Ev_dim<-age_limits$lump_age[age_limits$species=="ELVI"]+2
Fs_dim<-age_limits$lump_age[age_limits$species=="FESU"]+2
Pa_dim<-age_limits$lump_age[age_limits$species=="POAL"]+2
Pu_dim<-age_limits$lump_age[age_limits$species=="POAU"]+2
Ps_dim<-age_limits$lump_age[age_limits$species=="POSY"]+2

## source life history function
source("analysis/Lifehistorytrait_function.R")

## set up data frame to store life history metrics
LH_dat<-data.frame(R0_em=rep(NA,n_post),R0_ep=rep(NA,n_post),
                               G_em=rep(NA,n_post),G_ep=rep(NA,n_post),
                               lambda_em=rep(NA,n_post),lambda_ep=rep(NA,n_post),
                               pRep_em=rep(NA,n_post),pRep_ep=rep(NA,n_post),
                               La_em=rep(NA,n_post),La_ep=rep(NA,n_post),
                               matlifexp_em=rep(NA,n_post),matlifexp_ep=rep(NA,n_post),
                               meanelexp_em=rep(NA,n_post),meanelexp_ep=rep(NA,n_post),
                               entropyd_em=rep(NA,n_post),entropyd_ep=rep(NA,n_post),
                               entropyk_em=rep(NA,n_post),entropyk_ep=rep(NA,n_post),
                               gini_em=rep(NA,n_post),gini_ep=rep(NA,n_post))
Ap_lifehistorypost<-Er_lifehistorypost<-Ev_lifehistorypost<-Fs_lifehistorypost<-Pa_lifehistorypost<-Pu_lifehistorypost<-Ps_lifehistorypost<-LH_dat

for(i in 1:n_post){
  ## grab recruitment estimates for this posterior draw -- gives a 1x7 vector
  r1_em<-exp(recruit_params$alpha1[recruit_i[i],])
  r1_ep<-exp(recruit_params$alpha1[recruit_i[i],]+recruit_params$beta1[recruit_i[i],])
  r2_em<-exp(recruit_params$alpha2[recruit_i[i],])
  r2_ep<-exp(recruit_params$alpha2[recruit_i[i],]+recruit_params$beta2[recruit_i[i],])
  
  ## Agrostis perennans
  Ap_em_surv <- invlogit(c(Ap_par$beta_Ap[survfert_i[i],"(Intercept)"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5"]))
  Ap_ep_surv <- invlogit(c(Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                           Ap_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5"]+Ap_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"]))
  Ap_em_fert <- exp(c(Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5"]))
  Ap_ep_fert <- exp(c(Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                      Ap_fert_par$beta_Ap[survfert_i[i],"(Intercept)"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(endo_01)1"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5"]+Ap_fert_par$beta_Ap[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"]))
  ## assemble matrices
  Ap_em_U<-Ap_em_F<-Ap_ep_U<-Ap_ep_F<-matrix(0,Ap_dim,Ap_dim)
  ## inf production goes in top row of F (except first element)
  Ap_em_F[1,2:Ap_dim]<-Ap_em_fert
  ## rest of the second row is recruitment from last year's infs
  Ap_em_F[2,1]<-r2_em[1]
  Ap_em_F[2,2:Ap_dim]<-r1_em[1]*Ap_em_fert
  ##subdiagonals are the inf bank recruitment rate and age-specific survival
  diag(Ap_em_U[-1,-ncol(Ap_em_U)])<-c(0,Ap_em_surv[1:5])
  ## survival of the terminal age class goes in the bottom right corner
  Ap_em_U[Ap_dim,Ap_dim]<-Ap_em_surv[6]
  ## Ap_eplus -- follows same structure
  Ap_ep_F[1,2:Ap_dim]<-Ap_ep_fert
  Ap_ep_F[2,1]<-r2_ep[1]
  Ap_ep_F[2,2:Ap_dim]<-r1_ep[1]*Ap_ep_fert
  diag(Ap_ep_U[-1,-ncol(Ap_ep_U)])<-c(0,Ap_ep_surv[1:5])
  Ap_ep_U[Ap_dim,Ap_dim]<-Ap_ep_surv[6]
  Ap_em<-lifeTimeRepEvents(Ap_em_U,Ap_em_F,2)
  Ap_ep<-lifeTimeRepEvents(Ap_ep_U,Ap_ep_F,2)
  ## store outputs
  Ap_lifehistorypost$R0_em[i]<-Ap_em$Ro
  Ap_lifehistorypost$G_em[i]<-Ap_em$generation.time
  Ap_lifehistorypost$lambda_em[i]<-Ap_em$lambda
  Ap_lifehistorypost$pRep_em[i]<-Ap_em$pRep
  Ap_lifehistorypost$La_em[i]<-Ap_em$La
  Ap_lifehistorypost$matlifexp_em[i]<-Ap_em$remainingMatureLifeExpectancy
  Ap_lifehistorypost$meanelexp_em[i]<-Ap_em$meanelexp
  Ap_lifehistorypost$entropyd_em[i]<-Ap_em$entropyd
  Ap_lifehistorypost$entropyk_em[i]<-Ap_em$entropyk
  Ap_lifehistorypost$gini_em[i]<-Ap_em$gini
  Ap_lifehistorypost$longevity_em[i]<-Ap_em$longevity
  
  Ap_lifehistorypost$R0_ep[i]<-Ap_ep$Ro
  Ap_lifehistorypost$G_ep[i]<-Ap_ep$generation.time
  Ap_lifehistorypost$lambda_ep[i]<-Ap_ep$lambda
  Ap_lifehistorypost$pRep_ep[i]<-Ap_ep$pRep
  Ap_lifehistorypost$La_ep[i]<-Ap_ep$La
  Ap_lifehistorypost$matlifexp_ep[i]<-Ap_ep$remainingMatureLifeExpectancy
  Ap_lifehistorypost$meanelexp_ep[i]<-Ap_ep$meanelexp
  Ap_lifehistorypost$entropyd_ep[i]<-Ap_ep$entropyd
  Ap_lifehistorypost$entropyk_ep[i]<-Ap_ep$entropyk
  Ap_lifehistorypost$gini_ep[i]<-Ap_ep$gini
  Ap_lifehistorypost$longevity_ep[i]<-Ap_ep$longevity
  
  ## Elymus villosus
  Er_em_surv <- invlogit(c(Er_par$beta_Er[survfert_i[i],"(Intercept)"],
                           Er_par$beta_Er[survfert_i[i],"(Intercept)"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)1"],
                           Er_par$beta_Er[survfert_i[i],"(Intercept)"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)2"]))
  Er_ep_surv <- invlogit(c(Er_par$beta_Er[survfert_i[i],"(Intercept)"]+Er_par$beta_Er[survfert_i[i],"as.factor(endo_01)1"],
                           Er_par$beta_Er[survfert_i[i],"(Intercept)"]+Er_par$beta_Er[survfert_i[i],"as.factor(endo_01)1"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)1"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Er_par$beta_Er[survfert_i[i],"(Intercept)"]+Er_par$beta_Er[survfert_i[i],"as.factor(endo_01)1"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)2"]+Er_par$beta_Er[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"]))
  Er_em_fert <- exp(c(Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"],
                      Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)1"],
                      Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)2"]))
  Er_ep_fert <- exp(c(Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(endo_01)1"],
                      Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(endo_01)1"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)1"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Er_par_fert$beta_Er[survfert_i[i],"(Intercept)"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(endo_01)1"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)2"]+Er_par_fert$beta_Er[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"]))
  Er_em_U<-Er_em_F<-Er_ep_U<-Er_ep_F<-matrix(0,Er_dim,Er_dim)
  Er_em_F[1,2:Er_dim]<-Er_em_fert
  Er_em_F[2,1]<-r2_em[2]
  Er_em_F[2,2:Er_dim]<-r1_em[2]*Er_em_fert
  diag(Er_em_U[-1,-ncol(Er_em_U)])<-c(0,Er_em_surv[1:2])
  Er_em_U[Er_dim,Er_dim]<-Er_em_surv[3]
  Er_ep_F[1,2:Er_dim]<-Er_ep_fert
  Er_ep_F[2,1]<-r2_ep[2]
  Er_ep_F[2,2:Er_dim]<-r1_ep[2]*Er_ep_fert
  diag(Er_ep_U[-1,-ncol(Er_ep_U)])<-c(0,Er_ep_surv[1:2])
  Er_ep_U[Er_dim,Er_dim]<-Er_ep_surv[3]
  Er_em<-lifeTimeRepEvents(Er_em_U,Er_em_F,2)
  Er_ep<-lifeTimeRepEvents(Er_ep_U,Er_ep_F,2)
  ## store outputs
  Er_lifehistorypost$R0_em[i]<-Er_em$Ro
  Er_lifehistorypost$G_em[i]<-Er_em$generation.time
  Er_lifehistorypost$lambda_em[i]<-Er_em$lambda
  Er_lifehistorypost$pRep_em[i]<-Er_em$pRep
  Er_lifehistorypost$La_em[i]<-Er_em$La
  Er_lifehistorypost$matlifexp_em[i]<-Er_em$remainingMatureLifeExpectancy
  Er_lifehistorypost$meanelexp_em[i]<-Er_em$meanelexp
  Er_lifehistorypost$entropyd_em[i]<-Er_em$entropyd
  Er_lifehistorypost$entropyk_em[i]<-Er_em$entropyk
  Er_lifehistorypost$gini_em[i]<-Er_em$gini
  Er_lifehistorypost$longevity_em[i]<-Er_em$longevity
  
  Er_lifehistorypost$R0_ep[i]<-Er_ep$Ro
  Er_lifehistorypost$G_ep[i]<-Er_ep$generation.time
  Er_lifehistorypost$lambda_ep[i]<-Er_ep$lambda
  Er_lifehistorypost$pRep_ep[i]<-Er_ep$pRep
  Er_lifehistorypost$La_ep[i]<-Er_ep$La
  Er_lifehistorypost$matlifexp_ep[i]<-Er_ep$remainingMatureLifeExpectancy
  Er_lifehistorypost$meanelexp_ep[i]<-Er_ep$meanelexp
  Er_lifehistorypost$entropyd_ep[i]<-Er_ep$entropyd
  Er_lifehistorypost$entropyk_ep[i]<-Er_ep$entropyk
  Er_lifehistorypost$gini_ep[i]<-Er_ep$gini
  Er_lifehistorypost$longevity_ep[i]<-Er_ep$longevity
  
  ##Elymus virginicus
  Ev_em_surv <- invlogit(c(Ev_par$beta_Ev[survfert_i[i],"(Intercept)"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)1"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)2"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)3"]))
  Ev_ep_surv <- invlogit(c(Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(endo_01)1"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)1"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)2"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                           Ev_par$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)3"]+Ev_par$beta_Ev[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"]))
  Ev_em_fert <- exp(c(Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)1"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)2"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)3"]))
  Ev_ep_fert <- exp(c(Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(endo_01)1"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)1"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)2"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                      Ev_par_fert$beta_Ev[survfert_i[i],"(Intercept)"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(endo_01)1"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)3"]+Ev_par_fert$beta_Ev[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"]))
  Ev_em_U<-Ev_em_F<-Ev_ep_U<-Ev_ep_F<-matrix(0,Ev_dim,Ev_dim)
  Ev_em_F[1,2:Ev_dim]<-Ev_em_fert
  Ev_em_F[2,1]<-r2_em[3]
  Ev_em_F[2,2:Ev_dim]<-r1_em[3]*Ev_em_fert
  diag(Ev_em_U[-1,-ncol(Ev_em_U)])<-c(0,Ev_em_surv[1:3])
  Ev_em_U[Ev_dim,Ev_dim]<-Ev_em_surv[4]
  Ev_ep_F[1,2:Ev_dim]<-Ev_ep_fert
  Ev_ep_F[2,1]<-r2_ep[3]
  Ev_ep_F[2,2:Ev_dim]<-r1_ep[3]*Ev_ep_fert
  diag(Ev_ep_U[-1,-ncol(Ev_ep_U)])<-c(0,Ev_ep_surv[1:3])
  Ev_ep_U[Ev_dim,Ev_dim]<-Ev_ep_surv[4]
  Ev_em<-lifeTimeRepEvents(Ev_em_U,Ev_em_F,2)
  Ev_ep<-lifeTimeRepEvents(Ev_ep_U,Ev_ep_F,2)
  ## store outputs
  Ev_lifehistorypost$R0_em[i]<-Ev_em$Ro
  Ev_lifehistorypost$G_em[i]<-Ev_em$generation.time
  Ev_lifehistorypost$lambda_em[i]<-Ev_em$lambda
  Ev_lifehistorypost$pRep_em[i]<-Ev_em$pRep
  Ev_lifehistorypost$La_em[i]<-Ev_em$La
  Ev_lifehistorypost$matlifexp_em[i]<-Ev_em$remainingMatureLifeExpectancy
  Ev_lifehistorypost$meanelexp_em[i]<-Ev_em$meanelexp
  Ev_lifehistorypost$entropyd_em[i]<-Ev_em$entropyd
  Ev_lifehistorypost$entropyk_em[i]<-Ev_em$entropyk
  Ev_lifehistorypost$gini_em[i]<-Ev_em$gini
  Ev_lifehistorypost$longevity_em[i]<-Ev_em$longevity
  
  Ev_lifehistorypost$R0_ep[i]<-Ev_ep$Ro
  Ev_lifehistorypost$G_ep[i]<-Ev_ep$generation.time
  Ev_lifehistorypost$lambda_ep[i]<-Ev_ep$lambda
  Ev_lifehistorypost$pRep_ep[i]<-Ev_ep$pRep
  Ev_lifehistorypost$La_ep[i]<-Ev_ep$La
  Ev_lifehistorypost$matlifexp_ep[i]<-Ev_ep$remainingMatureLifeExpectancy
  Ev_lifehistorypost$meanelexp_ep[i]<-Ev_ep$meanelexp
  Ev_lifehistorypost$entropyd_ep[i]<-Ev_ep$entropyd
  Ev_lifehistorypost$entropyk_ep[i]<-Ev_ep$entropyk
  Ev_lifehistorypost$gini_ep[i]<-Ev_ep$gini
  Ev_lifehistorypost$longevity_ep[i]<-Ev_ep$longevity
  
  ## Festuca subverticillata
  Fs_em_surv <- invlogit(c(Fs_par$beta_Fs[survfert_i[i],"(Intercept)"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)2"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)3"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)4"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)5"]))
  Fs_ep_surv <- invlogit(c(Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)2"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)3"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)4"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                           Fs_par$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)5"]+Fs_par$beta_Fs[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"]))
  Fs_em_fert <- exp(c(Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)2"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)3"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)4"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)5"]))
  Fs_ep_fert <- exp(c(Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)2"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)3"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)4"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                      Fs_par_fert$beta_Fs[survfert_i[i],"(Intercept)"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(endo_01)1"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)5"]+Fs_par_fert$beta_Fs[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"]))
  Fs_em_U<-Fs_em_F<-Fs_ep_U<-Fs_ep_F<-matrix(0,Fs_dim,Fs_dim)
  Fs_em_F[1,2:Fs_dim]<-Fs_em_fert
  Fs_em_F[2,1]<-r2_em[4]
  Fs_em_F[2,2:Fs_dim]<-r1_em[4]*Fs_em_fert
  diag(Fs_em_U[-1,-ncol(Fs_em_U)])<-c(0,Fs_em_surv[1:5])
  Fs_em_U[Fs_dim,Fs_dim]<-Fs_em_surv[6]
  Fs_ep_F[1,2:Fs_dim]<-Fs_ep_fert
  Fs_ep_F[2,1]<-r2_ep[4]
  Fs_ep_F[2,2:Fs_dim]<-r1_ep[4]*Fs_ep_fert
  diag(Fs_ep_U[-1,-ncol(Fs_ep_U)])<-c(0,Fs_ep_surv[1:5])
  Fs_ep_U[Fs_dim,Fs_dim]<-Fs_ep_surv[6]
  Fs_em<-lifeTimeRepEvents(Fs_em_U,Fs_em_F,2)
  Fs_ep<-lifeTimeRepEvents(Fs_ep_U,Fs_ep_F,2)
  ##store outputs
  Fs_lifehistorypost$R0_em[i]<-Fs_em$Ro
  Fs_lifehistorypost$G_em[i]<-Fs_em$generation.time
  Fs_lifehistorypost$lambda_em[i]<-Fs_em$lambda
  Fs_lifehistorypost$pRep_em[i]<-Fs_em$pRep
  Fs_lifehistorypost$La_em[i]<-Fs_em$La
  Fs_lifehistorypost$matlifexp_em[i]<-Fs_em$remainingMatureLifeExpectancy
  Fs_lifehistorypost$meanelexp_em[i]<-Fs_em$meanelexp
  Fs_lifehistorypost$entropyd_em[i]<-Fs_em$entropyd
  Fs_lifehistorypost$entropyk_em[i]<-Fs_em$entropyk
  Fs_lifehistorypost$gini_em[i]<-Fs_em$gini
  Fs_lifehistorypost$longevity_em[i]<-Fs_em$longevity
  
  Fs_lifehistorypost$R0_ep[i]<-Fs_ep$Ro
  Fs_lifehistorypost$G_ep[i]<-Fs_ep$generation.time
  Fs_lifehistorypost$lambda_ep[i]<-Fs_ep$lambda
  Fs_lifehistorypost$pRep_ep[i]<-Fs_ep$pRep
  Fs_lifehistorypost$La_ep[i]<-Fs_ep$La
  Fs_lifehistorypost$matlifexp_ep[i]<-Fs_ep$remainingMatureLifeExpectancy
  Fs_lifehistorypost$meanelexp_ep[i]<-Fs_ep$meanelexp
  Fs_lifehistorypost$entropyd_ep[i]<-Fs_ep$entropyd
  Fs_lifehistorypost$entropyk_ep[i]<-Fs_ep$entropyk
  Fs_lifehistorypost$gini_ep[i]<-Fs_ep$gini
  Fs_lifehistorypost$longevity_ep[i]<-Fs_ep$longevity
  
  
  ##Poa alsodes
  Pa_em_surv <- invlogit(c(Pa_par$beta_Pa[survfert_i[i],"(Intercept)"],
                           Pa_par$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)1"],
                           Pa_par$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)2"]))
  Pa_ep_surv <- invlogit(c(Pa_par$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(endo_01)1"],
                           Pa_par$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(endo_01)1"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)1"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Pa_par$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(endo_01)1"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)2"]+Pa_par$beta_Pa[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"]))
  Pa_em_fert <- exp(c(Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"],
                      Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)1"],
                      Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)2"]))
  Pa_ep_fert <- exp(c(Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(endo_01)1"],
                      Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)1"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Pa_par_fert$beta_Pa[survfert_i[i],"(Intercept)"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(endo_01)1"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)2"]+Pa_par_fert$beta_Pa[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"]))
  Pa_em_U<-Pa_em_F<-Pa_ep_U<-Pa_ep_F<-matrix(0,Pa_dim,Pa_dim)
  Pa_em_F[1,2:Pa_dim]<-Pa_em_fert
  Pa_em_F[2,1]<-r2_em[5]
  Pa_em_F[2,2:Pa_dim]<-r1_em[5]*Pa_em_fert
  diag(Pa_em_U[-1,-ncol(Pa_em_U)])<-c(0,Pa_em_surv[1:2])
  Pa_em_U[Pa_dim,Pa_dim]<-Pa_em_surv[3]
  Pa_ep_F[1,2:Pa_dim]<-Pa_ep_fert
  Pa_ep_F[2,1]<-r2_ep[5]
  Pa_ep_F[2,2:Pa_dim]<-r1_ep[5]*Pa_ep_fert
  diag(Pa_ep_U[-1,-ncol(Pa_ep_U)])<-c(0,Pa_ep_surv[1:2])
  Pa_ep_U[Pa_dim,Pa_dim]<-Pa_ep_surv[3]
  Pa_em<-lifeTimeRepEvents(Pa_em_U,Pa_em_F,2)
  Pa_ep<-lifeTimeRepEvents(Pa_ep_U,Pa_ep_F,2)
  ##store outputs
  Pa_lifehistorypost$R0_em[i]<-Pa_em$Ro
  Pa_lifehistorypost$G_em[i]<-Pa_em$generation.time
  Pa_lifehistorypost$lambda_em[i]<-Pa_em$lambda
  Pa_lifehistorypost$pRep_em[i]<-Pa_em$pRep
  Pa_lifehistorypost$La_em[i]<-Pa_em$La
  Pa_lifehistorypost$matlifexp_em[i]<-Pa_em$remainingMatureLifeExpectancy
  Pa_lifehistorypost$meanelexp_em[i]<-Pa_em$meanelexp
  Pa_lifehistorypost$entropyd_em[i]<-Pa_em$entropyd
  Pa_lifehistorypost$entropyk_em[i]<-Pa_em$entropyk
  Pa_lifehistorypost$gini_em[i]<-Pa_em$gini
  Pa_lifehistorypost$longevity_em[i]<-Pa_em$longevity
  
  Pa_lifehistorypost$R0_ep[i]<-Pa_ep$Ro
  Pa_lifehistorypost$G_ep[i]<-Pa_ep$generation.time
  Pa_lifehistorypost$lambda_ep[i]<-Pa_ep$lambda
  Pa_lifehistorypost$pRep_ep[i]<-Pa_ep$pRep
  Pa_lifehistorypost$La_ep[i]<-Pa_ep$La
  Pa_lifehistorypost$matlifexp_ep[i]<-Pa_ep$remainingMatureLifeExpectancy
  Pa_lifehistorypost$meanelexp_ep[i]<-Pa_ep$meanelexp
  Pa_lifehistorypost$entropyd_ep[i]<-Pa_ep$entropyd
  Pa_lifehistorypost$entropyk_ep[i]<-Pa_ep$entropyk
  Pa_lifehistorypost$gini_ep[i]<-Pa_ep$gini
  Pa_lifehistorypost$longevity_ep[i]<-Pa_ep$longevity
  
  ##Poa autumnalis
  Pu_em_surv <- invlogit(c(Pu_par$beta_Pu[survfert_i[i],"(Intercept)"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)1"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)2"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)3"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)4"]))
  Pu_ep_surv <- invlogit(c(Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(endo_01)1"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)1"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)2"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)3"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                           Pu_par$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)4"]+Pu_par$beta_Pu[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"]))
  Pu_em_fert <- exp(c(Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)1"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)2"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)3"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)4"]))
  Pu_ep_fert <- exp(c(Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(endo_01)1"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)1"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)2"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)3"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                      Pu_par_fert$beta_Pu[survfert_i[i],"(Intercept)"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(endo_01)1"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)4"]+Pu_par_fert$beta_Pu[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"]))
  Pu_em_U<-Pu_em_F<-Pu_ep_U<-Pu_ep_F<-matrix(0,Pu_dim,Pu_dim)
  Pu_em_F[1,2:Pu_dim]<-Pu_em_fert
  Pu_em_F[2,1]<-r2_em[6]
  Pu_em_F[2,2:Pu_dim]<-r1_em[6]*Pu_em_fert
  diag(Pu_em_U[-1,-ncol(Pu_em_U)])<-c(0,Pu_em_surv[1:4])
  Pu_em_U[Pu_dim,Pu_dim]<-Pu_em_surv[5]
  Pu_ep_F[1,2:Pu_dim]<-Pu_ep_fert
  Pu_ep_F[2,1]<-r2_ep[6]
  Pu_ep_F[2,2:Pu_dim]<-r1_ep[6]*Pu_ep_fert
  diag(Pu_ep_U[-1,-ncol(Pu_ep_U)])<-c(0,Pu_ep_surv[1:4])
  Pu_ep_U[Pu_dim,Pu_dim]<-Pu_ep_surv[5]
  ##store outputs
  Pu_em<-lifeTimeRepEvents(Pu_em_U,Pu_em_F,2)
  Pu_ep<-lifeTimeRepEvents(Pu_ep_U,Pu_ep_F,2)
  Pu_lifehistorypost$R0_em[i]<-Pu_em$Ro
  Pu_lifehistorypost$G_em[i]<-Pu_em$generation.time
  Pu_lifehistorypost$lambda_em[i]<-Pu_em$lambda
  Pu_lifehistorypost$pRep_em[i]<-Pu_em$pRep
  Pu_lifehistorypost$La_em[i]<-Pu_em$La
  Pu_lifehistorypost$matlifexp_em[i]<-Pu_em$remainingMatureLifeExpectancy
  Pu_lifehistorypost$meanelexp_em[i]<-Pu_em$meanelexp
  Pu_lifehistorypost$entropyd_em[i]<-Pu_em$entropyd
  Pu_lifehistorypost$entropyk_em[i]<-Pu_em$entropyk
  Pu_lifehistorypost$gini_em[i]<-Pu_em$gini
  Pu_lifehistorypost$longevity_em[i]<-Pu_em$longevity
  
  Pu_lifehistorypost$R0_ep[i]<-Pu_ep$Ro
  Pu_lifehistorypost$G_ep[i]<-Pu_ep$generation.time
  Pu_lifehistorypost$lambda_ep[i]<-Pu_ep$lambda
  Pu_lifehistorypost$pRep_ep[i]<-Pu_ep$pRep
  Pu_lifehistorypost$La_ep[i]<-Pu_ep$La
  Pu_lifehistorypost$matlifexp_ep[i]<-Pu_ep$remainingMatureLifeExpectancy
  Pu_lifehistorypost$meanelexp_ep[i]<-Pu_ep$meanelexp
  Pu_lifehistorypost$entropyd_ep[i]<-Pu_ep$entropyd
  Pu_lifehistorypost$entropyk_ep[i]<-Pu_ep$entropyk
  Pu_lifehistorypost$gini_ep[i]<-Pu_ep$gini
  Pu_lifehistorypost$longevity_ep[i]<-Pu_ep$longevity
  
  ##Poa sylvestris
  Ps_em_surv <- invlogit(c(Ps_par$beta_Ps[survfert_i[i],"(Intercept)"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)2"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)3"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)4"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)5"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)6"]))
  Ps_ep_surv <- invlogit(c(Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)2"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)3"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)4"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)5"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"],
                           Ps_par$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)6"]+Ps_par$beta_Ps[survfert_i[i],"as.factor(age_lump)6:as.factor(endo_01)1"]))
  Ps_em_fert <- exp(c(Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)2"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)3"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)4"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)5"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)6"]))
  Ps_ep_fert <- exp(c(Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)1:as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)2"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)2:as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)3"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)3:as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)4"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)4:as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)5"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)5:as.factor(endo_01)1"],
                      Ps_par_fert$beta_Ps[survfert_i[i],"(Intercept)"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(endo_01)1"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)6"]+Ps_par_fert$beta_Ps[survfert_i[i],"as.factor(age_lump)6:as.factor(endo_01)1"]))
  Ps_em_U<-Ps_em_F<-Ps_ep_U<-Ps_ep_F<-matrix(0,Ps_dim,Ps_dim)
  Ps_em_F[1,2:Ps_dim]<-Ps_em_fert
  Ps_em_F[2,1]<-r2_em[7]
  Ps_em_F[2,2:Ps_dim]<-r1_em[7]*Ps_em_fert
  diag(Ps_em_U[-1,-ncol(Ps_em_U)])<-c(0,Ps_em_surv[1:6])
  Ps_em_U[Ps_dim,Ps_dim]<-Ps_em_surv[7]
  Ps_ep_F[1,2:Ps_dim]<-Ps_ep_fert
  Ps_ep_F[2,1]<-r2_ep[7]
  Ps_ep_F[2,2:Ps_dim]<-r1_ep[7]*Ps_ep_fert
  diag(Ps_ep_U[-1,-ncol(Ps_ep_U)])<-c(0,Ps_ep_surv[1:6])
  Ps_ep_U[Ps_dim,Ps_dim]<-Ps_ep_surv[7]
  Ps_em<-lifeTimeRepEvents(Ps_em_U,Ps_em_F,2)
  Ps_ep<-lifeTimeRepEvents(Ps_ep_U,Ps_ep_F,2)
  ##store outputs
  Ps_lifehistorypost$R0_em[i]<-Ps_em$Ro
  Ps_lifehistorypost$G_em[i]<-Ps_em$generation.time
  Ps_lifehistorypost$lambda_em[i]<-Ps_em$lambda
  Ps_lifehistorypost$pRep_em[i]<-Ps_em$pRep
  Ps_lifehistorypost$La_em[i]<-Ps_em$La
  Ps_lifehistorypost$matlifexp_em[i]<-Ps_em$remainingMatureLifeExpectancy
  Ps_lifehistorypost$meanelexp_em[i]<-Ps_em$meanelexp
  Ps_lifehistorypost$entropyd_em[i]<-Ps_em$entropyd
  Ps_lifehistorypost$entropyk_em[i]<-Ps_em$entropyk
  Ps_lifehistorypost$gini_em[i]<-Ps_em$gini
  Ps_lifehistorypost$longevity_em[i]<-Ps_em$longevity
  
  Ps_lifehistorypost$R0_ep[i]<-Ps_ep$Ro
  Ps_lifehistorypost$G_ep[i]<-Ps_ep$generation.time
  Ps_lifehistorypost$lambda_ep[i]<-Ps_ep$lambda
  Ps_lifehistorypost$pRep_ep[i]<-Ps_ep$pRep
  Ps_lifehistorypost$La_ep[i]<-Ps_ep$La
  Ps_lifehistorypost$matlifexp_ep[i]<-Ps_ep$remainingMatureLifeExpectancy
  Ps_lifehistorypost$meanelexp_ep[i]<-Ps_ep$meanelexp
  Ps_lifehistorypost$entropyd_ep[i]<-Ps_ep$entropyd
  Ps_lifehistorypost$entropyk_ep[i]<-Ps_ep$entropyk
  Ps_lifehistorypost$gini_ep[i]<-Ps_ep$gini
  Ps_lifehistorypost$longevity_ep[i]<-Ps_ep$longevity
  
  
}##end posterior sample loop

## combine species in to single data frame
Ap_lifehistorypost$species<-"AGPE"
Er_lifehistorypost$species<-"ELRI"
Ev_lifehistorypost$species<-"ELVI"
Fs_lifehistorypost$species<-"FESU"
Pa_lifehistorypost$species<-"POAL"
Pu_lifehistorypost$species<-"POAU"
Ps_lifehistorypost$species<-"POSY"

lifehistorypost<-bind_rows(Ap_lifehistorypost,
                           Er_lifehistorypost,
                           Ev_lifehistorypost,
                           Fs_lifehistorypost,
                           Pa_lifehistorypost,
                           Pu_lifehistorypost,
                           Ps_lifehistorypost) %>% 
  mutate(species = factor(species))

write.csv(lifehistorypost,"analysis/lifehistorypost.csv")


# Assemble matrices from posterior medians --------------------------------------------------
##Ap eminus
Ap_dim<-age_limits$lump_age[age_limits$species=="AGPE"]+2
Ap_em_mat<-matrix(0,Ap_dim,Ap_dim)
## inf production goes in top row (except first element)
Ap_em_mat[1,2:Ap_dim]<-Ap_em_fert[3,]
## rest of the second row is recruitment from last year's infs
Ap_em_mat[2,2:Ap_dim]<-r1_em[3,1]*Ap_em_fert[3,]
##subdiagonals are the inf bank recruitment rate and age-specific survival
diag(Ap_em_mat[-1,-ncol(Ap_em_mat)])<-c(r2_em[3,1],Ap_em_surv[3,1:6])
## survival of the terminal age class goes in the bottom right corner
Ap_em_mat[Ap_dim,Ap_dim]<-Ap_em_surv[3,7]
## Ap_eplus
Ap_ep_mat<-matrix(0,Ap_dim,Ap_dim)
Ap_ep_mat[1,2:Ap_dim]<-Ap_ep_fert[3,]
Ap_ep_mat[2,2:Ap_dim]<-r1_ep[3,1]*Ap_ep_fert[3,]
diag(Ap_ep_mat[-1,-ncol(Ap_ep_mat)])<-c(r2_ep[3,1],Ap_ep_surv[3,1:6])
Ap_ep_mat[Ap_dim,Ap_dim]<-Ap_ep_surv[3,7]

#Er eminus
Er_dim<-age_limits$lump_age[age_limits$species=="ELRI"]+2
Er_em_mat<-matrix(0,Er_dim,Er_dim)
Er_em_mat[1,2:Er_dim]<-Er_em_fert[3,]
Er_em_mat[2,2:Er_dim]<-r1_em[3,2]*Er_em_fert[3,]
diag(Er_em_mat[-1,-ncol(Er_em_mat)])<-c(r2_em[3,2],Er_em_surv[3,1:4])
Er_em_mat[Er_dim,Er_dim]<-Er_em_surv[3,5]
#Er eplus
Er_ep_mat<-matrix(0,Er_dim,Er_dim)
Er_ep_mat[1,2:Er_dim]<-Er_ep_fert[3,]
Er_ep_mat[2,2:Er_dim]<-r1_ep[3,2]*Er_ep_fert[3,]
diag(Er_ep_mat[-1,-ncol(Er_ep_mat)])<-c(r2_ep[3,2],Er_ep_surv[3,1:4])
Er_ep_mat[Er_dim,Er_dim]<-Er_ep_surv[3,5]

#Ev eminus
Ev_dim<-age_limits$lump_age[age_limits$species=="ELVI"]+2
Ev_em_mat<-matrix(0,Ev_dim,Ev_dim)
Ev_em_mat[1,2:Ev_dim]<-Ev_em_fert[3,]
Ev_em_mat[2,2:Ev_dim]<-r1_em[3,3]*Ev_em_fert[3,]
diag(Ev_em_mat[-1,-ncol(Ev_em_mat)])<-c(r2_em[3,3],Ev_em_surv[3,1:4])
Ev_em_mat[Ev_dim,Ev_dim]<-Ev_em_surv[3,5]
#Ev eplus
Ev_ep_mat<-matrix(0,Ev_dim,Ev_dim)
Ev_ep_mat[1,2:Ev_dim]<-Ev_ep_fert[3,]
Ev_ep_mat[2,2:Ev_dim]<-r1_ep[3,3]*Ev_ep_fert[3,]
diag(Ev_ep_mat[-1,-ncol(Ev_ep_mat)])<-c(r2_ep[3,3],Ev_ep_surv[3,1:4])
Ev_ep_mat[Ev_dim,Ev_dim]<-Ev_ep_surv[3,5]

#Fs eminus
Fs_dim<-age_limits$lump_age[age_limits$species=="FESU"]+2
Fs_em_mat<-matrix(0,Fs_dim,Fs_dim)
Fs_em_mat[1,2:Fs_dim]<-Fs_em_fert[3,]
Fs_em_mat[2,2:Fs_dim]<-r1_em[3,4]*Fs_em_fert[3,]
diag(Fs_em_mat[-1,-ncol(Fs_em_mat)])<-c(r2_em[3,4],Fs_em_surv[3,1:6])
Fs_em_mat[Fs_dim,Fs_dim]<-Fs_em_surv[3,7]
#Fs eplus
Fs_ep_mat<-matrix(0,Fs_dim,Fs_dim)
Fs_ep_mat[1,2:Fs_dim]<-Fs_ep_fert[3,]
Fs_ep_mat[2,2:Fs_dim]<-r1_ep[3,4]*Fs_ep_fert[3,]
diag(Fs_ep_mat[-1,-ncol(Fs_ep_mat)])<-c(r2_ep[3,4],Fs_ep_surv[3,1:6])
Fs_ep_mat[Fs_dim,Fs_dim]<-Fs_ep_surv[3,7]

#Pa eminus
Pa_dim<-age_limits$lump_age[age_limits$species=="POAL"]+2
Pa_em_mat<-matrix(0,Pa_dim,Pa_dim)
Pa_em_mat[1,2:Pa_dim]<-Pa_em_fert[3,]
Pa_em_mat[2,2:Pa_dim]<-r1_em[3,5]*Pa_em_fert[3,]
diag(Pa_em_mat[-1,-ncol(Pa_em_mat)])<-c(r2_em[3,5],Pa_em_surv[3,1:2])
Pa_em_mat[Pa_dim,Pa_dim]<-Pa_em_surv[3,3]
#Pa eplus
Pa_ep_mat<-matrix(0,Pa_dim,Pa_dim)
Pa_ep_mat[1,2:Pa_dim]<-Pa_ep_fert[3,]
Pa_ep_mat[2,2:Pa_dim]<-r1_ep[3,5]*Pa_ep_fert[3,]
diag(Pa_ep_mat[-1,-ncol(Pa_ep_mat)])<-c(r2_ep[3,5],Pa_ep_surv[3,1:2])
Pa_ep_mat[Pa_dim,Pa_dim]<-Pa_ep_surv[3,3]

#Pu eminus
Pu_dim<-age_limits$lump_age[age_limits$species=="POAU"]+2
Pu_em_mat<-matrix(0,Pu_dim,Pu_dim)
Pu_em_mat[1,2:Pu_dim]<-Pu_em_fert[3,]
Pu_em_mat[2,2:Pu_dim]<-r1_em[3,6]*Pu_em_fert[3,]
diag(Pu_em_mat[-1,-ncol(Pu_em_mat)])<-c(r2_em[3,6],Pu_em_surv[3,1:5])
Pu_em_mat[Pu_dim,Pu_dim]<-Pu_em_surv[3,6]
#Pu eplus
Pu_ep_mat<-matrix(0,Pu_dim,Pu_dim)
Pu_ep_mat[1,2:Pu_dim]<-Pu_ep_fert[3,]
Pu_ep_mat[2,2:Pu_dim]<-r1_ep[3,6]*Pu_ep_fert[3,]
diag(Pu_ep_mat[-1,-ncol(Pu_ep_mat)])<-c(r2_ep[3,6],Pu_ep_surv[3,1:5])
Pu_ep_mat[Pu_dim,Pu_dim]<-Pu_ep_surv[3,6]

#Ps eminus
Ps_dim<-age_limits$lump_age[age_limits$species=="POSY"]+2
Ps_em_mat<-matrix(0,Ps_dim,Ps_dim)
Ps_em_mat[1,2:Ps_dim]<-Ps_em_fert[3,]
Ps_em_mat[2,2:Ps_dim]<-r1_em[3,7]*Ps_em_fert[3,]
diag(Ps_em_mat[-1,-ncol(Ps_em_mat)])<-c(r2_em[3,7],Ps_em_surv[3,1:7])
Ps_em_mat[Ps_dim,Ps_dim]<-Ps_em_surv[3,8]
#Ps eplus
Ps_ep_mat<-matrix(0,Ps_dim,Ps_dim)
Ps_ep_mat[1,2:Ps_dim]<-Ps_ep_fert[3,]
Ps_ep_mat[2,2:Ps_dim]<-r1_ep[3,7]*Ps_ep_fert[3,]
diag(Ps_ep_mat[-1,-ncol(Ps_ep_mat)])<-c(r2_ep[3,7],Ps_ep_surv[3,1:7])
Ps_ep_mat[Ps_dim,Ps_dim]<-Ps_ep_surv[3,8]

## bundle matrices into list
ltreb_matrix_list<-list(list(Ap_em_mat,Ap_ep_mat),
                        list(Er_em_mat,Er_ep_mat),
                        list(Ev_em_mat,Ev_ep_mat),
                        list(Fs_em_mat,Fs_ep_mat),
                        list(Pa_em_mat,Pa_ep_mat),
                        list(Pu_em_mat,Pu_ep_mat),
                        list(Ps_em_mat,Ps_ep_mat))
names(ltreb_matrix_list)<-spp_list
names(ltreb_matrix_list$A.p.)<-names(ltreb_matrix_list$F.s.)<-names(ltreb_matrix_list$E.vil.)<-names(ltreb_matrix_list$E.vir.)<-names(ltreb_matrix_list$P.al.)<-names(ltreb_matrix_list$P.au.)<-names(ltreb_matrix_list$P.s.)<-c("Em","Ep")
lapply(ltreb_matrix_list$A.p., lambda)
lapply(ltreb_matrix_list$F.s., lambda)
lapply(ltreb_matrix_list$E.vil., lambda)
lapply(ltreb_matrix_list$E.vir., lambda)
lapply(ltreb_matrix_list$P.al., lambda)
lapply(ltreb_matrix_list$P.au., lambda)
lapply(ltreb_matrix_list$P.s., lambda)


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

