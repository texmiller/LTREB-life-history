## Purpose: fit models for age- and endo-specific survival and fertility
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(scales)

## functions
invlogit<-function(x){exp(x)/(1+exp(x))}

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
X_Ap<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ap_surv)

Er_surv <- surv_data %>% filter(species=="ELRI") %>% droplevels()
X_Er<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Er_surv)

Ev_surv <- surv_data %>% filter(species=="ELVI") %>% droplevels()
X_Ev<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ev_surv)

Fs_surv <- surv_data %>% filter(species=="FESU") %>% droplevels()
X_Fs<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Fs_surv)

Pa_surv <- surv_data %>% filter(species=="POAL") %>% droplevels()
X_Pa<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pa_surv)

Pu_surv <- surv_data %>% filter(species=="POAU") %>% droplevels()
X_Pu<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Pu_surv)

Ps_surv <- surv_data %>% filter(species=="POSY") %>% droplevels()
X_Ps<-model.matrix(~as.factor(age_lump) * as.factor(endo_01),data=Ps_surv)

stan_dat_surv <- list(n_spp=7,
                 n_years=max(surv_data$year_index),
                 n_plots=max(surv_data$plot),
                 ##AGPE
                 y_Ap=Ap_surv$surv_t1, 
                 n_Ap=length(Ap_surv$surv_t1),
                 beta_Ap_dim=ncol(X_Ap),
                 X_Ap=X_Ap,
                 year_Ap=Ap_surv$year_index,
                 plot_Ap=Ap_surv$plot,
                 ##ELRI
                 y_Er=Er_surv$surv_t1, 
                 n_Er=length(Er_surv$surv_t1),
                 beta_Er_dim=ncol(X_Er),
                 X_Er=X_Er,
                 year_Er=Er_surv$year_index,
                 plot_Er=Er_surv$plot,
                 ##ELVI
                 y_Ev=Ev_surv$surv_t1, 
                 n_Ev=length(Ev_surv$surv_t1),
                 beta_Ev_dim=ncol(X_Ev),
                 X_Ev=X_Ev,
                 year_Ev=Ev_surv$year_index,
                 plot_Ev=Ev_surv$plot,
                 ##FESU
                 y_Fs=Fs_surv$surv_t1, 
                 n_Fs=length(Fs_surv$surv_t1),
                 beta_Fs_dim=ncol(X_Fs),
                 X_Fs=X_Fs,
                 year_Fs=Fs_surv$year_index,
                 plot_Fs=Fs_surv$plot,
                 ##POAL
                 y_Pa=Pa_surv$surv_t1, 
                 n_Pa=length(Pa_surv$surv_t1),
                 beta_Pa_dim=ncol(X_Pa),
                 X_Pa=X_Pa,
                 year_Pa=Pa_surv$year_index,
                 plot_Pa=Pa_surv$plot,
                 ##POAU
                 y_Pu=Pu_surv$surv_t1, 
                 n_Pu=length(Pu_surv$surv_t1),
                 beta_Pu_dim=ncol(X_Pu),
                 X_Pu=X_Pu,
                 year_Pu=Pu_surv$year_index,
                 plot_Pu=Pu_surv$plot,
                 ##POAU
                 y_Ps=Ps_surv$surv_t1, 
                 n_Ps=length(Ps_surv$surv_t1),
                 beta_Ps_dim=ncol(X_Ps),
                 X_Ps=X_Ps,
                 year_Ps=Ps_surv$year_index,
                 plot_Ps=Ps_surv$plot)

survival_model <- stan_model("analysis/Stan/ltreb_age_survival.stan")
surv_fit<-sampling(survival_model,data = stan_dat_surv,
                       chains=3,
                       #control = list(adapt_delta=0.99,stepsize=0.1),
                       iter=10000,thin=2,
                       pars = c("beta_Ap","beta_Er",
                                "beta_Ev","beta_Fs",
                                "beta_Pa","beta_Pu",
                                "beta_Ps","sigma_year","sigma_plot"), 
                       save_warmup=F)
write_rds(surv_fit,"analysis/Stan/surv_fit.rds")

## check a few trace plots
bayesplot::mcmc_trace(surv_fit,pars = c("sigma_year","sigma_plot"))
bayesplot::mcmc_trace(surv_fit,pars = c("beta_Ap[1]","beta_Er[1]"))

## wrangle parameter indices to get age- and endo-specific survival
## Agrostis perennans
Ap_par <- rstan::extract(surv_fit,pars="beta_Ap")
colnames(Ap_par$beta_Ap)<-colnames(X_Ap)
quantile_probs<-c(0.1,0.25,0.5,0.75,0.9)
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
colnames(Er_par$beta_Er)<-colnames(X_Er)
quantile_probs<-c(0.1,0.25,0.5,0.75,0.9)
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

## nice figure
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(Ap_surv$age_lump,Ap_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,6.5))
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

plot(Er_surv$age_lump,Er_surv$surv_t1,type="n",xlab="Age group",ylab="Survival",
     xlim=c(-0.5,4.5))
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

