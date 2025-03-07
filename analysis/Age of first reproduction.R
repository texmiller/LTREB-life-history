library(tidyverse)
library(scales)
library(xtable)
library(lme4)
library(optimx)
library(dplyr)
library(readr)


## set working directory
bell = "/Users/bell/Documents/GitHub/LTREB-life-history/data prep"
setwd(bell)

## read in QAQC'd data
## as of 4.18.2024 there are still some data issues (see comments) but this is 
## far enough along to proceed
ltreb_allplants<-read.csv("ltreb_allspp_qaqc.csv") %>% 
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

ltreb_age_lump %>% 
  group_by(species, id, endo_01) %>% 
  filter(flw_count_t >0) %>% 
  mutate(min_age = min(age_lump)) %>% 
  select (species, id, min_age, endo_01, year_t1, plot) -> min_age_repro

#Perlim anlysis (mean of age of first repro)
min_age_repro %>% 
  group_by(species, endo_01) %>% 
  summarise(mean_min_age = mean(min_age)) %>% 
  view()

#Possion modle of age of first reproduction with year and plot as random effects
mod_min_age_repro <- glmer(min_age ~ endo_01 * species + (1 | year_t1) + (1 | plot), 
                   data = min_age_repro, 
                   family = "poisson")
#collecing estamated coefficients
age_repro = summary(mod_min_age_repro)$coefficients[,1]

min_repro_agpe_em = exp(age_repro[1])
min_repro_agpe_ep = exp(age_repro[1] + age_repro[2])

min_repro_elri_em = exp(age_repro[1] + age_repro[3])
min_repro_elri_ep = exp(age_repro[1] + age_repro[3] + age_repro[2] +  age_repro[9])

min_repro_elvi_em = exp(age_repro[1] + age_repro[4])
min_repro_elvi_ep = exp(age_repro[4] + age_repro[2] +  age_repro[10])

min_repro_fesu_em = exp(age_repro[1] + age_repro[5])
min_repro_fesu_ep = exp(age_repro[1] + age_repro[5] + age_repro[2] +  age_repro[11])

min_repro_poal_em = exp(age_repro[1] + age_repro[6])
min_repro_poal_ep = exp(age_repro[1] + age_repro[6] + age_repro[2] +  age_repro[12])

min_repro_poau_em = exp(age_repro[1] + age_repro[7])
min_repro_poau_ep = exp(age_repro[1] + age_repro[7] + age_repro[2] + age_repro[13])

min_repro_posy_em = exp(age_repro[1] + age_repro[8])
min_repro_posy_ep = exp(age_repro[1] + age_repro[8] + age_repro[2] + age_repro[14])


#Creating datafram with all the information
min_age_repro_modle_results = data.frame(
  species= rep(c("agpe", "elri", "elvi", "fesu", "poal", "poau", "posy"), each= 2),
                                        endo = rep(c(0,1), times = 7),
                                        min_repro =c(min_repro_agpe_em, min_repro_agpe_ep,
                                                     min_repro_elri_em, min_repro_elri_ep,
                                                     min_repro_elvi_em, min_repro_elvi_ep,
                                                     min_repro_fesu_em, min_repro_fesu_ep,
                                                     min_repro_poal_em, min_repro_poal_ep,
                                                     min_repro_poau_em, min_repro_poau_ep,
                                                     min_repro_posy_em, min_repro_posy_ep),
  color = rep(c("tomato","cornflowerblue"), times = 7)
  )



ggplot(data = min_age_repro_modle_results, aes(x = species, y = min_repro, fill = color)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("cornflowerblue","tomato"), labels = c("E+", "E-")) + # Customize labels as needed
  labs(
    fill = "Endophyte Status",
    x = "Species",
    y = "Expected Age of First Reproduction"
  ) +
  theme_minimal()

bell = "/Users/bell/Documents/GitHub/LTREB-life-history"
setwd(bell)
  

lifehistorypost<-read.csv("analysis/lifehistorypost.csv")

# PCA ---------------------------------------------------------------------
library(factoextra)
## first reduce data frame down to averages
lifehistorypostmean<-lifehistorypost %>% 
  group_by(species) %>% 
  summarise_all(mean) %>% select(-X)

## select traits that will go into PCA
lifehistorypostmean %>% 
  select(species,R0_em,R0_ep) %>% 
  pivot_longer(R0_em:R0_ep,names_to="endo",values_to="R0")->R0
lifehistorypostmean %>% 
  select(species,G_em,G_ep) %>% 
  pivot_longer(G_em:G_ep,names_to="endo",values_to="G")->G
lifehistorypostmean %>% 
  select(species,meanelexp_em,meanelexp_ep) %>% 
  pivot_longer(meanelexp_em:meanelexp_ep,names_to="endo",values_to="meanelexp")->meanelexp
lifehistorypostmean %>% 
  select(species,matlifexp_em,matlifexp_ep) %>% 
  pivot_longer(matlifexp_em:matlifexp_ep,names_to="endo",values_to="matlifexp")->matlifexp
lifehistorypostmean %>% 
  select(species,longevity_em,longevity_ep) %>% 
  pivot_longer(longevity_em:longevity_ep,names_to="endo",values_to="longevity")->longevity
lifehistorypostmean %>% 
  select(species,entropyd_em,entropyd_ep) %>% 
  pivot_longer(entropyd_em:entropyd_ep,names_to="endo",values_to="entropyd")->entropyd

pca.dat<-bind_cols(R0$R0,G$G,meanelexp$meanelexp,longevity$longevity,entropyd$entropyd, min_age_repro_modle_results$min_repro)
names(pca.dat)<-c("R0","Gen time","Life expect","Longevity","EntropyD", "min_repro")
row.names(pca.dat)<-c("AGPE-","AGPE+",
                      "ELRI-","ELRI+",
                      "ELVI-","ELVI+",
                      "FESU-","FESU+",
                      "POAL-","POAL+",
                      "POAU-","POAU+",
                      "POSY-","POSY+")

lifehistory_pca<-prcomp(pca.dat,scale=T)
plot(lifehistory_pca$x,pch=c(1,16))
arrows(lifehistory_pca$x[c(1,3,5,7,9,11,13),1],
       lifehistory_pca$x[c(1,3,5,7,9,11,13),2],
       lifehistory_pca$x[c(2,4,6,8,10,12,14),1],
       lifehistory_pca$x[c(2,4,6,8,10,12,14),2])

fviz_pca_var(lifehistory_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

groups<-as.factor(R0$endo)
fviz_pca_biplot(lifehistory_pca, repel = TRUE,
                col.var = "gray", # Variables color
                col.ind = groups, # color by groups
                palette = c("tomato",  "cornflowerblue")  # Individuals color
)+
  annotate("segment",x=lifehistory_pca$x[c(1,3,5,7,9,11,13),1], 
           y=lifehistory_pca$x[c(1,3,5,7,9,11,13),2], 
           xend=lifehistory_pca$x[c(2,4,6,8,10,12,14),1], 
           yend=lifehistory_pca$x[c(2,4,6,8,10,12,14),2],
           arrow = arrow(length=unit(.4, 'cm')),col=alpha("cornflowerblue",0.25))+
  theme(legend.position = "none")

# posterior density plots of single traits ---------------------------------
## generate E+/E- contrasts
lifehistorypost$R0_diff<-lifehistorypost$R0_ep-lifehistorypost$R0_em
lifehistorypost$G_diff<-lifehistorypost$G_ep-lifehistorypost$G_em
lifehistorypost$lambda_diff<-lifehistorypost$lambda_ep-lifehistorypost$lambda_em
lifehistorypost$pRep_diff<-lifehistorypost$pRep_ep-lifehistorypost$pRep_em
lifehistorypost$La_diff<-lifehistorypost$La_ep-lifehistorypost$La_em
lifehistorypost$matlifexp_diff<-lifehistorypost$matlifexp_ep-lifehistorypost$matlifexp_em
lifehistorypost$meanelexp_diff<-lifehistorypost$meanelexp_ep-lifehistorypost$meanelexp_em
lifehistorypost$entropyd_diff<-lifehistorypost$entropyd_ep-lifehistorypost$entropyd_em
lifehistorypost$entropyk_diff<-lifehistorypost$entropyk_ep-lifehistorypost$entropyk_em
lifehistorypost$gini_diff<-lifehistorypost$gini_ep-lifehistorypost$gini_em
lifehistorypost$longevity_diff<-lifehistorypost$longevity_ep-lifehistorypost$longevity_em

## some quick and dirty plots of life history effects
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(pr_epos_lambda = mean(lambda_diff>0),
            pr_epos_R0 = mean(R0_diff>0),
            pr_epos_G = mean(G_diff>0),
            pr_epos_meanelexp = mean(meanelexp_diff>0),
            pr_epos_longevity = mean(longevity_diff>0))

##lambda
ggplot(lifehistorypost,aes(lambda_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##R0
ggplot(lifehistorypost,aes(R0_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")


##G
ggplot(lifehistorypost,aes(G_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##pRep-- there is no variation here - not sure if this is calculated correctly
ggplot(lifehistorypost,aes(pRep_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##La--same
ggplot(lifehistorypost,aes(La_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##mean life expectancy
ggplot(lifehistorypost,aes(meanelexp_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##life expectancy from maturity
ggplot(lifehistorypost,aes(matlifexp_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##entropy D
ggplot(lifehistorypost,aes(entropyd_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##entropy K
ggplot(lifehistorypost,aes(entropyk_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

## Gini -- no variation here
ggplot(lifehistorypost,aes(gini_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

## longevity
ggplot(lifehistorypost,aes(longevity_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")


         
       