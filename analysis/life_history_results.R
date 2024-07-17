setwd("C:/Users/tm9/Dropbox/github/LTREB-life-history")
library(tidyverse)
library(factoextra)
##read in life history outputs
lifehistorypost<-read.csv("analysis/lifehistorypost.csv")

# PCA ---------------------------------------------------------------------

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

pca.dat<-bind_cols(R0$R0,G$G,meanelexp$meanelexp,longevity$longevity,entropyd$entropyd)
names(pca.dat)<-c("R0","Gen time","Life expect","Longevity","EntropyD")
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


                