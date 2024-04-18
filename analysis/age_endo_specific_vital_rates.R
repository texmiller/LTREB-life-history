## Purpose: fit models for age- and endo-specific survival and fertility
library(tidyverse)

## read in QAQC'd data
ltreb<-read.csv("./data prep/ltreb_allspp_qaqc.csv") %>% 
  ## drop original plants
  filter(original==0) %>% 
  ## drop LOAR (so few recruits)
  filter(species!="LOAR") %>% 
  ## drop rows with -1 ages (this is first appearance in)


## what are the age limits we can use for each species
table(ltreb$age,ltreb$species,ltreb$endo_01)

## three vital rates: survival, fertility, and recruitment