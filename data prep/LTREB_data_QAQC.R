## Authors: Bell Scherick, Tom Miller
## Purpose: prep grass-endo demographic data for analysis
## LDW data are already QAQC'd through EDI package -- that just needs to be loaded
## POAU data need QAQC, done here

library(readxl)
library(tidyverse)
## Note on data read-in: read_excel was coercing some data to NA
## I saved the most complete and up to date spreadsheet as .csv and that is read in here
## read_excel still used for plot read-in because that is small and simple

# Poa autumnalis in TX ----------------------------------------------------

## read in raw POAU data from Miller Lab drive folder
tom_dir<-"G:/Shared drives/Miller Lab/LTREB/POAU/"
##bell_dir<-
use_dir<-tom_dir
poau_demog <- read.csv(paste0(use_dir,"POAU_SFAEF_demography_20082023.csv"))
poau_plots <- read_excel(paste0(use_dir,"POAU SFAEF demography complete2.xlsx"),
                               sheet = "POAU_plot_assignment")
poau <- left_join(poau_demog,poau_plots,by="Plot")
## change name of org/rec
poau %>% rename(org_rec = Org.rec) %>% 
  ## add row and individual ID
  ## originals have ID numbers 1-20 that overlap with recruit numbers
  ## so need to include both tag and ID
  mutate(rowID = row_number(),
         indID = paste(Plot,Tag,ID))->poau

## do we have birth years and obs years for everything
sum(!is.na(poau$year_t)) ## how many year t entries
sum(!is.na(poau$year_t1)) ## how many year t1 entries -- more year t1's is bc of recruits
sum(!is.na(poau$year_recruit)) ## how many year_recruit entries
## this tells us that we are missing lots of year t entries for "real" data
#View(poau %>% filter(is.na(year_t) & !is.na(Seedling_t1)))
## we will fix this by subtracting 1 from year t1 -- we checked and these are all real obs
poau[is.na(poau$year_t),"year_t"]<-poau$year_t1[is.na(poau$year_t)]-1
## we stumbled on this one, garbage entry that we figured out should be zero
poau[which(poau$Seedling_t1==1 & poau$tiller_number_t>0),"Seedling_t1"]<-0

## next issue: new plants that were not seedlings
poau %>% 
  filter(Seedling_t1==1 & (tiller_number_t1>1 | inf_number_spring_t1>0)) -> inspect
inspect$spring_survival_notes_t1
change_to_seedling0 <- c(30,31,32,35,36,37,39,42,50,56,58,59,60,61,62,63,64,65,66,67,68)
## we do not believe that these are seedlings, and most of the notes indicate this, so change it
not_seedlings<-inspect[change_to_seedling0,]
## find all occurrences in the database
poau[poau$indID %in% not_seedlings$indID,"year_recruit"]<-NA

## now everything left in inspect that is flowering we will change to NA birth year bc we don't believe they are recruits
na_these_birth_years<-as.vector(inspect[inspect$inf_number_spring_t1>0,"indID"])
poau[poau$indID %in% na_these_birth_years,"year_recruit"]<-NA

## make sure that each individual carries its birth year every time it appears
poau %>% 
  group_by(indID) %>% select(year_recruit) %>% 
  summarize(birth_year=mean(year_recruit,na.rm=T))->birthyears
## find any non-integer birth years (means they have different birth years in different rows)
unique(birthyears$birth_year)
problem_births<-c(
birthyears[which(birthyears$birth_year>2011 & birthyears$birth_year<2012),"indID"]$indID,
birthyears[which(birthyears$birth_year>2013 & birthyears$birth_year<2014),"indID"]$indID,
birthyears[which(birthyears$birth_year>2016 & birthyears$birth_year<2017),"indID"]$indID,
birthyears[which(birthyears$birth_year>2018 & birthyears$birth_year<2019),"indID"]$indID,
birthyears[which(birthyears$birth_year>2021 & birthyears$birth_year<2022),"indID"]$indID,
birthyears[which(birthyears$birth_year>2022 & birthyears$birth_year<2023),"indID"]$indID)
poau %>% filter(indID %in% problem_births) %>% group_by(indID) %>% View
## There was a list of problem plants that I corrected directly in the raw data
## There is one remaining that I cannot figure out (notes suggest some re-labeling)
## drop this one problem plant, then merge birth years back into full demog data, and replace the year_recruit variable
poau %>% 
  filter(indID != "79 NA 228") %>% 
  left_join(.,birthyears,by="indID")->poau

## does every unique ID have a unique org/rec value?
poau %>% select(indID,org_rec) %>% distinct() -> poau_ind
## no, there are lots of NAs but I can tell that there are recruit plants based on ID numbers
poau$org_rec[is.na(poau$org_rec)]<-0

## is every new individual given seedling_t1==1 upon its first appearance?
poau %>% filter(year_recruit==year_t1) %>% 
  filter(is.na(Seedling_t1)) %>% View
# yes. And if we are not in the year of recruitment seedling t1 should never =1
poau %>% filter(year_recruit!=year_t1) %>% 
  filter(Seedling_t1==1) %>% View
# good

## add age as the difference of year_t and year_recruit
## note that newborns get age -1 when they appear in year t1
## but age zero when they first appear in year t (we'll drop the -1's)
poau$age<-poau$year_t-poau$birth_year
## assign age as NA for original plants
poau$age[poau$org_rec==1]<-NA
## there should be as many -1's as there are 0's
table(poau$age) #good

## if there is a tiller number count but inf count is NA, assume inf count is zero
poau$inf_number_spring_t[!is.na(poau$tiller_number_t) & is.na(poau$inf_number_spring_t)]<-0
poau$inf_number_spring_t1[!is.na(poau$tiller_number_t1) & is.na(poau$inf_number_spring_t1)]<-0

## weird stuff with spikelet counts
## there are two columns for spike_b_t but one of them has only two observations
subset(poau,!is.na(poau$spike_b_t.1)) %>% View()
## so just proceed with the main spike_b_t

## the c count has non-numeric entries
unique(poau$spike_c_t)
poau$spike_c_t[which(poau$spike_c_t=="too young")]<-NA
poau$spike_c_t[which(poau$spike_c_t=="-")]<-NA
unique(poau$spike_c_t1)
poau$spike_c_t1[which(poau$spike_c_t1=="too young")]<-NA
poau$spike_c_t1[which(poau$spike_c_t1=="-")]<-NA
## now coerce spikelet counts to numeric
poau$spike_c_t <- as.numeric(poau$spike_c_t)
poau$spike_c_t1 <- as.numeric(poau$spike_c_t1)

## for plants with spikelet counts, calculate their average - this takes a few seconds
poau %>% 
  rowwise %>% 
  mutate(mean_spike_t = mean(c_across(c(spike_a_t,spike_b_t,spike_c_t)),na.rm=T),
         mean_spike_t1 = mean(c_across(c(spike_a_t1,spike_b_t1,spike_c_t1)),na.rm=T))->poau

## finally, drop the 2019-2021 "transition year"
## if there are non-NA entries in year t 2019, do not use their year t1 data
poau %>% mutate(trans_year = paste(year_t,year_t1,sep=".")) %>% 
  filter(trans_year != "2019.2021") -> poau
## note that 2020-2021 is still here because these are rows for plants first
## tagged in 2021

## other QA/QC -- check ranges of variables and class types
str(poau) # two problems
unique(poau$spring_survival_t1)
unique(poau$inf_number_spring_t1)
## reassign NA and coerce to numeric
poau$spring_survival_t1[poau$spring_survival_t1=="TNF"]<-NA
poau$spring_survival_t1<-as.numeric(poau$spring_survival_t1)
poau$inf_number_spring_t1[poau$inf_number_spring_t1=="."]<-NA
poau$inf_number_spring_t1<-as.numeric(poau$inf_number_spring_t1)

## check that size and flower count are always integer -- all good
which(poau$tiller_number_t-floor(poau$tiller_number_t)>0)
which(poau$tiller_number_t1-floor(poau$tiller_number_t1)>0)
which(poau$inf_number_spring_t-floor(poau$inf_number_spring_t)>0)
which(poau$inf_number_spring_t1-floor(poau$inf_number_spring_t1)>0)

# Indiana LDW -------------------------------------------------------------

## These data are published on EDI here
## https://doi.org/10.6073/pasta/ea7db07a578fb030a173f37f76596b62 (Accessed 2024-03-26).
## read in data and apply some of the data transformations used above
indiana<-read.csv("data prep/LDW_LTREB_20072022.csv")
str(indiana)
indiana %>% 
  rowwise %>% 
  mutate(mean_spike_t = mean(c_across(c(spike_a_t,spike_b_t,spike_c_t)),na.rm=T),
         mean_spike_t1 = mean(c_across(c(spike_a_t1,spike_b_t1,spike_c_t1)),na.rm=T))->indiana

## check that size and flower count are always integer
which(indiana$flw_count_t-floor(indiana$flw_count_t)>0)#fine
which(indiana$flw_count_t1-floor(indiana$flw_count_t1)>0)#fine
indiana[which(indiana$size_t-floor(indiana$size_t)>0),] %>% View
indiana[which(indiana$size_t1-floor(indiana$size_t1)>0),] %>% View
## not sure what's up with these so I will reassign as NA
indiana[which(indiana$size_t-floor(indiana$size_t)>0),"size_t"]<-NA
indiana[which(indiana$size_t1-floor(indiana$size_t1)>0),"size_t1"]<-NA

## add age as the difference of year_t and year_recruit
indiana$age<-indiana$year_t-indiana$birth
## assign age as NA for original plants
indiana$age[indiana$origin_01==0]<-NA
## are there as many -1's as there are 0's?
table(indiana$age) #no there are way more 0s
## this must be because the recruit data were managed differently in the early years
## all the -1s (first appearance in year t1) started in 2016 or later

## check that -1 ages are always cases with no data in year_t - FUCK
indiana[which(indiana$age==-1 & !is.na(indiana$size_t)),]->bad_plants
indiana %>% filter(id %in% bad_plants$id) %>% View

## many of these have NA id, why?
indiana[which(is.na(indiana$id)),]->no_id

## apply rule that if size is non-NA and inf count is NA, then inf count should be zero
indiana$flw_count_t[!is.na(indiana$size_t) & is.na(indiana$flw_count_t)]<-0
indiana$flw_count_t1[!is.na(indiana$size_t1) & is.na(indiana$flw_count_t1)]<-0

## check that each time an individual appears it carries the same birth year
indiana %>% 
  group_by(id) %>% select(birth) %>% 
  summarize(birth_year=mean(birth,na.rm=T))->birthyears
## find any non-integer birth years (means they have different birth years in different rows)
problems<-which(birthyears$birth_year - floor(birthyears$birth_year) != 0)
##there are this many individuals with more than one birth year
length(problems) #ugh
## who are they?
indiana %>% filter(id %in% birthyears$id[problems]) %>% View
## are these usually only off by one year?
indiana %>% filter(id %in% birthyears$id[problems]) %>% 
  group_by(id) %>% 
  summarise(min_birth = min(birth),
            max_birth = max(birth),
            range = max_birth-min_birth) %>% View
## jesus -- one plant has a range of 10 years

## some of the problems are likely numbers that we reused, which I could know
## if the plant died but the number reappeared later. In other cases it is likely 
## a data entry or copying error, and it is probably safe to assume the earlier year
## is the correct one, as long as everything else checks out
## I will go through these manually and hand-pick cases where the earliest
## year is not the correct birth year

## For now I want to keep moving, so I am just going to drop "problems".
## Losing 596 rows. Tom will return to this (-TM 4/9/2024)
indiana_no_problems <- indiana %>% filter(!(id %in% birthyears$id[problems]))

##prepare indiana data for row bind with poau data
indiana_no_problems %>% 
  ##change org_rec so that original=1 / recruit=0
  mutate(original=abs(origin_01-1)) %>% 
  select(species,plot,endo_01,id,original,endo_status_from_check,birth,year_t,age,
         size_t,flw_count_t,mean_spike_t,year_t1,surv_t1,
         size_t1,flw_count_t1,mean_spike_t1) -> indiana_for_merge

##prepare POAU data for row bind with indiana data
poau %>% 
  mutate(species="POAU") %>% 
  ## drop 2022-23 just to have the same years as indiana
  filter(year_t<2022) %>% 
  select(species,Plot,Plot_endo_status,indID,org_rec,Endo,birth_year,year_t,age,
         tiller_number_t,inf_number_spring_t,mean_spike_t,year_t1,spring_survival_t1,
         tiller_number_t1,inf_number_spring_t1,mean_spike_t1) %>% 
  rename(species=species,plot=Plot,endo_01=Plot_endo_status,id=indID,original=org_rec,
         endo_status_from_check=Endo,birth=birth_year,year_t=year_t,age=age,
         size_t=tiller_number_t,flw_count_t=inf_number_spring_t,mean_spike_t=mean_spike_t,
         year_t1=year_t1,surv_t1=spring_survival_t1,
         size_t1=tiller_number_t1,flw_count_t1=inf_number_spring_t1,mean_spike_t1=mean_spike_t1)-> poau_for_merge

## should have same number of columns
dim(poau_for_merge);dim(indiana_for_merge)


## write out a cleaned, derived data frame
ltreb<-bind_rows(indiana_for_merge,poau_for_merge)
names(ltreb)
write.csv(ltreb,"data prep/ltreb_allspp_qaqc.csv")
