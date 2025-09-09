#This script 
#Calculates community weighted mean (CWM) traits by plot and treatment
#tests whether CWM traits vary by treatment and habitat (harsh vs. lush serpentine)
rm(list = ls()) # clears console

traitdat = read.csv("Output data/Traits_PCscores_all.csv") #see 5_Analyses_TraitPCA.R for script generating this file
#note traitdat is using the transformed trait values as they were input into the PCA in script 5
summary(traitdat)
str(traitdat)
dim(traitdat)

abundat = read.csv("Cleaned data/SBRAWbyspecies_tomatch2017_05112023_cleancode.csv") %>% #see 3_Final_namecleaning.R for script generating this file, this matches raw seedling counts from Eskelinen et al. 2021, seedbankgrowoutallyrs_wide_totnumseedlings_tomatch2017.csv
          mutate(species = clean_code, rawcount = relcover) %>% #despite naming in file, it's actually raw counts not relative cover
          select(-clean_code, -family, -annual, -MarinaSp, -nat.inv, -Functional.group, -relcover) #these are updated in traitdat, can remove from here
#JRG notes, all abundances using the RAW data file are raw abundances, not relative, and has adapted scripts accordingly

summary(abundat)
str(abundat)
dim(abundat)

abund_trait_dat = left_join(abundat, traitdat, by= "species") %>% #
                  mutate(Line = as.factor(Line), habitat=as.factor(habitat), watering=as.factor(watering), 
                         fertilization=as.factor(fertilization), WFtreatment=as.factor(WFtreatment), Plot=as.factor(Plot), 
                         species = as.factor(species), origin = as.factor(origin))
dim(abund_trait_dat)
summary(abund_trait_dat)

#get relative abundances
rawab_plots=abund_trait_dat%>%
  group_by(Plot,WFtreatment)%>%
  dplyr::summarise(totplotincmissing=sum(rawcount)) #previously was done as a presence/absence (sum(as.numeric(as.character()...)))

#below calculates the total abundance for species in those plots for which we had trait data
rawab_wtraits =abund_trait_dat%>%
  filter(is.na(Mass)==F & is.na(CN)==F)%>% #removes rows that are missing seed mass or cn value
  group_by(Plot,WFtreatment)%>%
  dplyr::summarise(totplottraits=sum(rawcount))

rawab = full_join(rawab_plots, rawab_wtraits)
summary(rawab)

rawab = full_join(rawab_plots, rawab_wtraits)


all_abtraits = left_join(abund_trait_dat, rawab) %>%
               mutate(relab = rawcount/totplotincmissing) %>%
               mutate(relab_traitsonly = rawcount/totplottraits) #shouldn't use this

#double check that relative abundance is correct (should add to 1 for each plot)
summary(all_abtraits %>% 
  group_by(Plot) %>%
  summarize(sumrelab = sum(relab, na.rm=T)))
    
#### calculate CWMs ####
#Note that JRG did this slightly differently than ECE dissertation 
#first, get trait value * relative abundance for each species in each plot 
cwmdat_a = all_abtraits %>%
         mutate(across(c(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                         Texture, Compact), ~.x*relab, .names = "{.col}_CWM")) %>%
         mutate(across(starts_with("PC"), ~.x*relab,.names = "{.col}_CWM") ) 

#next, sum for each plot to get a plot level CWM
cwmdat = cwmdat_a %>%
         group_by(Line, habitat, watering, fertilization, WFtreatment, Plot) %>%
         summarize(across(ends_with("_CWM"), ~sum(.x, na.rm=T)))%>%
         ungroup()

summary(cwmdat)        
dim(cwmdat) #90 plots
str(cwmdat)
#write.csv(cwmdat, "Output data/CWMs_traits_PCs.csv")