#this file creates a species list for the supplement, as needed.
rm(list = ls()) # clears everything

#load libraries
library(tidyverse)

#load data
traitdat = read.csv("Data/transformed_traitdata.csv",header = T)  #created in Data prep_trait_transformations.R
          
summary(traitdat)
dim(traitdat)

spmeta = read.csv("Data/species_metadata.csv") %>%
         rename(Fullname = Species)
summary(spmeta)
dim(spmeta)

abundat = read.csv("Data/Seedbank_raw_abundances.csv") 

#merge data frames and create species lists
sptraitdat = traitdat %>%
             distinct(species)

dim(sptraitdat) #67 species with trait data
listsptraitdat = as.list(sptraitdat)

spabundat = abundat%>%
            group_by(species) %>% #if add treatment, list is very long, and those data are in Eskelinen et al. 2021
            summarize(minabund = min(rawcount,  na.rm=T), maxabund = max(rawcount, na.rm=T),
                      meanabund = mean(rawcount, na.rm=T)) %>%
            mutate(Abundance_range = paste(minabund, maxabund, sep = "-")) %>%
            dplyr::select(species, Abundance_range, meanabund)
            
dim(spabundat) #108 species in seedbank 
summary(spabundat)

#add column to seed bank species list for whether we have trait data

spdat = spabundat %>%
        mutate(InSeedbank = as.factor("Yes")) %>% 
        mutate(TraitData= as.factor(if_else(species %in% c(listsptraitdat$species), "Yes", "No"))) %>%
        rename(Species = species)
summary(spdat) #61 species in seed bank data that have traits

onlytraits = tibble(Species = setdiff(sptraitdat$species, spabundat$species)) %>% #6 species that we have traits but weren't in seedbank
             mutate(InSeedbank = as.factor("No")) %>%
             mutate(TraitData = as.factor("Yes")) %>%
             mutate(Abundance_range = NA, meanabund = NA)

allsbsp = rbind.data.frame(spdat, onlytraits) %>%
          rename(Abbreviation = Species)
summary(allsbsp)

#merge with full names and families
allsbsp2 = left_join(allsbsp, spmeta, by = "Abbreviation") %>%
           dplyr::select(Fullname, Abbreviation, Family = family, "Life History" = annual, Functional.group, Origin = origin,
                         meanabund, InSeedbank, TraitData) %>% #keep meanabund, range makes table too big
           dplyr::filter(Fullname != "NA") %>%
           arrange(Fullname)


summary(allsbsp2)
dim(allsbsp2) 
#write.csv(allsbsp2, "Data/Species list_with seed bank mean abund.csv")
