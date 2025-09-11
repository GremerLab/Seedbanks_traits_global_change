#this file creates a species list for the supplement, as needed.
rm(list = ls()) # clears everything

#load libraries
library(tidyverse)

#load data
traitdat = read.csv("Cleaned data/transformed_traitdata.csv",header = T) %>% #created in 4_trait_transformations.R
           rename(Species = clean_code) 
summary(traitdat)
dim(traitdat)

spmeta = read.csv("Raw data/species_metadata.csv") %>%
         rename(Fullname = Species, Species = code) %>%
          mutate(origin = as.factor(ifelse(nat.inv=="native","Native","Non-Native"))) %>%
         select( -MarinaSp, -nat.inv, -tolumpname, -falsename, -species_code) %>%
         distinct()
summary(spmeta)
dim(spmeta)

abundat = read.csv("Cleaned data/SBRAWbyspecies_tomatch2017_05112023_cleancode.csv") %>%
          rename(Species = clean_code) %>%
          select( -MarinaSp, -nat.inv)

#merge data frames and create species lists
sptraitdat = traitdat %>%
             distinct(Species)

dim(sptraitdat) #67 species with trait data
listsptraitdat = as.list(sptraitdat)

spabundat = abundat%>%
            distinct(Species)
dim(spabundat) #108 species in seedbank


#add column to seed bank species list for whether we have trait data

spdat = spabundat %>%
        mutate(InSeedbank = as.factor("Yes")) %>% 
        mutate(TraitData= as.factor(if_else(Species %in% c(listsptraitdat$Species), "Yes", "No")))
summary(spdat) #61 species in seed bank data that have traits

onlytraits = tibble(Species = setdiff(sptraitdat$Species, spabundat$Species)) %>% #6 species that we have traits but weren't in seedbank
             mutate(InSeedbank = as.factor("No")) %>%
             mutate(TraitData = as.factor("Yes"))

allsbsp = rbind.data.frame(spdat, onlytraits)
summary(allsbsp)

#merge with full names and families
allsbsp2 = left_join(allsbsp, spmeta, by = "Species") %>%
           select(Fullname, SpeciesCode = Species, Annual = annual, Functional.group, Origin = origin,
                  InSeedbank, TraitData)


summary(allsbsp2)
dim(allsbsp2) #why did it add 12 species? 
#write.csv(allsbsp2, "Output data/Species list.csv")
