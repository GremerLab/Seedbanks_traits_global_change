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
            distinct(species)
dim(spabundat) #108 species in seedbank


#add column to seed bank species list for whether we have trait data

spdat = spabundat %>%
        mutate(InSeedbank = as.factor("Yes")) %>% 
        mutate(TraitData= as.factor(if_else(species %in% c(listsptraitdat$species), "Yes", "No"))) %>%
        rename(Species = species)
summary(spdat) #61 species in seed bank data that have traits

onlytraits = tibble(Species = setdiff(sptraitdat$species, spabundat$species)) %>% #6 species that we have traits but weren't in seedbank
             mutate(InSeedbank = as.factor("No")) %>%
             mutate(TraitData = as.factor("Yes"))

allsbsp = rbind.data.frame(spdat, onlytraits) %>%
          rename(Abbreviation = Species)
summary(allsbsp)

#merge with full names and families
allsbsp2 = left_join(allsbsp, spmeta, by = "Abbreviation") %>%
           dplyr::select(Fullname, Abbreviation, Family = family, "Life History" = annual, Functional.group, Origin = origin,
                  InSeedbank, TraitData)


summary(allsbsp2)
dim(allsbsp2) 
#write.csv(allsbsp2, "Data/Species list.csv")
