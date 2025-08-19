## Things to run in referencing script
# rm(list = ls()) # clears everything
# library(tidyverse)

# dat2fixnames=read.csv("Seed Trait Paper/Analysis_data/SBRA_traitdat_Oct11.csv",header = T)
# dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
# dat2fixnames$species_code=NULL

## old species lists
# spnames=read.csv("Seed Trait Paper/Master_Species_name_key.csv", header=T)
# spnames=read.csv("Seed Trait Paper/species_metadata_Nov22.csv",header=T, strip.white = T)
# spnames=read.csv("Seed Trait Paper/updated_name_issues1_updated.csv",header=T, strip.white = T)

##Autorun from here#####

# Species name data
spnames=read.csv("Cleaned data/updated_name_issues2_withIDsforall_Dec2022.csv",header=T, strip.white = T)

if("IDnum" %in% names(dat2fixnames)==T){
  spnames=spnames%>%
    mutate(clean_code=toupper(as.character(species_code)))%>%
    select(-Species_name,-falsename,-tolumpname, species_code)%>%
    distinct()
  
  #Easy join by ID number, left_join discards unneeded names
  cleannamedat=dat2fixnames%>%
    left_join(spnames)%>%
    distinct()
  }else{
  print("ID number missing")
}

## FILGal is lumping with filcal but is invasive

rm(spnames,dat2fixnames)

