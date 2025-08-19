#TO DO on this script: need missing files 

rm(list = ls()) # clears everything
library(tidyverse)

#filelog
## IN
# mf_sm_summary_120222.csv, cn_summary_120222.csv,SCP_summary_120222.csv
# SBRAbyspecies_tomatch2017_120222.csv
# Germ_expt_Nov2022.csv
# abovegroundcommunitydata2010_relabundance_28June19.csv
# abovegroundcommunitydata2017_relabundance_remove0cols_28June19

## Run within this script
# source("Seed Trait Paper/Scripts/name_cleaning_tosource_spcode.R")
# source("Seed Trait Paper/Scripts/name_cleaning_tosource_ID.R")

## out
# trait_data_Dec22.csv
# Germ_Dec22.csv
# SBRA_trait_long_Dec2022.csv
# AGRA_2010_long_Dec2022.csv,AGRA_2017_long_Dec2022.csv
# AGRA_summarydat_Dec2022.csv

#### read in and name clean for germination and relative abundance data ####
# dat2fixnames=read.csv("Seed Trait Paper/clean_final_subdata/SBRAbyspecies_tomatch2017_120222.csv",header = T,strip.white = T)
dat2fixnames=read.csv("Cleaned data/SBRAWbyspecies_tomatch2017_05112023.csv",header = T,strip.white = T)
dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
dat2fixnames$species_code=NULL
dat2fixnames$obsID=c(1:length(dat2fixnames$Line))

# input dat2fixnames with column prevname for the code
 source("Scripts/2_name_cleaning_tosource_spcode.R")
 # output = cleannamedat with column clean_code
SBRAdat=cleannamedat

## Fix name issues that arrise from grouping species. Ex. Agoseris and Agohet are grouped, but while names are all changed to Agohet, values from each are kept as they are separate observations from the seed bank abundance dataframe from Anu (even though one of them is a zero.)
SBRAdat=SBRAdat%>%
  filter(is.na(clean_code)==F)
SBRAdat=SBRAdat%>%
  group_by(Plot,IDnum)%>%
  summarise(truecover=sum(relcover))%>%
  inner_join(SBRAdat)%>%
  select(-obsID,-relcover)%>%
  distinct()

datGerm=read.csv("Raw data/Germ_expt_Nov2022.csv",header = T)
datGerm=datGerm%>%mutate(IDnum = Species..)
names(datGerm)
head(datGerm)
dat2fixnames=datGerm%>%
  select(-Species.Name,-Species..)
source("Scripts/2_name_cleaning_tosource_ID.R")
GermExt=cleannamedat

### read in SCP, CN, and MFsm, VM, starch, and clean names before combining #####
SCP=read.csv("Cleaned data/SCP_summary_120222.csv",header=T,strip.white = T)
CN=read.csv("Cleaned data/cn_summary_120222.csv",header=T,strip.white=T)
MFsm=read.csv("Cleaned data/mf_sm_summary_120222.csv",header = T,strip.white = T)

head(SCP)
head(CN)
dat2fixnames=CN
dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
dat2fixnames$species_code=NULL
dat2fixnames$obsID=c(1:length(dat2fixnames$prevname))
source("Scripts/2_name_cleaning_tosource_spcode.R")
head(cleannamedat)
CN=cleannamedat

## input dat2fixnames with column prevname for the code
dat2fixnames=MFsm
dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
dat2fixnames$species_code=NULL
dat2fixnames$obsID=c(1:length(dat2fixnames$prevname))
source("Scripts/2_name_cleaning_tosource_spcode.R")
head(cleannamedat)
MFsm=cleannamedat

###
head(SCP)
dat2fixnames=SCP
source("Scripts/2_name_cleaning_tosource_ID.R")
head(cleannamedat)
SCP=cleannamedat

traitdat=MFsm%>%
  select(-obsID,-species)%>%
  right_join(SCP)

traitdat=traitdat%>%
  right_join(CN)%>%
  select(-species_code,-obsID)

#### MISSING FILES ####
## Starch
Starch=read.csv("Seed Trait Paper/current_csv files/Starch_20221222.csv", header=T, strip.white = T)

Starch=Starch%>%select(IDnum="Idnum",starch,mucilage)
traitdat2=traitdat%>%
  left_join(Starch)


## Morphological categorization
Morph=read.csv("Seed Trait Paper/current_csv files/seedtraitphotos_categories_Dec2022_v1.csv", header=T, strip.white = T)
names(Morph)
Morph1=Morph[c(1,4:6,8)]
# head(Morph2)
##Issue avefat has both types many short 1 long

# How many categories (currently 11, change to 7)
Morph1%>%
  group_by(category, subcategory)%>%
  summarise(num=n())

Morph2=Morph1%>%
  mutate(lump_dispcat=case_when(sp..== 6 ~ 'elong_long',
                                category =='balloon'~'balloon',
                                category == 'flat'~ 'flat',
                                subcategory =='many long'~'elong_long',
                                subcategory =='one long'~'elong_long',
                                subcategory =='many short'~'elong_short',
                                subcategory =='one short'~'elong_short',
                                category =='other'~'other',
                                TRUE ~ subcategory))

Morph2%>%
  group_by(lump_dispcat)%>%
  summarise(num=n())
Morph2%>%
  group_by(fruit_type,germ_is_disp)%>%
  summarise(num=n())
Morph2%>%
  group_by(lump_dispcat,germ_is_disp)%>%
  summarise(num=n())
  

Morph3=Morph2%>%
  select(IDnum="sp..",germ_is_disp,lump_dispcat)%>%
  distinct()

traitdat3=traitdat2%>%
  full_join(Morph3)

## Videometer data
# VM=read.csv("Seed Trait Paper/clean_final_subdata/VM_byspecies.csv",header=T,strip.white = T) # the full data is a lot of variables. For now using sub data.
VM=read.csv("Seed Trait Paper/clean_final_subdata/VM_byspecies_Jan20.csv",header=T,strip.white = T) # the full data is a lot of variables. For now using sub data.

traitdat=traitdat3%>%
  left_join(VM)


traitwithgerm=traitdat2%>%
  right_join(GermExt)%>%
  select(-species_code)

names(traitwithgerm)
names(traitdat)

### SCT data
SCT=read.csv("Seed Trait Paper/clean_final_subdata/SCT_micrometer.csv",header=T)

traitdat=traitdat%>%
  left_join(SCT)

# write.csv(traitdat,"Seed Trait Paper/clean_final_subdata/trait_data_Apr11.csv",row.names = F)
# write.csv(traitdat,"Seed Trait Paper/clean_final_subdata/trait_data_Jan20.csv",row.names = F)
# write.csv(traitdat,"Seed Trait Paper/clean_final_subdata/trait_data_Jan17.csv",row.names = F)
  # write.csv(traitdat,"Seed Trait Paper/clean_final_subdata/trait_data_Dec23.csv",row.names = F)
 # write.csv(traitwithgerm,"Seed Trait Paper/clean_final_subdata/Germ_Dec22.csv",row.names = F)


SBRA_traitdat=SBRAdat%>%
  full_join(traitdat, by =c("IDnum"="IDnum","clean_code"="clean_code"))%>%
  distinct()


temp=SBRA_traitdat%>%
  anti_join(SBRAdat) #12 species that have some trait data but not any ind in the seedbank.

names(SBRAdat)
"AGOHET" %in% SBRAdat$clean_code
"BROELE" %in% SBRAdat$clean_code
"CLAGRA" %in% SBRAdat$clean_code
"CLARKIASP" %in% SBRAdat$clean_code
"MEDPOL" %in% SBRAdat$clean_code

# write.csv(SBRA_traitdat,"Seed Trait Paper/clean_final_subdata/SBRAW_trait_long_April2023.csv",row.names = F)
# write.csv(SBRA_traitdat,"Seed Trait Paper/clean_final_subdata/SBRA_trait_long_Jan2023v2.csv",row.names = F)
# write.csv(SBRA_traitdat,"Seed Trait Paper/clean_final_subdata/SBRA_trait_long_Jan23.csv",row.names = F)
# write.csv(SBRA_traitdat,"Seed Trait Paper/clean_final_subdata/SBRA_trait_long_Dec2022v3.csv",row.names = F)

traitdat2%>%filter(clean_code=="MEDPOL")
SBRA_traitdat%>%filter(clean_code=="MEDPOL")
splist1=unique(SBRA_traitdat$IDnum)
splist2=unique(traitdat2$IDnum)
a=splist2 %in% splist1
data.frame(a,splist2)%>%
  filter(a==0)

######### Aboveground #####
### Look at species distributions from 2010. Fix names. 
#2010
dat=read.csv("Data/abovegroundcommunitydata2010_relabundance_28June19.csv", header=T)
names(dat)
dat=dat[10:length(dat)]
dat1=dat%>%
  pivot_longer( names_to = "Species", cols = -Plot)%>%
  filter(value>0)%>%
  group_by(Species)%>%
  summarise(nplots2010=n(),maxcov2010=max(value))

head(dat1)

dat2fixnames=dat%>%
  pivot_longer( names_to = "prevname", cols = -Plot, values_to = "AGRA_2010")%>%
  filter(AGRA_2010>0)

dat2fixnames$obsID=c(1:length(dat2fixnames$Plot))
head(dat2fixnames)

# input dat2fixnames with column prevname for the code
source("Seed Trait Paper/Scripts/2_name_cleaning_tosource_spcode.R")
# output = cleannamedat with column clean_code
AGRA_2010dat=cleannamedat
head(cleannamedat)
AGRA_2010dat$obsID=NULL

#2017 species and plots aboveground cover
# dat=read.csv("Data/abovegroundcommunitydata2017_relabundance_remove0cols_28June19.csv", header=T) # relative cover is 127 species in 90 plots

dat=read.csv("Anu Paper/datafiles/aboveground community 2017_raw abund_for Jenny_zerocolumns removed from Anu_27Feb2020_lumpspecies.csv",header=T) # raw is 91 rows by 107 columns: row 91 is species totals
names(dat)
dat=dat[13:length(dat)]

## OVERVIEW OF COMMON SPECIES
dat2=dat%>%
  pivot_longer( names_to = "Species", cols = -Plot)%>%
  filter(value>0)%>%
  group_by(Species)%>%
  summarise(nplots2017=n(),maxcov2017=max(value))

head(dat2)

dat2fixnames=dat%>%
  pivot_longer( names_to = "prevname", cols = -Plot, values_to = "AGRA_2017")%>%
  filter(AGRA_2017>0)

dat2fixnames$obsID=c(1:length(dat2fixnames$Plot))
head(dat2fixnames)

# input dat2fixnames with column prevname for the code
source("Seed Trait Paper/Scripts/2_name_cleaning_tosource_spcode.R")
# output = cleannamedat with column clean_code
AGRA_2017dat=cleannamedat
head(cleannamedat)
AGRA_2017dat$obsID=NULL
# 
# # 
# write.csv(AGRA_2010dat,"Seed Trait Paper/clean_final_subdata/AGRA_2010_long_Dec2022.csv",row.names = F)
# write.csv(AGRA_2017dat,"Seed Trait Paper/clean_final_subdata/AGRA_2017_long_Dec2022.csv",row.names = F)

## AGRA_summarydat
dat2fixnames=dat1%>%
  full_join(dat2)%>%
  mutate(prevname=Species)%>%
  select(-Species)
dat2fixnames$obsID=c(1:length(dat2fixnames$prevname))
source("Seed Trait Paper/Scripts/2_name_cleaning_tosource_spcode.R")
head(cleannamedat)
AGRA_summarydat=cleannamedat%>%
  select(-obsID)%>%
  distinct()
head(AGRA_summarydat)
# write.csv(AGRA_summarydat,"Seed Trait Paper/current_csv files/AGRA_summarydat_Dec2022.csv",row.names = F)


## Combined AGRAdata
comboAGRA=AGRA_2010dat%>%
  full_join(AGRA_2017dat)

comboAGBG=comboAGRA%>%
  full_join(SBRA_traitdat)

## quick check for duplicates
temp=comboAGBG%>%
  group_by(Plot,clean_code,IDnum)%>%
  summarise(numobs=n(),nvalues=sum(is.na(truecover)==F))%>%
  distinct()

temp=temp%>%filter(nvalues>1)
unique(temp$clean_code) 
unique(temp$Plot)

temp2=temp%>%
  inner_join(comboAGBG)%>%
  select(-AGRA_2010)%>%
  distinct()

temp3=temp2%>%anti_join(temp)

dat%>%
  filter(Plot=="107")

## issue is that in 2010 data there is Tribif, Trigra, and Tribifgra. We combine these into one category= Trigra
# one in 2017

# Fix
temp=comboAGBG%>%
  group_by(Plot,clean_code,IDnum)%>%
  summarise(numobs=n(),nvalues=sum(is.na(truecover)==F))%>%
  filter(nvalues>1)%>%
  distinct()

OGdat=comboAGBG%>%
  anti_join(temp)

tofixdat=comboAGBG%>%
  inner_join(temp)%>%
  group_by(all_of(c(-starts_with("ABRA"))))%>%
  summarise(nrows=n())
