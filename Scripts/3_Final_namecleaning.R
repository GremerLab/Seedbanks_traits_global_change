rm(list = ls()) # clears everything
library(tidyverse)

#filelog
## IN
# mf_sm_summary_120222.csv, cn_summary_120222.csv,SCP_summary_120222.csv
# SBRAWbyspecies_tomatch2017_120222.csv
# Germ_expt_Nov2022.csv

## Run within this script
# source("Seed Trait Paper/Scripts/name_cleaning_tosource_spcode.R")
# source("Seed Trait Paper/Scripts/name_cleaning_tosource_ID.R")

## out
# trait_data_Dec22.csv
# Germ_Dec22.csv
# SBRA_trait_long_Dec2022.csv

#### read in and name clean for germination and relative abundance data ####
# dat2fixnames=read.csv("Seed Trait Paper/clean_final_subdata/SBRAbyspecies_tomatch2017_120222.csv",header = T,strip.white = T)
dat2fixnames=read.csv("Raw data/SBRAWbyspecies_tomatch2017_05112023.csv",header = T,strip.white = T)
dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
dat2fixnames$species_code=NULL
dat2fixnames$obsID=c(1:length(dat2fixnames$Line))

# input dat2fixnames with column prevname for the code
 source("Scripts/2_name_cleaning_tosource_spcode.R")
 # output = cleannamedat with column clean_code
SBRAdat=cleannamedat

#write.csv(cleannamedat, "Cleaned data/SBRAWbyspecies_tomatch2017_05112023_cleancode.csv")
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
SCP=read.csv("Cleaned data/SCP_summary.csv",header=T,strip.white = T)
CN=read.csv("Cleaned data/cn_summary.csv",header=T,strip.white=T)
MFsm=read.csv("Cleaned data/mf_sm_summary.csv",header = T,strip.white = T)

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
Starch=read.csv("Cleaned data/Starch.csv", header=T, strip.white = T)

Starch=Starch%>%select(IDnum="IDnum",starch,mucilage)
traitdat2=traitdat%>%
  left_join(Starch)


## Morphological categorization
Morph=read.csv("Raw data/seedtraitphotos_categories_Dec2022_v1.csv", header=T, strip.white = T)
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
VM=read.csv("Cleaned data/VM_byspecies.csv",header=T,strip.white = T) # the full data is a lot of variables. For now using sub data.

traitdat=traitdat3%>%
  left_join(VM)


traitwithgerm=traitdat2%>%
  right_join(GermExt)%>%
  select(-species_code)

names(traitwithgerm)
names(traitdat)

### SCT data
SCT=read.csv("Cleaned data/SCT_micrometer.csv",header=T)

traitdat=traitdat%>%
  left_join(SCT)

# write.csv(traitdat,"Seed Trait Paper/clean_final_subdata/trait_data.csv",row.names = F)
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

# write.csv(SBRA_traitdat,"Cleaned data/SBRAW_trait_long.csv",row.names = F)

traitdat2%>%filter(clean_code=="MEDPOL")
SBRA_traitdat%>%filter(clean_code=="MEDPOL")
splist1=unique(SBRA_traitdat$IDnum)
splist2=unique(traitdat2$IDnum)
a=splist2 %in% splist1
data.frame(a,splist2)%>%
  filter(a==0)


