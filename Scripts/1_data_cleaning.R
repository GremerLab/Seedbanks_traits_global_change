rm(list = ls()) # clears everything
library(tidyverse)

#### File tracking ####
## IN 
# dry_seed_coat_permeabilityALL.csv
# list.files(pattern = "Wet seed coat", path = "Seed Trait Paper/current_csv files/"
# seedmass_20200228.csv
# CN.csv
# Seed_CN_Marina.csv
# 20210610_Marina_SpeciesList.csv
# seedbankgrowoutallyrs_wide_relabun_LUMPunknowns_tomatch2017.csv

## OUT
# SCP : sead coat permiability
# SCP_summary
# MF_sm: maternal family seed mass
# MF_sm_summary
# cn_summary_full: carbon/nitrogen ratio ELISE + Marina's samples (MCL species only)
# SBRA_tomatch2017 : seed bank relative abundance - including 0's i.e. to match aboveground 2017 data

#### read in, clean, and combine csv for seed coat permiability ####
SCP_d=read.csv("Raw data/dry_seed_coat_permeabilityALL.csv",header=T, strip.white = T)

#remove extra columns
SCP_d=SCP_d[1:5]
names(SCP_d)[1]<-"IDnum"
names(SCP_d)
SCP_d=SCP_d%>%select(-Notes)

#read in multiple wet SCP files
library(readxl)
files = list.files(pattern = "Wet seed coat", path = "Raw data/", full.names = TRUE); files

SCP_w =   map_df(files, ~ read_csv(files))%>%
  distinct()

str(SCP_w)
SCP_w1=SCP_w%>%select(IDnum="sp #",number_seeds,wet_mass,Rep="Letter_group",Notes)

SCP_w1=SCP_w1%>%
  filter(is.na(wet_mass)==F)

SCP=SCP_d%>%
  full_join(SCP_w1)%>%
  mutate(repdryperseed=dry.mass/number_seeds,
         repwetperseed=wet_mass/number_seeds,
         rep_pctch=wet_mass/dry.mass,
         fulldata=ifelse(number_seeds>0,1,0))%>%
  distinct()

SCP_summary=SCP%>%
  group_by(IDnum)%>%
  summarise(          spdryperseed=mean(repdryperseed,na.rm=T),
                        spwetperseed=mean(repwetperseed,na.rm=T),
                        spPerm_pctch=mean(rep_pctch,na.rm=T),
                      sd_spwetperseed=sd(repwetperseed,na.rm = T),
                      sd_spdryperseed=sd(repdryperseed,na.rm = T),
                      sd_spPerm_pctch=sd(rep_pctch,na.rm = T),
                      SS_SCP=sum(fulldata,na.rm=T))
head(SCP_summary)

### Seed mass ####
dat=read.csv("Raw data/seedmass_20200228.csv",header=T)

names(dat)
unique(dat$species_code)

# fix NA's
dat=dat%>%
  mutate(species_code=ifelse(is.na(species_code),species,species_code))

dat1=dat%>%
  select(species_code,total_seed_mass_g,nseeds=Total._seeds_weighed,species)%>%
  mutate(seed_mass=total_seed_mass_g/nseeds)

table(dat1$species_code)

MF_sm_summary=dat1%>%
  group_by(species_code,species)%>%
  summarise(nseedssp=sum(nseeds),nfam=n(),
            spdrymf=mean(seed_mass,na.rm=T),
            sd_spdrymf=sd(seed_mass,na.rm=T))

MF_sm=dat1%>%
  select(species_code,seed_mass)%>%
  full_join(MF_sm_summary)

#rm(dat,dat1)

## CN data cleaning (Elise's and Marina's samples) ####
cndat=read.csv("Raw data/CN.csv",header=T)
names(cndat)

cndat1=cndat%>%
  select(C=Total.C..µg.,N=Total.N..µg.,ID=Sample.ID, cn_samplemass=Sample.Weight..mg..from.Sample.List)%>%
  filter(is.na(C)==F)%>%
  mutate(newID=gsub("-1","",ID))%>%
  mutate(newID=gsub("-light","",newID))%>%
  mutate(newID=gsub("-2","",newID))%>%
  mutate(newID=gsub("-dark","",newID))

cn_summary=cndat1%>%
  group_by(species_code=newID)%>%
  summarise(C=mean(C),N=mean(N),cn_samplemass=mean(cn_samplemass,na.rm=T),nsamples=n())%>%
  mutate(cn=C/N)

## Marina's seed CN information
datMcn=read.csv("Raw data/Seed_CN_Marina.csv", header=T, strip.white = T)

datMcn=datMcn%>%
  select(ID=Sample.ID,C=C.ug, N=N.ug,cn_samplemass=sample.wgt.mg)%>%
  mutate(cn=C/N)

#Marina's metadata
datM=read.csv("Raw data/20210610_Marina_SpeciesList.csv",header = T,strip.white = T)

#select seeds only from Mclaughlin
datM=datM%>%
  select(Species,species.alt,code,accept.code,site,ID,type)%>%
  filter(type!="cultivated" & site=="MCL-ML")%>%
  distinct()

names(datMcn)
sum(datMcn$ID %in% datM$ID) #10 species from Marina's CN values are from McLaughlin
datM$ID

#pull out cn values from MCL (samples Marina ran for Elise)
datMcn2=datMcn%>%
  filter(str_detect(ID, "MCL") == TRUE |str_detect(ID, "EE"))

# create code without MCL for combining
datMcn2=datMcn2%>%
  left_join(datM)%>%
  select(-type,-site,-accept.code)%>%
  mutate(code=ifelse(is.na(code),substr(ID,4,9),code))

#FIX NAMES
datMcn2$code2 <- recode(datMcn2$code, "LOGFIL" = "FILCAL",   
                            "ZELTRI" = "CENTRI", 
                            "ACMWRA" = "LOTWRA") 

# Elise ran 2 samples from previous CN batch see just below
a=toupper(cn_summary$species_code)
datMcn2$code2 %in% a 

# Join Elise's and Marina's CN data, the 2 overlapping species are in both. Don't need to include these.
# must adjust names to be in the same format - first letter upper, rest lower
cn_summary_full=datMcn2%>%
  mutate(nsamples=1)%>%
  filter(code != "CHLPOM" & code !="LINDIC")%>%
  select(code2, C,N,cn,cn_samplemass,nsamples)%>%
  mutate(species_code=paste(substr(code2,1,1),substr(tolower(code2),2,length(code2)),sep="")) %>%
  full_join(cn_summary)%>%
  select(-code2)
 
head(cn_summary_full)

#rm(cndat,datM,datMcn,datMcn2, cn_summary)
## quick CN comparison for accuracy ####
# 2 overlaping with marina LINDIC CHLPOM
cndat1%>%filter(ID !=newID)%>%mutate(cn=C/N)%>%arrange(newID)
## the precision varies with species but generally +/- 3
# colspa ~0 diff
# Avefat, Galmur ~0.5 diff
# Broele, Hemcon ~ 3+

#CHLPOM = 2 in original (15.5,13.5) vs Marina 10.5 ( 4 pts lower)
#LINDIC = 1 sample in original (17.2) vs Marina 17.5  (SAME)
## These are close enough to ignore - I don't have enough data to really quantify any differences.

# NOTE CENSOL light and dark are very different!dark is much higher (far LESS nitrogen) than light seed

#### Relative Abundance of seedbank dataprep ####

# dat1=read.csv("Anu Paper/datafiles/seedbankgrowoutallyrs_wide_relabun_LUMPunknowns_tomatch2017.csv") # 90 by 129 lumped unknowns
# dat1=read.csv("Anu Paper/datafiles/seedbankgrowoutallyrs_wide__totnumseedlings_tomatch2017_withlumpedunknowns.csv") # 90 by 129 lumped unknowns # This is same as relative abundance data. switch to raw abundance values April 2023
# dat1=read.csv("Anu Paper/datafiles/seedbankgrowoutallyrs_wide_abovefirst_totnumseedlings.csv") # 90 by 147 unknowns, has NAS
dat1=read.csv("Raw data/seedbankgrowoutallyrs_wide_totnumseedlings_tomatch2017.csv") # 90 by 147 unknowns
# has zeros
names(dat1)


## using dat1: seed bank rel abundances with unknowns lumped - change to long format.
names(dat1)
relcover3=dat1[,11:length(dat1)]%>%
  # select(!contains("UNK")) %>%
  pivot_longer(!Plot, values_to = "relcover",names_to ="species_code",values_drop_na = T)

SBRA_tomatch2017=dat1%>%
  select("Line","habitat","watering","fertilization","WFtreatment","Plot")%>%
  full_join(relcover3)#%>%
  # filter(relcover>0)  For now keep all rows - even 0's to preserve empty species. 

# write.csv(SBRA_tomatch2017,"Seed Trait Paper/current_csv files/cleandat_SBbyspecies_tomatch2017.csv",row.names = F)
head(SBRA_tomatch2017)
unique(SBRA_tomatch2017$species_code)

rm(dat1,relcover3)

### write out files to clean_final_subdata folder ####

# write.csv(SCP,"Seed Trait Paper/clean_final_subdata/SCP_120222.csv",row.names = F)
# write.csv(SCP_summary,"Seed Trait Paper/clean_final_subdata/SCP_summary_120222.csv",row.names = F)
# write.csv(MF_sm,"Seed Trait Paper/clean_final_subdata/mf_sm_120222.csv",row.names = F)
# write.csv(MF_sm_summary,"Seed Trait Paper/clean_final_subdata/mf_sm_summary_120222.csv",row.names = F)
# write.csv(cn_summary_full,"Seed Trait Paper/clean_final_subdata/cn_summary_120222.csv",row.names = F)
# write.csv(SBRA_tomatch2017,"Seed Trait Paper/clean_final_subdata/SBRAbyspecies_tomatch2017_120222.csv",row.names = F)

# write.csv(SBRA_tomatch2017,"Seed Trait Paper/clean_final_subdata/SBRAWbyspecies_tomatch2017_05112023.csv",row.names = F)

