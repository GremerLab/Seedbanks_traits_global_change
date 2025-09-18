#edit origin = native vs non-native once have all files

rm(list = ls()) # clears everything
library(tidyverse)

#### File tracking ####
## IN (in order of Table 1)
# SCT_microscopeandzoom.csv
# dry_seed_coat_permeabilityALL.csv
# list.files(pattern = "Wet seed coat", path = "Seed Trait Paper/current_csv files/"
# Starch_20221222.csv  ###MISSING FILE ####
# CN.csv, Seed_CN_Marina.csv, 20210610_Marina_SpeciesList.csv
# seedmass_20200228.csv
# Files from videometer: list.files(pattern = "_blobs_", path = "Raw data/VMExports/", full.names = TRUE)
# seedtraitphotos_categories_Dec2022_v1.csv #### missing file ####
# seedbankgrowoutallyrs_wide_relabun_LUMPunknowns_tomatch2017.csv

## OUT
# Seed coat thickness: SCT_micrometer.csv
#Seed coat permeability: SCP_120222.csv, SCP_summary_120222.csv
# Starch: Starch_20221222.csv 
# Seed mass: MF_sm: maternal family seed mass - mf_sm_120222.csv
# Seed mass summary: mf_sm_summary_120222.csv
# cn_summary_full: carbon/nitrogen ratio ELISE + Marina's samples (MCL species only): cn_summary_120222.csv
# SBRA_tomatch2017 : seed bank relative abundance - including 0's i.e. to match aboveground 2017 data

#### read in and clean SCT data ####
SCTdat=read.csv("Raw data/SCT_microscopeandzoom.csv", header=T)
names(SCTdat)
head(SCTdat)
unique(SCTdat$clean_code)
doubles=c("double_TRIALB", "double_PETPRO","double_EPIBRA")
## remove doubled species and zoom columns
SCTdat=SCTdat%>%
  select(-starts_with("Photo"),-Cross_zoom,-Ext_Zoom)%>%
  filter(clean_code %in% doubles ==F)

names(SCTdat)
# don't need to fix names
str(SCTdat)

SCTdat1=SCTdat%>%
  mutate(mulongonly=rowMeans(subset(SCTdat,select = c(long_1,long_2,long_3,long_4,long_5)), na.rm = TRUE),
         muSCT=rowMeans(subset(SCTdat,select=c(cross_1,long_1,long_2,long_3,long_4,long_5)),na.rm = T))%>%
  group_by(IDnum,clean_code)%>%
  dplyr::summarise(muSCT=mean(muSCT,na.rm = T),SCTlongonly=mean(mulongonly,na.rm = T))

SCTdat2=SCTdat1%>%select(-muSCT)
# write.csv(SCTdat2,"Cleaned data/SCT_micrometer.csv",row.names = F)


#### read in, clean, and combine csv for seed coat permeability ####
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

# write.csv(SCP,"Cleaned data/SCP.csv",row.names = F)
# write.csv(SCP_summary,"Cleaned data/SCP_summary.csv",row.names = F)

#### Starch content ####
Starch=read.csv("Raw data/Starch.csv", header=T, strip.white = T) #missing file, need to get from Elise

Starch=Starch%>%select(IDnum="Idnum",starch,mucilage)

#write.csv(Starch, "Cleaned data/Starch.csv")

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

# write.csv(cn_summary_full,"Cleaned data/cn_summary.csv",row.names = F)


### Seed mass ####
massdat=read.csv("Raw data/seedmass_20200228.csv",header=T)

names(massdat)
unique(massdat$species_code)

# fix NA's
massdat=massdat%>%
  mutate(species_code=ifelse(is.na(species_code),species,species_code))

dat1=massdat%>%
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

#rm(massdat,massdat1)
# write.csv(MF_sm,"Cleaned data/mf_sm.csv",row.names = F)
# write.csv(MF_sm_summary,"Cleaned data/mf_sm_summary.csv",row.names = F)


#### Length, 3D shape, compactness, seed texture from Videometer ####
#read in multiple files
library(readxl)
files = list.files(pattern = "_blobs_", path = "Raw data/VMExports/", full.names = TRUE); files

# sizes = map_df(files, ~read_xlsx(., na = "NA", col_types = c("numeric", "text", "text", "text", "text", "text", "numeric", "numeric", "text"))) 
VM_1 =   map_df(files, ~ read_xlsx(.,na="NA"))
names(VM_1)

## Select species name and apply correct name information
dat2fixnames=VM_1%>%
  select(Org.Filename)%>%
  mutate(prevname=paste(substr(Org.Filename,1,1),tolower(substr(Org.Filename,2,6)),sep=""))%>%
  distinct()

length(unique(dat2fixnames$prevname))
dat2fixnames

source("Scripts/2_name_cleaning_tosource_spcode.R")

# combine data with clean names
VM_1=cleannamedat%>%
  select(IDnum,clean_code,Org.Filename)%>%
  right_join(VM_1)

# combine data with clean names
VM_1=cleannamedat%>%
  select(IDnum,clean_code,Org.Filename)%>%
  right_join(VM_1)

## remove unnecessary columns
names(VM_1)
VMdat1=VM_1%>%
  mutate(texture=VM_1$`AutoCorrelationEnergy1 [9]`)%>%
  select(-starts_with("AutoC"), -starts_with("Multi"),
         -starts_with("CIEL"),-starts_with("Vert"),
         -starts_with("Perp"))%>%
  distinct()
names(VMdat1)
head(VMdat1)
length(unique(VMdat1$IDnum))

## Remove repeats
# some doubles for lumpers Clarkiasp, Galiumsp
# some true doubles
VMdat1=VMdat1%>%
  filter(IDnum<72)

VMdat1=VMdat1%>%
  group_by(Org.Filename,`Area (mm2)`,Saturation)%>%
  summarise(newBlobId=mean(BlobId))%>%
  distinct()%>%
  left_join(VMdat1)%>%
  select(-BlobId)%>%
  distinct()
VMdat1$BlobId=VMdat1$newBlobId

## specifically deal with known thycur and mimdou issue where photos were with pods.
to_remove=VMdat1%>%
  filter(str_detect(Org.Filename,"MIMDOU") == T | 
           str_detect(Org.Filename,"THYCUR") == T )%>%
  filter(str_detect(Org.Filename, "seed") == F)

VMdat1=VMdat1%>%
  anti_join(to_remove)

## determine blob orderof 2 photos
VMdat1$localblob=c(1:length(VMdat1$BlobId))
blobstomatch=VMdat1%>%
  select(localblob,IDnum,Org.Filename,`Area (mm2)`,BlobId,Saturation)%>%
  group_by(Org.Filename,IDnum)%>%
  summarise(localblobMIN=min(localblob))%>%
  full_join(VMdat1)%>%
  mutate(localblobrelative=localblob-localblobMIN)%>%
  select(-localblob,-BlobId)

blobstomatch2=blobstomatch%>%
  group_by(IDnum,clean_code,localblobrelative)%>%
  summarise(avgarea=mean(`Area (mm2)`,rm.na=T),
            num=n())%>%
  inner_join(blobstomatch)%>%
  mutate(LargestArea=ifelse(`Area (mm2)`>=avgarea,1,0))

firstblobs=blobstomatch2%>%
  filter(LargestArea==1)%>%
  select(-num,-newBlobId, -LargestArea)

secondblobs=blobstomatch2%>%
  filter(LargestArea==0)%>%
  select(IDnum,clean_code,localblobrelative,Org.Filename2=Org.Filename,LH_area=`Area (mm2)`,LH_L=`Length (mm)`,H=`Width (mm)`,LH_Circle=`Compactness Circle`)

blobs=secondblobs%>%
  full_join(firstblobs)
#View(blobs)

## secondblob will have all of the second photos.
# I am defining the Length as the longest axis of the side with the most area. The width as the second axis of the side with the most area. The height is the shorter axis from the other dimension (LH). I am keeping the values from the photo with the largest area as the main data source. A few variables are kept for the height axis (area, Circleness, filename).

VMdatblobs=blobs%>%
  select(IDnum,clean_code,
         H,
         L="Length (mm)",
         LH_L,
         W="Width (mm)",
         area_LW="Area (mm2)",
         Perimeter1,Hue,PCof19WL=IHSIntensityMean1,texture,
         Rectangularity1, LH_Circle, LW_Circle=`Compactness Circle`,`Compactness Ellipse`,BetaShape_a,BetaShape_b,Saturation,LH_area,Org.Filename,Org.Filename2)%>%
  mutate(H=ifelse(is.na(H),W,H),
         Wadj=W/L,
         Ladj=1,
         Hadj=H/L,
         shape=var(c(Wadj,Hadj,Ladj)))%>% #third dimension.
  select(-Ladj,-Hadj,-Wadj)
names(VMdatblobs)

VM_summary=VMdatblobs%>%
  mutate(H=ifelse(is.na(H),W,H),
         Wadj=W/L,
         Ladj=1,
         Hadj=H/L,
         shape=var(c(Wadj,Hadj,Ladj)))%>%
  select(-Ladj,-Hadj,-Wadj)%>%
  select(where(is.numeric))%>%
  group_by(IDnum,clean_code)%>%
  summarise(across(everything(),c(mean,sd)))

# write.csv(VM_summary,"Cleaned data/VM_byspecies.csv",row.names = F)

# write.csv(VMdatblobs,"Cleaned data/VM_byblobs.csv",row.names = F)


head(VM_summary)
names(VM_summary)
quickdat=VM_summary%>%
  select(IDnum,clean_code,
         shape_1,Perimeter1_1,L_1,area_LW_1,
         sd_shape="shape_2",sd_perim="Perimeter1_2",sd_L="L_2",sdarea="area_LW_2")

# data_without_na %>%
#   select(where(is.numeric)) %>%
#   summarise(across(everything(), mean))
quickdat
# write.csv(quickdat,"Cleaned data/VM_byspecies_mainvariables.csv",row.names = F)

#### Dispersal appendages ####
Morph=read.csv("Raw data/seedtraitphotos_categories_Dec2022_v1.csv", header=T, strip.white = T) #### missing file ####
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

#write out morphological data here # 

#### Relative Abundance of seedbank dataprep ####
abundat=read.csv("Raw data/seedbankgrowoutallyrs_wide_totnumseedlings_tomatch2017.csv") # 90 by 147 unknowns
# has zeros
names(abundat)


## using abundat: seed bank rel abundances with unknowns lumped - change to long format.
names(abundat)
relcover3=abundat[,11:length(abundat)]%>%
  # select(!contains("UNK")) %>%
  pivot_longer(!Plot, values_to = "relcover",names_to ="species_code",values_drop_na = T)

SBRA_tomatch2017=abundat%>%
  select("Line","habitat","watering","fertilization","WFtreatment","Plot")%>%
  full_join(relcover3)#%>%
  # filter(relcover>0)  For now keep all rows - even 0's to preserve empty species. 

head(SBRA_tomatch2017)
unique(SBRA_tomatch2017$species_code)


# write.csv(SBRA_tomatch2017,"Cleaned data/SBRAWbyspecies_tomatch2017.csv",row.names = F)




