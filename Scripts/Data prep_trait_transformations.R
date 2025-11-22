### 
### this script transforms and formats traits for PCA and CWM analyses 
rm(list = ls()) # clears everything

#load libraries
library(tidyverse)
library(modelsummary)


###load data
traitdat=read.csv("Data/traitdata.csv",header=T,strip.white = T) 
summary(traitdat)
dim(traitdat)

#### impute texture and Starch data for species missing measurements ####
a=median(traitdat$Texture,na.rm = T)
b=median(traitdat$Intensity,na.rm = T)

sum(is.na(traitdat$Texture)) # 7 sp without texture

missingsp = traitdat%>%
  filter(is.na(Texture)==T)
missingsp
unique(missingsp$species)

traitdat%>%
  filter(is.na(Starch))# last 2 values to uptraitdate (naspul and micdou)

traitdat=traitdat%>%
  mutate(Texture=ifelse(is.na(Texture) == T,a,Texture),Intensity=ifelse(is.na(Intensity)==T,b,Intensity),
         Starch=ifelse(species =="NASPUL",1,ifelse(species=="MICDOU",0,Starch)))%>%
  filter(is.na(Starch)==F) 

summary(traitdat)


##Transform trait data ##
traitdat_trans=traitdat%>%
  mutate(Compact=exp(Compact),
         Texture=log(Texture),
         Shape=log(Shape),
         Perimeter=log(Perimeter), 
         area=log(area),
         Length=log(Length), 
         SCP=log(SCP),
         Mass=log(Mass), 
         SCT=sqrt(SCT)) %>%
  select(-nseedssp, -nfam)


datasummary_skim(traitdat_trans)

names(traitdat_trans)

#save transformed trait data for analyses
#write.csv(traitdat_trans, "Data/transformed_traitdata2.csv")


