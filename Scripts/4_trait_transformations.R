### 
rm(list = ls()) # clears everything
library(vegan)
library(tidyverse)
library(ggbiplot)
library(modelsummary)
library(ggplot2)

###NOTE: this uses 
traitdat=read.csv("Cleaned data/trait_data_Apr11.csv",header=T,strip.white = T) #created in 3_Final_namecleaning.R
# traitdat=read.csv("Seed Trait Paper/clean_final_subtraitdata/trait_traitdata_Jan20.csv",header=T,strip.white = T)

# removes standard deviation of all variables
traitdat2=traitdat%>% dplyr::select(-ends_with("_2"))

names(traitdat2)
# select only columns to transform
traitdat3=traitdat2%>%
     dplyr::select(clean_code,
                SCT=SCTlongonly,
                mass=spdryperseed,
                SCP=spPerm_pctch,C,N,cn,cn_samplemass,
                L=L_1,
                W= W_1,
                area=area_LW_1,
                Fit_shape=shape_1,
                Intensity=PCof19WL_1,
                P=Perimeter1_1,
                texture_1,
                Fit_rect=Rectangularity1_1,	
                Fit_compact=Compactness.Ellipse_1,
                germ_is_disp,
                lump_dispcat,
                mucilage,
                starch,
                areaLH=LH_area_1,annual,
                family,
                Functional.group,
                nat.inv,
                spwetperseed)%>%
  mutate(C=C/cn_samplemass, N=N/cn_samplemass,
         how3d=1-ifelse(is.na(areaLH),1,areaLH/area),
         SCTmm=SCT/80,
         SCTperL=(SCTmm/L),
         SCPpc=(spwetperseed-mass)/mass,
         SCTperW = (SCTmm/W),
         SCTperMass = (SCTmm)/mass)%>%
  #adjusting dispersal trait, where 0 = none, 1 = ephemeral, and 2 = persistent structure
  mutate(germ_is_disp=ifelse(germ_is_disp=="no",1,0),
         C=C/1000) %>%
    mutate(disp2 = if_else(germ_is_disp == 0 & lump_dispcat %in% c("elong_long", "elong_short", "flat", "balloon"), 2, germ_is_disp))%>%
  dplyr::select(-cn_samplemass, - areaLH,-spwetperseed)%>%
  distinct()
datasummary_skim(traitdat3)
#

traitdat_trans=traitdat3%>%
  mutate(SA_2D=log(P/area),Fit_compact=exp(Fit_compact),
         texture_1=log(texture_1),
         Fit_shape=log(Fit_shape),
         P=log(P), area=log(area),
         L=log(L), SCP=log(SCP),SCPpc=log(SCPpc),
         mass=log(mass), SCT=sqrt(SCT),SCTperL=sqrt(SCTperL),
         SCTperW = sqrt(SCTperW), SCTperMass = sqrt(SCTperMass))
datasummary_skim(traitdat_trans)

names(traitdat_trans)
traitdat5=traitdat_trans
row.names(traitdat5)=traitdat_trans$clean_code

#### impute texture and starch data for species missing measurements ####
a=median(traitdat5$texture_1,na.rm = T)
b=median(traitdat5$Intensity,na.rm = T)

traitdat5=traitdat5%>%
  filter(is.na(L)==F)
sum(is.na(traitdat5$texture_1)) # 7 sp without texture

missingsp = traitdat5%>%
  filter(is.na(texture_1)==T)
missingsp
unique(missingsp$clean_code)

traitdat5%>%
  filter(is.na(starch))# last 2 values to uptraitdate (naspul and micdou)

traitdat5=traitdat5%>%
  mutate(texture_1=ifelse(is.na(texture_1) == T,a,texture_1),Intensity=ifelse(is.na(Intensity)==T,b,Intensity),
         starch=ifelse(clean_code =="NASPUL",1,ifelse(clean_code=="MICDOU",0,starch)))%>%
  filter(is.na(starch)==F) #%>%
 # mutate(germ_is_disp=ifelse(germ_is_disp=="no",1,0),
 #        C=C/1000)
summary(traitdat5)
#write.csv(traitdat5,"Cleaned data/transformed_traitdata.csv")


#do the same for traitdat3, untransformed traitdata 
traitdat3=traitdat3%>%
  filter(is.na(L)==F)
sum(is.na(traitdat3$texture_1)) # 7 sp without texture

traitdat3%>%
  filter(is.na(starch))# last 2 values to uptraitdate (naspul and micdou)

traitdat3=traitdat3%>%
  mutate(texture_1=ifelse(is.na(texture_1) == T,a,texture_1),Intensity=ifelse(is.na(Intensity)==T,b,Intensity),
         starch=ifelse(clean_code =="NASPUL",1,ifelse(clean_code=="MICDOU",0,starch)))%>%
  filter(is.na(starch)==F) 

#also save raw trait values (not transformed)
row.names(traitdat3) = traitdat3$clean_code
summary(traitdat3)

#write.csv(traitdat3,"Cleaned data/untransformed_traitdata.csv")
