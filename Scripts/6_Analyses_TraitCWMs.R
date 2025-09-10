#This script 
#Calculates community weighted mean (CWM) traits by plot and treatment
#tests whether CWM traits vary by treatment and habitat (harsh vs. lush serpentine)
rm(list = ls()) # clears console

#load libraries
library(vegan)
library(tidyverse)
library(modelsummary)
library(nlme)
library(broom.mixed) #to output pvalues of lme models to table
library(MuMIn)
library(modelsummary)
library(kableExtra)
library(ggplot2)

traitdat = read.csv("Output data/Traits_PCscores_all.csv") #see 5_Analyses_TraitPCA.R for script generating this file
#note traitdat is using the transformed trait values as they were input into the PCA in script 5
summary(traitdat)
str(traitdat)
dim(traitdat)

abundat = read.csv("Cleaned data/SBRAWbyspecies_tomatch2017_05112023_cleancode.csv") %>% #see 3_Final_namecleaning.R for script generating this file, this matches raw seedling counts from Eskelinen et al. 2021, seedbankgrowoutallyrs_wide_totnumseedlings_tomatch2017.csv
          mutate(species = clean_code, rawcount = relcover) %>% #despite naming in file, it's actually raw counts not relative cover
          select(-clean_code, -family, -annual, -MarinaSp, -nat.inv, -Functional.group, -relcover) #these are updated in traitdat, can remove from here
#JRG notes, all abundances using the RAW data file are raw abundances, not relative, and has adapted scripts accordingly

summary(abundat)
str(abundat)
dim(abundat)

abund_trait_dat = left_join(abundat, traitdat, by= "species") %>% #
                  mutate(Line = as.factor(Line), habitat=as.factor(habitat), watering=as.factor(watering), 
                         fertilization=as.factor(fertilization), WFtreatment=as.factor(WFtreatment), Plot=as.factor(Plot), 
                         species = as.factor(species), origin = as.factor(origin))
dim(abund_trait_dat)
summary(abund_trait_dat)

#get relative abundances
rawab_plots=abund_trait_dat%>%
  group_by(Plot,WFtreatment)%>%
  dplyr::summarise(totplotincmissing=sum(rawcount)) #previously was done as a presence/absence (sum(as.numeric(as.character()...)))

#below calculates the total abundance for species in those plots for which we had trait data
rawab_wtraits =abund_trait_dat%>%
  filter(is.na(Mass)==F & is.na(CN)==F)%>% #removes rows that are missing seed mass or cn value
  group_by(Plot,WFtreatment)%>%
  dplyr::summarise(totplottraits=sum(rawcount))

rawab = full_join(rawab_plots, rawab_wtraits)
summary(rawab)

rawab = full_join(rawab_plots, rawab_wtraits)


all_abtraits = left_join(abund_trait_dat, rawab) %>%
               mutate(relab = rawcount/totplotincmissing) %>%
               mutate(relab_traitsonly = rawcount/totplottraits) 
#all_abtraits will have some relab_traitsonly greater than 1, this is species that don't have trait values
#so, their CWM trait values will be NA, so no need to filter them out here
summary(all_abtraits)
all_abtraits %>% filter(relab_traitsonly > 1)
all_abtraits %>% filter(is.na(Mass) == TRUE)

#double check that relative abundance is correct (should add to 1 for each plot)
summary(all_abtraits %>% 
  group_by(Plot) %>%
  summarize(sumrelab = sum(relab, na.rm=T)))

#check what coverage we have for species that we do have traits
summary(all_abtraits %>% 
          filter(is.na(Mass)==F & is.na(CN)==F)  %>% 
          group_by(Plot) %>%
          summarize(sumrelab = sum(relab, na.rm=T))) #median = 90%, mean = 79.9%
    
#### calculate CWMs ####
#Note that JRG did this slightly differently than ECE dissertation 
#first, get trait value * relative abundance for each species in each plot
#NOTE: this uses the relative abundance for the species for which we have traits.  If want "true" relative abundance, change to "relab"
cwmdat_a = all_abtraits %>%
         mutate(across(c(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                         Texture, Compact), ~.x*relab_traitsonly, .names = "{.col}_CWM")) %>%
         mutate(across(starts_with("PC"), ~.x*relab_traitsonly,.names = "{.col}_CWM") ) 

#next, sum for each plot to get a plot level CWM
cwmdat = cwmdat_a %>%
         group_by(Line, habitat, watering, fertilization, WFtreatment, Plot) %>%
         summarize(across(ends_with("_CWM"), ~sum(.x, na.rm=T)))%>%
         ungroup()%>%
         mutate(watering=ifelse(watering=="unwatered","none",watering),
         fertilization=ifelse(fertilization=="unfertilized","none",fertilization))

summary(cwmdat)        
dim(cwmdat) #90 plots
str(cwmdat)
#write.csv(cwmdat, "Output data/CWMs_traits_PCs.csv")

#### Model outputs for individual traits ####  
indtraitlist = c("SCT_CWM","Mass_CWM" ,"SCP_CWM", "CN_CWM" ,"Length_CWM", "Starch_CWM", "Shape_CWM")
length(indtraitlist)

for (i in 1:length(indtraitlist)){
  modeldat=cwmdat
  dataversion="all plots" # alternative is to use testdatPC which only includes plots with a minimum of 70% cover
  response= indtraitlist[i] 
  COL=which((response == names(modeldat))==T)
  temp=as.data.frame(modeldat[COL])
  colnames(temp)[1]="responsevar"
  modeldat=cbind(modeldat,temp)
  names(modeldat)
  
  ##  with 3 way interaction 
  m0=lme((responsevar)~ habitat*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
  #m0ml=lme((responsevar)~ habitat*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
  #JRG: we want something like this:  a= anova(m0ml)$'p-value'
 anova_output0 = as.data.frame(anova(m0)) %>%
               mutate(factor = rownames(.))
 anova_output0_table = anova_output0 %>%
          select(-factor) %>%
          kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )

  if(round(anova_output0$'p-value'[anova_output0$factor == "habitat:watering:fertilization"],3) >0.05){ #condition 1
    m1=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
    #m1ml=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
    # Create the kable table
    anova_output1 = as.data.frame(anova(m1)) %>%
      mutate(factor = rownames(.))
    anova_output1_table = anova_output1%>%
        select(-factor) %>%
          kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
          if(round(anova_output1$'p-value'[5],3) >0.05 & round(anova_output1$'p-value'[6],3) >0.05 & round(anova_output1$'p-value'[7],3) >0.05){ #condition 2
            m2=lme((responsevar)~ habitat+ watering+ fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
            anova_output2 = as.data.frame(anova(m2))%>%
              mutate(factor = rownames(.))
            anova_output2_table = anova_output2 %>%
              select(-factor) %>%
              kable(
                digits = 3,
                caption = "ANOVA for Fixed Effects",
                col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
              ) %>%
              kable_styling(
                bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE
              )
            #output results for main factors only 
            kableExtra::save_kable(anova_output2_table, file = paste("Output tables/",response,"_main_both.html",sep=""))      
            } else {
            #output results for 2 ways if they are significant 
            kableExtra::save_kable(anova_output1_table, file = paste("Output tables/",response,"_2way_both.html",sep="")) }     
  } else { #output results for 3 way if it is significant 
    kableExtra::save_kable(anova_output0_table, file = paste("Output tables/",response,"_3way_both.html",sep=""))
  }
   
 }


#### Model outputs for CWM traits calculated from PC scores (multivariate traits) ####
cwmtraitlist = c("PC1_CWM","PC2_CWM" ,"PC3_CWM", "PC4_CWM")
length(cwmtraitlist)

for (i in 1:length(cwmtraitlist)){
  modeldat=cwmdat
  dataversion="all plots" # alternative is to use testdatPC which only includes plots with a minimum of 70% cover
  response= cwmtraitlist[i] 
  COL=which((response == names(modeldat))==T)
  temp=as.data.frame(modeldat[COL])
  colnames(temp)[1]="responsevar"
  modeldat=cbind(modeldat,temp)
  names(modeldat)
  
  ##  with 3 way interaction 
  m0=lme((responsevar)~ habitat*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
  #m0ml=lme((responsevar)~ habitat*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
  #JRG: we want something like this:  a= anova(m0ml)$'p-value'
  anova_output0 = as.data.frame(anova(m0))  %>%
    mutate(factor = rownames(.))
  anova_output0_table = anova_output0 %>%
    select(-factor) %>%
    kable(
      digits = 3,
      caption = "ANOVA for Fixed Effects",
      col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
  
  if(round(anova_output0$'p-value'[anova_output0$factor == "habitat:watering:fertilization"],3) >0.05){ #condition 1
    m1=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
    #m1ml=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
    # Create the kable table
    anova_output1 = as.data.frame(anova(m1))  %>%
      mutate(factor = rownames(.))
    anova_output1_table = anova_output1%>%
      select(-factor) %>%
      kable(
        digits = 3,
        caption = "ANOVA for Fixed Effects",
        col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover", "condensed"),
        full_width = FALSE
      )
    if(round(anova_output1$'p-value'[5],3) >0.05 & round(anova_output1$'p-value'[6],3) >0.05 & round(anova_output1$'p-value'[7],3) >0.05){ #condition 2
      m2=lme((responsevar)~ habitat+ watering+ fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
      anova_output2 = as.data.frame(anova(m2))  %>%
        mutate(factor = rownames(.))
      anova_output2_table = anova_output2 %>%
        select(-factor) %>%
        kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
      #output results for main factors only 
      kableExtra::save_kable(anova_output2_table, file = paste("Output tables/",response,"_main_both.html",sep=""))      
    } else {
      #output results for 2 ways if they are significant 
      kableExtra::save_kable(anova_output1_table, file = paste("Output tables/",response,"_2way_both.html",sep="")) }     
  } else { #output results for 3 way if it is significant 
    kableExtra::save_kable(anova_output0_table, file = paste("Output tables/",response,"_3way_both.html",sep=""))
  }
  
}


#### selected post-hoc contrasts ####


#### Figures ####
