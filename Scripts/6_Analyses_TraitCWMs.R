#This script 
#Calculates community weighted mean (CWM) traits by plot and treatment
#tests whether CWM traits vary by treatment and habitat (harsh vs. lush serpentine)
rm(list = ls()) # clears console

#load libraries
library(vegan)
library(tidyverse)
library(modelsummary)
library(nlme)
library(MuMIn)
library(modelsummary)
library(kableExtra) #tables of model output
library(emmeans) #post hoc contrasts
library(ggplot2)
library(cowplot) #for arranging panels in figures
library(multcomp) #for getting letters for post hoc contrasts on figures

traitdat = read.csv("Output data/Traits_PCscores_all.csv") #see 5_Analyses_TraitPCA.R for script generating this file
#note traitdat is using the transformed trait values as they were input into the PCA in script 5
summary(traitdat)
str(traitdat)
dim(traitdat)

abundat = read.csv("Cleaned data/SBRAWbyspecies_tomatch2017_cleancode.csv") %>% #see 3_Final_namecleaning.R for script generating this file, this matches raw seedling counts from Eskelinen et al. 2021, seedbankgrowoutallyrs_wide_totnumseedlings_tomatch2017.csv
          mutate(species = clean_code, rawcount = relcover) #despite naming in file, it's actually raw counts not relative cover
         
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
               mutate(relab_traitsonly = rawcount/totplottraits) %>%
               mutate(habitat = fct_recode(habitat, "Harsh" = "Harshserp", "Lush" = "Lushserp")) %>%
               mutate(WFtreatment = fct_recode(WFtreatment, "W" = "W", "N" = "F", "WN" = "FW"))
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
relabsummary = all_abtraits %>% 
          filter(is.na(Mass)==F & is.na(CN)==F)  %>% 
          group_by(Plot) %>%
          summarize(sumrelab = sum(relab, na.rm=T)) #median = 90%, mean = 79.9%
sd(relabsummary$sumrelab)

#save trait and abundance data 
#write.csv(all_abtraits, "Cleaned data/Transformed trait data_with abundance.csv")

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
         mutate(watering=as.factor(ifelse(watering=="unwatered","none","watered")),
         fertilization=as.factor(ifelse(fertilization=="unfertilized","none","fertilized"))) %>%
         mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

summary(cwmdat)        
dim(cwmdat) #90 plots
str(cwmdat)
#write.csv(cwmdat, "Output data/CWMs_traits_PCs.csv")

#### Model outputs for individual traits ####  
indtraitlist = c("SCT_CWM","Mass_CWM" ,"SCP_CWM", "CN_CWM" ,"Length_CWM", "Starch_CWM", "Shape_CWM")
length(indtraitlist)
output_ind = list() 


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
               mutate(factor = rownames(.), trait = indtraitlist[i])
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
      mutate(factor = rownames(.), trait = indtraitlist[i])
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
              mutate(factor = rownames(.), trait = indtraitlist[i])
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
            output_ind[[i]]= anova_output2
            write.csv(anova_output2, file = paste("Output tables/",response,"_main_both.csv",sep="")) 
            kableExtra::save_kable(anova_output2_table, file = paste("Output tables/",response,"_main_both.html",sep=""))      
            } else {
            #output results for 2 ways if they are significant 
              output_ind[[i]]= anova_output1
               write.csv(anova_output1, file = paste("Output tables/",response,"_2way_both.csv",sep=""))
            kableExtra::save_kable(anova_output1_table, file = paste("Output tables/",response,"_2way_both.html",sep="")) }     
  } else { #output results for 3 way if it is significant 
    output_ind[[i]]= anova_output0
    write.csv(anova_output0, file = paste("Output tables/",response,"_3way_both.csv",sep=""))
    kableExtra::save_kable(anova_output0_table, file = paste("Output tables/",response,"_3way_both.html",sep=""))
  }
   
 }

#pull anovas together for table
all_anova_indtraits = do.call(rbind.data.frame, output_ind) %>%
                      rename(Fval = 'F-value', p= 'p-value') %>%
                      mutate(p = ifelse(round(p,3) == 0, "<0.001",round(p, 3))) %>%
                      mutate(Fval = round(Fval, 3))
#write.csv(all_anova_indtraits, file = "Output tables/CWM_indtraits_ANOVA_all.csv")

#Better formatting for table 2
all_anova_indtraits_wider = all_anova_indtraits %>%
              mutate(DF = paste(numDF, denDF, sep = ","), 
                     FP= paste(Fval,p, sep = ", ")) %>%
              select(-numDF, -denDF, -Fval, -p) %>%
              pivot_wider(
                id_cols = c(trait),
                names_from = factor, 
                values_from = c(DF, FP)
              ) %>%
            rename(DF = DF_habitat) %>%
            select(!starts_with("DF_")) %>%
            rename_with(~str_remove(.x, "FP_")) 
#Nicer tohave sub columns for F and p, so pivot longer again
all_anova_indtraits_wider2 = all_anova_indtraits %>%
        mutate(Fval = as.character(Fval)) %>%
        mutate(DF = paste(numDF, denDF, sep = ",")) %>%
        select(-numDF, -denDF) %>%
        pivot_wider(
        id_cols = c(trait),
        names_from = factor, 
        values_from = c(DF, Fval, p)
      ) %>%
      rename(DF = DF_habitat) %>%
      select(!starts_with("DF_")) %>%
     pivot_longer(
       cols = c(starts_with("Fval_"), starts_with("p_")),
       names_to= "statfact",
       values_to = "value"
     ) %>%
   separate(statfact, sep = "_", into = c("stat", "factor")) %>% #then go wider again, this seems silly
  pivot_wider(
   id_cols = c(trait, DF, stat),
   names_from= factor,
    values_from = value
  ) 
 
          
#write.csv(all_anova_indtraits_wider2, file = "Output tables/CWM_indtraits_ANOVA_all_wide.csv")


#### Model outputs for CWM traits calculated from PC scores (multivariate traits) ####
pctraitlist = c("PC1_CWM","PC2_CWM" ,"PC3_CWM", "PC4_CWM")
length(pctraitlist)
output_pc = list() 

for (i in 1:length(pctraitlist)){
  modeldat=cwmdat
  dataversion="all plots" # alternative is to use testdatPC which only includes plots with a minimum of 70% cover
  response= pctraitlist[i] 
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
    mutate(factor = rownames(.), trait = pctraitlist[i])
  anova_output0_table = anova_output0 %>%
    select(-factor) %>%
    kable(
      digits = 3,
      caption = "ANOVA for Fixed Effects",
      col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "Fval-value", "P-value")
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
      mutate(factor = rownames(.), trait = pctraitlist[i])
    anova_output1_table = anova_output1%>%
      select(-factor) %>%
      kable(
        digits = 3,
        caption = "ANOVA for Fixed Effects",
        col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "Fval-value", "P-value")
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover", "condensed"),
        full_width = FALSE
      )
    if(round(anova_output1$'p-value'[5],3) >0.05 & round(anova_output1$'p-value'[6],3) >0.05 & round(anova_output1$'p-value'[7],3) >0.05){ #condition 2
      m2=lme((responsevar)~ habitat+ watering+ fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
      anova_output2 = as.data.frame(anova(m2))  %>%
        mutate(factor = rownames(.), trait = pctraitlist[i])
      anova_output2_table = anova_output2 %>%
        select(-factor) %>%
        kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "Fval-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
      #output results for main factors only 
      output_pc[[i]]= anova_output2
      write.csv(anova_output2, file = paste("Output tables/",response,"_main_both.csv",sep="")) 
      kableExtra::save_kable(anova_output2_table, file = paste("Output tables/",response,"_main_both.html",sep=""))      
    } else {
      #output results for 2 ways if they are significant 
      output_pc[[i]]= anova_output1
      write.csv(anova_output1, file = paste("Output tables/",response,"_2way_both.csv",sep="")) 
      kableExtra::save_kable(anova_output1_table, file = paste("Output tables/",response,"_2way_both.html",sep="")) }     
  } else { #output results for 3 way if it is significant 
    output_pc[[i]]= anova_output0
    write.csv(anova_output0, file = paste("Output tables/",response,"_2way_both.csv",sep="")) 
    kableExtra::save_kable(anova_output0_table, file = paste("Output tables/",response,"_3way_both.html",sep=""))
  }
  
}

#pull anovas together for table
all_anova_pctraits = do.call(rbind.data.frame, output_pc) %>%
  rename(Fval = 'F-value', p= 'p-value') %>%
  mutate(p = ifelse(round(p,3) == 0, "<0.001",round(p, 3))) %>%
  mutate(Fval = round(Fval, 3))
#write.csv(all_anova_pctraits, file = "Output tables/CWM_pctraits_ANOVA_all.csv")

#Better formatting for alternate table 2
all_anova_pctraits_wider = all_anova_pctraits %>%
  mutate(DF = paste(numDF, denDF, sep = ","), 
         FP= paste(Fval,p, sep = ", ")) %>%
  select(-numDF, -denDF, -Fval, -p) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, FP)
  ) %>%
  rename(DF = DF_habitat) %>%
  select(!starts_with("DF_")) %>%
  rename_with(~str_remove(.x, "FP_")) 
#Nicer to have sub columns for F and p, so pivot longer again
all_anova_pctraits_wider2 = all_anova_pctraits %>%
  mutate(Fval = as.character(Fval)) %>%
  mutate(DF = paste(numDF, denDF, sep = ",")) %>%
  select(-numDF, -denDF) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, Fval, p)
  ) %>%
  rename(DF = DF_habitat) %>%
  select(!starts_with("DF_")) %>%
  pivot_longer(
    cols = c(starts_with("Fval_"), starts_with("p_")),
    names_to= "statfact",
    values_to = "value"
  ) %>%
  separate(statfact, sep = "_", into = c("stat", "factor")) %>% #then go wider again, this seems silly
  pivot_wider(
    id_cols = c(trait, DF, stat),
    names_from= factor,
    values_from = value
  ) 


#write.csv(all_anova_pctraits_wider2, file = "Output tables/CWM_pctraits_ANOVA_all_wide.csv")

#### Post-hoc contrasts ####

##SCT ##
#3 way is sig
SCT_lm =lme(SCT_CWM ~ habitat*watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(SCT_lm)

SCT_emm = emmeans(SCT_lm,  ~ habitat*watering*fertilization) 

cld_SCT = as.data.frame(cld(SCT_emm, 
                  adjust = "Tukey",     # p-value adjustment
                  Letters = letters,    # Specify letters to use
                  alpha = 0.05,         # Significance level
                  reversed = TRUE) )  %>%    # Sort means in decreasing order
          rename(SCT_letters = .group)  %>% # Rename the letters column
          mutate(WFtreatment = as.factor(case_when(
                watering == "none" & fertilization == "none" ~ "C",
                watering == "watered" & fertilization == "none" ~ "W",
                watering == "none" & fertilization == "fertilized" ~ "N",
                watering == "watered" & fertilization == "fertilized" ~ "WN"
          ))) %>%
      mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))
          
## Mass ##
#2 way: hab x fert is sig #
Mass_lm=lme(Mass_CWM ~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(Mass_lm)

Mass_emm = emmeans(Mass_lm,  ~ habitat*watering*fertilization) #keep all factors in here to generate letters

cld_Mass = as.data.frame(cld(Mass_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(Mass_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

## SCP ##
#2 way: hab x watering is sig #
SCP_lm=lme(SCP_CWM ~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(SCP_lm) 

SCP_emm = emmeans(SCP_lm,  ~ habitat*watering*fertilization) 

cld_SCP = as.data.frame(cld(SCP_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(SCP_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

## CN ##
#Only main effects of habitat and fertilization are sig #
CN_lm =lme(CN_CWM~ habitat+ watering+ fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(CN_lm) 

CN_emm = emmeans(CN_lm,  ~ habitat*watering*fertilization) 

cld_CN = as.data.frame(cld(CN_emm, 
                               adjust = "Tukey",     # p-value adjustment
                               Letters = letters,    # Specify letters to use
                               alpha = 0.05,         # Significance level
                               reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(CN_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

## Length ## 
#2 way: hab x fert is sig #
Length_lm=lme(Length_CWM~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(Length_lm) 

Length_emm = emmeans(Length_lm,  ~ habitat*watering*fertilization) 

cld_Length = as.data.frame(cld(Length_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(Length_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

## Starch ## 
#2 ways: hab x fert and hab x watering are sig #
Starch_lm =lme(Starch_CWM ~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(Starch_lm) 

Starch_emm = emmeans(Starch_lm,  ~ habitat*watering*fertilization) 

cld_Starch = as.data.frame(cld(Starch_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(Starch_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))
## Shape ##
#habitat x fert is sig, habitat x watering is marginally sig#
Shape_lm =lme(Shape_CWM ~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(Shape_lm) 

Shape_emm = emmeans(Shape_lm,  ~ habitat*watering*fertilization) 

cld_Shape = as.data.frame(cld(Shape_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(Shape_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))


#### Figures ####
SCTplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=SCT_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM SCT", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
         labs(fill = "Habitat")

a= SCTplot + geom_text(data = cld_SCT, aes( x = WFtreatment_order, y = 2.1, group = habitat, label = SCT_letters ),
                    position = position_dodge(width = 0.9)) 
a
a_alt = a + facet_grid(~habitat)
a_alt

Massplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=Mass_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Mass", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

b= Massplot + geom_text(data = cld_Mass, aes( x = WFtreatment_order, y = 2, group = habitat, label = Mass_letters ),
                    position = position_dodge(width = 0.9))
b_alt = b + facet_grid(~habitat)
b_alt

SCPplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=SCP_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM SCP", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

c= SCPplot + geom_text(data = cld_SCP, aes( x = WFtreatment_order, y = 2.1, group = habitat, label = SCP_letters ),
                    position = position_dodge(width = 0.9))
c_alt = c + facet_grid(~habitat)
c_alt

CNplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=CN_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM CN", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

d= CNplot + geom_text(data = cld_CN, aes( x = WFtreatment_order, y = 26, group = habitat, label = CN_letters ),
                    position = position_dodge(width = 0.9))
d_alt = d + facet_grid(~habitat)
d_alt

Lengthplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=Length_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Length", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

e =Lengthplot + geom_text(data = cld_Length, aes( x = WFtreatment_order, y = 2.5, group = habitat, label = Length_letters ),
                    position = position_dodge(width = 0.9))
e
e_alt = e + facet_grid(~habitat)
e_alt

Starchplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=Starch_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Starch", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

f= Starchplot + geom_text(data = cld_Starch, aes( x = WFtreatment_order, y = 1.05, group = habitat, label = Starch_letters ),
                    position = position_dodge(width = 0.9))
f
f_alt = f + facet_grid(~habitat)
f_alt

Shapeplot = ggplot(data=cwmdat,aes(x=WFtreatment_order, y=Shape_CWM, fill = habitat))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Shape", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(habitat),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "Habitat")

g= Shapeplot + geom_text(data = cld_Shape, aes( x = WFtreatment_order, y = -1.65, group = habitat, label = Shape_letters ),
                    position = position_dodge(width = 0.9))
g
g_alt = g + facet_grid(~habitat)
g_alt
### save figure 2: CWM traits ###
#need to drop one of the graphs to make it even...
plot_grid(a + theme(legend.position = "none"),
          b+ theme(legend.position = "none"),
          c+ theme(legend.position = "none"),
          d+ theme(legend.position = "none"),
          e+ theme(legend.position = "none"), #ECE didn't have length in her graph, could add/remove here
          f, 
          g,
          ncol = 2, nrow = 4,
          labels = c("A.", "B.", "C.", "D.", "E.", "F.", "G."), label_size=14)

#ggsave("Plots/Fig2_CWMresponses.jpg", height = 10, width = 15)
e
#ggsave("Plots/Fig2_CWMresponses_length.jpg", height = 5, width = 8)

#next make alternate faceted panel
plot_grid(a_alt + theme(legend.position = "none"),
          b_alt+ theme(legend.position = "none"),
          c_alt+ theme(legend.position = "none"),
          d_alt+ theme(legend.position = "none"),
          e_alt + theme(legend.position = "none"), #ECE didn't have length in her graph, could add/remove here
          f_alt, 
          g_alt,
          ncol = 2,
          labels = c("A.", "B.", "C.", "D.", "E.", "F.", "G."), label_size=14)
#ggsave("Plots/Fig2_CWMresponses_faceted.jpg", height = 10, width = 10)
