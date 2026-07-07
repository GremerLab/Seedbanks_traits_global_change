#This script 
#Calculates community weighted mean (CWM) traits by plot and treatment
#tests whether CWM traits vary by treatment and habitat (harsh vs. lush serpentine)
#rm(list = ls()) # clears console

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

traitdat = read.csv("Data/Traits_PCscores_all.csv") #see Analyses_TraitPCA.R which generates this file
#note traitdat uses the transformed trait values 
summary(traitdat)
str(traitdat)
dim(traitdat)

abundat = read.csv("Data/Seedbank_raw_abundances.csv")  
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
  dplyr::summarise(totplotincmissing=sum(rawcount)) #

#below calculates the total abundance for species in those plots for which we had trait data
rawab_wtraits =abund_trait_dat%>%
  filter(is.na(Mass)==F & is.na(PC2)==F)%>% #removes rows that are missing seed mass or PC2 value
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
          filter(is.na(Mass)==F & is.na(PC2)==F)  %>% 
          group_by(Plot) %>%
          summarize(sumrelab = sum(relab, na.rm=T)) #median = 90%, mean = 79.9%
summary(relabsummary)
sd(relabsummary$sumrelab)

#save trait and abundance data 
#write.csv(all_abtraits, "Data/Transformed trait data_with abundance.csv")

#### calculate CWMs ####
#NOTE: this uses the relative abundance for the species for which we have traits.  
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
#write.csv(cwmdat, "Data/CWMs_traits_PCs.csv")

#### Model outputs for individual traits ####  
indtraitlist = c("SCT_CWM","Mass_CWM" ,"SCP_CWM", "CN_CWM" ,"Length_CWM", "Starch_CWM", "Shape_CWM")
length(indtraitlist)
output_ind = list() 
mod_est_ind = list()

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
  
 anova_output0 = as.data.frame(anova(m0)) %>%
               mutate(factor = rownames(.), trait = indtraitlist[i])
 anova_output0_table = anova_output0 %>%
          dplyr:: select(-factor) %>%
          kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
 
   model_output0 = as.data.frame(m0$coefficients$fixed) %>%
     mutate(factor = rownames(.), trait = indtraitlist[i]) %>%
     rename(Est = "m0$coefficients$fixed")
   
  if(round(anova_output0$'p-value'[anova_output0$factor == "habitat:watering:fertilization"],3) >0.05){ #condition 1
    m1=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
    #m1ml=lme((responsevar)~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
    # Create the kable table
    anova_output1 = as.data.frame(anova(m1)) %>%
      mutate(factor = rownames(.), trait = indtraitlist[i])
    
    anova_output1_table = anova_output1%>%
      dplyr::select(-factor) %>%
          kable(
          digits = 3,
          caption = "ANOVA for Fixed Effects",
          col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
    
    model_output1 = as.data.frame(m1$coefficients$fixed) %>%
      mutate(factor = rownames(.), trait = indtraitlist[i]) %>%
      rename(Est = "m1$coefficients$fixed")
    
          if(round(anova_output1$'p-value'[5],3) >0.05 & round(anova_output1$'p-value'[6],3) >0.05 & round(anova_output1$'p-value'[7],3) >0.05){ #condition 2
            m2=lme((responsevar)~ habitat+ watering+ fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
            anova_output2 = as.data.frame(anova(m2))%>%
              mutate(factor = rownames(.), trait = indtraitlist[i])
            
            anova_output2_table = anova_output2 %>%
              dplyr::select(-factor) %>%
              kable(
                digits = 3,
                caption = "ANOVA for Fixed Effects",
                col.names = c("Term", "Sum of Squares", "Mean Squares", "Num DF", "Den DF", "F-value", "P-value")
              ) %>%
              kable_styling(
                bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE
              )
            
            model_output2 = as.data.frame(m2$coefficients$fixed) %>%
              mutate(factor = rownames(.), trait = indtraitlist[i]) %>%
              rename(Est = "m2$coefficients$fixed")
            
            #output results for main factors only 
            output_ind[[i]]= anova_output2
            mod_est_ind[[i]]= model_output2
            write.csv(anova_output2, file = paste("Output tables/",response,"_main_both.csv",sep="")) 
            write.csv(model_output2, file = paste("Output tables/",response,"_main_both_est.csv",sep="")) 
            kableExtra::save_kable(anova_output2_table, file = paste("Output tables/",response,"_main_both.html",sep=""))      
            } else {
            #output results for 2 ways if they are significant 
              output_ind[[i]]= anova_output1
              mod_est_ind[[i]]= model_output1
               write.csv(anova_output1, file = paste("Output tables/",response,"_2way_both.csv",sep=""))
               write.csv(model_output1, file = paste("Output tables/",response,"_2way_both_est.csv",sep="")) 
            kableExtra::save_kable(anova_output1_table, file = paste("Output tables/",response,"_2way_both.html",sep="")) }     
  } else { #output results for 3 way if it is significant 
    output_ind[[i]]= anova_output0
    mod_est_ind[[i]]= model_output0
    write.csv(anova_output0, file = paste("Output tables/",response,"_3way_both.csv",sep=""))
    write.csv(model_output0, file = paste("Output tables/",response,"_3way_both_est.csv",sep="")) 
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
              dplyr::select(-numDF, -denDF, -Fval, -p) %>%
              pivot_wider(
                id_cols = c(trait),
                names_from = factor, 
                values_from = c(DF, FP)
              ) %>%
            rename(DF = DF_habitat) %>%
            dplyr::select(!starts_with("DF_")) %>%
            rename_with(~str_remove(.x, "FP_")) 
#Nicer to have sub columns for F and p, so pivot longer again
all_anova_indtraits_wider2 = all_anova_indtraits %>%
        mutate(Fval = as.character(Fval)) %>%
        mutate(DF = paste(numDF, denDF, sep = ",")) %>%
        dplyr::select(-numDF, -denDF) %>%
        pivot_wider(
        id_cols = c(trait),
        names_from = factor, 
        values_from = c(DF, Fval, p)
      ) %>%
      rename(DF = DF_habitat) %>%
      dplyr::select(!starts_with("DF_")) %>%
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


#put model estimates together and match with anovas
all_est_indtraits = do.call(rbind.data.frame, mod_est_ind)

all_anova_est = left_join(all_anova_indtraits,all_est_indtraits)

#format estimates
all_anova_indtraits_wider = all_est_indtraits %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = Est
  )  

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
    dplyr::select(-factor) %>%
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
      dplyr::select(-factor) %>%
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
        dplyr::select(-factor) %>%
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
  dplyr::select(-numDF, -denDF, -Fval, -p) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, FP)
  ) %>%
  rename(DF = DF_habitat) %>%
  dplyr::select(!starts_with("DF_")) %>%
  rename_with(~str_remove(.x, "FP_")) 
#Nicer to have sub columns for F and p, so pivot longer again
all_anova_pctraits_wider2 = all_anova_pctraits %>%
  mutate(Fval = as.character(Fval)) %>%
  mutate(DF = paste(numDF, denDF, sep = ",")) %>%
  dplyr::select(-numDF, -denDF) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, Fval, p)
  ) %>%
  rename(DF = DF_habitat) %>%
  dplyr::select(!starts_with("DF_")) %>%
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
#calculate trait means and standard errors
# standard error function 
std_error <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

cwmsummaries = cwmdat %>%
               group_by(habitat, watering, fertilization, WFtreatment, WFtreatment_order) %>%
               summarise(
                 n = n(),
                 across(
                   .cols = c(contains("_CWM")),
                   .fns = list(mean = ~mean(., na.rm= TRUE), sd = ~sd(., na.rm = TRUE), se = ~std_error(.x, na.rm= TRUE)),
                   .names = "{.col}_{.fn}"
                 )) %>%
                ungroup()
#pivot longer to have means, sds, and se by habitat, treatment, and trait               
cwmsummaries_long = cwmsummaries %>%
                    pivot_longer(
                      cols = c(contains("CWM")),
                      names_to = "variable",
                      values_to = "value"
                    ) %>%
                   separate(variable, sep = "_", into= c("trait", "CWM", "stat")) %>%
                  #pivot wider again to have stats in separate columns
                  pivot_wider(
                    id_cols = c(habitat, watering, fertilization, WFtreatment, WFtreatment_order, n, trait, CWM),
                    names_from = stat, 
                    values_from = value
                  ) %>%
                   dplyr::select(- CWM)

#SCT
SCTplot = ggplot(data = subset(cwmsummaries_long, trait == "SCT"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+ 
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "SCT"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM SCT", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
SCTplot  

#SCP
SCPplot = ggplot(data = subset(cwmsummaries_long, trait == "SCP"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "SCP"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM SCP", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
SCPplot  

#Length
Lengthplot = ggplot(data = subset(cwmsummaries_long, trait == "Length"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "Length"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM Length", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
Lengthplot  

#Mass
Massplot = ggplot(data = subset(cwmsummaries_long, trait == "Mass"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "Mass"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM Mass", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
Massplot  

#CN
CNplot = ggplot(data = subset(cwmsummaries_long, trait == "CN"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "CN"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM CN", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
CNplot  

#Shape
Shapeplot = ggplot(data = subset(cwmsummaries_long, trait == "Shape"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "Shape"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM Shape", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
Shapeplot  

#### Figure 2: CWM individual traits ####
#changed order of traits to align with table 1
plot_grid(SCTplot + theme(legend.position = "none") +
            theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")) + 
            labs(title = "Figure 2") + theme(plot.title = element_text(vjust = 5, hjust = -.05))+
            theme(axis.title.x = element_blank()),
          Massplot + theme(legend.position = "none")+
            theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")) + 
            labs(title = "    ") + theme(plot.title = element_text(vjust = 5, hjust = -.05))+
            theme(axis.title.x = element_blank()),
          SCPplot + theme(legend.position = "none")+
            theme(axis.title.x = element_blank()),
          CNplot + theme(legend.position = "none")+
            theme(axis.title.x = element_blank()),
          Lengthplot,
          Shapeplot + theme(
            legend.title = element_blank(),
            legend.text = element_blank(),
            legend.key = element_blank(),
            legend.background = element_blank())  + 
            guides(shape=guide_legend(override.aes=list(shape=NA))),
          ncol = 2, byrow= T,
          labels = c("A.", "B.", "C.", "D.", "E.", "F.", "G."), label_size=14, vjust = 4)
#ggsave("Plots/Fig2_CWMresponses_scatterplot.jpg", height = 10, width = 10)
#ggsave("Plots/Fig2_CWMresponses_scatterplot.pdf", height = 10, width = 10)



#### Figure S3: starch CWM ####
#Starch
Starchplot = ggplot(data = subset(cwmsummaries_long, trait == "Starch"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "Starch"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM Starch", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
Starchplot  
#ggsave("Plots/FigS3_CWMresponses_starch_scatterplot.jpg", height = 5, width = 5)
#ggsave("Plots/FigS3_CWMresponses_starch_scatterplot.pdf", height = 5, width = 5)


#### PC CWM contrasts and figures ####
## PC1 ##
#2 way is sig#
PC1_lm=lme(PC1_CWM~ habitat*watering + habitat*fertilization + watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(PC1_lm) 

PC1_emm = emmeans(PC1_lm,  ~ habitat*watering*fertilization) 

cld_PC1 = as.data.frame(cld(PC1_emm, 
                               adjust = "Tukey",     # p-value adjustment
                               Letters = letters,    # Specify letters to use
                               alpha = 0.05,         # Significance level
                               reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(PC1_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

#PC1
PC1plot = ggplot(data = subset(cwmsummaries_long, trait == "PC1"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "PC1"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM PC1", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
PC1plot  


## PC2 ##
#Only main effects are sig #
PC2_lm =lme(PC2_CWM~ habitat+ watering+ fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(PC2_lm) 

PC2_emm = emmeans(PC2_lm,  ~ habitat*watering*fertilization) 

cld_PC2 = as.data.frame(cld(PC2_emm, 
                           adjust = "Tukey",     # p-value adjustment
                           Letters = letters,    # Specify letters to use
                           alpha = 0.05,         # Significance level
                           reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(PC2_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

#PC2
PC2plot = ggplot(data = subset(cwmsummaries_long, trait == "PC2"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "PC2"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM PC2", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
PC2plot 

## PC3 ##
#2 way is sig#
PC3_lm=lme(PC3_CWM~ habitat*watering + habitat*fertilization + watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(PC3_lm) 

PC3_emm = emmeans(PC3_lm,  ~ habitat*watering*fertilization) 

cld_PC3 = as.data.frame(cld(PC3_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(PC3_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

#PC3
PC3plot = ggplot(data = subset(cwmsummaries_long, trait == "PC3"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "PC3"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM PC3", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
PC3plot 


## PC4 ##
#2 way is sig#
PC4_lm=lme(PC4_CWM~ habitat*watering + habitat*fertilization + watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")

anova(PC4_lm) 

PC4_emm = emmeans(PC4_lm,  ~ habitat*watering*fertilization) 

cld_PC4 = as.data.frame(cld(PC4_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(PC4_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "W",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "WN"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

#PC4
PC4plot = ggplot(data = subset(cwmsummaries_long, trait == "PC4"),aes(x=WFtreatment_order, y= mean, group = habitat, shape = habitat))+
  geom_point( size = 4)+
  scale_shape_manual(values = c(15,5)) +
  geom_errorbar(data = subset(cwmsummaries_long, trait == "PC4"), aes(ymin = mean - se, ymax = mean + se), width = 0.2) + 
  theme_bw() + 
  labs(x = "Treatment", y = "CWM PC4", shape = "Habitat") + 
  theme(legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
PC4plot 
 
#### Fig S2 ####
plot_grid(PC1plot + theme(legend.position = "none", axis.title.x = element_blank()),
          PC2plot+ theme(legend.position = "none",  axis.title.x = element_blank()),
          PC3plot, 
          PC4plot + theme(
            legend.title = element_blank(),
            legend.text = element_blank(),
            legend.key = element_blank(),
            legend.background = element_blank())  + 
            guides(shape=guide_legend(override.aes=list(shape=NA))),
          ncol = 2,
          labels = c("A.", "B.", "C.", "D."), label_size=14)
#ggsave("Plots/FigS2_CWMresponses_PCs_scatterplot.jpg", height = 10, width = 10)
#ggsave("Plots/FigS2_CWMresponses_PCs_scatterplot.pdf", height = 10, width = 10)
