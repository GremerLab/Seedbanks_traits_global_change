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
               mutate(relab_traitsonly = rawcount/totplottraits) %>%
               mutate(habitat = fct_recode(habitat, "Harsh" = "Harshserp", "Lush" = "Lushserp")) %>%
               mutate(WFtreatment = fct_recode(WFtreatment, "R" = "W", "N" = "F", "NR" = "FW"))
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
         mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR"))) #to match Eskelinen et al. 2021

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


#### Post-hoc contrasts ####
library(multcomp)
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
                watering == "watered" & fertilization == "none" ~ "R",
                watering == "none" & fertilization == "fertilized" ~ "N",
                watering == "watered" & fertilization == "fertilized" ~ "NR"
          ))) %>%
      mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))
          
## Mass ##
#2 way: hab x fert is sig #
Mass_lm=lme(Mass_CWM ~ habitat*watering+habitat*fertilization+watering*fertilization,random=~1|Line, data=cwmdat,na.action=na.exclude,method = "REML")
anova(mass_lm)

Mass_emm = emmeans(Mass_lm,  ~ habitat*watering*fertilization) #keep all factors in here to generate letters

cld_Mass = as.data.frame(cld(Mass_emm, 
                            adjust = "Tukey",     # p-value adjustment
                            Letters = letters,    # Specify letters to use
                            alpha = 0.05,         # Significance level
                            reversed = TRUE) )  %>%    # Sort means in decreasing order
  rename(Mass_letters = .group)  %>% # Rename the letters column
  mutate(WFtreatment = as.factor(case_when(
    watering == "none" & fertilization == "none" ~ "C",
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))

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
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))

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
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))

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
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))

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
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))
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
    watering == "watered" & fertilization == "none" ~ "R",
    watering == "none" & fertilization == "fertilized" ~ "N",
    watering == "watered" & fertilization == "fertilized" ~ "NR"
  ))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","R", "N", "NR")))


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
### save figure 2: CWM traits ###
plot_grid(a + theme(legend.position = "none"),
          b+ theme(legend.position = "none"),
          c+ theme(legend.position = "none"),
          d+ theme(legend.position = "none"),
         # e+ theme(legend.position = "none"), #ECE didn't have length in her graph, could add/remove here
          f, 
          g,
          ncol = 2, nrow = 3,
          labels = c("A.", "B.", "C.", "D.", "E.", "F.", "G."), label_size=14)
#note, can't get formatting to look good with really long panel

#ggsave("Plots/Fig2_CWMresponses.jpg", height = 10, width = 15)
e
#ggsave("Plots/Fig2_CWMresponses_length.jpg", height = 5, width = 8)

