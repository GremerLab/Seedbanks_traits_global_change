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

#load trait and abundance data
all_abtraits = read.csv("Data/Transformed trait data_with abundance.csv") %>% #
  mutate(Line = as.factor(Line), habitat=as.factor(habitat), watering=as.factor(watering), 
         fertilization=as.factor(fertilization), WFtreatment=as.factor(WFtreatment), Plot=as.factor(Plot), 
         species = as.factor(species), origin = as.factor(origin)) 
summary(all_abtraits)


#drop species for which origin is NA
nas_origins = all_abtraits %>%
              filter(is.na(origin) == T)
unique(nas_origins$species)
summary(nas_origins) #most are species for which we do not have traits and those that were grouped by genus

all_abtraits = all_abtraits %>%
               filter(is.na(origin) == F)
summary(all_abtraits)

#### calculate CWMs ####
#NOTE: this uses the relative abundance for the species for which we have traits.  
cwmdat_a = all_abtraits %>%
  mutate(across(c(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                  Texture, Compact), ~.x*relab_traitsonly, .names = "{.col}_CWM")) %>%
  mutate(across(starts_with("PC"), ~.x*relab_traitsonly,.names = "{.col}_CWM") ) 

  
#next, sum for each plot to get a plot level CWM
#here adding origin to get native vs non-native patterns
cwmdat = cwmdat_a %>%
  group_by(origin, Line, habitat, watering, fertilization, WFtreatment, Plot) %>%
  summarize(across(ends_with("_CWM"), ~sum(.x, na.rm=T)))%>%
  ungroup()%>%
  mutate(watering=as.factor(ifelse(watering=="unwatered","none","watered")),
         fertilization=as.factor(ifelse(fertilization=="unfertilized","none","fertilized"))) %>%
  mutate(WFtreatment_order = factor(WFtreatment, levels = c("C","W", "N", "WN")))

summary(cwmdat)        
dim(cwmdat) #90 plots x 2 origins
str(cwmdat)


##filter to harsh soils which had stronger responses to treatments ##
cwmdat_harsh = cwmdat %>%
               filter(habitat == "Harsh")
summary(cwmdat_harsh)

#### test for CWM responses to treatments for targeted functional traits ####
traitlist = c("SCT_CWM","SCP_CWM", "Length_CWM", "Shape_CWM")
length(traitlist)
output_ind = list() 
mod_est_ind = list()

for (i in 1:length(traitlist)){
  modeldat=cwmdat_harsh
  response= traitlist[i] 
  COL=which((response == names(modeldat))==T)
  temp=as.data.frame(modeldat[COL])
  colnames(temp)[1]="responsevar"
  modeldat=cbind(modeldat,temp)
  names(modeldat)
  
  ##  with 3 way interaction 
  m0=lme((responsevar)~ origin*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
  #m0ml=lme((responsevar)~ origin*watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
  
  anova_output0 = as.data.frame(anova(m0)) %>%
    mutate(factor = rownames(.), trait = traitlist[i])
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
    mutate(factor = rownames(.), trait = traitlist[i]) %>%
    rename(Est = "m0$coefficients$fixed")
  
  if(round(anova_output0$'p-value'[anova_output0$factor == "origin:watering:fertilization"],3) >0.05){ #condition 1
    m1=lme((responsevar)~ origin*watering+origin*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
    #m1ml=lme((responsevar)~ origin*watering+origin*fertilization+watering*fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "ML")
    # Create the kable table
    anova_output1 = as.data.frame(anova(m1)) %>%
      mutate(factor = rownames(.), trait = traitlist[i])
    
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
      mutate(factor = rownames(.), trait = traitlist[i]) %>%
      rename(Est = "m1$coefficients$fixed")
    
    if(round(anova_output1$'p-value'[5],3) >0.05 & round(anova_output1$'p-value'[6],3) >0.05 & round(anova_output1$'p-value'[7],3) >0.05){ #condition 2
      m2=lme((responsevar)~ origin+ watering+ fertilization,random=~1|Line, data=modeldat,na.action=na.exclude,method = "REML")
      anova_output2 = as.data.frame(anova(m2))%>%
        mutate(factor = rownames(.), trait = traitlist[i])
      
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
        mutate(factor = rownames(.), trait = traitlist[i]) %>%
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
anova_cwm_origin_traits = do.call(rbind.data.frame, output_ind) %>%
  rename(Fval = 'F-value', p= 'p-value') %>%
  mutate(p = ifelse(round(p,3) == 0, "<0.001",round(p, 3))) %>%
  mutate(Fval = round(Fval, 3))
#write.csv(anova_cwm_origin_traits, file = "Output tables/CWM_origintraits_ANOVA_all.csv")

#Better formatting for table 2
anova_cwm_origin_traits_wider = anova_cwm_origin_traits %>%
  mutate(DF = paste(numDF, denDF, sep = ","), 
         FP= paste(Fval,p, sep = ", ")) %>%
  dplyr::select(-numDF, -denDF, -Fval, -p) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, FP)
  ) %>%
  rename(DF = DF_origin) %>%
  dplyr::select(!starts_with("DF_")) %>%
  rename_with(~str_remove(.x, "FP_")) 

#Nicer to have sub columns for F and p, so pivot longer again
anova_cwm_origin_traits_wider2 = anova_cwm_origin_traits %>%
  mutate(Fval = as.character(Fval)) %>%
  mutate(DF = paste(numDF, denDF, sep = ",")) %>%
  dplyr::select(-numDF, -denDF) %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = c(DF, Fval, p)
  ) %>%
  rename(DF = DF_origin) %>%
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


#write.csv(anova_cwm_origin_traits_wider2, file = "Output tables/CWM_origin_harsh_ANOVA_all_wide.csv")


#put model estimates together and match with anovas
origin_harsh_est_indtraits = do.call(rbind.data.frame, mod_est_ind)

origin_harsh_anova_est = left_join(anova_cwm_origin_traits,origin_harsh_est_indtraits)

#format estimates
origin_harsh_est_indtraits_wider = origin_harsh_est_indtraits %>%
  pivot_wider(
    id_cols = c(trait),
    names_from = factor, 
    values_from = Est
  )  

#### Post-hoc contrasts ####

##SCT ##
#2 way interactions between origin x watering and origin x fertilization are significant
SCT_lm =lme(SCT_CWM ~ origin*watering + origin*fertilization,random=~1|Line, data=cwmdat_harsh,na.action=na.exclude,method = "REML")
anova(SCT_lm)

SCT_emm = emmeans(SCT_lm,  ~ origin*watering*fertilization) #keep all factors to generate letters

cld_SCT_harshorigin  = as.data.frame(cld(SCT_emm, 
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

##SCP ##
#Only main effects of origin and watering are significant
SCP_lm =lme(SCP_CWM ~ origin + watering + fertilization,random=~1|Line, data=cwmdat_harsh,na.action=na.exclude,method = "REML")
anova(SCP_lm)

SCP_emm = emmeans(SCP_lm,  ~ origin*watering*fertilization) #keep all factors to generate letters

cld_SCP_harshorigin  = as.data.frame(cld(SCP_emm, 
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

##length ##
#Origin x fertilization is sig
Length_lm =lme(Length_CWM ~  watering + origin*fertilization,random=~1|Line, data=cwmdat_harsh,na.action=na.exclude,method = "REML")
anova(Length_lm)

Length_emm = emmeans(Length_lm,  ~ origin*watering*fertilization) #keep all factors to generate letters

cld_Length_harshorigin  = as.data.frame(cld(Length_emm, 
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

##Shape ##
#3 way interaction is sig
Shape_lm =lme(Shape_CWM ~  watering*origin*fertilization,random=~1|Line, data=cwmdat_harsh,na.action=na.exclude,method = "REML")
anova(Shape_lm)

Shape_emm = emmeans(Shape_lm,  ~ origin*watering*fertilization) #keep all factors to generate letters

cld_Shape_harshorigin = as.data.frame(cld(Shape_emm, 
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
SCTplot_harshorigin = ggplot(data=cwmdat_harsh,aes(x=WFtreatment_order, y=SCT_CWM, fill = origin))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM SCT", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(origin),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "origin")

a_harshorigin= SCTplot_harshorigin + geom_text(data = cld_SCT_harshorigin, aes( x = WFtreatment_order, y = 2.1, group = origin, label = SCT_letters ),
                       position = position_dodge(width = 0.9))  + facet_grid(~origin)
a_harshorigin

##SCP ##
SCPplot_harshorigin = ggplot(data=cwmdat_harsh,aes(x=WFtreatment_order, y=SCP_CWM, fill = origin))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM SCP", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(origin),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "origin")

b_harshorigin= SCPplot_harshorigin + geom_text(data = cld_SCP_harshorigin, aes( x = WFtreatment_order, y = 2.1, group = origin, label = SCP_letters ),
                                               position = position_dodge(width = 0.9))  + facet_grid(~origin)
b_harshorigin

##Length ##
Lengthplot_harshorigin = ggplot(data=cwmdat_harsh,aes(x=WFtreatment_order, y=Length_CWM, fill = origin))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Length", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(origin),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "origin")

c_harshorigin= Lengthplot_harshorigin + geom_text(data = cld_Length_harshorigin, aes( x = WFtreatment_order, y = 2.4, group = origin, label = Length_letters ),
                                               position = position_dodge(width = 0.9))  + facet_grid(~origin)
c_harshorigin


##Shape ##
Shapeplot_harshorigin = ggplot(data=cwmdat_harsh,aes(x=WFtreatment_order, y=Shape_CWM, fill = origin))+
  geom_boxplot(position = position_dodge(.9))+
  labs(title="",subtitle = "",y="CWM Shape", x="")+
  scale_fill_grey(start=0.9, end=0.4) +
  scale_color_grey(start=0.4, end=0.7) +
  #facet_grid(rows = vars(type), cols=vars(origin),switch="y",scales = "free")+
  theme_bw()+
  theme(panel.spacing = unit(0, units = "cm"), 
        legend.position = "bottom",legend.title = element_blank(),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  +
  labs(fill = "origin")

d_harshorigin= Shapeplot_harshorigin + geom_text(data = cld_Shape_harshorigin, aes( x = WFtreatment_order, y = 2.1, group = origin, label = Shape_letters ),
                                                  position = position_dodge(width = 0.9))  + facet_grid(~origin)
d_harshorigin

## panel figure ##
#### Fig. 3: CWM traits in harsh soils by origin ####
plot_grid(a_harshorigin + theme(legend.position = "none"),
          b_harshorigin+ theme(legend.position = "none"),
          c_harshorigin + theme(legend.position = "none"),
          d_harshorigin + theme(legend.position = "none"),
          ncol = 2, byrow= T,
          labels = c("A.", "B.", "C.", "D."), label_size=14)
#ggsave("Plots/Fig3_CWMresponses_harsh_origin_faceted.jpg", height = 10, width = 10)
