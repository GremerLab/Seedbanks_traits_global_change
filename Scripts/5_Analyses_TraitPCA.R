#This script:
#performs principle components analysis on seed functional traits
#tests for differences in seed functional among functional groups and between origins (native vs. non-native)
#generates the trait PCA figures and associated supplemental figures

rm(list = ls()) # clears everything
library(vegan)
library(tidyverse)
library(ggbiplot)
library(modelsummary)
library(ggplot2)
library(cowplot)
library(ggcorrplot)
### data prep ####
# input transformed trait data from 4_trait_transformations
traitdat=read.csv("Cleaned data/transformed_traitdata.csv",header = T) %>%
    #NOTE, will need to put code below into data_cleaning/prep scripts once have all missing files
  mutate(species = as.factor(clean_code), annual = as.factor(annual), Functional.group = as.factor(Functional.group), 
         lump_dispcat = as.factor(lump_dispcat)) %>%
  mutate(origin = as.factor(ifelse(nat.inv=="native","Native","Non-Native"))) %>%
  select(- clean_code, -nat.inv, -SCP) %>% #need clarification on difference between SCP and SCPpc.  SCPps used for Elwood dissertation
  rename(Mass = mass,SCP=SCPpc,Carbon = C, CN =cn,Length = L,Starch = starch, Shape="Fit_shape", Disp = disp2,
         Texture="texture_1", Compact="Fit_compact")

names(traitdat)
summary(traitdat)
dim(traitdat)
rownames(traitdat) = traitdat$species
length(unique(traitdat$species))
#split up dataframe into values vs metadata to run PCA
traitdat_4PCA =traitdat%>%
  dplyr::select(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                Texture, Compact)#

traitdat_info =traitdat%>%
  dplyr::select(species, annual, family, Functional.group, origin, lump_dispcat, mucilage)

names(traitdat_4PCA)
names(traitdat_info)

####Trait PCA ####
trait.pca<-prcomp(traitdat_4PCA, scale=TRUE, center = T)
summary(trait.pca) 
round(trait.pca$rotation,2)
plot(trait.pca, type = "lines")

# plot with species names
ggbiplot::ggbiplot(trait.pca, labels = traitdat_info$species)+
  theme_classic()
# ggsave("Plots/PCA12_specieslabels.jpg", height = 5, width = 5, dpi = 300)
ggbiplot::ggbiplot(trait.pca, labels = traitdat_info$species, choices = c(3,4))+
  theme_classic()
# ggsave("Plots/PCA34_species labels.jpg", height = 5, width = 5, dpi = 300)

## extract loadings for each species
ind.coord=trait.pca$x
speciesvals=cbind(species=rownames(traitdat), ind.coord)
speciesvals=as.data.frame(speciesvals)
loadingvals=trait.pca$rotation
loadingvals=(cbind(trait =paste(rownames(loadingvals)), loadingvals))
loadingvals=as.data.frame(loadingvals)

# export PCA trait loadings and species values
# write.csv(loadingvals,"Output data/Trait_PCA_traitloadings.csv")
# write.csv(speciesvals,"Output data/Trait_PCA_speciesscores.csv")

### Graph PCAs with functional type and status ####
all2= speciesvals %>% #joining species PCA scores with the main trait dataframe
  inner_join(traitdat)%>%
  mutate_at(vars(contains('PC')), list(as.numeric)) %>%
  select(-X)

names(all2)
dim(all2)
summary(all2)

#write.csv(all2, "Output data/Traits_PCscores_all.csv")
library(ggfortify)

funtype_mean = all2 %>%
               group_by(Functional.group)%>%
               dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4),n())
#add trait category to PCA output
loadingvals = loadingvals%>%
              mutate_at(vars(contains('PC')), list(as.numeric)) %>%
              mutate(category = as.factor(c("Barrier", "Morphological", "Barrier", "Chemical", 
                                            "Chemical", "Morphological", "Chemical",
                       "Morphological", "Morphological", "Morphological","Morphological"))) %>% 
             mutate(linetype_group = case_when(
               category == "Morphological" ~ "solid",
               category == "Barrier" ~ "longdash", 
               category == "Chemical" ~ "dotted")) %>%
            mutate(color_group = case_when(
              category == "Morphological" ~ "black",
              category == "Barrier" ~ "black", 
              category == "Chemical" ~ "gray35"))

pc12 <- autoplot(trait.pca, data = all2, colour = 'origin', shape =  "Functional.group", loadings = F, size =2, scale = 0,
                 x=1, y=2) + 
        scale_color_grey(start = 0.4, end = 0.7)+ 
        theme_bw()+
       
        geom_segment(data = loadingvals, aes(x=0, y=0, xend = PC1*5, yend = PC2*5, linetype = category), #0.8 just scales loading arrows to fit graph, autoplot does this automatically when plotting loading = T
                     arrow = arrow(length = unit(0.2, "cm"), type= "closed"),  linewidth = 0.75, color = "black")  + #, color = loadingvals$color_group
        #geom_text(data = loadingvals, mapping = aes(label = trait, x = PC1*5.5, y=PC2*5.5), size = 6 ) +  #multiply by scores to move a little away from arrow
        scale_linetype_manual(values = c("dashed", "dotted", "solid") )+ 
        geom_text_repel(data = loadingvals, aes(x =  PC1*5.5, y = PC2*5.6, label= trait), size=6 )
pc12
  #add means for functional groups
fig1a <- pc12 +  geom_point(x=funtype_mean$`mean(PC1)`[1],y=funtype_mean$`mean(PC2)`[1], size=6, shape = 1, stroke =1.5)+
  geom_point(x=funtype_mean$`mean(PC1)`[2],y=funtype_mean$`mean(PC2)`[2], size=6, shape = 2, stroke =1.5) + #stroke controls outline width
  geom_point(x=funtype_mean$`mean(PC1)`[3],y=funtype_mean$`mean(PC2)`[3], size=6, shape = 0, stroke =1.5) +
  theme(legend.direction ="horizontal", legend.position = "bottom", legend.title = element_blank(), 
        text = element_text(size = 18), legend.key.width = unit(2, "line")) + 
  guides(linetype = "none")
fig1a 

pc34 <- autoplot(trait.pca, data = all2, colour = 'origin', shape =  "Functional.group", loadings = F, size =2, scale = 0,
                 x=3, y=4) + 
  scale_color_grey(start = 0.4, end = 0.7)+ 
  theme_bw()+
  geom_segment(data = loadingvals, aes(x=0, y=0, xend = PC3*5, yend = PC4*5, linetype = category), #0.8 just scales loading arrows to fit graph, autoplot does this automatically when plotting loading = T
               arrow = arrow(length = unit(0.2, "cm"), type= "closed"),  size = 0.75, color = "black")  + #, color = loadingvals$color_group
  geom_text(data = loadingvals, mapping = aes(label = trait, x = PC3*6, y=PC4*5.3), size = 6) +  #multiply by 5.4 to move a little away from arrow
  scale_linetype_manual(values = c("dashed", "dotted", "solid") ) #+ 
 # geom_text_repel(data = loadingvals, aes(x =  PC1*5, y = PC2*5, label= trait), size=6 )
pc34
#add means for functional groups
fig1b <- pc34 +  geom_point(x=funtype_mean$`mean(PC3)`[1],y=funtype_mean$`mean(PC4)`[1], size=6, shape = 1, stroke =1.5)+
  geom_point(x=funtype_mean$`mean(PC3)`[2],y=funtype_mean$`mean(PC4)`[2], size=6, shape = 2, stroke =1.5) +
  geom_point(x=funtype_mean$`mean(PC3)`[3],y=funtype_mean$`mean(PC4)`[3], size=6, shape = 0, stroke =1.5) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.title = element_blank(), 
        text = element_text(size = 16),  legend.key.width = unit(2, "line")) + 
  guides(shape = "none", color = "none")
fig1b 


####Figure 1: Trait PCA ####
plot_grid(fig1a , fig1b  , labels = c("A.", "B."), label_size=14)

#ggsave("Plots/Fig1_TraitPCA.jpg", height = 10, width = 20)

####Figure SX: Trait PCA with convex hulls ####
#add convex hulls for forb, grass, Nfixer
pca_hull_fxnalgroup12 <- 
  all2 %>% 
  group_by(Functional.group) %>% 
  slice(chull(PC1, PC2)) %>%
  mutate(fillgroup = case_when(
        Functional.group == "Forb" ~ "#D95F02FF",
        Functional.group == "Grass" ~ "#1B9E77FF",
        Functional.group == "N-fixer" ~ "#7570B3FF"
  ))
  

figSXa <- 
  fig1a +
  geom_polygon(data = pca_hull_fxnalgroup12,
               aes(fill = Functional.group,
                   color = Functional.group),
               alpha = 0.3,
               show.legend = T) +
  scale_fill_grey(start = 1, end = 0.2) +
  guides( color = "none") #, shape = guide_legend(override.aes = list(fill = NA))

figSXa 

#add convex hulls for forb, grass, Nfixer
pca_hull_fxnalgroup34 <- 
  all2 %>% 
  group_by(Functional.group) %>% 
  slice(chull(PC3, PC4))%>%
  mutate(fillgroup = case_when(
    Functional.group == "Forb" ~ "#D95F02FF",
    Functional.group == "Grass" ~ "#1B9E77FF",
    Functional.group == "N-fixer" ~ "#7570B3FF"
  ))

figSXb <- fig1b +
  geom_polygon(data = pca_hull_fxnalgroup34,
               aes(fill = Functional.group,
                   color = Functional.group),
               alpha = 0.3,
               show.legend = T) +
  scale_fill_grey(start = 1, end = 0.2) +
  guides( color = "none",fill = "none", linetype = guide_legend(override.aes = list(fill = NA)))  #
  

figSXb 

####Fig 1alt: Trait PCA with convex hulls for functional groups ####
plot_grid(figSXa , figSXb  , labels = c("A.", "B."), label_size=14)
#ggsave("Plots/Fig1alt_TraitPCA_withconvexhulls.jpg", height = 10, width = 20)


#### test significance of differences in PC values among functional groups and origin ####
#interaction was never significant for origin x functional group (p>0.09) in any anovas or PERManovas so dropped interaction
permanova_pc1 = adonis2(all2$PC1 ~ origin+Functional.group, method = "euc", data = all2, by= "terms")
permanova_pc1

m1=aov(data=all2,PC1~origin+Functional.group)
anova(m1)
TukeyHSD(m1)
#permanova and anova are concordant, and PC1 varies by functional group and origin 
#Grasses have lower PC1 than forbs, N-fixers have higher PC1 than grasses
#native have higher PC1 than non-native 

permanova_pc2 = adonis2(all2$PC2 ~ origin+Functional.group, method = "euc", data = all2, by= "terms")
permanova_pc2
m2=aov(data=all2,PC2~origin+Functional.group)
anova(m2)
TukeyHSD(m2)
#only sig difference is for functional group
#N fixers have higher PC2 than forbs and grasses

permanova_pc3 = adonis2(all2$PC3 ~ origin+Functional.group, method = "euc", data = all2, by= "terms")
permanova_pc3
m3=aov(data=all2,PC3~origin+Functional.group)
TukeyHSD(m3)
anova(m3)
#origin is sig in permanova
#permanova is more sig than anova, which is marginal p-0.053
#non-native have higher PC3 than native

permanova_pc4 = adonis2(all2$PC4 ~ origin+Functional.group, method = "euc", data = all2, by= "terms")
permanova_pc4
m4=aov(data=all2,PC4~origin+Functional.group)
TukeyHSD(m4)
anova(m4)
#functional group is sig, same in both permanova and anova, but contrasts aren't sig.  
#n-fixers have lower PC4 than natives (p=0.07)

#get means #not sure if we need this 
all2%>%
  group_by(origin)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

all2%>%
  group_by(Functional.group)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

all2=all2%>%
  mutate(bi_annual=ifelse(annual=="annual", "A","P"))%>%
  mutate(groupings=paste(origin,Functional.group,sep=""))%>%
  mutate(groupings=ifelse(bi_annual=="P"& Functional.group=="Forb",paste(bi_annual,groupings,sep = ""),groupings))#

all2$groupings

all2%>%
  group_by(groupings)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

#### PCA without the 7 species missing texture ####
#from 4_trait_transfromations.R: "CHAGLA" "HYPGLA" "LASCAL" "MICDOU" "SISBEL" "TRIGRA" "RIGLEP"
traitdat2 = traitdat %>%
            filter(!species %in% c("CHAGLA", "HYPGLA", "LASCAL", "MICDOU", "SISBEL", "TRIGRA", "RIGLEP"))
dim(traitdat)
dim(traitdat2)

traitdat_4PCA2 =traitdat2%>%
  dplyr::select(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                Texture, Compact)#

traitdat_info2 =traitdat2%>%
  dplyr::select(species, annual, family, Functional.group, origin, lump_dispcat, mucilage)

## trait PCA without 7 species missing texture
trait.pca2<-prcomp(traitdat_4PCA2, scale=TRUE, center = T)
summary(trait.pca2) 
round(trait.pca2$rotation,2)
plot(trait.pca2, type = "lines")

## extract loadings for each species
ind.coord2=trait.pca2$x
speciesvals2=cbind(species=rownames(traitdat2), ind.coord2)
speciesvals2=as.data.frame(speciesvals2)
loadingvals2=trait.pca2$rotation
loadingvals2=(cbind(trait =paste(rownames(loadingvals2)), loadingvals2))
loadingvals2=as.data.frame(loadingvals2)

### Graph PCAs with functional type and status ####
all22= speciesvals2 %>% #joining species PCA scores with the main trait dataframe
  inner_join(traitdat2)%>%
  mutate_at(vars(contains('PC')), list(as.numeric)) %>%
  select(-X)

names(all22)
dim(all22)
summary(all22)

#write.csv(all22, "Output data/Traits_PCscores_wo_7species.csv")

funtype_mean2 = all22 %>%
  group_by(Functional.group)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4),n())

#add trait category to PCA output
loadingvals2 = loadingvals2%>%
  mutate_at(vars(contains('PC')), list(as.numeric)) %>%
  mutate(category = as.factor(c("Barrier", "Morphological", "Barrier", "Chemical", 
                                "Chemical", "Morphological", "Chemical",
                                "Morphological", "Morphological", "Morphological","Morphological"))) %>% 
  mutate(linetype_group = case_when(
    category == "Morphological" ~ "solid",
    category == "Barrier" ~ "longdash", 
    category == "Chemical" ~ "dotted")) %>%
  mutate(color_group = case_when(
    category == "Morphological" ~ "black",
    category == "Barrier" ~ "black", 
    category == "Chemical" ~ "gray35"))

pc122 <- autoplot(trait.pca2, data = all22, colour = 'origin', shape =  "Functional.group", loadings = F, size =2, scale = 0,
                 x=1, y=2) + 
  scale_color_grey(start = 0.4, end = 0.7)+ 
  theme_bw()+
  
  geom_segment(data = loadingvals2, aes(x=0, y=0, xend = PC1*5, yend = PC2*5, linetype = category), #0.8 just scales loading arrows to fit graph, autoplot does this automatically when plotting loading = T
               arrow = arrow(length = unit(0.2, "cm"), type= "closed"),  size = 0.75, color = "black")  + #, color = loadingvals$color_group
  geom_text(data = loadingvals2, mapping = aes(label = trait, x = PC1*5.4, y=PC2*5.4)) +  #multiply by 5.4 to move a little away from arrow
  scale_linetype_manual(values = c("dashed", "dotted", "solid") )

#add means for functional groups
fig1a_alt <- pc122 +  geom_point(x=funtype_mean2$`mean(PC1)`[1],y=funtype_mean2$`mean(PC2)`[1], size=6, shape = 1, stroke=1.5)+
  geom_point(x=funtype_mean2$`mean(PC1)`[2],y=funtype_mean2$`mean(PC2)`[2], size=6, shape = 2, stroke=1.5) +
  geom_point(x=funtype_mean2$`mean(PC1)`[3],y=funtype_mean2$`mean(PC2)`[3], size=6, shape = 0, stroke=1.5) +
  theme(legend.direction ="horizontal", legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 14)) + 
  guides(linetype = "none")
fig1a_alt

pc342 <- autoplot(trait.pca2, data = all22, colour = 'origin', shape =  "Functional.group", loadings = F, size =2, scale = 0,
                 x=3, y=4) + 
  scale_color_grey(start = 0.4, end = 0.7)+ 
  theme_bw()+
  geom_segment(data = loadingvals2, aes(x=0, y=0, xend = PC3*5, yend = PC4*5, linetype = category), #0.8 just scales loading arrows to fit graph, autoplot does this automatically when plotting loading = T
               arrow = arrow(length = unit(0.2, "cm"), type= "closed"),  size = 0.75, color = "black")  + #, color = loadingvals$color_group
  geom_text(data = loadingvals2, mapping = aes(label = trait, x = PC3*5.4, y=PC4*5)) +  #multiply by 5.4 to move a little away from arrow
  scale_linetype_manual(values = c("dashed", "dotted", "solid") )

#add means for functional groups
fig1b_alt <- pc342 +  geom_point(x=funtype_mean2$`mean(PC3)`[1],y=funtype_mean2$`mean(PC4)`[1], size=6, shape = 1, stroke=1.5)+
  geom_point(x=funtype_mean2$`mean(PC3)`[2],y=funtype_mean2$`mean(PC4)`[2], size=6, shape = 2, stroke=1.5) +
  geom_point(x=funtype_mean2$`mean(PC3)`[3],y=funtype_mean2$`mean(PC4)`[3], size=6, shape = 0, stroke=1.5) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 14)) + 
  guides(shape = "none", color = "none")
fig1b_alt 


####Alt Trait PCA figure without 7 species missing texture ####
plot_grid(fig1a_alt , fig1b_alt  , labels = c("A.", "B."), label_size=14)

#### Test for multivariate trait differences without 7 species ####
#interaction was never significant for origin x functional group (p>0.09) in any anovas or PERManovas so dropped interaction
permanova_pc12 = adonis2(all22$PC1 ~ origin+Functional.group, method = "euc", data = all22, by= "terms")
permanova_pc12

m12=aov(data=all22,PC1~origin+Functional.group)
anova(m12)
TukeyHSD(m12)
#permanova and anova are concordant

permanova_pc22 = adonis2(all22$PC2 ~ origin+Functional.group, method = "euc", data = all22, by= "terms")
permanova_pc22
m22=aov(data=all22,PC2~origin+Functional.group)
anova(m22)
TukeyHSD(m22)


permanova_pc32 = adonis2(all22$PC3 ~ origin+Functional.group, method = "euc", data = all22, by= "terms")
permanova_pc32
m32=aov(data=all22,PC3~origin+Functional.group)
TukeyHSD(m32)
anova(m32)

permanova_pc42 = adonis2(all22$PC4 ~ origin+Functional.group, method = "euc", data = all22, by= "terms")
permanova_pc42
m42=aov(data=all22,PC4~origin+Functional.group)
TukeyHSD(m42)
summary(m42)

####traitdat#### Trait correlation matrix, Fig S2 ####
traits_all = traitdat %>%
             select(SCT, Mass, SCP, Carbon, CN, Length, Starch, Shape, Disp, Texture, Compact, mucilage, Intensity, Perimeter =P)
traitcor_all =as.data.frame(cor(traits_all))
pmat=ggcorrplot::cor_pmat(x = traits_all)
FigS2 = ggcorrplot::ggcorrplot(traitcor_all,type="upper",p.mat = pmat, lab=T,insig = "blank")
FigS2
#ggsave("Plots/FigS2_Traits_correlationmatrix.jpg", height = 10, width = 10)
