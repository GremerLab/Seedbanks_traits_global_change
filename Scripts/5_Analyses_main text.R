rm(list = ls()) # clears everything
library(vegan)
library(tidyverse)
library(ggbiplot)
library(modelsummary)
library(ggplot2)

### data prep ####
# input transformed trait data from 4_traittransformations
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

#split up dataframe into values vs metadata to run PCA
traitdat_4PCA =traitdat%>%
  dplyr::select(SCT,Mass ,SCP,Carbon, CN ,Length, Starch, Shape, Disp,
                Texture, Compact)#

traitdat_info =traitdat%>%
  dplyr::select(species, annual, family, Functional.group, origin, lump_dispcat, mucilage)

names(traitdat_4PCA)
names(traitdat_info)

#### run PCA ####
trait.pca<-prcomp(traitdat_4PCA, scale=TRUE, center = T)
summary(trait.pca) 
round(trait.pca$rotation,2)
plot(trait.pca, type = "lines")

# plot with species names
ggbiplot::ggbiplot(trait.pca, labels = traitdat_info$species)+
  theme_classic()
# ggsave("Plots/PCA12_specieslabels.png", height = 5, width = 5, dpi = 300)
ggbiplot::ggbiplot(trait.pca, labels = traitdat_info$species, choices = c(3,4))+
  theme_classic()
# ggsave("Plots/PCA34_species labels.png", height = 5, width = 5, dpi = 300)

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
               category == "Barrier" ~ "twodash", 
               category == "Chemical" ~ "dotdash")) %>%
            mutate(color_group = case_when(
              category == "Morphological" ~ "black",
              category == "Barrier" ~ "black", 
              category == "Chemical" ~ "gray35"))

#see how to manually set linetypes and colors for the geom_segment - may need to separate out that call.  
pc12 <- autoplot(trait.pca, data = all2, colour = 'origin', shape =  "Functional.group", loadings = F, size =2, scale = 0) + 
        geom_segment(data = loadingvals, aes(x=0, y=0, xend = PC1*5, yend = PC2*5, linetype = category, color = category), #0.8 just scales loading arrows to fit graph, autoplot does this automatically when plotting loading = T
                     arrow = arrow(length = unit(0.2, "cm")),  size = 0.6)  + #color = "black",
        geom_text(data = loadingvals, mapping = aes(label = trait, x = PC1*5.5, y=PC2*5.5)) + #multiply by .53 to move a little away from arrow
  scale_color_grey(start = 0.4, end = 0.7)+ 
  theme_bw()+
  theme(legend.direction ='horizontal', 
        legend.position = 'bottom',
        legend.title = element_blank())  

  #add means for functional groups
fig1a <- pc12 +  geom_point(x=funtype_mean$`mean(PC1)`[1],y=funtype_mean$`mean(PC2)`[1], size=6, shape = 1)+
  geom_point(x=funtype_mean$`mean(PC1)`[2],y=funtype_mean$`mean(PC2)`[2], size=6, shape = 2) +
  geom_point(x=funtype_mean$`mean(PC1)`[3],y=funtype_mean$`mean(PC2)`[3], size=6, shape = 0)
fig1a

##start here, make the same plot for fig 2b

p12 <- autoplot(trait.pca , data = all2, colour = 'origin',
                loadings = TRUE, loadings.colour = "black", #loadings.colour = c("green", "purple", "green", "orange", "orange", "purple", "orange", "purple", "purple", "purple", "purple"), #  "black",
                loadings.linetype = c(trait.pca$linetype),
                size=2,scale = 0,
                shape="Functional.group",loadings.label=T,loadings.label.size = 3.5,x = 1,y=2, 
                loadings.label.colour = "black",
                #loadings.label.colour= c("green", "purple", "green", "orange", "orange", "purple", "orange", "purple", "purple", "purple", "purple"), #  "black",
                loadings.label.vjust=1.5, loadings.label.hjust=0.5)+
  
  scale_color_grey(start = 0.4, end = 0.7)+ 
  theme_bw()+
  theme(legend.direction ='horizontal', 
        legend.position = 'bottom',
        legend.title = element_blank())

p12 + geom_point(x=funtype_mean$`mean(PC1)`[1],y=funtype_mean$`mean(PC2)`[1], size=6, shape = 1)+
      geom_point(x=funtype_mean$`mean(PC1)`[2],y=funtype_mean$`mean(PC2)`[2], size=6, shape = 2) +
      geom_point(x=funtype_mean$`mean(PC1)`[3],y=funtype_mean$`mean(PC2)`[3], size=6, shape = 0)

#add convex hulls for forb, grass, Nfixer
pca_hull_fxnalgroup12 <- 
  all2 %>% 
  group_by(Functional.group) %>% 
  slice(chull(PC1, PC2))

chull_functionalgroup12 <- 
  p12 +
  geom_polygon(data = pca_hull_fxnalgroup12,
               aes(fill = Functional.group,
                   colour = Functional.group),
               alpha = 0.3,
               show.legend = T) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") 

chull_functionalgroup12

pca_hull_nativenon12 <- 
  all2 %>% 
  group_by(origin) %>% 
  slice(chull(PC1, PC2))

chull_nativenon12 <-  p12 + geom_polygon(data = pca_hull_nativenon12,
                                         aes(fill = origin,
                                             colour = origin),
                                         alpha = 0.3,
                                         show.legend = T) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") 
chull_nativenon12

#do this for PC3 and PC4 too.
p34 <- autoplot(trait.pca_prcomp , data = all2, #colour = 'origin',
                loadings = TRUE, loadings.colour = c("green", "purple", "green", "orange", "orange", "purple", "orange", "purple", "purple", "purple", "purple"), #  "black",
                size=2,scale = 0,
                shape="Functional.group",loadings.label=T,loadings.label.size = 3.5,x = 3,y=4, 
                loadings.label.colour= c("green", "purple", "green", "orange", "orange", "purple", "orange", "purple", "purple", "purple", "purple"), #  "black",
                loadings.label.vjust=1.5, loadings.label.hjust=0.5)+
  
  scale_color_grey(start = 0.4, end = 0.7)+
  theme_bw()+
  theme(legend.direction ='horizontal', 
        legend.position = 'bottom',
        legend.title = element_blank())
p34

#add convex hulls for forb, grass, Nfixer
pca_hull_fxnalgroup34 <- 
  all2 %>% 
  group_by(Functional.group) %>% 
  slice(chull(PC3, PC4))

chull_functionalgroup34 <- 
  p34 +
  geom_polygon(data = pca_hull_fxnalgroup34,
               aes(fill = Functional.group,
                   colour = Functional.group),
               alpha = 0.3,
               show.legend = T) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") 

chull_functionalgroup34

pca_hull_nativenon34 <- 
  all2 %>% 
  group_by(origin) %>% 
  slice(chull(PC3, PC4))

chull_nativenon34 <-  p34 + geom_polygon(data = pca_hull_nativenon34,
                                         aes(fill = origin,
                                             colour = origin),
                                         alpha = 0.3,
                                         show.legend = T) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") 
chull_nativenon34


#### test significance of differences between PC values ####
m=aov(data=all2,PC1~nat.inv+Functional.group)
anova(m)
TukeyHSD(m)
m=aov(data=all2,PC2~nat.inv+Functional.group)
TukeyHSD(m)
anova(m)
m=aov(data=all2,PC3~nat.inv+Functional.group)
TukeyHSD(m)
anova(m)
m=aov(data=all2,PC4~nat.inv+Functional.group)
TukeyHSD(m)
summary(m)

all2%>%
  group_by(nat.inv)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

all2%>%
  group_by(Functional.group)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

all2=all2%>%
  mutate(bi_annual=ifelse(annual=="annual", "A","P"))%>%
  mutate(groupings=paste(nat.inv,Functional.group,sep=""))%>%
  mutate(groupings=ifelse(bi_annual=="P"& Functional.group=="Forb",paste(bi_annual,groupings,sep = ""),groupings))#

all2$groupings

all2%>%
  group_by(groupings)%>%
  dplyr::summarise(mean(PC1),mean(PC2),mean(PC3),mean(PC4))

### Run PC without the 7 species missing texture, very similar output - see script 5_PCA_final

### correlation
names(traitdat)
dat_trans=traitdat%>%
  dplyr::select(SCT,mass,SCP=SCPpc,C,"CN_ratio"=cn,"size"=L,starch, shape="Fit_shape",disp ="germ_is_disp",texture="texture_1",compact="Fit_compact",mucilage,"intensity"=Intensity,P)#

temp=as.data.frame(cor(dat_trans))
pmat=ggcorrplot::cor_pmat(x = dat_trans)
ggcorrplot::ggcorrplot(temp,type="upper",p.mat = pmat, lab=T,insig = "blank")

