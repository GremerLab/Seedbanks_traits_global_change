

## Things to run in referencing script
# rm(list = ls()) # clears everything
# library(tidyverse)

# dat2fixnames=read.csv("Seed Trait Paper/Analysis_data/SBRA_traitdat_Oct11.csv",header = T)
# dat2fixnames$prevname=dat2fixnames$species_code  ##here matching with species_code
# dat2fixnames$species_code=NULL

# old species list
# spnames=read.csv("Seed Trait Paper/Master_Species_name_key.csv", header=T)
# spnames=read.csv("Seed Trait Paper/species_metadata_Nov22.csv",header=T, strip.white = T)
# spnames=read.csv("Seed Trait Paper/updated_name_issues1_updated.csv",header=T, strip.white = T)
# spnames=read.csv("Seed Trait Paper/current_csv files/updated_name_issues2_withIDsforall_Dec2022_v2.csv",header=T, strip.white = T)

##Autorun from here#####
          
# Species name data
spnames=read.csv("Cleaned data/updated_name_issues2_withIDsforall_Dec2022.csv",header=T, strip.white = T)
spnames$code=NULL
spnames$Species_name=NULL

names(spnames)
names(dat2fixnames)


#first easy working codes: semi to keep only rows where they match up /exist in both datasets
# right not in left - dropped
A=dat2fixnames%>%
  inner_join(spnames, by=c("prevname"="species_code"))%>%
  select(-falsename,-tolumpname,species_code=prevname)%>%
  distinct()%>%
  mutate(prevname=species_code)


#now to check for any identified name issues

AA=dat2fixnames%>% # check if any known issue names are in dataset
  inner_join(spnames,by=c("prevname"="falsename"))%>%
  select(-tolumpname)%>%
  distinct()

AAA=dat2fixnames%>% # check if any known issue names are in dataset
  inner_join(spnames,by=c("prevname"="tolumpname"))%>%
  select(-falsename)%>%
  distinct()


if("obsID" %in% names(A)){
Anums=A$obsID
AAnums=AA$obsID
AA=AA%>%
  filter(obsID %in% Anums ==F)

AAAnums=AAA$obsID
AAA=AAA%>%
  filter(obsID %in% Anums ==F)
}

A4=rbind(A,AA,AAA)%>%
  distinct()%>%
  filter(is.na(prevname)==F)


## deal with too many?
# temp=as.data.frame(table(A4$obsID))%>%
#   filter(Freq>1)
# 
# temp$obsID=as.integer(temp$Var1)
# temp2=A4%>%
#   semi_join(temp)%>%
#   distinct()%>%
#   filter(relcover>0)
# 
# unique(temp2$species_code)
# names(temp)
# A4$obsID

# now to identify nonmatchers
OGnewdatnames=unique(dat2fixnames$prevname)
yesnameissue=OGnewdatnames %in% unique(A4$prevname)
B=as.data.frame(yesnameissue,OGnewdatnames)%>%
  filter(yesnameissue==FALSE)

issues_out=row.names(B)

cleannamedat=A4%>%
  mutate(clean_code=toupper(as.character(species_code)))%>%
  select(-prevname,-species_code)%>%
  distinct()

outnames=as.matrix(spnames)
temp=matrix(NA,length(issues_out),length(spnames))
temp[,8]=issues_out
outnames=rbind(outnames,temp)

write.csv(outnames,"Cleaned data/updated_name_issues2.csv", row.names=F)

if(length(issues_out)==0){
  print("No name issues")
}

rm(A,AA,AAA,A4,B,issues_out,temp,dat2fixnames,spnames)

