## Heterozygosity 

library(data.table)
library(tidyverse)
library(FSA)

#new = readRDS(file = "../../../../../Desktop/new.RDS")


glf <- read.table(file = 'both_species_hetero.beagle.gz', header=TRUE)[,-c(1:3)]
glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))
EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])
maxiter <- 20000
tol <- 1e-8

for(sample in 1:dim(glf)[2])
{
  
  for (iter in 1:maxiter)
  {
    upd <- EMstep (SFS[,sample], glf[,sample,])
    if (sqrt(sum((upd - SFS[,sample])^2)) < tol)
      break;
    SFS[,sample] <- upd
  }
  if (iter == maxiter) warning("increase maximum number of iterations")
}

#By site
#read in metadata, and bam list
expDesign=data.frame(read.csv("../../expDesign.csv", header=T))
bams=data.frame(read.csv("../both_species_bams", header = FALSE)) # list of bam files for RDA outlier loci
colnames(bams) = c("bam_name")

#matching adapter names between metadata and bams file
expDesign = expDesign %>%
  mutate(bam_name = paste(Sample, ".trim.bt2.bam"),
         bam_name = str_replace(bam_name," ", ""),
         pop = str_replace(Site, " ", "_"),
         species = str_replace(Species, " ", "_")) %>%
  select(bam_name, species, pop)

fornicata = expDesign %>%
  filter(species == "Crepidula_fornicata")

plana = expDesign %>%
  filter(species == "Crepidula_plana")


flist=list()
#for loop to match SFS output for each sample and calculate population level heterozygosities
for(i in 1:length(levels(as.factor(fornicata$pop)))) {
  varname=levels(as.factor(fornicata$pop))[i]
  varsub=subset(fornicata, pop==varname)
  varnums=as.character(which(bams$bam_name %in% varsub$bam_name))
  sfs=SFS[2,c(as.numeric(varnums))]
  Site=replicate(length(varnums), varname)
  df=data.frame(sfs, Site)
  flist[[i]]=df
}

fornicata_list= rbindlist(flist)

plist=list()
#for loop to match SFS output for each sample and calculate population level heterozygosities
for(i in 1:length(levels(as.factor(plana$pop)))) {
  varname=levels(as.factor(plana$pop))[i]
  varsub=subset(plana, pop==varname)
  varnums=as.character(which(bams$bam_name %in% varsub$bam_name))
  sfs=SFS[2,c(as.numeric(varnums))]
  Site=replicate(length(varnums), varname)
  df=data.frame(sfs, Site)
  plist[[i]]=df
}

plana_list= rbindlist(plist)

#NAMES=levels(as.factor(All$Site)) 
#colors=c('#9ad790ff','#fa9fb5ff','#d8efd1ff','#2b8cbeff','#ae017eff',
#         '#4eb3d3ff','#f768a1ff','#d40000ff','#ff5555ff','#00441bff',
#         '#ff8500ff','#dd3497ff','#ccebc5ff','#810f7cff','#238b45ff','#7bccc4ff','#0c2c84ff')
#names(colors)=NAMES

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

fornicata_list$Site=factor(fornicata_list$Site,levels = c("Cape_May", "Newport", "Beverly", "Kettle_Cove", "Robbinston"))
plana_list$Site=factor(plana_list$Site,levels = c("Cape_May", "Newport", "Beverly", "Kettle_Cove"))

colours = c("Robbinston" = "#264D59", "Kettle_Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape_May" = "#D46C4E")


f_plot = 
  ggplot(fornicata_list, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colours, guide="none")+
  theme_classic(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 35, vjust = , hjust=1)) +
  stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")

p_plot = 
  ggplot(plana_list, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colours, guide="none")+
  theme_classic(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 35, vjust = , hjust=1)) +
  stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")


dunnres=dunnTest(sfs ~ Site,
                 data=All,
                 method="bh")
df=dunnres[["res"]]