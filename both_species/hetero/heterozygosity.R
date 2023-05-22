## Heterozygosity 

library(data.table)
library(tidyverse)
library(FSA)

new = readRDS(file = "../../../../../Desktop/new.RDS")


glf <- read.table(file = '../both_species_results.beagle.gz', header=TRUE)[,-c(1:3)]

glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))
EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])
maxiter <- 10000
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
#read in metadata
samples=data.frame(read.csv("../../expDesign.csv", header=T))#[,3:4])

#matching adapter names between metadata and bams file
samples = samples %>%
  mutate(bam_name = paste(Sample, ".trim.bt2.bam")) %>%
  mutate(bam_name = str_replace(bam_name," ", ""))

#samples <- data.frame(lapply(samples, function(x) {
#  gsub("trim$", "trim.bt2.bam", x)
#}))
row.names(samples)=samples$bam_name

bams=data.frame(read.csv("../both_species_bams", header = FALSE)) # list of bam files for RDA outlier loci
#bams <- data.frame(lapply(bams, function(x) {
#  gsub(".*AstrangiaAndOculina/", "", x)
#}))

colnames(bams)<- "Adapter"
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]
new$AstOc=ifelse(new$pop %in% c("NC","NC_Radio","WH","VA","Texas_PortA","Texas_Packery","RI","Panacea"), "Astrangia","Oculina")
new$AstOc=ifelse(new$Sample %in% c("A1T2","A3T2"),"Oculina", new$AstOc)
rep_str=c("^NC$"='Triangle Wrecks', 'NC_Radio'='Radio Island (Astrangia)',
          'VA'='J.B. Eskridge Wreck','WH'='Woods Hole', 'RI'='Fort Wetherill',
          'Texas_PortA'='Port Aransas','Texas_Packery'='Packery', 'Oculina_Panacea'='Panama City',
          'Oculina_NC'='Radio Island (Oculina)','Georgia_R2'='R2 Tower', 'Georgia_J'='J Reef', 
          'Fort_Pierce'='Fort Pierce','Cape_Florida'='Cape Florida')
new$pop <- str_replace_all(new$pop, rep_str)


test = samples %>%
  mutate(pop = str_replace(Site, " ", "_"))

dflist=list()
#for loop to match SFS output for each sample and calculate population level heterozygosities
for(i in 1:length(levels(as.factor(test$pop)))) {
  varname=levels(as.factor(test$pop))[i]
  varsub=subset(test, pop==varname)
  varnums=as.character(which(bams$Sample %in% varsub$bam_name))
  sfs=SFS[2,c(as.numeric(varnums))]
  Site=replicate(length(varnums), varname)
  df=data.frame(sfs, Site)
  dflist[[i]]=df
}
All= rbindlist(dflist)

NAMES=levels(as.factor(All$Site)) 
colors=c('#9ad790ff','#fa9fb5ff','#d8efd1ff','#2b8cbeff','#ae017eff',
         '#4eb3d3ff','#f768a1ff','#d40000ff','#ff5555ff','#00441bff',
         '#ff8500ff','#dd3497ff','#ccebc5ff','#810f7cff','#238b45ff','#7bccc4ff','#0c2c84ff')
names(colors)=NAMES

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
All$Site=factor(All$Site,levels = c("Cape_May", "Newport", "Beverly", "Kettle_Cove", "Robbinston"))

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")


ggplot(All, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colours, guide="none")+
  theme_classic(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 35, vjust = , hjust=1))
p+stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")

dunnres=dunnTest(sfs ~ Site,
                 data=All,
                 method="bh")
df=dunnres[["res"]]