# both species pcangsd

# pcangsd fornicata

library(tidyverse)
library(RcppCNPy)

expDesign = read.csv("../../expDesign.csv") 

write.table(expDesign, "../i2p", sep = "\t", row.names = F)

# assembling the input table
dir="../" # path to input files
inName="both_species_pcangsd.admix.2.Q" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="i2p" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

#------------

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=T)

npops=ncol(tbl)
i2p=read.table(paste(dir,pops,sep=""),header=T)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
tbl$pop=as.factor(tbl$pop)

levels(tbl$pop) = c("Robbinston", "Kettle Cove", "Beverly", "Newport", "Cape May")

ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.05,angle=45,cex=0.8)
