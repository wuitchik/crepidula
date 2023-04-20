### This script builds a dendrogram for looking at how related samples are based on SNPs

library(tidyverse)

### Crepidula fornicata ###

#First we start by reading in the bams files.
bams=read.table("../fornicata_bams")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("../fornicata_results.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")

pdf("Fornicata_dendrogram.pdf",
    width = 8,
    height = 4)

plot(hc,cex=0.7) 

dev.off()

### Oculina arbuscula ###

bams=read.table("Oculina/bam_list")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("Oculina/coral_ang.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")

pdf("Oculina/Oculina_dendrogram.pdf",
    width = 8,
    height = 4)

plot(hc,cex=0.7) 

dev.off()
