## PCoA based on identity by state (IBS) based on single read resampling

#First confirm clones, or weirdo's
ma = as.matrix(read.table("../myresult.withclones.ibsMat"))
bams=data.frame(read.csv("../bams", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
samples=data.frame(read.csv("../Adapter2Sample.csv"))
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]