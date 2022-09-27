#First look for clones
setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Astrangia&Oculina/")
ma = as.matrix(read.table("myresult_withclones.ibsMat"))
bams=data.frame(read.csv("Allbams", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"

#read in meta data
samples=data.frame(read.csv("Adapter2Sample.csv", header=T))#[,3:4])
#change column names
colnames(samples)<-c("Sample","pop","Adapter")

#if necessary change bam extensions to match meta data names
samples <- data.frame(lapply(samples, function(x) {
                    gsub("trim$", "trim.bt2.bam", x)
                }))

#merge
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]


dimnames(ma)=list(new$Sample)   
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.35)

#checking after removing clones
ma = as.matrix(read.table("myresult_withoutclones.ibsMat"))
bams=data.frame(read.csv("bamscl", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
samples=data.frame(read.csv("Adapter2Sample.csv", header=T))#[,3:4])
colnames(samples)<-c("Sample","pop","Adapter")

samples <- data.frame(lapply(samples, function(x) {
  gsub("trim$", "trim.bt2.bam", x)
}))

new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]

dimnames(ma)=list(new$Sample)
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7)
