getwd()

setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Astrangia&Oculina/Admixture")
#read in metadata
samples=data.frame(read.csv("../Adapter2Sample.csv", header=T))#[,3:4])
#rename columns for samples
colnames(samples)<-c("Sample","pop","Adapter","apo.sym","sex")

samples <- data.frame(lapply(samples, function(x) {
  gsub("trim$", "trim.bt2.bam", x)
}))

row.names(samples)=samples$Adapter

bams=data.frame(read.csv("../bamscl", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]

#write pop table for raxml
#write.table(new[-1,3],file="phyallpops",quote = FALSE, row.names = FALSE)

#Nsites=Nnew$Site
site=new$pop



# settign up colors for plotting
# palette(rainbow(length(unique(Nsites))))
# colors=as.numeric(as.factor(Nsites))
# colpops=as.numeric(as.factor(sort(unique(Nsites))))



palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))
#--------------------
#run above code 'loading invidividual.." first
# covariance / PCA 

library(vegan)
# choose either of the following two covarince matrices:
co = as.matrix(read.table("myresult_withoutclones.covMat")) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
dimnames(co)=list(new$Sample)
# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(site)
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
kmeans.axes
#2
# This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
# Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering

axes2plot=c(1,2)  
#quartz()
cc=pp0
plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19,col=colors)
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cc,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops)

# unscaled, to identify outliers
n2identify=2
cmd=pp0
cmd$CA$u
# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=rownames(co),n=n2identify,cex=0.7)
#
library(dplyr)
library(ggplot2)

pca_s <- as.data.frame(cmd$CA$u[,axes2plot])

pca_s$site=as.factor(new$pop)
pca_s$site=as.factor(pca_s$site)
pca_s$label=as.factor(new$Sample)
levels(pca_s$site) <- list(NC_TriangleWrecks  = "NC", NC_RadioIsland = "NC_Radio",Virginia="VA",Massachusetts="WH",Rhode_Island="RI",
                           TX_PortA="Texas_PortA",TX_Packery="Texas_Packery", NC_RadioIsland_Oculina="Oculina_NC",FL_Panacea="Panacea")

colors=c('#e3dbdbff','#bae4bc','#00ccc4ff',"#2b8cbe", "#0C2C84","#66a61e","#d40000ff","#ff8500ff","#ff8080ff")
pca_s$site=factor(pca_s$site,levels=c("NC_RadioIsland","NC_TriangleWrecks","Virginia","Rhode_Island","Massachusetts","NC_RadioIsland_Oculina","TX_Packery","TX_PortA","FL_Panacea"))

eigenvals=pp0$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)

#ggplot(pca_s, aes(MDS1, MDS2, label=pca_s$label)) +
ggplot(pca_s, aes(MDS2, MDS3)) +
  geom_point(size=3.5, pch=21,aes(fill=as.factor(pca_s$site)),colour="black")  +
  theme_classic(base_size = 20) +
  #geom_text(hjust=0, vjust=0)+
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(pca_s$site)),show.legend = NULL)+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS3 (",formatC((varexplained[3]*100), digits = 2, format = "f"),"%)"))# +
  #scale_colour_manual(values=colors, name="", labels = c("lineage1", "lineage2", "lineage3", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  #scale_fill_manual(values=colors, name="") #
#all, 23
axes2plot=c(2,3)  

pca_s <- as.data.frame(cmd$CA$u[,axes2plot])

pca_s$site=as.factor(new$pop)
pca_s$site=as.factor(pca_s$site)
levels(pca_s$site) <- list(NC_TriangleWrecks  = "NC", NC_RadioIsland = "NC_Radio",Virginia="VA",Massachusetts="WH",Rhode_Island="RI",
                           TX_PortA="Texas_PortA",TX_Packery="Texas_Packery", NC_RadioIsland_Oculina="Oculina_NC",FL_Panacea="Panacea")

colors=c('#e3dbdbff','#bae4bc','#00ccc4ff',"#2b8cbe", "#0C2C84","#66a61e","#d40000ff","#ff8500ff","#ff8080ff")
pca_s$site=factor(pca_s$site,levels=c("NC_RadioIsland","NC_TriangleWrecks","Virginia","Rhode_Island","Massachusetts","NC_RadioIsland_Oculina","TX_Packery","TX_PortA","FL_Panacea"))

eigenvals=pp0$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)

ggplot(pca_s, aes(MDS2, MDS3)) +
  geom_point(size=3.5, pch=21,aes(fill=as.factor(pca_s$site)),colour="black")  +
  theme_classic(base_size = 20) +
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(pca_s$site)),show.legend = NULL)+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS3 (",formatC((varexplained[3]*100), digits = 2, format = "f"),"%)")) +
  #scale_colour_manual(values=colors, name="", labels = c("lineage1", "lineage2", "lineage3", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  scale_fill_manual(values=colors, name="") #


