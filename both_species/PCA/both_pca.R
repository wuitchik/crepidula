# Principal component analysis

library(tidyverse)
library(vegan)

# metadata
expDesign = read.csv("../../expDesign.csv")

# choose either of the following two covarince matrices:
co = as.matrix(read.table("../both_species_results.covMat")) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
dimnames(co)=list(expDesign$Sample)

# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(expDesign$Site)
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~expDesign$Site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~expDesign$Site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
kmeans.axes

axes2plot=c(1,2)  

cc=pp0

colors=as.numeric(as.factor(expDesign$Species))
colpops=as.numeric(as.factor(sort(expDesign$Site)))

plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19, col = colors)
ordispider(cc,choices= axes2plot,groups=expDesign$Species,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= expDesign$Species,draw="polygon",col=colors)


ma <- as.matrix(read.table("../both_species_results.ibsMat"))
e <- eigen(ma)
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2")
ordiellipse(e$vectors[,1:2],groups=expDesign$Species,draw="polygon", col=colpops, label = T)

scores <- e$vectors[,1:2]
df <- cbind(scores,expDesign)
df <- unname(df)
colnames(df) <- c("PC1","PC2","Sample","Species","Site","Lane")

ggplot(df,aes(x=PC1,y=PC2,color=Site))+
  geom_point(size=2)+
  theme_cowplot()

