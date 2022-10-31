# Principal component analysis

library(tidyverse)
library(vegan)

# metadata
expDesign = read.csv("../expDesign.csv") %>%
  filter(Species == "Crepidula fornicata") %>%
  mutate(Site = factor(Site, levels = c("Robbinston", "Kettle Cove", "Beverly", "Newport", "Cape May")))

# covariance matric
co = as.matrix(read.table("fornicata_results.covMat"))

#  PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(expDesign$Site)
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~expDesign$Quality,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~expDesign$Site,conds)

#  eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
kmeans.axes #[1] only one is really doing much
axes2plot=c(1,2)  

pca = as.data.frame(pp$CA$u[,axes2plot])

eigenvals=pp$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")

ggplot(pca, aes(MDS1, MDS2)) +
  geom_point(size=3.5, pch=21,aes(fill=expDesign$Site))  +
  theme_classic(base_size = 20) +
  xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  scale_fill_manual(values = colours) +
  guides(fill = guide_legend(title = "Site"))




## IBS

ma = as.matrix(read.table("fornicata_results.ibsMat"))

# partial PCoA, taking out the variation due to coverage
pp1=capscale(ma~1)
pp2=capscale(ma~1+Condition(expDesign$Quality))

pca1 = as.data.frame(pp1$CA$u[,axes2plot])
pca2 = as.data.frame(pp2$CA$u[,axes2plot])

eigenvals1=pp1$CA$eig
eigenvals2=pp2$CA$eig
#calc variance explained
varexplained1=eigenvals1/sum(eigenvals1)
varexplained2=eigenvals2/sum(eigenvals2)


ggplot(pca1, aes(MDS1, MDS2)) +
  geom_point(size=3.5, pch=21,aes(fill=expDesign$Site))  +
  theme_classic(base_size = 20) +
  xlab(paste0("MDS1 (",formatC((varexplained1[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained1[2]*100), digits = 2, format = "f"),"%)")) +
  scale_fill_manual(values = colours) +
  guides(fill = guide_legend(title = "Site"))

ggplot(pca2, aes(MDS1, MDS2)) +
  geom_point(size=3.5, pch=21,aes(fill=expDesign$Site))  +
  theme_classic(base_size = 20) +
  xlab(paste0("MDS1 (",formatC((varexplained2[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained2[2]*100), digits = 2, format = "f"),"%)")) +
  scale_fill_manual(values = colours) +
  guides(fill = guide_legend(title = "Site"))






