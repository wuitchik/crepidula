 ## fornicata outlier and neutral loci selection using PCANGSD

library(RcppCNPy)
library(bigutilsr)
library(qvalue)
library(GenomicRanges)
library(tidyverse)

filename="fornicata_outlier_out.pcadapt.zscores.npy"
zscores <- npyLoad(filename)

#number of PCs you care about
PCs=1
dflist2=list()
for (K in 1:PCs){
  NAME=paste("PC",K, sep="")
  print(NAME)
  dist=(zscores[,K] - median(zscores[,K]))^2
  #colnames(dist)=NAME
  #dflist[K]=dist
  pvals=data.frame(pchisq(dist, df=K, lower.tail=F))
  print(pvals)
  dflist2[K]=pvals
  names(dflist2)[K]=NAME
  
}

pval=bind_rows(dflist2, .id = "column_label")
p<-read.table("fornicata_Sites4PCAngsd.sites",colC=c("factor","integer"))
dflist=list()
for (i in 1:PCs){
  df=cbind(p,pval[,i])
  colnames(df)=c("contig","pos","pval")
  df$qval <- qvalue(df$pval)$qvalues
  alpha <- 0.005
  df_outliers=subset(df, qval < alpha)
  print(nrow(df_outliers))
  df_outliers=select(df_outliers, -c("pval","qval"))#removing p and qvals for fishers
  dflist[[i]]=df_outliers
  names(dflist)[i]=colnames(pval)[i]
}

fornicata_outliers = dflist



write.table(fornicata_outliers, "f_outliers.txt", sep = "\t", row.names = F)
