 ## fornicata outlier and neutral loci selection using PCANGSD

library(RcppCNPy)
library(bigutilsr)
library(qvalue)
library(GenomicRanges)
library(tidyverse)

filename="outlier_loci/plana_outlier_out.pcadapt.zscores.npy"
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
p<-read.table("outlier_loci/plana_Sites4PCAngsd.sites",colC=c("factor","integer"))
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


write.table(df_outliers, "p_outliers.txt", sep = "\t", row.names = FALSE)


# Neutral sites, choosing p-value > 0.05


neutral_list=list()
for (i in 1:PCs){
  df=cbind(p,pval[,i])
  colnames(df)=c("contig","pos","pval")
  df$qval <- qvalue(df$pval)$qvalues
  alpha <- 0.05
  df_neutral=subset(df, qval > alpha)
  print(nrow(df_neutral))
  df_neutral=select(df_neutral, -c("pval","qval"))#removing p and qvals for fishers
  dflist[[i]]=df_neutral
  names(dflist)[i]=colnames(pval)[i]
}

write.table(df_neutral, "p_neutral.txt", sep = "\t", row.names = FALSE)
