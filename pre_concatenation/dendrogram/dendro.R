# dendro

library(tidyverse)

# filter 1
data1 = read.table("f1_results.ibsMat") %>%
  as.matrix()

expDesign = read.delim("filtered_expDesign.csv", sep = ",", header = T)
dimnames(data1)=list(expDesign$Sample) 

hc=hclust(as.dist(data1),"ave")
plot(hc,cex=0.7)

#filter 2

data2 = read.table("f2_results.ibsMat") %>%
  as.matrix()

dimnames(data2)=list(expDesign$Sample) 
  
#data2[is.na(data2)] = 1

hc=hclust(as.dist(data2),"ave")
plot(hc,cex=0.7)

# filter 3
data3 = read.table("f3_results.ibsMat") %>%
  as.matrix()

#data3[is.na(data3)] = 1

dimnames(data3)=list(expDesign$Sample) 

hc3=hclust(as.dist(data3),"ave")
plot(hc3,cex=0.7)

# filter 4
data4 = read.table("f4_results.ibsMat") %>%
  as.matrix()

#data4[is.na(data4)] = 1

dimnames(data4)=list(expDesign$Sample) 

hc=hclust(as.dist(data4),"ave")
plot(hc,cex=0.7)

# filter 5
data5 = read.table("filter_5_results.ibsMat") %>%
  as.matrix()

#data5[is.na(data5)] = 1

#expDesign = read.delim("expDesign.csv", sep = ",", header = T)
dimnames(data5)=list(expDesign$Sample) 

hc=hclust(as.dist(data5),"ave")
plot(hc,cex=0.7)

# filter 6
data6 = read.table("f6_results.ibsMat") %>%
  as.matrix()

#data6[is.na(data6)] = 1

#expDesign = read.delim("expDesign.csv", sep = ",", header = T)
dimnames(data6)=list(expDesign$Sample) 

hc=hclust(as.dist(data6),"ave")
plot(hc,cex=0.7)


### Fornicata 

data1 = read.table("fornicata_f1.ibsMat") %>%
  as.matrix()

expDesign = read.delim("expDesign.csv", sep = ",", header = T) %>%
  filter(Species == 'Crepidula fornicata')
dimnames(data1)=list(expDesign$Sample) 

hc=hclust(as.dist(data1),"ave")
plot(hc,cex=0.7)

# 2

data2 = read.table("fornicata_f2.ibsMat") %>%
  as.matrix()

dimnames(data2)=list(expDesign$Sample) 

hc=hclust(as.dist(data2),"ave")
plot(hc,cex=0.7)

# 3

data3 = read.table("fornicata_f3.ibsMat") %>%
  as.matrix()

dimnames(data3)=list(expDesign$Sample) 

hc=hclust(as.dist(data3),"ave")
plot(hc,cex=0.7)





