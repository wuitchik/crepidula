
setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Astrangia&Oculina/Admixture")

# assembling the input table
dir=("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Astrangia&Oculina/Admixture/") # path to input files
#ngsAdmix
#inName="mydata_k4_r6.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
#inName="myresult2.6.Q"
inName="myresult_dep1.4.Q"fin
k=4

samples=data.frame(read.csv("../Adapter2Sample.csv", header=T))#[,3:4])
colnames(samples)<-c("Sample","pop","Adapter")

samples <- data.frame(lapply(samples, function(x) {
  gsub("trim$", "trim.bt2.bam", x)
}))

row.names(samples)=samples$Adapter

bams=data.frame(read.csv("../bamscl", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"


#------------

#npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))

tbl=read.table(paste(dir,inName,sep=""),header=F)
row.names(tbl)= bams$Adapter

#new<-merge(tbl,i2p,by="row.names")
new1<-transform(merge(tbl,samples,by=0), row.names=Row.names, Row.names=NULL)

#create a bspsops that is ordered for demo analysis
# ordered_bspops=new1[,6:7]
# ordered_bspops <- ordered_bspops[order(ordered_bspops$pop),]
# write.table(ordered_bspops, file="bspops_ordered", sep = '\t', row.names = F, 
#             col.names = F, 
#             quote = F)

#create a bspops that is atlantic only
#atlantic=subset(new1,pop!="Texas_Packery" & pop!="Texas_PortA" & pop!="Panacea" & pop!="Oculina_NC")[,6:7]
#atlantic <- atlantic[order(atlantic$pop),]
#write.table(atlantic, file="bspops_atlantic", sep = '\t', row.names = F, 
#            col.names = F, 
#            quote = F)


tbl=new1[1:(2+k)]


colnames(tbl)=c(paste("V",c(1:(k)),  sep=""), "ind","pop")
#tbl$ind=tbl$adapter
head(tbl)
tbl$pop=as.factor(tbl$pop)


# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
levels(tbl$pop) <- list(NC_TriangleWrecks  = "NC", NC_RadioIsland = "NC_Radio",Virginia="VA",Massachusetts="WH",Rhode_Island="RI",
                        TX_PortA="Texas_PortA",TX_Packery="Texas_Packery", NC_RadioIsland_Oculina="Oculina_NC",FL_Panacea="Panacea")

tbl$pop=factor(tbl$pop,levels=c("TX_PortA","TX_Packery","FL_Panacea","NC_RadioIsland","NC_RadioIsland_Oculina","NC_TriangleWrecks","Virginia","Rhode_Island","Massachusetts"))
#colors=c("#7fabd3ff","#00bfffff", "#6af2ffff","#325fa2ff","#00ff00ff","#006400FF","#325fa2ff")
colors=c("#66a61e","#d7191c","#2b83ba","yellow", "grey","pink")

#3,2,4,1
#tbl$pop=factor(tbl$pop, levels=c("Sekisei","Oura Bay","Kushima","Kochi","Shirahama","Kushimoto","Amakusa"))
#Download below from https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
source("C:/Users/james/My Drive/Documents/BOSTON/Davies/Range expansion/Code/plot_admixture_v5_function.R")
ords=plotAdmixture(data=tbl,npops=k,grouping.method="distance", vshift=-.65,colors = colors)


# recording cluster affiliations
cluster.admix=apply(tbl[,1:k],1,function(x) {return(paste(which(x>0.1),collapse=".")) })
cluster.assign=as.data.frame(cluster.admix)
write.csv(cluster.assign, file="cluster.assign.0.1.ngsADMIX.csv")
#Gulf Astrangia =1, Atlantic Astrangia=4, Atlantic Oculina=3, Gulf Oculina=2
write.table(row.names(subset(cluster.assign, cluster.admix=='1')), file=file("bamscl_GulfAst"), quote=FALSE,col.names = F,row.names = F)
write.table(row.names(subset(cluster.assign, cluster.admix=='4')), file=file("bamscl_AtlAst"), quote=FALSE,col.names = F,row.names = F)
write.table(row.names(subset(cluster.assign, cluster.admix=='2')), file=file("bamscl_GulfOc"), quote=FALSE,col.names = F,row.names = F)
write.table(row.names(subset(cluster.assign, cluster.admix=='3')), file=file("bamscl_AtlOC"), quote=FALSE,col.names = F,row.names = F)

save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))

###################
#### Atlantic Only ###########
setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Code/Admixture/")

# assembling the input table
dir=("C:/Users/james/My Drive/Documents/BOSTON/Davies/Astrangia/Code/Admixture/") # path to input files
#ngsAdmix
#inName="mydata_k4_r6.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
inName="myresult_Atlantic.vcf.3.Q"
inName="myresult_Atlantic.vcf.5.Q"


samples=data.frame(read.csv("../Adapter2Sample.csv", header=T))#[,3:4])
colnames(samples)<-c("Sample","pop","Adapter")

samples <- data.frame(lapply(samples, function(x) {
  gsub("trim", "trim.bt2.bam", x)
}))

row.names(samples)=samples$Adapter

bams=data.frame(read.csv("./bamscl_atlantic", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
bams <- data.frame(lapply(bams, function(x) {
  gsub(".*/", "", x)
}))

#------------

tbl=read.table(paste(dir,inName,sep=""),header=F)
row.names(tbl)= bams$Adapter
new1<-transform(merge(tbl,samples,by=0), row.names=Row.names, Row.names=NULL)



#k=5
tbl=new1[1:7]
colnames(tbl)=c("V1","V2","V3","V4","V5","ind","pop")

#k=4
tbl=new1[1:6]
#k=3
tbl=new1[1:5]
head(tbl)
colnames(tbl)=c("V1","V2","V3","ind","pop")

#k=4
# colnames(tbl)=c("V1","V2","V3","V4","ind","pop")
# #tbl$ind=tbl$adapter
# head(tbl)

tbl$pop=as.factor(tbl$pop)
levels(tbl$pop) <- list(NC_TriangleWrecks  = "NC", NC_RadioIsland = "NC_Radio",Virginia="VA",MA="WH",RI="RI")
tbl$pop <- factor(tbl$pop, levels = c("NC_RadioIsland","NC_TriangleWrecks","Virginia", "RI", "MA"))


#K number equals pops number
npops=3

colors=c("#0C2C84","#e3dbdbff","#00ccc4ff")
#switch out the yellow for #66a61e

#3,2,4,1
#tbl$pop=factor(tbl$pop, levels=c("Sekisei","Oura Bay","Kushima","Kochi","Shirahama","Kushimoto","Amakusa"))
#Download below from https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
source("C:/Users/james/My Drive/Documents/BOSTON/Davies/Range expansion/Code/plot_admixture_v6_function.R")
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance", vshift=-.12, colors=colors)


# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.5),collapse=".")) })
cluster.assign=as.data.frame(cluster.admix)
write.csv(cluster.assign, file="cluster.assign.csv")
write.csv(cluster.assign, file="cluster.assign.HALF.ngsADMIX.csv")
save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))
