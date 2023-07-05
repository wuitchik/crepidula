FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 14 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 3"
# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b bams_pcyl -GL 1 $FILTERS $TODO -P 1 -out myresult_pcyl_withclones4relate

zcat myresult_pcyl_withclones4relate.mafs.gz | cut -f5 |sed 1d >freq
NIND=`cat bams_pcyl | wc -l`
/projectnb/davieslab/jfifer/BIN/ngsRelate/ngsRelate -f freq -g myresult_pcyl_withclones4relate.glf.gz -n $NIND -z bams_pcyl -O output.res

#get kinship cooefs
cut -f1,2,3,4,18 output.res