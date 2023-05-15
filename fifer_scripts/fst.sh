##Analysis of genetic divergence## 
#Calculating SAF
#158 inds
# note pop1= RED, pop2=BLUE, pop3=TAN ,pop4=GREEN

module load angsd/0.923
# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 126 "
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"

angsd -b ../../bamscl_3pop -GL 1 $FILTERS $TODO -P 1 -out AllSites
angsd sites index AllSites 
export GENOME_REF=/projectnb/davieslab/jfifer/Japan_rad/ref_genome/Amil_v2.01/Amil.v2.01.chrs.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"


# In the following lines, set minInd to 75-90% of each pop's sample size
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites AllSites -b ../bamscl_pop2 -GL 1 -P 1 -minInd 16 $TODO -out pop2.out # 21
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites AllSites -b ../bamscl_pop4 -GL 1 -P 1 -minInd 34 $TODO -out pop4.out # 43
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites AllSites -b ../bamscl_northpop -GL 1 -P 1 -minInd 75 $TODO -out popN.out #94

# generating per-population SFS
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS pop2.out.saf.idx >pop2.sfs
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS popN.out.saf.idx >popN.sfs
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS pop4.out.saf.idx >pop4.sfs

# global Fst between populations
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS  fst index ../pop2.out.saf.idx ../popN.out.saf.idx -sfs ../pop2.popN.sfs -fstout pop2popN.fst
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS  fst index ../pop2.out.saf.idx ../pop4.out.saf.idx -sfs ../pop2.pop4.sfs -fstout pop2pop4.fst
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS  fst index ../pop4.out.saf.idx ../popN.out.saf.idx -sfs ../pop4.popN.sfs -fstout pop4popN.fst

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Heterozygosity##
FILTERS="-maxHetFreq 0.5 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 126"
TODO="-ref $GENOME_REF -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -bam /projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/bamscl_3pop -GL 1 -P 1 $TODO $FILTERS -out lineage_ref
#scp beagle.gz and use script XXXX.R to calculate heterozygosity 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~