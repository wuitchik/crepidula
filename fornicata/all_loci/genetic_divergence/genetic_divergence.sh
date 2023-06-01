#======= ANDSD on all well-genotyped sites: genetic diversity stats, SFS for demographic analysis

# made separate files listing bams for each population (without clones and replicates)
# populations include: beverly, cape_may, kettle, newport and robbinston


# generating list of filtered SNP sites for SFS production (note: no filters that distort allele frequency!):
# sb - strand bias filter; only use for 2bRAD, GBS or WGS (not for ddRAD or RADseq)
# hetbias - detects weird heterozygotes because they have unequal representation of alleles
# maxHetFreq - filters out lumped paralogs 
# set minInd to 80% of all your individuals (depending on the results from quality control step)

# beverly -minInd 13
# cape_may -minInd 16
# kettle -minInd 14
# newport -minInd 15
# robbinston -mind 13

# If the goal is genome-wide diversity statistics, consider running this on a fraction of genome (a few Mb) - use angsd options -r or -rf
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 69 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2"

# ANGSD command:
angsd -b fornicata_bams -GL 1 -P 1 $FILTERS $TODO -out sfilt

# individual heterozygosities (proportion of heterozygotes across SNPs that pass filters)
R heterozygosity_beagle.R sfilt.beagle.gz
# this script (by Nathaniel Pope) outputs an R data bundle containing AFS (rows) for each individual (columns). The proportion of heterozygotes is the second row.

# collecting and indexing filter-passing sites
zcat sfilt.mafs.gz | cut -f 1,2 | tail -n +2 >allSites
angsd sites index allSites

# estimating site frequency likelihoods for each population, also saving allele frequencies (for genome scan) 
export GENOME_REF=cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"

angsd -sites allSites -b f_cape_may_bams -GL 1 -P 1 $TODO -out FC_all
angsd -sites allSites -b f_beverly_bams -GL 1 -P 1 $TODO -out FB_all
angsd -sites allSites -b FK_all_bams -GL 1 -P 1 $TODO -out FK_all
angsd -sites allSites -b f_newport_bams -GL 1 -P 1 $TODO -out FN_all
angsd -sites allSites -b f_robbinston_bams -GL 1 -P 1 $TODO -out FR_all

# generating per-population SFS
realSFS FC_all.saf.idx > FC_all.sfs
realSFS FB_all.saf.idx > FB_all.sfs
realSFS FK_all.saf.idx > FK_all.sfs
realSFS FN_all.saf.idx > FN_all.sfs
realSFS FR_all.saf.idx > FR_all.sfs

# writing down 2d-SFS priors
realSFS FC_all.saf.idx FR_all.saf.idx -P 24 > FC_all_2_FR_all.sfs ; realSFS fst index FC_all.saf.idx FR_all.saf.idx -sfs FC_all_2_FR_all.sfs -fstout FC_all_2_FR_all 
realSFS FC_all.saf.idx FN_all.saf.idx -P 24 > FC_all_2_FN_all.sfs ; realSFS fst index FC_all.saf.idx FN_all.saf.idx -sfs FC_all_2_FN_all.sfs -fstout FC_all_2_FN_all
realSFS FC_all.saf.idx FK_all.saf.idx -P 24 > FC_all_2_FK_all.sfs ; realSFS fst index FC_all.saf.idx FK_all.saf.idx -sfs FC_all_2_FK_all.sfs -fstout FC_all_2_FK_all
realSFS FC_all.saf.idx FB_all.saf.idx -P 24 > FC_all_2_FB_all.sfs ; realSFS fst index FC_all.saf.idx FB_all.saf.idx -sfs FC_all_2_FB_all.sfs -fstout FC_all_2_FB_all

realSFS FB_all.saf.idx FR_all.saf.idx -P 24 > FB_all_2_FR_all.sfs ; realSFS fst index FB_all.saf.idx FR_all.saf.idx -sfs FB_all_2_FR_all.sfs -fstout FB_all_2_FR_all 
realSFS FB_all.saf.idx FN_all.saf.idx -P 24 > FB_all_2_FN_all.sfs ; realSFS fst index FB_all.saf.idx FN_all.saf.idx -sfs FB_all_2_FN_all.sfs -fstout FB_all_2_FN_all
realSFS FB_all.saf.idx FK_all.saf.idx -P 24 > FB_all_2_FK_all.sfs ; realSFS fst index FB_all.saf.idx FK_all.saf.idx -sfs FB_all_2_FK_all.sfs -fstout FB_all_2_FK_all

realSFS FK_all.saf.idx FR_all.saf.idx -P 24 > FK_all_2_FR_all.sfs ; realSFS fst index FK_all.saf.idx FR_all.saf.idx -sfs FK_all_2_FR_all.sfs -fstout FK_all_2_FR_all
realSFS FK_all.saf.idx FN_all.saf.idx -P 24 > FK_all_2_FN_all.sfs ; realSFS fst index FK_all.saf.idx FN_all.saf.idx -sfs FK_all_2_FN_all.sfs -fstout FK_all_2_FN_all

realSFS FR_all.saf.idx FN_all.saf.idx -P 24 > FR_all_2_FN_all.sfs ; realSFS fst index FR_all.saf.idx FN_all.saf.idx -sfs FR_all_2_FN_all.sfs -fstout FR_all_2_FN_all


# global Fst between populations
realSFS fst stats FC_all_2_FR_all.fst.idx # 0.018727	0.028121
realSFS fst stats FC_all_2_FN_all.fst.idx # 0.014067	0.020421
realSFS fst stats FC_all_2_FK_all.fst.idx # 0.015018	0.028155
realSFS fst stats FC_all_2_FB_all.fst.idx # 0.016777	0.015620
realSFS fst stats FB_all_2_FR_all.fst.idx # 0.017519	0.018382
realSFS fst stats FB_all_2_FN_all.fst.idx # 0.016606	0.023183
realSFS fst stats FB_all_2_FK_all.fst.idx # 0.016747	0.020613
realSFS fst stats FK_all_2_FR_all.fst.idx # 0.018477	0.019868
realSFS fst stats FK_all_2_FN_all.fst.idx # 0.014750	0.022944
realSFS fst stats FR_all_2_FN_all.fst.idx # 0.017396	0.023154


