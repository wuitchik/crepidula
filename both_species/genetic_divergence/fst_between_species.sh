#======= ANDSD on all well-genotyped sites: genetic diversity stats, SFS for demographic analysis

# made separate files listing bams for each population (without clones and replicates)
# populations include: beverly, cape_may, kettle, newport and robbinston


# generating list of filtered SNP sites for SFS production (note: no filters that distort allele frequency!):
# sb - strand bias filter; only use for 2bRAD, GBS or WGS (not for ddRAD or RADseq)
# hetbias - detects weird heterozygotes because they have unequal representation of alleles
# maxHetFreq - filters out lumped paralogs 
# set minInd to 80% of all your individuals (depending on the results from quality control step)

# beverly -minInd 15
# cape_may -minInd 16
# kettle -minInd 13
# newport -minInd 16


# If the goal is genome-wide diversity statistics, consider running this on a fraction of genome (a few Mb) - use angsd options -r or -rf
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 130 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2"

# ANGSD command:
angsd -b both_species_bams -GL 1 -P 1 $FILTERS $TODO -out both_sfilt

# individual heterozygosities (proportion of heterozygotes across SNPs that pass filters)
R heterozygosity_beagle.R sfilt.beagle.gz
# this script (by Nathaniel Pope) outputs an R data bundle containing AFS (rows) for each individual (columns). The proportion of heterozygotes is the second row.

# collecting and indexing filter-passing sites
zcat both_sfilt.mafs.gz | cut -f 1,2 | tail -n +2 > both_allSites
angsd sites index both_allSites

# estimating site frequency likelihoods for each population, also saving allele frequencies (for genome scan) 
export GENOME_REF=cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"

# writing down 2d-SFS priors

# within pops
realSFS FC_all.saf.idx PC_all.saf.idx -P 24 > FC_all_2_PC_all.sfs ; realSFS fst index FC_all.saf.idx PC_all.saf.idx -sfs FC_all_2_PC_all.sfs -fstout FC_all_2_PC_all
realSFS FN_all.saf.idx PN_all.saf.idx -P 24 > FN_all_2_PN_all.sfs ; realSFS fst index FN_all.saf.idx PN_all.saf.idx -sfs FN_all_2_PN_all.sfs -fstout FN_all_2_PN_all
realSFS FB_all.saf.idx PB_all.saf.idx -P 24 > FB_all_2_PB_all.sfs ; realSFS fst index FB_all.saf.idx PB_all.saf.idx -sfs FB_all_2_PB_all.sfs -fstout FB_all_2_PB_all
realSFS FK_all.saf.idx PK_all.saf.idx -P 24 > FK_all_2_PK_all.sfs ; realSFS fst index FK_all.saf.idx PK_all.saf.idx -sfs FK_all_2_PK_all.sfs -fstout FK_all_2_PK_all

# between pops

realSFS FC_all.saf.idx PN_all.saf.idx -P 24 > FC_all_2_PN_all.sfs ; realSFS fst index FC_all.saf.idx PN_all.saf.idx -sfs FC_all_2_PN_all.sfs -fstout FC_all_2_PN_all
realSFS FC_all.saf.idx PB_all.saf.idx -P 24 > FC_all_2_PB_all.sfs ; realSFS fst index FC_all.saf.idx PB_all.saf.idx -sfs FC_all_2_PB_all.sfs -fstout FC_all_2_PB_all
realSFS FC_all.saf.idx PK_all.saf.idx -P 24 > FC_all_2_PK_all.sfs ; realSFS fst index FC_all.saf.idx PK_all.saf.idx -sfs FC_all_2_PK_all.sfs -fstout FC_all_2_PK_all

realSFS FN_all.saf.idx PC_all.saf.idx -P 24 > FN_all_2_PC_all.sfs ; realSFS fst index FN_all.saf.idx PC_all.saf.idx -sfs FN_all_2_PC_all.sfs -fstout FN_all_2_PC_all
realSFS FN_all.saf.idx PB_all.saf.idx -P 24 > FN_all_2_PB_all.sfs ; realSFS fst index FN_all.saf.idx PB_all.saf.idx -sfs FN_all_2_PB_all.sfs -fstout FN_all_2_PB_all
realSFS FN_all.saf.idx PK_all.saf.idx -P 24 > FN_all_2_PK_all.sfs ; realSFS fst index FN_all.saf.idx PK_all.saf.idx -sfs FN_all_2_PK_all.sfs -fstout FN_all_2_PK_all

realSFS FB_all.saf.idx PC_all.saf.idx -P 24 > FB_all_2_PC_all.sfs ; realSFS fst index FB_all.saf.idx PC_all.saf.idx -sfs FB_all_2_PC_all.sfs -fstout FB_all_2_PC_all
realSFS FB_all.saf.idx PN_all.saf.idx -P 24 > FB_all_2_PN_all.sfs ; realSFS fst index FB_all.saf.idx PN_all.saf.idx -sfs FB_all_2_PN_all.sfs -fstout FB_all_2_PN_all
realSFS FB_all.saf.idx PK_all.saf.idx -P 24 > FB_all_2_PK_all.sfs ; realSFS fst index FB_all.saf.idx PK_all.saf.idx -sfs FB_all_2_PK_all.sfs -fstout FB_all_2_PK_all

realSFS FK_all.saf.idx PC_all.saf.idx -P 24 > FK_all_2_PC_all.sfs ; realSFS fst index FK_all.saf.idx PC_all.saf.idx -sfs FK_all_2_PC_all.sfs -fstout FK_all_2_PC_all
realSFS FK_all.saf.idx PN_all.saf.idx -P 24 > FK_all_2_PN_all.sfs ; realSFS fst index FK_all.saf.idx PN_all.saf.idx -sfs FK_all_2_PN_all.sfs -fstout FK_all_2_PN_all
realSFS FK_all.saf.idx PB_all.saf.idx -P 24 > FK_all_2_PB_all.sfs ; realSFS fst index FK_all.saf.idx PB_all.saf.idx -sfs FK_all_2_PB_all.sfs -fstout FK_all_2_PB_all


# global Fst between populations
realSFS fst stats FC_all_2_PC_all.fst.idx # 0.022151	0.567758
realSFS fst stats FN_all_2_PN_all.fst.idx # 0.029775	0.601601
realSFS fst stats FB_all_2_PB_all.fst.idx # 0.022066	0.355131
realSFS fst stats FK_all_2_PK_all.fst.idx # 0.031526	0.607662

realSFS fst stats FC_all_2_PN_all.fst.idx # 0.028166	0.598099
realSFS fst stats FC_all_2_PB_all.fst.idx # 0.017951	0.358750
realSFS fst stats FC_all_2_PK_all.fst.idx # 0.030262	0.607609

realSFS fst stats FN_all_2_PC_all.fst.idx # 0.028936	0.565859
realSFS fst stats FN_all_2_PB_all.fst.idx # 0.027976	0.379096
realSFS fst stats FN_all_2_PK_all.fst.idx # 0.030041	0.608804

realSFS fst stats FB_all_2_PC_all.fst.idx # 0.026548	0.569375
realSFS fst stats FB_all_2_PN_all.fst.idx # 0.036898	0.596257
realSFS fst stats FB_all_2_PK_all.fst.idx # 0.034837	0.604402

realSFS fst stats FK_all_2_PC_all.fst.idx # 0.028053	0.550662
realSFS fst stats FK_all_2_PN_all.fst.idx # 0.031345	0.598368
realSFS fst stats FK_all_2_PB_all.fst.idx # 0.027940	0.358385



