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


################ to here

angsd -sites allSites -b p_cape_may_bams -GL 1 -P 1 -minInd 16 $TODO -out PC_all
angsd -sites allSites -b p_beverly_bams -GL 1 -P 1 -minInd 15 $TODO -out PB_all
angsd -sites allSites -b p_kettle_cove_bams -GL 1 -P 1 -minInd 13 $TODO -out PK_all
angsd -sites allSites -b p_newport_bams -GL 1 -P 1 -minInd 16 $TODO -out PN_all

# generating per-population SFS
realSFS PC_all.saf.idx > PC_all.sfs
realSFS PB_all.saf.idx > PB_all.sfs
realSFS PK_all.saf.idx > PK_all.sfs
realSFS PN_all.saf.idx > PN_all.sfs

# writing down 2d-SFS priors
realSFS PC_all.saf.idx PN_all.saf.idx -P 24 > PC_all_2_PN_all.sfs ; realSFS fst index PC_all.saf.idx PN_all.saf.idx -sfs PC_all_2_PN_all.sfs -fstout PC_all_2_PN_all
realSFS PC_all.saf.idx PK_all.saf.idx -P 24 > PC_all_2_PK_all.sfs ; realSFS fst index PC_all.saf.idx PK_all.saf.idx -sfs PC_all_2_PK_all.sfs -fstout PC_all_2_PK_all
realSFS PC_all.saf.idx PB_all.saf.idx -P 24 > PC_all_2_PB_all.sfs ; realSFS fst index PC_all.saf.idx PB_all.saf.idx -sfs PC_all_2_PB_all.sfs -fstout PC_all_2_PB_all

realSFS PB_all.saf.idx PN_all.saf.idx -P 24 > PB_all_2_PN_all.sfs ; realSFS fst index PB_all.saf.idx PN_all.saf.idx -sfs PB_all_2_PN_all.sfs -fstout PB_all_2_PN_all
realSFS PB_all.saf.idx PK_all.saf.idx -P 24 > PB_all_2_PK_all.sfs ; realSFS fst index PB_all.saf.idx PK_all.saf.idx -sfs PB_all_2_PK_all.sfs -fstout PB_all_2_PK_all

realSFS PK_all.saf.idx PN_all.saf.idx -P 24 > PK_all_2_PN_all.sfs ; realSFS fst index PK_all.saf.idx PN_all.saf.idx -sfs PK_all_2_PN_all.sfs -fstout PK_all_2_PN_all



# global Fst between populations
realSFS fst stats PC_all_2_PN_all.fst.idx # 0.013560	0.007467
realSFS fst stats PC_all_2_PK_all.fst.idx # 0.012498	0.011726
realSFS fst stats PC_all_2_PB_all.fst.idx # 0.013381	0.008348
realSFS fst stats PB_all_2_PN_all.fst.idx # 0.017689	0.006908
realSFS fst stats PB_all_2_PK_all.fst.idx # 0.013411	0.010094
realSFS fst stats PK_all_2_PN_all.fst.idx # 0.016482	0.016119


