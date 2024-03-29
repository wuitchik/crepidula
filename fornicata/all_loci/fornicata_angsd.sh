#!/bin/bash -l

#$ -V # inherit the submission environment
#$ -cwd #start job in submission directory
#$ -N fornicata  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be

ls F*.bam > fornicata_bams

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 69 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd 1"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b fornicata_bams -GL 1 $FILTERS $TODO -P 1 -out fornicata_results