#--------------- assessing site depth
#re-did angsd without clones (removed the .bam names from 'bams' & named it 'bamscl')
#remember to change 'minInd' to new 80% (mine is now minInd 99 for 80% from 124 samples)
nano angsd_cl.sh
#added this to top:
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd #start job in submission directory
#$ -N angsd_cl.sh # job name, anything you want
#$ -l h_rt=24:00:00
#$ -M thenicolakriefall@gmail.com
#$ -m be
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 99 -snp_pval 1e-5 -minMaf 0.1 -minIndDepth 8 -hwe_pval 1e-5"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2 -dumpCounts 2"
angsd -b bamscl -GL 1 $FILTERS $TODO -P 1 -out clresult
#exit & save
# Assessing average site depth coverage per sample using vcftools:
gunzip clresult.vcf.gz
# You have to do this strange thing where you remove "(angsd version)" from the header of your vcf file, or else vcftools won't work
nano clresult.vcf
#first line reads:
##fileformat=VCFv4.2(angsd version)
# manually delete "(angsd version)" from this line so it just reads "VCFv4.2"
#exit & save
module load vcftools
vcftools --vcf clresult.vcf --depth --out clresult
#results (mean depth of sites per individual) are in a file ending in '.idepth'
#I transferred it to my desktop & arranged the column smallest -> largest in excel to check it out
#re-ran this step after removing 10 samples with an average of <7 reads, went from 2150 SNPs to 3474