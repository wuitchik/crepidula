#===================== A  N  G  S  D =====================
module load angsd

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1700 -minInd 85"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# listing all bam filenames 
ls *bam >bams

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd 

# summarizing results (using modified script by Matteo Fumagalli)
Rscript plotQC.R prefix=dd


# proportion of sites covered at >5x:
cat quality.txt

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs


#--------------- population structure

# Note: PCA and Admixture are not supposed to be run on data that contain clones or genotyping replicates. For PCA, these can be removed without rerunning ANGSD from the IBS distance matrix; but for ngsAdmix ANGSD must be rerun.

# manually removed the following duplicate samples
FN44B_2.trim.bt2.bam
PB22A_2.trim.bt2.bam
PK6_2.trim.bt2.bam

# then separated bam list into C. fornicata and C. plana

######## C. fornicata (n = 92) #######

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 920 -minInd 45"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
angsd -b fornicata_bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd 

# summarizing results (using modified script by Matteo Fumagalli)
Rscript plotQC.R prefix = fornicata_dd


# Generating genotype likelihoods from highly confident (non-sequencing-error) SNPs
# set minInd to 75-80% of your total number of bams
# if you expect very highly differentiated populations with nearly fixed alternative alleles, remove '-hwe_pval 1e-5' form FILTERS
# -doGeno 8 : genotype likelihood format setting for ngsLD; if you want to run PCA, use -doGeno 32 (but I recommend using ibsMat for all ordination work)
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 69 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
angsd -b fornicata_bams -GL 1 $FILTERS $TODO -P 1 -out fornicata_result

# how many SNPs?
NSITES=`zcat fornicata_result.mafs.gz | wc -l`
echo $NSITES


# ------- PCAngsd: estimating admixture, kinship, and SNP covariance (note: inbreeding is better estimated as individual heterozygosities later, based on all sites, to avoid ascertainment bias)

module load python2
python /projectnb/coral/crepidula/bin/pcangsd/pcangsd/pcangsd.py -beagle fornicata_result.beagle.gz -admix -o fornicata_pcangsd -inbreed 2 -kinship -selection -threads 12



## C. plana (n = 78)


