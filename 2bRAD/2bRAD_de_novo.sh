###### 2bRAD de novo ####

# raw backup
/projectnb/coral/raw_backup/crepidula

# working directory
/projectnb/coral/crepidula

# copy fastq's into working directory, my files were on lanes 5-6
cp *.gz /projectnb/coral/crepidula

# using Misha's pipeline, so we clone their scripts into our server
git clone https://github.com/z0on/2bRAD_denovo.git

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl 
chmod +x *.py
chmod +x *.R

# Let's unzip in parallel
for file in *.gz
	do echo "gunzip $file">> gunzip
done


scc6_qsub_launcher.py -N unzip -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile gunzip
qsub unzip_array.qsub


#==================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads


# run 2bRAD_trim_launch_dedup.pl
2bRAD_trim_launch_dedup.pl fastq > trims

# use that output file to run in parallel
scc6_qsub_launcher.py -N trim -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile trims
qsub trim_array.qsub

# this creates a bunch of files with specific adapter
# I renamed my files to correspond with their barcodes + adapters 
naming.sh

# quality filtering using cutadapt (see installation above)
module load cutadapt 

# for de novo analysis: removing reads with qualities at ends less than Q15
>trimse
for file in *.tr0; do
	echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;
done

scc6_qsub_launcher.py -N cutadapt_trim -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile trimse
qsub cutadapt_trim_array.qsub

# 'uniquing' ('stacking') individual trimmed fastq reads:

for file in *.trim 
	do echo "uniquerOne.pl $file > ${file}.uni" >> unii
done

scc6_qsub_launcher.py -N uniquing -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile unii
qsub uniquing_array.qsub

# collecting common tags (= major alleles)
# merging uniqued files (set minInd to >10, or >10% of total number of samples, whichever is greater)
mergeUniq.pl uni minInd=17 > all.uniq

# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab

# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta

# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
module load blast
module load cdhit/4.6.8

#!/bin/bash -l

#$ -V # inherit the submission environment
#$ -cwd #start job in submission directory
#$ -N cd_hit  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be 


cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

# -aL alignment coverage for the longer sequence,
# default 0.0 if set to 0.9, the alignment must covers 90% of the sequence

# -aS	alignment coverage for the shorter sequence, default 0.0
# if set to 0.9, the alignment must covers 90% of the sequence

# -g	1 or 0, default 0
# by cd-hit's default algorithm, a sequence is clustered to the first 
# cluster that meet the threshold (fast cluster). If set to 1, the program
# will cluster it into the most similar cluster that meet the threshold
# (accurate but slow mode)

# -c	sequence identity threshold, default 0.9
# this is the default cd-hit's "global sequence identity" calculated as:
# number of identical amino acids in alignment
# divided by the full length of the shorter sequence

# -M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
# -T	number of threads, default 1; with 0, all CPUs will be used



#------------
# making fake reference genome (of 17 chromosomes, based off of Crepidula unguiformis, Libertini et al. 2009) out of major-allele tags
# need bowtie2 and samtools for indexing

module load bowtie2
module load samtools

concatFasta.pl fasta=cdh_alltags.fas num=17 

# formatting fake genome
export GENOME_FASTA=cdh_alltags_cc.fasta
export GENOME_DICT=cdh_alltags_cc.dict 
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files 

# for denovo: map reads to fake genome: 
GENOME_FASTA=cdh_alltags_cc.fasta

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)

for file in *.trim
	do echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x cdh_alltags_cc.fasta -U $file -S ${file}.bt2.sam" >> maps
done 

scc6_qsub_launcher.py -N mapping -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile maps
qsub mapping_array.qsub

# execute all commands written to maps...

>alignmentRates
for F in `ls *.trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - mapping.o* -A 4 | tail -1 | perl -pe 's/\.o6994438.\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

scc6_qsub_launcher.py -N sam2bam -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile s2b
qsub sam2bam_array.qsub

# since I've already checked if technical replicates cluster in dendrograms, we are concatenating these files 
module load samtools

samtools merge FC12m.trim.bt2.bam FC12B.trim.bt2.bam FC12B_2.trim.bt2.bam
samtools merge FN44m.trim.bt2.bam FN44B.trim.bt2.bam FN44B_2.trim.bt2.bam 
samtools merge FR23m.trim.bt2.bam FR23C.trim.bt2.bam FR23C_2.trim.bt2.bam 
samtools merge PB22m.trim.bt2.bam PB22A.trim.bt2.bam PB22A_2.trim.bt2.bam
samtools merge PK6m.trim.bt2.bam PK6.trim.bt2.bam PK6_2.trim.bt2.bam 
samtools merge PN6m.trim.bt2.bam PN6.trim.bt2.bam PN6_2.trim.bt2.bam 


# move the extra files to a different directory

mv FC12B.trim.bt2.bam extra_bams
mv FC12B_2.trim.bt2.bam extra_bams
mv FN44B.trim.bt2.bam extra_bams
mv FN44B_2.trim.bt2.bam extra_bams
mv FR23C.trim.bt2.bam extra_bams
mv FR23C_2.trim.bt2.bam extra_bams
mv PB22A.trim.bt2.bam extra_bams
mv PB22A_2.trim.bt2.bam extra_bams
mv PB22A.trim.bt2.bam extra_bams
mv PB22A_2.trim.bt2.bam extra_bams
mv PK6.trim.bt2.bam extra_bams
mv PK6_2.trim.bt2.bam  extra_bams
mv PN6.trim.bt2.bam extra_bams
mv PN6_2.trim.bt2.bam extra_bams
mv FN33A.trim.bt2.bam extra_bams # removed because it sequenced poorly and was giving wonky results
mv PB13B.trim.bt2.bam extra_bams # removed because outlier on PCA


### Running pcangsd. 
# Had issues getting the python script to run properly on the cluster, so used a conda environment workaround. 

git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
conda env create -f environment.yml
conda activate pcangsd
conda install -c conda-forge cython
python setup.py build_ext --inplace
pip3 install -e .


# actually running pcangst

pcangsd -b fornicata_results.beagle.gz -o fornicata_pcangsd --admix 
pcangsd -b plana_results.beagle.gz -o plana_pcangsd --admix 
pcangsd -b both_species_results.beagle.gz -o both_species_pcangsd --admix 

