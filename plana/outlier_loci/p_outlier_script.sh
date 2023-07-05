
# run angsd
module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -hetbias_pval 1e-3 -skipTriallelic 1 -maxHetFreq 0.5 -minInd 69 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
angsd -b fornicata_bams -GL 1 $FILTERS $TODO -P 1 -out fornicata_Sites4PCAngsd

# had issues with compiling and installing pcangsd. Had to do the following steps
conda install -c conda-forge c-compiler
conda install -c conda-forge cxx-compiler
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
python setup.py build_ext --inplace
pip install --user -r requirements.txt
pip3 install -e .

## running pcangsd 
pcangsd --beagle fornicata_Sites4PCAngsd.beagle.gz --pcadapt --out fornicata_outlier_out

zcat fornicata_Sites4PCAngsd.mafs.gz| tail -n +2 | cut -f 1,2 > fornicata_Sites4PCAngsd.sites

### Run outlier_script.R to get a list of outlier loci and neutral sites

### Outliers are saved as f_outliers.txt, neutral loci are f_netural.txt

# convert bcf to vcf
module load htslib/1.16
module load bcftools/1.16
bcftools convert -O v -o fornicata_Sites4PCAngsd.vcf fornicata_Sites4PCAngsd.bcf

# index -sites parameter in angsd
angsd sites index f_outliers.txt

# Re run angsd on the outlier list
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -hetbias_pval 1e-3 -skipTriallelic 1 -maxHetFreq 0.5 -minInd 69 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
angsd -sites f_outliers.txt -b fornicata_bams -GL 1 $FILTERS $TODO -P 1 -out fornicata_outlier 

#### Neutral loci

# index -sites parameter in angsd
angsd sites index f_neutral.txt

# Re run angsd on the outlier list
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -hetbias_pval 1e-3 -skipTriallelic 1 -maxHetFreq 0.5 -minInd 69 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
angsd -sites f_neutral.txt -b fornicata_bams -GL 1 $FILTERS $TODO -P 1 -out fornicata_neutral
