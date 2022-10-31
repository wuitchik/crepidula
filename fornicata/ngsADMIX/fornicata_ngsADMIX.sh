#NgsAdmix for K from 2 to 5 : do not run if the dataset contains clones or genotyping replicates!
#for K in `seq 2 5` ;
#do
#NGSadmix -likes myresult.noLD.beagle.gz -K $K -P 10 -o mydata.noLD_k${K};
#done

# alternatively, to use real ADMIXTURE on called SNPs (requires plink and ADMIXTURE):
#gunzip myresult.noLD.vcf.gz
#module load admixture/1.3.0
#module load plink/1.90b6.4

#plink --vcf myresult.noLD.vcf --make-bed --allow-extra-chr 0 --out myresult.noLD
#for K in `seq 1 5`; \
#do admixture --cv myresult.noLD.bed $K | tee myresult.noLD_${K}.out; done

# which K is least CV error?
#grep -h CV myresult.noLD_*.out

#CV error (K=1): 0.53296
#CV error (K=2): 0.45103
#CV error (K=3): 0.43497
#CV error (K=4): 0.43349
#CV error (K=5): 0.44908