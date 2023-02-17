
module load ngsadmix
#NgsAdmix for K from 2 to 5 : do not run if the dataset contains clones or genotyping replicates!
for K in `seq 2 5` ;
do
NGSadmix -likes plana_results.beagle.gz -K $K -P 10 -o plana_k${K};
done

# which K is least CV error?
#grep -h CV myresult.noLD_*.out

#CV error (K=1): 0.53296
#CV error (K=2): 0.45103
#CV error (K=3): 0.43497
#CV error (K=4): 0.43349
#CV error (K=5): 0.44908