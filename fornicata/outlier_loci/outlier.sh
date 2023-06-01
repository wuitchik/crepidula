f=(*all*.saf.idx)
for ((i = 0; i < ${#f[@]}; i++)); do       for ((j = i + 1; j < ${#f[@]}; j++)); do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS ${f[i]} ${f[j]} > ${f[i]/%.out.saf.idx/}.${f[j]/%.out.saf.idx/}.sfs; done; done

f=(*all*.saf.idx)
for ((i = 0; i < ${#f[@]}; i++)); do       for ((j = i + 1; j < ${#f[@]}; j++)); do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS fst index ${f[i]} ${f[j]} -sfs ${f[i]/%.out.saf.idx/}.${f[j]/%.out.saf.idx/}.sfs -fstout ${f[i]/%.out.saf.idx/}.${f[j]/%.out.saf.idx/}; done; done

for i in *fst.idx; do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS  fst stats $i  >>fstoutputstats.tmp; echo $i >> fstoutputtmp; paste fstoutputtmp fstoutputstats.tmp >fstoutputstats.txt; rm -f fstoutputstats.tmp fsoutputtmp;done

