#!/bin/bash -l

smartpca=/usr/local/sw/EIG-8.0.0/bin/smartpca 
datasetdir=$1

#grep -v Ign ../HumanOrigins_v2/HumanOrigins2583.pops_uniq > ./HumanOrigins2583.pops.noignore


#par file
#genotypename:   $datasetdir/HOv2_hg38_highlowcov_wref_all.geno
#snpname:        $datasetdir/HOv2_hg38_highlowcov_wref_all.snp
#indivname:      $datasetdir/HOv2_hg38_highlowcov_wref_all.ind
#evecoutname:    $datasetdir/pca/HOv2_hg38_highlowcov_wref_all.evec
#evaloutname:    $datasetdir/pca/HOv2_hg38_highlowcov_wref_all.eval
#outliermode:    2
#numoutevec:     10
#poplistname:    $datasetdir/HumanOrigins2583.pops.noignore
#numthreads:     20

$smartpca -p par_pca_highlow_wref_all > pca_log_highlow_wref1
