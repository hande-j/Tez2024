#!/bin/bash -l

datasetdir=$1
# f4 with imputed and genotyped data

#par_f4_impvsgeno :
#genotypename:   $datasetdir/HOv2_hg38_highcov.geno
#snpname:        $datasetdir/HOv2_hg38_highcov.snp
#indivname:      $datasetdir/HOv2_hg38_highcov.ed.ind
#popfilename:    combsimpvsgeno_1709.txt
#f4mode:         YES


# comb file looks like:
#Mbuti English AKT16-ph AKT16-1x-gp80
#Mbuti Finnish AKT16-ph AKT16-1x-gp80
#Mbuti Italian_North AKT16-ph AKT16-1x-gp80

#run
/usr/local/sw/AdmixTools-7.0.2/bin/qpDstat -p par_f4_impvsgeno  > result_f4_impvsgeno


# f4 with refgenome-added dataset
# parfile with new dataset with refgenome
# par_f4_imp_wref
#genotypename:   $datasetdir/HOv2_hg38_highlowcov_wref.bed
#snpname:        $datasetdir/HOv2_hg38_highlowcov_wref.bim
#indivname:      $datasetdir/HOv2_hg38_highlowcov_wref.pedind
#popfilename:    combsimpvsgeno_wref.txt
#f4mode:         YES


# comb file looks like:
#Mbuti Href AKT16-ph AKT16-1x-gp80
#Mbuti Href AKT16-valff AKT16-1x-gp80
# run
/usr/local/sw/AdmixTools-7.0.2/bin/qpDstat -p par_f4_imp_wref > result_f4_impcov_wref
 