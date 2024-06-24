 #!/bin/bash -l

plink=/usr/local/sw/plink/1.9/plink
convertf=/usr/local/sw/EIG-8.0.0/bin/convertf
datasetdir=$1
paneldir=$2
glimpsedir=$3

# merging highcov imp&diplo&pshap to the HOv2_hg38 dataset in highcov_post_imputation script
# merging lowcov imputed data in glimpse2_pipeline_lowcov script

# make HOv2_hg38_highlowcov plink format, with parfile par_convert:
#genotypename:   HOv2_hg38_highlowcov.geno
#snpname:        HOv2_hg38_highlowcov.snp
#indivname:      HOv2_hg38_highlowcov.ind
#outputformat:   PACKEDPED
#genotypeoutname:        HOv2_hg38_highlowcov.bed
#snpoutname:     HOv2_hg38_highlowcov.bim
#indivoutname:   HOv2_hg38_highlowcov.fam
#outputgroup:    YES
#$convertf -p par_convert

# merge refgenome to the dataset (hov2 merged with imputed dowmsampled (original highcov) data )

##mv $paneldir/refgenome38_humanorigins.map $paneldir/refgenome38_humanorigins.ed1.map
##awk '{$2=$1":"$4; print }' OFS='\t' $paneldir/refgenome38_humanorigins.ed1.map $paneldir/refgenome38_humanorigins.map

#${plink} --file $paneldir/refgenome38_humanorigins --make-bed --out $paneldir/refgenome38_humanorigins

#$plink --bfile ./HOv2_hg38_highlowcov --bmerge $paneldir/refgenome38_humanorigins --out  ./HOv2_hg38_highlowcov_wref

# convert back to eig format, with parfile  par_convert_2:
#genotypename:   HOv2_hg38_highlowcov_wref.bed
#snpname:        HOv2_hg38_highlowcov_wref.bim
#indivname:      HOv2_hg38_highlowcov_wref.pedind
#outputformat:   EIGENSTRAT
#genotypeoutname:        HOv2_hg38_highlowcov_wref.geno
#snpoutname:     HOv2_hg38_highlowcov_wref.snp
#indivoutname:   HOv2_hg38_highlowcov_wref.ind
#outputgroup:    YES

#/usr/local/sw/EIG-8.0.0/bin/convertf -p par_convert_2


# merge lowcov pseudohaploid data to the dataset
# parfile  par_merge2:
#geno1:  $datasetdir/HOv2_hg38_highlowcov_wref.geno
#snp1:   $datasetdir/HOv2_hg38_highlowcov_wref.snp
#ind1:   $datasetdir/HOv2_hg38_highlowcov_wref.ind
#geno2:  $glimpsedir/lowcov120_hov2_38.geno
#snp2:   $glimpsedir/lowcov120_hov2_38.snp
#ind2:   $glimpsedir/lowcov120_hov2_38.ind
#genooutfilename: $datasetdir/HOv2_hg38_highlowcov_wref_all.geno
#snpoutfilename:  $datasetdir/HOv2_hg38_highlowcov_wref_all.snp
#indoutfilename:  $datasetdir/HOv2_hg38_highlowcov_wref_all.ind
#outputformat:   EIGENSTRAT
#docheck:        YES
#strandcheck:    NO
#hashcheck:      YES

#/usr/local/sw/EIG-8.0.0/bin/mergeit -p par_merge_2

# remove 10 inds
# wref dataset - no pshap lowcov merge
list=removelist.txt
#$plink  --remove $list --bfile HOv2_hg38_highlowcov --make-bed --out HOv2_hg38_highlowcov_rm10

#$plink --remove removelist2.txt --bfile $glimpsedir/lowcov120_hov2_38 --make-bed --out $glimpsedir/lowcov120_hov2_38_rm10

#$plink --bfile HOv2_hg38_highlowcov_rm10 --bmerge $paneldir/refgenome38_humanorigins --make-bed --out HOv2_hg38_highlowcov_rm10_wref

##cp HOv2_hg38_highlowcov_wref.pedind newrm.pedind

# manually edit out the ten from pedind
##mv newrm.pedind HOv2_hg38_highlowcov_rm10_wref.pedind

# convert to eig again 
#/usr/local/sw/EIG-8.0.0/bin/convertf -p par_convert_3

# convert lowcov rm data too
#$glimpsedir/lowcov120_hov2_38_rm10
#awk '{$2=sprintf($2"-ph"); $6=$2; print}' OFS='\t' $glimpsedir/lowcov120_hov2_38_rm10.fam > $glimpsedir/lowcov120_hov2_38_rm10.pedind

#/usr/local/sw/EIG-8.0.0/bin/convertf -p par_convert_4

#merge
#/usr/local/sw/EIG-8.0.0/bin/mergeit -p par_merge_3

