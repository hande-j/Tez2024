#!/bin/bash -l

#bcftools=/usr/local/sw/bcftools-1.9/bcftools
tabix=/usr/local/sw/bcftools-1.9/htslib-1.9/tabix
#plink=/usr/local/sw/plink-1.90/plink
plink=/usr/local/sw/plink/1.9/plink
#convertf=/usr/local/sw/EIG-7.2.1/bin/convertf
convertf=/usr/local/sw/EIG-8.0.0/bin/convertf
#mergeit=/usr/local/sw/EIG-7.2.1/bin/mergeit
mergeit=/usr/local/sw/EIG-8.0.0/bin/mergeit

GLIMPSE_chunk=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_chunk
GLIMPSE_phase=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_phase
GLIMPSE_ligate=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_ligate
GLIMPSE_split_reference=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_split_reference
GLIMPSE_concordance=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_concordance

##  	GLIMPSE v2 Pipeline for Low Coverage Genomes
refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
glimpsedir=$2
reference_panel=$glimpsedir/reference_panel
chrlist=$3
datasetdir=$4

#while read i; do
#bamfileall=$i
#name="$(awk '{print substr($0, 45)}' <<< ${bamfileall} | awk -F [_.-] '{print $1}')"
#name="$(awk -F '/' '{print $NF}'<<< $bamfileall | awk -F [_.-] '{print $1}')"
#echo "$name-v2"
#mkdir -p sixtyfour/samplename
#echo "$name-v2" > sixtyfour/samplename/${name}.v2.txt
#done < imputenewlist_64bamlist.txt

#mkdir -p $glimpsedir/sixtyfour/${name}
outdir=$glimpsedir/sixtyfour/${name}
maindir=$glimpsedir/sixtyfour

#mkdir -p ${outdir}/GLIMPSE_impute
#mkdir -p ${outdir}/GLIMPSE_ligate
#mkdir -p ${outdir}/GLIMPSE_concordance
#mkdir -p ${outdir}/GLIMPSE_validation

## Impute and Phase a Whole Chromosome
#while read i;
#do
#chrlist=$i
#BAM=${bamfileall}
#REF=${reference_panel}/split/1000GP

#while IFS="" read -r LINE || [ -n "$LINE" ];
#do
#	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
#	IRG=$(echo $LINE | cut -d" " -f3)
#	ORG=$(echo $LINE | cut -d" " -f4)
#	CHR=$(echo ${LINE} | cut -d" " -f2 )
#	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1) 
#	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
#	OUT=${outdir}/GLIMPSE_impute/${name}.${CHR}.${REGS}.${REGE}.bcf
#
#	#${GLIMPSE_phase} --bam-file ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}
#	
#done < ${reference_panel}/split/chunks.${chrlist}.txt

##  Ligate multiple chunks together 
#LST=${outdir}/GLIMPSE_ligate/list.${chrlist}.txt
#ls -1v ${outdir}/GLIMPSE_impute/${name}.${chrlist}.*.bcf > ${LST} 
#OUT=${outdir}/GLIMPSE_ligate/${name}.${chrlist}.ligated.bcf
#
#${GLIMPSE_ligate} --input ${LST} --output ${OUT}
#
#done < $glimpsedir/chrlist

##  Merge all samples by chromosomes
#mkdir -p $maindir/imputed_vcf
#mergelist=${maindir}/all.bcflist.$chrlist
#ls -v $maindir/*/GLIMPSE_ligate/*.${chrlist}.ligated.bcf > $mergelist
#${bcftools} merge -l ${mergelist} -Oz9 -o  ${maindir}/imputed_vcf/allsamples.sixtyfour.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/allsamples.sixtyfour.$chrlist.vcf.gz

#list=$maindir/vcflist1
#ls -v $maindir/imputed_vcf/allsamples.sixtyfour.*.vcf.gz > $maindir/vcflist1
#${bcftools} concat -f ${list} -Oz9 -o ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz

## Reheader for ligated bcf files
#while read i; do
#chrlist=$i;
##$bcftools query -l $maindir/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz

#samplename=sixtyfour_v2_namelist
#${bcftools} reheader -s ${samplename} ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz -o ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.nm.vcf.gz
#${bcftools} index $maindir/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.nm.vcf.gz

## GP filtering
#${bcftools} filter ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp80.tmp.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp80.tmp.vcf.gz

#${bcftools} filter ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.99' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp99.tmp.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp99.tmp.vcf.gz

## Stats
#mkdir -p $maindir/stats
#$bcftools stats -S $maindir/sixtyfour_v2_namelist ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp99.tmp.vcf.gz > ${maindir}/stats/allsamples.sixtyfour.v2.imputed.kgp.gp99.tmp.stats.txt
#$bcftools stats -S $maindir/sixtyfour_v2_namelist ${maindir}/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.gp80.tmp.vcf.gz > ${maindir}/stats/allsamples.sixtyfour.v2.imputed.kgp.gp80.tmp.stats.txt

# 	MERGE with prev. imputed other samples  - merging all imputed lowcov samples
#$bcftools merge $maindir/imputed_vcf/allsamples.sixtyfour.v2.imputed.kgp.vcf.gz $glimpsedir/fiftysix_v2/imputed_vcf/allsamples.v2.imputed.kgp.vcf.gz -Oz9 -o $maindir/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz
#$bcftools index $maindir/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz

## 	GP filtering for all samples together
#${bcftools} filter ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp80.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp80.vcf.gz
#
#${bcftools} filter ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.90' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp90.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp90.vcf.gz
#
#${bcftools} filter ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.95' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp95.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp95.vcf.gz
#
#${bcftools} filter ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.vcf.gz -i'FORMAT/GP>=0.99' -S. -Oz9 -o ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp99.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp99.vcf.gz

# stats
samplename=$glimpsedir/hundredtwentysamples_namelist

#mkdir -p $maindir/datasets

filename='allsamples.hundredtwenty.v2.imputed.kgp'
gpnums=("80" "90" "95" "99")
for i in "${gpnums[@]}"
do
number="$i"
#number="99"
#$bcftools stats -S $samplename $maindir/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp${number}.vcf.gz > $maindir/stats/$filename.gp${number}.stats.txt
# MAF5 - MAF1
regions=$glimpsedir/fiftysix/KGP_panel/1KGpanel.maf5.sites.tsv.gz 
regions1=$glimpsedir/fiftysix/KGP_panel/1KGpanel.maf1.sites.tsv.gz 
maf='maf5'
#${bcftools} view -R ${regions} ${maindir}/imputed_vcf/$filename.gp$number.vcf.gz -Oz9 -o ${maindir}/imputed_vcf/$filename.gp$number.$maf.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/$filename.gp$number.$maf.vcf.gz
#$bcftools stats -S $samplename $maindir/imputed_vcf/$filename.gp${number}.$maf.vcf.gz > $maindir/stats/$filename.gp${number}.$maf.stats.txt

maf='maf1'
#${bcftools} view -R ${regions1} ${maindir}/imputed_vcf/$filename.gp$number.vcf.gz -Oz9 -o ${maindir}/imputed_vcf/$filename.gp$number.$maf.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/$filename.gp$number.$maf.vcf.gz
#$bcftools stats -S $samplename $maindir/imputed_vcf/$filename.gp${number}.$maf.vcf.gz > $maindir/stats/$filename.gp${number}.$maf.stats.txt

# HO for MAF and NO MAF
humanoriginsv2=$datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38_pos.bed
mafnums=("maf5" "maf1")
for j in "${mafnums[@]}"
do
maf="$j"
# MAF5-1
#${bcftools} view -R ${humanoriginsv2} ${maindir}/imputed_vcf/$filename.gp$number.$maf.vcf.gz -Ou | $bcftools norm -m +any -Oz9 -o ${maindir}/imputed_vcf/$filename.gp$number.$maf.ho.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/$filename.gp$number.$maf.ho.vcf.gz

#$bcftools stats -S $samplename $maindir/imputed_vcf/$filename.gp${number}.$maf.ho.vcf.gz > $maindir/stats/$filename.gp${number}.$maf.ho.stats.txt
#$plink --vcf ${maindir}/imputed_vcf/$filename.gp$number.$maf.ho.vcf.gz --biallelic-only strict --recode --out ${maindir}/datasets/$filename.gp$number.$maf.ho
#awk -v num="$number" -v maf="$maf" '{$2=sprintf($1"-gp"num"-"maf); print}' OFS=" " ${maindir}/datasets/$filename.gp$number.$maf.ho.ped > ${maindir}/datasets/$filename.gp$number.$maf.ho.up.ped
done

# NOMAF
#${bcftools} view -R ${humanoriginsv2} ${maindir}/imputed_vcf/$filename.gp$number.vcf.gz -Ou | $bcftools norm -m +any -Oz9 -o ${maindir}/imputed_vcf/$filename.gp$number.ho.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/$filename.gp$number.ho.vcf.gz

#$bcftools stats -S $samplename $maindir/imputed_vcf/$filename.gp${number}.ho.vcf.gz > $maindir/stats/$filename.gp${number}.ho.stats.txt
#$plink --vcf ${maindir}/imputed_vcf/$filename.gp$number.ho.vcf.gz --recode --biallelic-only strict --out ${maindir}/datasets/$filename.gp$number.ho
#awk -v num="$number" '{$2=sprintf($1"-gp"num); print}' OFS=" " ${maindir}/datasets/$filename.gp$number.ho.ped > ${maindir}/datasets/$filename.gp$number.ho.up.ped
done 

# plink datasets of some HO subset vcfs - maf and no maf, gp80 gp99
#$plink --ped ${maindir}/datasets/$filename.gp80.ho.up.ped --map ${maindir}/datasets/$filename.gp80.ho.map --merge ${maindir}/datasets/$filename.gp80.maf5.ho.up.ped ${maindir}/datasets/$filename.gp80.maf5.ho.map --recode --out ${maindir}/datasets/$filename.gp80.mafno5.ho
#$plink --ped ${maindir}/datasets/$filename.gp90.ho.up.ped --map ${maindir}/datasets/$filename.gp90.ho.map --merge ${maindir}/datasets/$filename.gp90.maf5.ho.up.ped ${maindir}/datasets/$filename.gp90.maf5.ho.map --recode --out ${maindir}/datasets/$filename.gp90.mafno5.ho
#$plink --ped ${maindir}/datasets/$filename.gp99.ho.up.ped --map ${maindir}/datasets/$filename.gp99.ho.map --merge ${maindir}/datasets/$filename.gp99.maf5.ho.up.ped ${maindir}/datasets/$filename.gp99.maf5.ho.map --recode --out ${maindir}/datasets/$filename.gp99.mafno5.ho

#maf1
#$plink --file ${maindir}/datasets/$filename.gp80.mafno5.ho --merge ${maindir}/datasets/$filename.gp80.maf1.ho.up.ped ${maindir}/datasets/$filename.gp80.maf1.ho.map --recode --out ${maindir}/datasets/$filename.gp80.allmaf.ho
#$plink --file ${maindir}/datasets/$filename.gp90.mafno5.ho --merge ${maindir}/datasets/$filename.gp90.maf1.ho.up.ped ${maindir}/datasets/$filename.gp90.maf1.ho.map --recode --out ${maindir}/datasets/$filename.gp90.allmaf.ho
#$plink --file ${maindir}/datasets/$filename.gp99.mafno5.ho --merge ${maindir}/datasets/$filename.gp99.maf1.ho.up.ped ${maindir}/datasets/$filename.gp99.maf1.ho.map --recode --out ${maindir}/datasets/$filename.gp99.allmaf.ho

#ls $maindir/datasets/*allmaf*.ped > $maindir/datasets/datasetallmafped
#ls $maindir/datasets/*allmaf*.map > $maindir/datasets/datasetallmafmap
#paste $maindir/datasets/datasetallmafped $maindir/datasets/datasetallmafmap > $maindir/datasets/datasetallmaflist
#$plink --merge-list $maindir/datasets/datasetallmaflist --recode --out $maindir/datasets/$filename.allmaf.ho
#awk '{$2=substr($2,length($2)-length($2),length($2)-4); print $0}' OFS='\t' $maindir/datasets/$filename.allmaf.ho.map > $maindir/datasets/$filename.allmaf.ho.ed.map
#mv $maindir/datasets/$filename.allmaf.ho.map $maindir/datasets/$filename.allmaf.ho.og.map
#mv $maindir/datasets/$filename.allmaf.ho.ed.map $maindir/datasets/$filename.allmaf.ho.map
#$plink --file $maindir/datasets/$filename.allmaf.ho --make-bed --out $maindir/datasets/$filename.allmaf.ho

# Make pedind for the dataset
#awk '{$6 = $2; print}' OFS='\t' $maindir/datasets/$filename.allmaf.ho.fam > $maindir/datasets/$filename.allmaf.ho.pedind


## Make EIG dataset (HOv2)
#cat << EOF > par_cv_lc_eig
#genotypename:	datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.bed
#snpname:	datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.bim
#indivname:	datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.pedind
#outputformat:	EIGENSTRAT
#genooutfilename: datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.geno
#snpoutfilename: datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.snp
#indoutfilename:	datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.ind
#outputgroup: YES
#familynames: NO
#hashcheck: NO
#allowdups:  YES
#pordercheck: NO
#EOF
#
#$convertf -p par_cv_lc_eig

# ADD DATASET TO HOV2_HIGHCOV DATASET AND MAKE HOV2_HIGHCOV_LOWCOV DATASET (highcov has pshap called data, except that all are diploid)
#cat << EOF > par_m_lc_hohc
#geno1:	$datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.geno
#snp1:	$datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.snp
#ind1:	$datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.ed.ind
#geno2:	$glimpsedir/sixtyfour/datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.geno
#snp2:	$glimpsedir/sixtyfour/datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.snp
#ind2:	$glimpsedir/sixtyfour/datasets/allsamples.hundredtwenty.v2.imputed.kgp.allmaf.ho.ind
#genooutfilename: $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highlowcov.geno
#snpoutfilename:	 $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highlowcov.snp
#indoutfilename:	 $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highlowcov.ind
#outputformat:	EIGENSTRAT
#docheck:	YES
#strandcheck:	NO
#hashcheck:	YES
#EOF
#
#$mergeit -p par_m_lc_hohc

## called pseudohaploid genotypes (pileupcaller script)
## convert the (pileupcaller called) pseudohaploid lowcov dataset to packedped to merge with master dataset with refgen
#cat << EOF > par_cv_120_bed
#genotypename:	$glimpsedir/lowcov120_hov2_38.geno
#snpname:	$glimpsedir/lowcov120_hov2_38.snp
#indivname:	$glimpsedir/lowcov120_hov2_38.ind
#outputformat:	PACKEDPED
#genooutfilename: $glimpsedir/lowcov120_hov2_38.bed
#snpoutfilename: $glimpsedir/lowcov120_hov2_38.bim
#indoutfilename: $glimpsedir/lowcov120_hov2_38.pedind
#outputgroup: YES
#familynames: NO
#hashcheck: NO
#allowdups:  YES
#pordercheck: NO
#EOF
#
#$convertf -p par_cv_120_bed 

#edit ind, merge at masterdataset directory
#mv lowcov120_hov2_38.ind lowcov120_hov2_38.old.ind
#awk '{$1=sprintf($1"-ph"); $3=$1; print}' OFS='\t' lowcov120_hov2_38.old.ind > lowcov120_hov2_38.ind
