#!/bin/bash -l

#bcftools=/usr/local/sw/bcftools-1.9/bcftools
#bcftools=/usr/local/sw/bcftools-1.11/bcftools
bcftools=/usr/local/sw/bcftools-1.18/bcftools
plink=/usr/local/sw/plink-1.90/plink
convertf=/usr/local/sw/EIG-7.2.1/bin/convertf
mergeit=/usr/local/sw/EIG-7.2.1/bin/mergeit
glimpsedir=$1
maindir=$glimpsedir/highcov_samples
datasetdir=$2

# the files used for concordance- merge nine samples at four depths
#chrlist="chr22"
#arrdepth=("0.1x" "0.25x" "0.5x" "1x")
#for i in "${arrdepth[@]}"
#do
#depth="$i"
#depth="0.1x"
#input1=$maindir/GLIMPSE_ligate/cov_${depth}/imputed.$chrlist.$depth.ligated.old.bcf
#input2=$maindir/GLIMPSE_ligate/cov_${depth}_2/imputed_2.$chrlist.$depth.ligated.bcf
#$bcftools merge $input1 $input2 -Ob -o $maindir/GLIMPSE_ligate/imputedten.$chrlist.$depth.bcf
#$bcftools index $maindir/GLIMPSE_ligate/imputedten.$chrlist.$depth.bcf
#done

# separate WC1ds by sample
#chr=$1
#arrdepth=("0.1x" "0.25x" "0.5x" "1x")
#for i in "${arrdepth[@]}"
#do
#depth="$i"
#outfile=$maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.$depth.ligated.hdr.bcf
#infile=$maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.$depth.ligated.bcf
#$bcftools reheader -s samplename/WC1.txt -o $outfile $infile
#$bcftools index $outfile
#done

#$bcftools view -s WC1-0.1x $maindir/GLIMPSE_ligate/WC1ds/WC1ds.imputed.$chr.ligated.bcf -Ob -o $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.1x.ligated.bcf
#$bcftools view -s WC1-0.25x $maindir/GLIMPSE_ligate/WC1ds/WC1ds.imputed.$chr.ligated.bcf -Ob -o $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.25x.ligated.bcf
#$bcftools view -s WC1-0.5x $maindir/GLIMPSE_ligate/WC1ds/WC1ds.imputed.$chr.ligated.bcf -Ob -o $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.5x.ligated.bcf
#$bcftools view -s WC1-1x $maindir/GLIMPSE_ligate/WC1ds/WC1ds.imputed.$chr.ligated.bcf -Ob -o $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.1x.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.1x.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.1x.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.25x.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chr.0.5x.ligated.bcf

#input1=$maindir/GLIMPSE_ligate/imputedten.$chrlist.$depth.bcf
#input2=$maindir/GLIMPSE_ligate/WC1ds/WC1.imputed.$chrlist.$depth.ligated.hdr.bcf
#$bcftools merge $input1 $input2 -Ob -o $maindir/GLIMPSE_ligate/imputedallten.$chrlist.$depth.bcf
#$bcftools index $maindir/GLIMPSE_ligate/imputedallten.$chrlist.$depth.bcf

# Separete Each Sample and Depth for Concordance - SampleNames Do not Have Depth in Them
#arrdepth=("0.1x" "0.25x" "0.5x" "1x")
#for i in "${arrdepth[@]}"
#do
#depth="$i"
#infile=$maindir/GLIMPSE_ligate/imputedallten.chr22.$depth.bcf
#basefile=$maindir/GLIMPSE_ligate/imputedallten.chr22.$depth
#while read j;
#do
#name="$j"
#outfile=$maindir/$name/GLIMPSE_ligate/$name.$depth.v2.kgp.$chrlist.ligated.nm.bcf
#$bcftools view -S $maindir/samplename/$name.txt -Oz9 -o $outfile $infile
#$bcftools index $outfile
#done < tensamplevalidnamelist
#done

## Merged Imputed Files with Depth in SampleName
# GP filters for Merged Imputed File - Files Have Depth in SampleName
#infile=$maindir/imputed_vcf/downsampledten_imputed.vcf.gz
#basefile=$maindir/imputed_vcf/downsampledten_imputed

#number='80'
#${bcftools} filter $infile -i'FORMAT/GP>=0.80' -S. -Oz9 -o $basefile.gp$number.vcf.gz
#$bcftools index $basefile.gp$number.vcf.gz
#
#number='90'
#${bcftools} filter $infile -i'FORMAT/GP>=0.90' -S. -Oz9 -o $basefile.gp$number.vcf.gz
#$bcftools index $basefile.gp$number.vcf.gz
#
#number='95'
#${bcftools} filter $infile -i'FORMAT/GP>=0.95' -S. -Oz9 -o $basefile.gp$number.vcf.gz
#$bcftools index $basefile.gp$number.vcf.gz
#
#number='99'
#${bcftools} filter $infile -i'FORMAT/GP>=0.99' -S. -Oz9 -o $basefile.gp$number.vcf.gz
#$bcftools index $basefile.gp$number.vcf.gz

# Stats for GP filters
#$bcftools query -l $infile > imputed_vcf/downsamplelist.txt

#array=("80" "90" "95" "99")
#for i in "${array[@]}"
#do
#number="$i"
#$bcftools stats -S imputed_vcf/downsamplelist.txt  $basefile.gp$number.vcf.gz > $basefile.gp$number.stats.txt
#done

# Reheader gp99 data to not mix with gp80
#awk '{$2=sprintf($1"-gp80"); print}' imputed_vcf/downsamplelist.txt > imputed_vcf/names.gp80.txt
#awk '{$2=sprintf($1"-gp99"); print}' imputed_vcf/downsamplelist.txt > imputed_vcf/names.gp99.txt

#$bcftools reheader -s imputed_vcf/names.gp80.txt ${maindir}/imputed_vcf/downsampledten_imputed.gp80.vcf.gz -o ${maindir}/imputed_vcf/downsampledten_imputed.gp80.nm.vcf.gz
#$bcftools reheader -s imputed_vcf/names.gp99.txt ${maindir}/imputed_vcf/downsampledten_imputed.gp99.vcf.gz -o ${maindir}/imputed_vcf/downsampledten_imputed.gp99.nm.vcf.gz
#$bcftools index ${maindir}/imputed_vcf/downsampledten_imputed.gp80.nm.vcf.gz
#$bcftools index ${maindir}/imputed_vcf/downsampledten_imputed.gp99.nm.vcf.gz

#  MAF Filter
# Subset KGP-MAF5 Positions
regions=$glimpsedir/fiftysix/KGP_panel/1KGpanel_new.maf5.sites.tsv.gz
#${bcftools} view -R ${regions} ${maindir}/imputed_vcf/downsampledten_imputed.gp80.nm.vcf.gz -Oz9 -o ${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.vcf.gz
#${bcftools} view -R ${regions} ${maindir}/imputed_vcf/downsampledten_imputed.gp99.nm.vcf.gz -Oz9 -o ${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.vcf.gz


# Subset HumanOrigins_v2_hg38 positions from already gp(80 or 99) and maf5 filtered data 
humanoriginsv2=$datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38_pos.bed
#${bcftools} view ${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.vcf.gz -R ${humanoriginsv2} -Oz9 -o ${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.hov2.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.hov2.vcf.gz

#${bcftools} view ${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.vcf.gz -R ${humanoriginsv2} -Oz9 -o ${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.hov2.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.hov2.vcf.gz

# merge two files
infile1=${maindir}/imputed_vcf/downsampledten_imputed.gp80.maf5.hov2.vcf.gz
infile2=${maindir}/imputed_vcf/downsampledten_imputed.gp99.maf5.hov2.vcf.gz
#$bcftools merge $infile1 $infile2 -Oz9 -o ${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.vcf.gz
#$bcftools index ${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.vcf.gz
 
# merge separated (to biallelic) multiallelic snps to remove them later with plink 
#$bcftools norm -m +any  ${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.vcf.gz -Oz9 -o ${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.multi.vcf.gz
#$bcftools index ${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.multi.vcf.gz

# make dataset
file=${maindir}/imputed_vcf/highcovten_imputedgp8099.maf5.hov2.multi.vcf.gz
#$plink --vcf $file --recode --biallelic-only strict --out ${maindir}/imputed_vcf/highcovten.imputedgp8099.maf5.hov2

# edit map file
#mv imputed_vcf/highcovten.imputedgp8099.maf5.hov2.map imputed_vcf/highcovten.imputedgp8099.maf5.hov2.og.map
#awk '{$1="chr"$1; $2=substr($2,0,length($2)-4); print}' OFS='\t' imputed_vcf/highcovten.imputedgp8099.maf5.hov2.og.map  > imputed_vcf/highcovten.imputedgp8099.maf5.hov2.map
#$plink --file ${maindir}/imputed_vcf/highcovten.imputedgp8099.maf5.hov2 --make-bed --out ${maindir}/imputed_vcf/highcovten.imputedgp8099.maf5.hov2 

# Diploid Called and Filtered (Validation) Data, merge multiallelics, make plink dataset
file=$maindir/imputed_vcf/validation.allten.filtered2.hov2.id.vcf.gz
#$bcftools norm -m +any $file -Oz9 -o $maindir/imputed_vcf/validation.allten.filtered2.hov2.id.multi.vcf.gz
#$bcftools index $maindir/imputed_vcf/validation.allten.filtered2.hov2.id.multi.vcf.gz

#file=$maindir/imputed_vcf/validation.allten.filtered2.hov2.id.multi.vcf.gz
#$plink --vcf $file --recode --out $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2
#mv $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2.map $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2.og.map
#awk '{$1="chr"$1; print}' OFS='\t' $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2.og.map > $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2.map
#$plink --file $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2 --make-bed --out $maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2

# merge imputed and valid(diploid-filtered) 
file1=$maindir/imputed_vcf/highcovten.imputedgp8099.maf5.hov2
file2=$maindir/imputed_vcf/highcovten.valid.diplo.ft.hov2
file3=$maindir/imputed_vcf/highcovten.imputed.valid.diplo.hov2
#$plink --bfile $file1 --bmerge $file2 --out $file3

# convert to eigenstrat format
#$convertf -p par_convert_206

# merge this imputed and diploid dataset with pseudohaploid dataset (downsampled versions included)
#$mergeit -p par_merge_107
# geno1: datasets_hov2/highcovten.imputed.valid.geno
# snp1: datasets_hov2/highcovten.imputed.valid.snp
# ind1: datasets_hov2/highcovten.imputed.valid.ind
# geno2: datasets_hov2/highcov.fifty.haplo.geno.txt
# snp2: datasets_hov2/highcov.fifty.haplo.snp.txt
# ind2: datasets_hov2/highcov.fifty.haplo.ind.txt
# genooutfilename: datasets_hov2/highcov.merge.geno
# snpoutfilename: datasets_hov2/highcov.merge.snp
# indoutfilename: datasets_hov2/highcov.merge.ind
# outputformat: EIGENSTRAT
# docheck: YES
# strandcheck: NO
# hashcheck: YES

# merge the samples dataset with humanorigins_v2_hg38 dataset
#geno1: $datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38.geno
#snp1: $datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38.snp
#ind1: $datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38.ind
#geno2: $glimpsedir/highcov_samples/datasets_hov2/highcov.merge.geno
#snp2: $glimpsedir/highcov_samples/datasets_hov2/highcov.merge.snp
#ind2: $glimpsedir/highcov_samples/datasets_hov2/highcov.merge.ind
#genooutfilename: $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.geno
#snpoutfilename: $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.snp
#indoutfilename: $datasetdir/HumanOrigins_v2_hg38_masterdataset/HOv2_hg38_highcov.ind
#outputformat: EIGENSTRAT
#docheck: YES
#strandcheck: NO
#hashcheck: YES
#$mergeit -p par_merge_hov238_highcov


# bcftools stats for chr22
samplenames=$maindir/tensamplevalidnamelist

depths=("0.1x" "0.25x" "0.5x" "1x")
array=("80" "90" "95" "99")
for j in "${depths[@]}"
do
depth="$j"
	for i in "${array[@]}"
	do 
	num="$i"

#$bcftools stats -S $samplenames $maindir/GLIMPSE_ligate/imputedallten.chr22.${depth}.gp${num}.vcf.gz > $maindir/imputed_vcf/highten.chr22.$depth.gp$num.stats.txt

	done
done
