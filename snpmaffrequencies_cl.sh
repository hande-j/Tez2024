#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 1
#SBATCH -t 5-00:00:00
#SBATCH -J freq
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

bcftools=/usr/local/sw/bcftools-1.18/bcftools
plink=/usr/local/sw/plink/1.9/plink
outdir=/mnt/NEOGENE4/projects/medical_2020/glimpse
maindir=$outdir/sixtyfour

humanoriginspos=/mnt/NEOGENE4/projects/medical_2020/SNP_datasets/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38_pos.bed

samplesv2=$maindir/freqs/anatolianeo_24_over0.1xlist_v2.txt
file=$outdir/sixtyfour/imputed_vcf/allsamples.hundredtwenty.v2.imputed.kgp.gp80.vcf.gz

#$bcftools view $file -S $samplesv2 -Oz9 -o $outdir/sixtyfour/freqs/imputed_anatolianeo24_gp80.vcf.gz
#$bcftools index $outdir/sixtyfour/freqs/imputed_anatolianeo24_gp80.vcf.gz

#$bcftools view $file -S $samplesv2 -R $humanoriginspos -Oz9 -o $outdir/sixtyfour/freqs/imputed_anatolianeo24_gp99.hov2.vcf.gz
#$bcftools index $outdir/sixtyfour/freqs/imputed_anatolianeo24_gp99.hov2.vcf.gz


# merge w project neo samples > 0.1x  - april 2024

# first get the samples from project samples in a separate file
samplesneo=projectneo0.1xlist.txt
file=$outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp80.vcf.gz
#$bcftools view $file -S $samplesneo -Oz9 -o $outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp80.neo.vcf.gz
#$bcftools index $outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp80.neo.vcf.gz

file1=$outdir/sixtyfour/freqs/imputed_anatolianeo24_gp80.vcf.gz
file2=$outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp80.neo.vcf.gz
#$bcftools merge $file1 $file2 -Oz9 -o $maindir/freqs/anatolianeo41.v2.imputed.gp80.vcf.gz

# new subset vcf for gp99
file=$outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp99.vcf.gz
#$bcftools view $file -S $samplesneo -R $humanoriginspos -Oz9 -o $outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp99.neo.hov2.vcf.gz
#$bcftools index $outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp99.neo.hov2.vcf.gz

# merge gp99
file1=$outdir/sixtyfour/freqs/imputed_anatolianeo24_gp99.hov2.vcf.gz
file2=$outdir/project_samples_v2/imputed_vcf/projectallsamples.v2.imputed.kgp.gp99.neo.hov2.vcf.gz
#$bcftools merge $file1 $file2 -Oz9 -o $maindir/freqs/anatolianeo41.v2.imputed.gp99.hov2.vcf.gz



# 1 plink freq of imputed data

freqdir=$maindir/freqs
	
	# A- gp80

file=$freqdir/anatolianeo41.v2.imputed.gp80.vcf.gz
#$plink --vcf $file --recode --out $freqdir/imputed_anatolianeo41_gp80
#awk '{$2 = substr($2,1, length($2)-4); print}' OFS='\t' $freqdir/imputed_anatolianeo41_gp80.map > $freqdir/imputed_anatolianeo41_gp80.ed.map
#mv $freqdir/imputed_anatolianeo41_gp80.map $freqdir/imputed_anatolianeo41_gp80.old.map
#mv $freqdir/imputed_anatolianeo41_gp80.ed.map $freqdir/imputed_anatolianeo41_gp80.map

# 	remove the ped set

poslist=$maindir/freqs/kgpm5ho_snpchrpos.edit.txt

#$plink --file imputed_anatolianeo41_gp80 --extract $poslist --recode --out $freqdir/imputed_anatolianeo41_gp80_m5_extract
#$plink --file $freqdir/imputed_anatolianeo41_gp80_m5_extract --freq --out $freqdir/anatolianeo41gp80_imputedfrequencies_m5

poslist=$maindir/freqs/kgpm1_5ho_snpchrpos.edit.txt
#$plink --file imputed_anatolianeo41_gp80 --extract $poslist --recode --out $freqdir/imputed_anatolianeo41_gp80_m1_5_extract
#$plink --file $freqdir/imputed_anatolianeo41_gp80_m1_5_extract --freq --out $freqdir/anatolianeo41gp80_imputedfrequencies_m1_5

	# B - gp 99

file=$maindir/freqs/anatolianeo41.v2.imputed.gp99.hov2.vcf.gz
filebase="anatolianeo41.v2.imputed.gp99.hov2"

#$plink --vcf $file --recode --out $freqdir/$filebase
#awk '{$2 = substr($2,1, length($2)-4); print}' OFS='\t' $freqdir/$filebase.map > $freqdir/$filebase.ed.map
#mv $freqdir/$filebase.map $freqdir/$filebase.old.map
#mv $freqdir/$filebase.ed.map $freqdir/$filebase.map

poslist=$maindir/freqs/kgpm5ho_snpchrpos.edit.txt
$plink --file $freqdir/$filebase --extract $poslist --recode --out $freqdir/${filebase}_m5_extract
$plink --file $freqdir/${filebase}_m5_extract --freq --out $freqdir/${filebase}_imputedfrequencies_m5

poslist=$maindir/freqs/kgpm1_5ho_snpchrpos.edit.txt
$plink --file $freqdir/$filebase --extract $poslist --recode --out $freqdir/${filebase}_m1_5_extract
$plink --file $freqdir/${filebase}_m1_5_extract --freq --out $freqdir/${filebase}_imputedfrequencies_m1_5



# 2 - angsd freq of haploid called data

angsd=/usr/local/sw/angsd/angsd
bamlistananeo=$freqdir/anatolianeo41bamlist.txt

poslist=$maindir/freqs/kgpm5ho_snpchrpos.ed.txt
#$angsd -out $freqdir/anatolianeo41_angsdfreq_m5 -rf $poslist -doMajorMinor 1 -doMaf 3 -bam $bamlistananeo -GL 2

poslist=$maindir/freqs/kgpm1_5ho_snpchrpos.ed.txt
#$angsd -out $freqdir/anatolianeo41_angsdfreq_m1_5 -rf $poslist -doMajorMinor 1 -doMaf 3 -bam $bamlistananeo -GL 2


