#!/bin/bash -l

liftOver=/usr/local/sw/liftOver
bcftools=/usr/local/sw/bcftools-1.9/bcftools
plink=/usr/local/sw/plink-1.90/plink

refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta
inputdir=$2

#posbed=HumanOrigins.autosomal.bed
#file=$inputdir/HumanOrigins.autosomal
#awk '{$1 ="chr"$1; print $0}' $posbed > HOautosomal.wchr.bed

# lift the positions from hg37 to hg38
#${liftOver} HOautosomal.wchr.bed hg19ToHg38.over.chain.gz hoautosom.hg38lifted.bed hoautosom.hg38unlifted.bed

# unlifted pos edit - generating plink id
#awk '!/^#/' hoautosom.hg38unlifted.bed > hoautosom.unliftedpos.txt
#awk '{print $1=sprintf(substr($1, 4, length($1)))"\t"$3}' < hoautosom.unliftedpos.txt  > hoautosom.unliftedpos2.txt
#awk '{print $1":"$2}' hoautosom.unliftedpos2.txt > hoautosom.unliftedposids.txt

# creating new .map file 
#awk '{$2=sprintf( substr($1,4)":"$3); print $1,$2,'0',$3} ' OFS='\t' hoautosom.hg38lifted.bed > HOautosom_38.map

# copy the og dataset
#$plink --file $file --recode --out HO37autosom.tmp

## edit the map file id column
#mv HO37autosom.tmp.map HO37autosom.map
#awk '{$2 = sprintf($1":"$4) ; print}' OFS='\t' < HO37autosom.map > HO37autosom.tmp.map

## exclude unlifted positions
#$plink --file HO37autosom.tmp --exclude hoautosom.unliftedposids.txt --recode --out HO38autosom.tmp

## replace the 37-coord-map file with lifted-positions-map-file
#mv HO38autosom.tmp.map HO37autosom2.map
#cp HOautosom_38.map HO38autosom.tmp.map
#$plink --file HO38autosom.tmp --allow-extra-chr --chr chr1-chr22 --make-bed --recode --out HO38autosomal

# get bed file of positions
#awk '{$1=sprintf("chr"$1); print $1, $4-1, $4}' OFS='\t' HO38autosomal.map > HO38autosomalpos.bed
