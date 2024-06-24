#!/bin/bash -l

convertf=/usr/local/sw/EIG-7.2.1/bin/convertf
bcftools=/usr/local/sw/bcftools-1.9/bcftools
outdir=$1
## infile: chromosomes of KGP reference panel merged in a single file
infile=$outdir/1KGpanel_new.vcf.gz
outfile=$outdir/refgenome_humanoriginsposlist
datasetdir=$2
humanoriginspos=$datasetdir/HumanOrigins_v2_hg38_pos.bed

## get the chr, pos and REF genotype info of KGP panel at humanorigins positions in a file
## then get only the homozygous ref genotypes

#$bcftools view -H -T $humanoriginspos $infile | awk '{print $1,$2,$3,$4,$4}' OFS=" " > $outfile
#awk '{print $4, $4 }' OFS=" " < $outfile > $outfile.refalleles

## make .ped line of the alleles
#tr '\n' ' ' < $outfile.refalleles | sed 's/ $/\n/' > $outdir/refgenome_humanoriginspos.ped
#awk '{print "HG38","HG38","0","0","0","HRef38", $0 }' OFS=" " $outdir/refgenome_humanoriginspos.ped > $outdir/refgenome38_humanorigins.ped
## make .map line
#awk '{$1=substr($1,4,2); $3=substr($3,0,length($3)-3); print $1,$3, "0",$2}' OFS="\t" refgenome_humanoriginsposlist > refgenome38_humanorigins.map
## make .pedind line
#echo "HG38 HG38 0 0 0 HRef38" > $outdir/refgenome38_humanorigins.pedind
## convert the plink line to eig format
#$convertf -p par_convert_href38
