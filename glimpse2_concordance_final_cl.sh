#!/bin/bash

# Genotype Concordance on chromosome 22 with GLIMPSE2_concordance

#mkdir -p GLIMPSE_concordance
GLIMPSE_concordance=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_concordance
glimpsedir=$1
maindir=$glimpsedir/highcov_samples
chr="chr22"

#name=$1
#maindir=$glimpsedir/highcov_samples/$name

# A - KGP Panel 
#name="highten"
#arrdepth=("0.1x" "0.25x" "0.5x" "1x") 
#for i in "${arrdepth[@]}"
#do
#depth="$i"
#1 Panel: KGP, FreqFile: kgp.isec, AF-TAG: AF_NFE, filename: highten; GP: No and yes
#tag="nfe"
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.isec.bcf ${maindir}/GLIMPSE_validation/highcovten.$chr.validation.bcf ${maindir}/GLIMPSE_ligate/imputedallten.$chr.$depth.bcf" > $infile

# 2  Panel: KGP, FreqFile: kgp,hgdptgp isec, AFTAG: AF_NFE, filename: highten, GP:no and yes
#tag="isec1nfe"
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.$tag.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.hgdptgp.isec.vcf.gz ${maindir}/GLIMPSE_validation/highcovten.$chr.validation.bcf ${maindir}/GLIMPSE_ligate/imputedallten.$chr.$depth.bcf" > $infile

#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile


# 3 with GP
#gpnumber=("80" "90" "95" "99")
#	for j in "${gpnumber[@]}"
#	do
#	number="$j"
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.kgp.$tag.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.kgp.$tag.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.kgp.$tag.conc.log
# 1 with gp
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.isec.bcf ${maindir}/GLIMPSE_validation/highcovten.$chr.validation.bcf ${maindir}/GLIMPSE_ligate/imputedallten.$chr.$depth.gp$number.vcf.gz" > $infile
# 2 with GP
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.hgdptgp.isec.vcf.gz ${maindir}/GLIMPSE_validation/highcovten.$chr.validation.bcf ${maindir}/GLIMPSE_ligate/imputedallten.$chr.$depth.gp$number.vcf.gz" > $infile

#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile

#	done
#done


# B - HGDPTGP Panel
# Panel: HGDPTGP , freqfile: kgp-hgdptgp isec (isec1) , AF_TAG: AF_NFE, filename: highten.hgdptgp
name="highten.hgdptgp"
imputefilebase=$glimpsedir/HGDP_TGP/imputeddataset/GLIMPSE_ligate/highcov_ten.v2.hgdptgp
validfile=$maindir/GLIMPSE_validation/validation.allten.hgdptgp.chr22.nm.bcf
arrdepth=("0.1x" "0.25x" "0.5x" "1x") 
#for i in "${arrdepth[@]}"
#do
depth="$i"
 
tag="isec1nfe"
infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.conc.infile
outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.conc
logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.conc.log
# 1 - no gp
#highcov_ten.v2.hgdptgp.chr22.1x.vcf.gz
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.hgdptgp.isec.vcf.gz	$validfile  $imputefilebase.$chr.$depth.nm.vcf.gz " > $infile

#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile


# 2 - with GP
#gpnumber=("80" "90" "95" "99")
#	for j in "${gpnumber[@]}"
#	do
#	number="$j"
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.$tag.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.$tag.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$number.$tag.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.$chr.tmp.kgp.hgdptgp.isec.vcf.gz	$validfile  $imputefilebase.$chr.$depth.gp$number.nm.vcf.gz " > $infile

#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile

	#done
#done

