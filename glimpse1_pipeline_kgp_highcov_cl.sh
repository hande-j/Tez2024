#!/bin/bash -l

#bcftools=/usr/local/sw/bcftools-1.9/bcftools
bcftools=/usr/local/sw/bcftools-1.18/bcftools
tabix=/usr/local/sw/bcftools-1.9/htslib-1.9/tabix
GLIMPSE_chunk=/usr/local/sw/GLIMPSE-1.1.1/tutorial_hg19/bin/GLIMPSE_chunk
GLIMPSE_phase=/usr/local/sw/GLIMPSE-1.1.1/phase/bin/GLIMPSE_phase
GLIMPSE_ligate=/usr/local/sw/GLIMPSE-1.1.1/ligate/bin/GLIMPSE_ligate
GLIMPSE_sample=/usr/local/sw/GLIMPSE-1.1.1/sample/bin/GLIMPSE_sample
GLIMPSE_concordance=/usr/local/sw/GLIMPSE-1.1.1/tutorial/bin/GLIMPSE_concordance

##  	 GLIMPSE1 Pipeline  - (Downsampled Highcov 10)
refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
glimpsedir=$2
reference_panel=$glimpsedir/reference_panel
chrlist=$3
datasetdir=$4
bamfileall=$5
bamfile="$(awk -F [/] '{print $NF}' <<< $5)"
#echo "$bamfile"
name="$(awk '{sub(/.hg38./, "-", $0); print}' <<< $bamfile | awk '{sub(/.bam/,"", $0); print}' )"
#echo $name
# done before
#echo "$name" > samplename/$name.txt

dirname="highcov_v1"
mkdir -p $glimpsedir/highcov_samples/$dirname/${name}
outdir=$glimpsedir/highcov_samples/$dirname/${name}
maindir=$glimpsedir/highcov_samples/$dirname
olddir=$glimpsedir/highcov_samples/$name
samplename=$glimpsedir/highcov_samples/samplename/$name.txt
#cat $samplename

chrlist="chr22"

mkdir -p ${outdir}/${name}_vcf
mkdir -p ${outdir}/GLIMPSE_impute
mkdir -p ${outdir}/GLIMPSE_ligate
##mkdir -p ${outdir}/GLIMPSE_sample
mkdir -p ${outdir}/GLIMPSE_concordance
#mkdir -p ${outdir}/GLIMPSE_validation

# Split the Genome into Chunks
#${GLIMPSE_chunk} --input ${reference_panel}/1000GP.${chrlist}.sites.vcf.gz --region ${chrlist} --window-size 2000000 --buffer-size 200000 --output ${reference_panel}/chunks.${chrlist}.txt

#while read i;
#do
#chrlist=$i

# Computing GLs for a single individual
BAM=${bamfileall}
VCF=${reference_panel}/1000GP.${chrlist}.sites.vcf.gz
TSV=${reference_panel}/1000GP.${chrlist}.sites.tsv.gz
REFGEN=${ref}
OUT=${outdir}/${name}_vcf/${name}.hg38.${chrlist}.vcf.gz

#${bcftools} mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${chrlist} ${BAM} -Ou | ${bcftools} call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
#${bcftools} index -f ${OUT} 

# Impute and Phase a Whole Chromosome
#VCF=${outdir}/${name}_vcf/${name}.hg38.${chrlist}.vcf.gz
#REF=${reference_panel}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrlist}.filtered.shapeit2-duohmm-phased.bcf
#MAP=/usr/local/sw/GLIMPSE-1.1.1/maps/genetic_maps.b38/${chrlist}.b38.gmap.gz
#
#while IFS="" read -r LINE || [ -n "$LINE" ];
#do
#	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
#	IRG=$(echo $LINE | cut -d" " -f3)
#	ORG=$(echo $LINE | cut -d" " -f4)
#	OUT=${outdir}/GLIMPSE_impute/${name}.${chrlist}.imputed.${ID}.bcf
#
#	${GLIMPSE_phase} --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
#	${bcftools} index -f ${OUT}
#	
#done < ${reference_panel}/chunks.${chrlist}.txt
#
## Ligate multiple chunks together 
#LST=${outdir}/GLIMPSE_ligate/list.${chrlist}.txt
#ls ${outdir}/GLIMPSE_impute/${name}.${chrlist}.imputed.*.bcf > ${LST} 
#OUT=${outdir}/GLIMPSE_ligate/${name}.${chrlist}.merged.bcf
#
#${GLIMPSE_ligate} --input ${LST} --output ${OUT}
#${bcftools} index -f ${OUT} 

#done < $glimpsedir/chrlist


## Reheader - Concordance
chr="chr22"
sname=$name
name="$(awk '{sub(/-.*/,"",$0); print}' <<< $sname )"
echo $name
depth="$(awk '{sub(/.*-/,"",$0); print}' <<< $sname )"
echo $depth
olddir=$glimpsedir/highcov_samples/$name

samplename=$glimpsedir/highcov_samples/samplename/$name.txt
#FILE=${outdir}/GLIMPSE_ligate/${name}-$depth.${chr}.merged.bcf
#$bcftools reheader -s $samplename $FILE -o ${outdir}/GLIMPSE_ligate/${name}.$depth.${chr}.merged.hdr.bcf
#$bcftools index ${outdir}/GLIMPSE_ligate/${name}.$depth.${chr}.merged.hdr.bcf


# Merge all samples - chromosome22
# 	with depth in the names
#mkdir -p $maindir/imputed_vcf
#mergelist=${maindir}/all.bcflist.$chrlist
#ls $maindir/*/GLIMPSE_ligate/*-*x.${chrlist}.merged.hdr.bcf > $mergelist
#${bcftools} merge -l ${mergelist} -Oz9 -o  ${maindir}/imputed_vcf/highten.v1.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/highten.v1.$chrlist.vcf.gz

#ls $maindir/*/GLIMPSE_ligate/*.0.1x.${chrlist}.merged.hdr.bcf > $maindir/bcflist.v1.0.1x
#${bcftools} merge -l $maindir/bcflist.v1.0.1x -Oz9 -o  ${maindir}/imputed_vcf/highten.v1.0.1x.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/highten.v1.0.1x.$chrlist.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.0.1x.$chrlist.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.0.1x.$chrlist.gp80.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.0.1x.$chrlist.gp80.vcf.gz

#ls $maindir/*/GLIMPSE_ligate/*.0.25x.${chrlist}.merged.hdr.bcf > $maindir/bcflist.v1.0.25x
#${bcftools} merge -l $maindir/bcflist.v1.0.25x -Oz9 -o  ${maindir}/imputed_vcf/highten.v1.0.25x.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/highten.v1.0.25x.$chrlist.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.0.25x.$chrlist.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.0.25x.$chrlist.gp80.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.0.25x.$chrlist.gp80.vcf.gz

#ls $maindir/*/GLIMPSE_ligate/*.0.5x.${chrlist}.merged.hdr.bcf > $maindir/bcflist.v1.0.5x
#${bcftools} merge -l $maindir/bcflist.v1.0.5x -Oz9 -o  ${maindir}/imputed_vcf/highten.v1.0.5x.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/highten.v1.0.5x.$chrlist.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.0.5x.$chrlist.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.0.5x.$chrlist.gp80.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.0.5x.$chrlist.gp80.vcf.gz

# also adds 0.1x- manually edit the list
#ls $maindir/*/GLIMPSE_ligate/*.1x.${chrlist}.merged.hdr.bcf > $maindir/bcflist.v1.1x
#${bcftools} merge -l $maindir/bcflist.v1.1x -Oz9 -o  ${maindir}/imputed_vcf/highten.v1.1x.$chrlist.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/highten.v1.1x.$chrlist.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.1x.$chrlist.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.1x.$chrlist.gp80.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.1x.$chrlist.gp80.vcf.gz

##  GP filters and stats for comparison with v2
#arrdepth=("0.1x" "0.25x" "0.5x" "1x") 
#for i in "${arrdepth[@]}"
#do
#depth="$i"
##${bcftools} filter ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.vcf.gz -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp80.vcf.gz
##${bcftools} index  ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp80.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.vcf.gz -i'FORMAT/GP>=0.90' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp90.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp90.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.vcf.gz -i'FORMAT/GP>=0.95' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp95.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp95.vcf.gz
#${bcftools} filter ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.vcf.gz -i'FORMAT/GP>=0.99' -S. -Oz9 -o ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp99.vcf.gz
#${bcftools} index  ${maindir}/imputed_vcf/highten.v1.$depth.$chrlist.gp99.vcf.gz
#
##$bcftools stats -S $glimpsedir/highcov_samples/tensamplevalidnamelist ${maindir}/imputed_vcf/highten.v1.${depth}.$chrlist.gp80.vcf.gz > $maindir/imputed_vcf/v1.$depth.gp80.stats.txt
#$bcftools stats -S $glimpsedir/highcov_samples/tensamplevalidnamelist ${maindir}/imputed_vcf/highten.v1.${depth}.$chrlist.gp90.vcf.gz > $maindir/imputed_vcf/v1.$depth.gp90.stats.txt
#$bcftools stats -S $glimpsedir/highcov_samples/tensamplevalidnamelist ${maindir}/imputed_vcf/highten.v1.${depth}.$chrlist.gp95.vcf.gz > $maindir/imputed_vcf/v1.$depth.gp95.stats.txt
#$bcftools stats -S $glimpsedir/highcov_samples/tensamplevalidnamelist ${maindir}/imputed_vcf/highten.v1.${depth}.$chrlist.gp99.vcf.gz > $maindir/imputed_vcf/v1.$depth.gp99.stats.txt
#done


##      Concordance
GLIMPSE_concordance=/usr/local/sw/GLIMPSE-1.1.1/tutorial/bin/GLIMPSE_concordance

## all together
#chr="chr22"
#name="highten"
#depth="0.1x"
# no gp
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc.log
# gp
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp80.kgp.v1.isec1nfe.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp80.kgp.v1.isec1nfe.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp80.kgp.v1.isec1nfe.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz	$glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.gp80.vcf.gz " > $infile

#$GLIMPSE_concordance --input $infile --minDP 8 --output $outfile --minPROB 0.99 --bins 0.000 0.001 0.005 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag AF_nfe $logfile

## by samples
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.kgp.v1.isec1nfe.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz	$olddir/GLIMPSE_validation/$name.$chr.validation.hdr.bcf  $outdir/GLIMPSE_ligate/$name.$depth.$chr.merged.hdr.bcf " > $infile

#$GLIMPSE_concordance --input $infile --minDP 8 --output $outfile --minPROB 0.99 --bins 0.000 0.001 0.005 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag AF_nfe $logfile
 
## samples
#depth="1x"
#chr="chr22"
#name="highten"
#num="99"

# hgdp-kgp isec with glimpse 1 - 
# gp 90-95-99 files are identical
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.kgp.v1.isec1nfe.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.kgp.v1.isec1nfe.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.kgp.v1.isec1nfe.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz	$glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.gp${num}.vcf.gz " > $infile

## GLIMPSE1_conc
#$GLIMPSE_concordance --input $infile --minDP 8 --output $outfile --minPROB 0.99 --bins 0.000 0.001 0.005 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag AF_nfe $logfile

## 	KGP NFE with glimpse1_conc NO GP filter
#tag="kgpnfe"
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.concv1.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.concv1
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.concv1.log

#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf	$glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.vcf.gz " > $infile
#$GLIMPSE_concordance --input $infile --minDP 8 --output $outfile --minPROB 0.99 --bins 0.000 0.001 0.005 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag AF_nfe $logfile

chr="chr22"
name="highten"
depth="1x"
num="80"
tag="kgpnfe"
infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$num.$tag.v1.concv1.infile
outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$num.$tag.v1.concv1
logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp$num.$tag.v1.concv1.log

#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf	$glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.gp$num.vcf.gz" > $infile

#$GLIMPSE_concordance --input $infile --minDP 8 --output $outfile --minPROB 0.99 --bins 0.000 0.001 0.005 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag AF_nfe $logfile

#done



## Concordance with glimpse2_conc: (not used for plotting)

#GLIMPSE_concordance=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_concordance

# KGP sites, NFE tag
#tag="kgpnfe" # delete *kgp.v1.conc* files later
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.$tag.v1.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.$tag.v1.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.gp${num}.$tag.v1.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf $glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.gp${num}.vcf.gz " > $infile
# no gp
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf $glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.vcf.gz " > $infile

# glimpse2
#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile

# two - hgdptgp-kgp isec (isec1) , nfe tag
# isec1 sites nfe tag
#tag="isec1nfe"
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz $glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.gp${num}.vcf.gz " > $infile
# two- no gp
#infile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc.infile
#outfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc
#logfile=$maindir/GLIMPSE_concordance/$name.$depth.$chr.$tag.v1.conc.log
#echo "chr22 ${glimpsedir}/gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz $glimpsedir/highcov_samples/GLIMPSE_validation/highcovten.$chr.validation.bcf  $maindir/imputed_vcf/$name.v1.$depth.$chr.vcf.gz " > $infile

# glimpse 2 code
#$GLIMPSE_concordance --bins 0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 --min-val-dp 8 --min-val-gl 0.99 --af-tag AF_nfe --input $infile --output $outfile --log $logfile
#done
