#!/bin/bash -l

#bcftools=/usr/local/sw/bcftools-1.11/bcftools
#bcftools=/usr/local/sw/bcftools-1.9/bcftools
bcftools=/usr/local/sw/bcftools-1.18/bcftools
tabix=/usr/local/sw/bcftools-1.9/htslib-1.9/tabix
GLIMPSE_chunk=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_chunk
GLIMPSE_phase=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_phase
GLIMPSE_ligate=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_ligate
GLIMPSE_split_reference=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_split_reference
GLIMPSE_concordance=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_concordance

#	GLIMPSE v2 Pipeline  - HGDP_TGP chr 22

refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
glimpsedir=$2
bamfileall=$3

panelname="HGDP_TGP"
maindir=$glimpsedir/${panelname}
reference_panel=$glimpsedir/${panelname}/reference_panel
mkdir -p ${reference_panel}
outdir=$maindir/imputeddataset
chrlist="chr22"

mkdir -p ${outdir}/GLIMPSE_impute
mkdir -p ${outdir}/GLIMPSE_ligate
#mkdir -p ${outdir}/GLIMPSE_concordance
mkdir -p ${outdir}/GLIMPSE_validation

#	GLIMPSE v2 HGDPTGP Reference Panel Prep
#file - HGDP.TGP.merged.haplo.chr22.updated.vcf.gz
prefix='HGDP.TGP.merged.haplo'

#  	Split the Genome into Chunks
MAP=/usr/local/sw/GLIMPSE-1.1.1/maps/genetic_maps.b38/${chrlist}.b38.gmap.gz
SITES=${maindir}/${prefix}.${chrlist}.updated.sites.vcf.gz
#${GLIMPSE_chunk} --input ${SITES} --map ${MAP} --region ${chrlist} --sequential --output ${reference_panel}/chunks.${chrlist}.txt

# 	Make Reference Panel Binary GLIMPSE format
REF=${maindir}/${prefix}.${chrlist}.updated.bcf
MAP=/usr/local/sw/GLIMPSE-1.1.1/maps/genetic_maps.b38/${chrlist}.b38.gmap.gz
mkdir -p ${reference_panel}/split/

while IFS="" read -r LINE || [ -n "$LINE" ];
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)

	#${GLIMPSE_split_reference} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${reference_panel}/split/${panelname}

done <  ${reference_panel}/chunks.${chrlist}.txt


##  	Impute and Phase a Whole Chromosome
#while read i;
#do
#chrlist=$i

#BAM=${bamfileall}
#BAM=$glimpsedir/highcov_samples/wc1dsbamsamplename.txt
#BAM=$glimpsedir/highcov_samples/highcov_ds_2.txt
#name="highcov_2"

REF=${reference_panel}/split/HGDP_TGP

#while IFS="" read -r LINE || [ -n "$LINE" ];
#do
#	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2 )
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1) 
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	
	#OUT=${outdir}/GLIMPSE_impute/${name}.${CHR}.${REGS}.${REGE}.bcf
	#OUT=${maindir}/GLIMPSE_impute/wc1ds.v2.hgdptgp.${CHR}.${REGS}.${REGE}.bcf
	#OUT=${outdir}/GLIMPSE_impute/$name.v2.hgdptgp.${CHR}.${REGS}.${REGE}.bcf


#	${GLIMPSE_phase} --bam-list ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}
	
#done < ${reference_panel}/chunks.${chrlist}.txt

##  Ligate multiple chunks together 
LST=${outdir}/GLIMPSE_ligate/$name.list.${chrlist}.txt
#ls -1v ${outdir}/GLIMPSE_impute/$name.v2.hgdptgp.${chrlist}.*.bcf > ${LST} 
OUT=${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.bcf
#${GLIMPSE_ligate} --input ${LST} --output ${OUT}
#done < $glimpsedir/chrlist


## 		Stats - in highcov_samples/
#filebase=(wc1ds.v2.hgdptgp.chr22.ligated.bcf .bcf)
#mkdir -p $maindir/stats
#$bcftools stats -S $maindir/wc1newsamplenames $maindir/GLIMPSE_ligate/wc1ds.v2.hgdptgp.chr22.ligated.bcf > $maindir/stats/wc1ds.v2.hgdptgp.chr22.stats.txt

#file=$outdir/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.nm.bcf 
#filebase=$outdir/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}
#filename=$name.v2.hgdptgp.$chrlist
file=$outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.${chrlist}.ligated.nm.bcf 
filebase=$outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.${chrlist}
filename=highcov_ten.v2.hgdptgp.$chrlist
#number='80'
#${bcftools} filter $file -i'FORMAT/GP>=0.80' -S. -Oz9 -o ${filebase}.gp$number.vcf.gz
#${bcftools} index  ${filebase}.gp$number.vcf.gz
#
#number='90'
#${bcftools} filter $file -i'FORMAT/GP>=0.90' -S. -Oz9 -o ${filebase}.gp$number.vcf.gz
#${bcftools} index  ${filebase}.gp$number.vcf.gz
#
#number='95'
#${bcftools} filter $file -i'FORMAT/GP>=0.95' -S. -Oz9 -o ${filebase}.gp$number.vcf.gz
#${bcftools} index  ${filebase}.gp$number.vcf.gz
#
#number='99'
#${bcftools} filter $file -i'FORMAT/GP>=0.99' -S. -Oz9 -o ${filebase}.gp$number.vcf.gz
#${bcftools} index  ${filebase}.gp$number.vcf.gz

### $bcftools query -l $file > $outdir/highcovds_samplelist

#gpnums=("80" "90" "95" "99")
#for i in "${gpnums[@]}"
#do
#number="$i"
##$bcftools stats -S $outdir/forty_namelist ${filebase}.gp${number}.vcf.gz > $outdir/stats/$filename.gp${number}.stats.txt
#done 

#samples=$glimpsedir/highcov_samples/twenty_dsnamelist
#$bcftools view -S $samples ${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.nm.bcf -Ob -o $outdir/GLIMPSE_ligate/highcov_1.v2.hgdptgp.chr22.ligated.bcf
#$bcftools index $outdir/GLIMPSE_ligate/highcov_1.v2.hgdptgp.chr22.ligated.bcf
#
#infile1=${outdir}/GLIMPSE_ligate/highcov_1.v2.hgdptgp.${chrlist}.ligated.bcf
#infile2=${outdir}/GLIMPSE_ligate/highcov_2.v2.hgdptgp.${chrlist}.ligated.bcf
#$bcftools merge $infile1 $infile2 -Ob -o $outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.chr22.ligated.bcf
#$bcftools index $outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.chr22.ligated.bcf

## 		Reheader - separete by depth - take by samplenames
#$bcftools query -l  $outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.chr22.ligated.bcf > $outdir/forty_dsnamelist
# gsub , vi edited file
#samplename=$outdir/forty_namelist
#name="highcov_ten"
#${bcftools} reheader -s ${samplename} ${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.bcf -o ${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.nm.bcf
#${bcftools} index ${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.nm.bcf

## 		Separate by Sample and each Depth
dir1=$glimpsedir/highcov_samples
#name="highcov_ten"
#infile=${outdir}/GLIMPSE_ligate/$name.v2.hgdptgp.${chrlist}.ligated.nm.bcf
#while read i;
#do 
#name="$i"
#name_list=$dir1/$name.samplelist
#outfile=${dir1}/$name/GLIMPSE_ligate/$name.all.v2.hgdptgp.${chrlist}.ligated.bcf
#$bcftools view -S $name_list -Oz9 -o $outfile $infile
#$bcftools index $outfile

#while read j;
#do sample="$j"
#	infile2=$outfile
#	outfile2=$dir1/$name/GLIMPSE_ligate/$sample.v2.hgdptgp.$chrlist.ligated.bcf
#	outfile3=$dir1/$name/GLIMPSE_ligate/$sample.v2.hgdptgp.$chrlist.ligated.nm.bcf
	#$bcftools view -S $dir1/samplename/$sample.txt -Oz9 -o $outfile2 $infile2
	#$bcftools index $outfile2
	#$bcftools reheader -s $dir1/samplename/$name.txt -o $outfile3 $outfile2
	#$bcftools index $outfile3
#done < $name_list 
#$dir1/$name.samplelist

#done < $dir1/tensamplevalidnamelist

file=$outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.${chrlist}.ligated.nm.bcf 
filebase=$outdir/GLIMPSE_ligate/highcov_ten.v2.hgdptgp.${chrlist}
filename=highcov_ten.v2.hgdptgp.$chrlist

##  Separate by Depth
#grep -- "25x" imputeddataset/forty_namelist > $outdir/highcov_0.25x_names
#grep -- "0.5x" imputeddataset/forty_namelist > $outdir/highcov_0.5x_names
#grep -- "0.1x" imputeddataset/forty_namelist > $outdir/highcov_0.1x_names
#grep -- "-1x" imputeddataset/forty_namelist > $outdir/highcov_1x_names
#mkdir -p $outdir/GLIMPSE_concordance
#gpnums=("80" "90" "95" "99")
#for i in "${gpnums[@]}"
#do
#number="$i"

depthlist=("0.1x" "0.25x" "0.5x" "1x")
for j in "${depthlist[@]}"
do 
depth="$j"
#infile=${filebase}.gp$number.vcf.gz
#outfile=${filebase}.$depth.gp$number.vcf.gz
#infile=$file
#outfile=$filebase.$depth.vcf.gz
#list=$outdir/highcov_${depth}_names
#$bcftools view -S $list -Oz9 -o $outfile $infile
#$bcftools index $outfile
#infile=${filebase}.$depth.gp$number.vcf.gz
#outfile=${filebase}.$depth.gp$number.nm.vcf.gz

infile=$filebase.$depth.vcf.gz
outfile=$filebase.$depth.nm.vcf.gz
$bcftools reheader -s $outdir/highcov_onlynames -o $outfile $infile
$bcftools index $outfile
done 
#done

