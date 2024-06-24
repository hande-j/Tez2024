#!/bin/bash -l

bcftools=/usr/local/sw/bcftools-1.9/bcftools
tabix=/usr/local/sw/bcftools-1.9/htslib-1.9/tabix
GLIMPSE_chunk=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_chunk
GLIMPSE_phase=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_phase
GLIMPSE_ligate=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_ligate
GLIMPSE_split_reference=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_split_reference
GLIMPSE_concordance=/usr/local/sw/GLIMPSE/tutorial/bin/GLIMPSE2_concordance

##      	GLIMPSE v2 Pipeline  -  KGP  -  Downsampled 10 HighCov Samples
refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
glimpsedir=$2
reference_panel=$glimpsedir/reference_panel
chrlist=$3
datasetdir=$4

#bamfileall=$5
#name="$(awk -F [/] '{print $1}' <<< ${bamfileall} )"

arrdepth=( "0.1x" "0.25x" "0.5x" "1x" )
#while read j
#do 
#bamfileall=$j
#name="$(awk '{print substr($0, 45)}' <<< ${bamfileall} | awk -F [_.-] '{print $1}')"
#for i in "${arrdepth[@]}"
#do
#depth="$i"	
#echo "$name-$depth" > samplename/$name-$depth.txt
#done < tenbamlist

#depth="$(awk '{print substr($0, length($0)-7, 4 ) }'  <<< ${bamfileall} )"
##echo "$depth"

#mkdir -p $glimpsedir/highcov_samples/${name}
#outdir=$glimpsedir/highcov_samples/${name}
#maindir=$glimpsedir/highcov_samples

#chrlist="chr22"

#mkdir -p ${outdir}/GLIMPSE_impute
#mkdir -p ${outdir}/GLIMPSE_ligate
#mkdir -p ${outdir}/GLIMPSE_concordance
#mkdir -p ${outdir}/GLIMPSE_validation

##  	Split the Genome into Chunks
#MAP=/usr/local/sw/GLIMPSE-1.1.1/maps/genetic_maps.b38/${chrlist}.b38.gmap.gz
#${GLIMPSE_chunk} --input ${reference_panel}/1000GP.${chrlist}.sites.vcf.gz --map ${MAP} --region ${chrlist} --sequential --output ${maindir}/chunks.${chrlist}.txt

## 	Make Reference Panel Binary GLIMPSE format
#REF=${reference_panel}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrlist}.filtered.shapeit2-duohmm-phased.bcf
#MAP=/usr/local/sw/GLIMPSE-1.1.1/maps/genetic_maps.b38/${chrlist}.b38.gmap.gz

#while IFS="" read -r LINE || [ -n "$LINE" ];
#do
#	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
#	IRG=$(echo $LINE | cut -d" " -f3)
#	ORG=$(echo $LINE | cut -d" " -f4)
#
#	${GLIMPSE_split_reference} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${reference_panel}/split/1000GP
#
#done <  ${maindir}/chunks.${chrlist}.txt

##  	Impute and Phase a Whole Chromosome
#while read i;
#do
#chrlist=$i

#BAM=${bamfileall}
#BAM=bamnamelist_0.1x.txt
#BAM=wc1dsbamsamplename.txt
#BAM=bamlist0.1x.txt
#depth='0.1x'
#REF=${reference_panel}/split/1000GP

#mkdir -p ${maindir}/GLIMPSE_impute/cov_$depth
#mkdir -p ${maindir}/GLIMPSE_impute/cov_${depth}_2
#mkdir -p ${maindir}/GLIMPSE_impute/WC1ds
#
#while IFS="" read -r LINE || [ -n "$LINE" ];
#do
#	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
#	IRG=$(echo $LINE | cut -d" " -f3)
#	ORG=$(echo $LINE | cut -d" " -f4)
#	CHR=$(echo ${LINE} | cut -d" " -f2 )
#	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1) 
#	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
#	
    #OUT=${outdir}/GLIMPSE_impute/${name}.${depth}.${CHR}.${REGS}.${REGE}.bcf
#	#OUT=${maindir}/GLIMPSE_impute/cov_$depth/imputed.${chrlist}.${depth}.${REGS}.${REGE}.bcf
#	#OUT=${maindir}/GLIMPSE_impute/cov_${depth}_2/imputed_2.${chrlist}.${depth}.${REGS}.${REGE}.bcf
#	#OUT=${maindir}/GLIMPSE_impute/WC1ds/WC1ds.imputed.${chrlist}.${REGS}.${REGE}.bcf
#
#	#${GLIMPSE_phase} --bam-file ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}
#	${GLIMPSE_phase} --bam-list ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}
#	
#done < ${reference_panel}/split/chunks.${chrlist}.txt

#mkdir -p ${maindir}/GLIMPSE_ligate/cov_${depth}_2
#mkdir -p ${maindir}/GLIMPSE_ligate/WC1ds

##  Ligate multiple chunks together 
#LST=${maindir}/GLIMPSE_ligate/cov_${depth}_2/list.${chrlist}.$depth.txt
#LST=${outdir}/GLIMPSE_ligate/list.${chrlist}.$depth.txt
#ls -1v ${outdir}/GLIMPSE_impute/${name}.${depth}.${chrlist}.*.bcf > ${LST} 
#ls -1v ${maindir}/GLIMPSE_impute/cov_${depth}_2/imputed_2.${chrlist}.$depth.*.bcf > ${LST} 
#OUT=${outdir}/GLIMPSE_ligate/${name}.${depth}.${chrlist}.ligated.bcf
#OUT=${maindir}/GLIMPSE_ligate/cov_${depth}_2/imputed_2.${chrlist}.${depth}.ligated.bcf

#LST=${maindir}/GLIMPSE_ligate/WC1ds/list.${chrlist}.txt
#ls -1v ${maindir}/GLIMPSE_impute/WC1ds/WC1ds.imputed.${chrlist}.*.bcf > ${LST} 
#OUT=${maindir}/GLIMPSE_ligate/WC1ds/WC1ds.imputed.${chrlist}.ligated.bcf
#
#${GLIMPSE_ligate} --input ${LST} --output ${OUT}

#done < $glimpsedir/chrlistrev


## Reheader for ligated bcf files
#while read i; do
#chrlist=$i;
#name="imputed_2"
#depth="1x"
#samplename=./newfive1
#${bcftools} reheader -s ${samplename} ${maindir}/GLIMPSE_ligate/cov_${depth}_2/${name}.${chrlist}.$depth.ligated.bcf -o ${maindir}/GLIMPSE_ligate/cov_${depth}_2/${name}.${chrlist}.$depth.ligated.hdr.bcf
#${bcftools} index  ${maindir}/GLIMPSE_ligate/cov_${depth}_2/${name}.${chrlist}.$depth.ligated.hdr.bcf
#done < ../chrlist

#while read i; do
#chrlist=$i
#ls -v ${maindir}/GLIMPSE_ligate/*_2/*.$chrlist.*ligated.hdr.bcf > ${maindir}/bcflist2.$chrlist
#done < ../chrlist

#chrlist=$1
#$bcftools merge -l $maindir/bcflist2.$chrlist -Ob -o $maindir/GLIMPSE_ligate/downsampled_newfive.$chrlist.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/downsampled_newfive.$chrlist.ligated.bcf

#chrlist=$1
#$bcftools merge $maindir/GLIMPSE_ligate/downsampled_five.$chrlist.ligated.bcf $maindir/GLIMPSE_ligate/downsampled_newfive.$chrlist.ligated.bcf -Ob -o $maindir/GLIMPSE_ligate/downsampled_ten.$chrlist.ligated.bcf
#$bcftools index $maindir/GLIMPSE_ligate/downsampled_ten.$chrlist.ligated.bcf

samplename=samplename/AKT16.txt
#$bcftools reheader -s $samplename  AKT16/GLIMPSE_ligate/AKT16.0.1x.chr22.ligated.bcf -o AKT16/GLIMPSE_ligate/AKT16.0.1x.chr22.ligated.hdr.bcf
#$bcftools index AKT16/GLIMPSE_ligate/AKT16.0.1x.chr22.ligated.hdr.bcf

# concat chrs
#list=list103
#ls -v $maindir/GLIMPSE_ligate/downsampled_ten.*.ligated.bcf > $list
#$bcftools concat -f $list -Oz9 -o $maindir/imputed_vcf/downsampledten_imputed.vcf.gz
#$bcftools index $maindir/imputed_vcf/downsampledten_imputed.vcf.gz
