#!/bin/bash -l
refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
glimpsedir=$2
sample=$3
samplename="$(awk '{print substr($0, 45)}' <<< ${sample} | awk -F [_.-] '{print $1}')" 
depth=$4

mkdir -p $glimpsedir/highcov_samples/${samplename}
outdir=$glimpsedir/highcov_samples/${samplename}

new="0.1"
ratio=$(awk "BEGIN {print $new/$depth}")
echo "ratio: ${ratio}"

samtools view -T ${ref} -s ${ratio} -bo ${outdir}/${samplename}.hg38.${new}x.bam ${sample}
samtools index ${outdir}/${samplename}.hg38.${new}x.bam

new="0.25"
ratio=$(awk "BEGIN {print $new/$depth}")
echo "ratio: ${ratio}"

samtools view -T ${ref} -s ${ratio} -bo ${outdir}/${samplename}.hg38.${new}x.bam ${sample}
samtools index ${outdir}/${samplename}.hg38.${new}x.bam

new="0.5"
ratio=$(awk "BEGIN {print $new/$depth}")
echo "ratio: ${ratio}"

samtools view -T ${ref} -s ${ratio} -bo ${outdir}/${samplename}.hg38.${new}x.bam ${sample}
samtools index ${outdir}/${samplename}.hg38.${new}x.bam

new="1"
ratio=$(awk "BEGIN {print $new/$depth}")
echo "ratio: ${ratio}"

samtools view -T ${ref} -s ${ratio} -bo ${outdir}/${samplename}.hg38.${new}x.bam ${sample}
samtools index ${outdir}/${samplename}.hg38.${new}x.bam
