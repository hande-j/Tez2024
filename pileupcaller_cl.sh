#!/bin/bash -l

pileupCaller=/usr/local/sw/sequenceTools-1.5.3.1/bin/pileupCaller
samtools=samtools
refdir=$1
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
datasetdir=$2
bedfile=$datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38_pos.bed
eigensnpfile=$datasetdir/HumanOrigins_v2_hg38/HumanOrigins_v2_hg38.snp

bamlist=$3
dataname=$4
datadir=$5

${samtools} mpileup -R -B -q30 -Q30 -l ${bedfile} -f ${ref} $(cat ${bamlist} | xargs) > $datadir/${dataname}.pileup.txt
sampleNames=$(grep -o '[^/]*$' ${bamlist} | cut -d "." -f1 | cut -d "_" -f1 | cut -d "-" -f1 | cut -f1,2 -d "_" | tr '\n' ',' | sed 's/,$//')
$pileupCaller --randomHaploid --sampleNames ${sampleNames} -f ${eigensnpfile} -e $datadir/${dataname}. < $datadir/${dataname}.pileup.txt
