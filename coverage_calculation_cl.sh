#!/bin/bash -l
mkdir -p coverages
refdir=$1
refgenome=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta
genomeCov=/usr/local/sw/bedtools2-2.25.0/bin/genomeCoverageBed 

# Calculate whole genome coverage with quality filter

while read i; do 
bamfile=$i
name="$(awk -F [/] '{print $NF}' <<< ${bamfile} )"
echo $name

if [ -f coverages/${name}.cov.txt ]
then
	echo "file already created"
else
	echo "file not present"

	samtools view -b -q 30 -F 4 ${bamfile} | ${genomeCov} -ibam - -g ${refgenome} > coverages/${name}.cov.txt
	coverage=$(grep genome coverages/${name}.cov.txt | awk '{NUM+=$2*$3; DEN+=$3} END {print NUM/DEN}')
	echo "${name}	${coverage}" > coverages/${name}.genomecoverage.txt

fi

done < bamlist

cat coverages/*.genomecoverage.txt > downsample_new_genomecoverages.txt
 