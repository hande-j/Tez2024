#!/bin/bash -l

# Downloading the KGP Phase 3 Reference Panel

# step 1
glimpsedir=$1
# list of vcf.gz files
refpanelfile=$2
#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${refpanelfile} $glimpsedir/reference_panel/

# step 2
while read i j;
do
filename=$i
MD5_code=$j

echo -e "Checking MD5 sums of downloaded files..."
vcf_md5_to_test=${MD5_code}
vcf_md5_from_file=$(md5sum ${filename} | cut -d " " -f1)

if [[ $vcf_md5_to_test == $vcf_md5_from_file ]]; then
    #    echo -e "\n\e[92mSUCCESS\e[39m\nMD5 sum of downloaded vcf.gz file matches the expected hash!"
else
    #    echo -e "\n\e[91mFAILURE\e[39m\nMD5 sum of downloaded vcf.gz file does not match the expected hash. Please try downloading the file again!"
    #    exit 1
fi
done < ../file_md5_list > md5checkresult2.txt

# step 3
bcftools=/usr/local/sw/bcftools-1.9/bcftools
filename=$1
filebase=$(basename ${filename} .vcf.gz)

#${bcftools} norm -m -any ${filename} -Ou --threads 4 | ${bcftools} view -m 2 -M 2 -v snps --threads 4 -Ob -o ${filebase}.bcf
#${bcftools} index -f ${filebase}.bcf --threads 4

