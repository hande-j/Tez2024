#!/bin/bash -l

bcftools=/usr/local/sw/bcftools-1.9/bcftools
tabix=/usr/local/sw/bcftools-1.9/htslib-1.9/tabix
outdir=$1

# Creating KGP_Panel MAF Filtered Sites VCF

# merge the chromosomes with a list of reference_panel vcf files
list=${outdir}/onekgvcflist_true
${bcftools} concat -f ${list} -Oz9 -o ${outdir}/imputed_vcf/1KGpanel_new.vcf.gz
${bcftools} index ${outdir}/imputed_vcf/1KGpanel_new.vcf.gz

# MAF>5 Positions
FILE=${outdir}/1KGpanel_new.vcf.gz
${bcftools} view ${FILE} -q 0.05:minor -Oz9 -o 1KGpanel_new.maf5.vcf.gz
${bcftools} index 1KGpanel_new.maf5.vcf.gz

${bcftools} query -f '%CHROM\t%POS\n' ${outdir}/1KGpanel_new.maf5.vcf.gz | bgzip -c > ${outdir}/1KGpanel.maf5.sites.new.tsv.gz
${tabix} -s1 -b2 -e2 ${outdir}/1KGpanel.maf5.sites.new.tsv.gz

# MAF>1 Positions
FILE=${outdir}/1KGpanel_new.vcf.gz
${bcftools} view ${FILE} -q 0.01:minor -Oz9 -o ${outdir}/1KGpanel_new.maf1.vcf.gz
${bcftools} index ${outdir}/1KGpanel_new.maf1.vcf.gz

${bcftools} query -f '%CHROM\t%POS\n' ${outdir}/1KGpanel_new.maf1.vcf.gz | bgzip -c > ${outdir}/1KGpanel.maf1.sites.tsv.gz
${tabix} -s1 -b2 -e2 ${outdir}/1KGpanel.maf1.sites.tsv.gz
