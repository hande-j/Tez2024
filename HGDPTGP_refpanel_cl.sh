#!/bin/bash -l

# Using KGP+HGDP merged data as Reference Panel from Phased Haplotypes - chr22

# download the data 
#wget -c https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr22.merged.bcf

filename=hgdp.tgp.gwaspy.merged.chr22.merged.bcf
#filebase=$(basename $filename .bcf)
#chrlist="chr22"
#${bcftools} norm -m -any ${filename} -Ou --threads 4 | ${bcftools} view -m 2 -M 2 -v snps --threads 4 -Ob -o HGDP.TGP.merged.haplo.chr22.bcf
#${bcftools} index -f HGDP.TGP.merged.haplo.chr22.bcf --threads 4

#${bcftools} view -G HGDP.TGP.merged.haplo.${chrlist}.bcf -Oz -o HGDP.TGP.merged.haplo.${chrlist}.sites.vcf.gz
#${bcftools} index -f HGDP.TGP.merged.haplo.${chrlist}.sites.vcf.gz

#${bcftools} query -f'%CHROM\t%POS\t%REF,%ALT\n' HGDP.TGP.merged.haplo.${chrlist}.sites.vcf.gz | bgzip -c > HGDP.TGP.merged.haplo.${chrlist}.sites.tsv.gz
#${tabix} -s1 -b2 -e2 HGDP.TGP.merged.haplo.${chrlist}.sites.tsv.gz

# error stating AN is needed ->

# update vcf file header to include AN (AN=#samples*2 )
#head updatevcf.txt
#chr22   10510380        8198
#chr22   10510424        8198
#chr22   10510438        8198
#chr22   10510450        8198
#chr22   10510530        8198
#chr22   10510545        8198
#chr22   10510569        8198
#chr22   10510597        8198
#chr22   10513324        8198
#chr22   10514011        8198

#bgzip updatevcf.txt -c > updatevcf.tsv.gz
#tabix -s1 -b2 -e2 updatevcf.tsv.gz
#echo '##INFO=<ID=AN,Number=1, Type=Integer, Description="Allele Number">'> anheader.txt

#$bcftools annotate -h anheader.txt -a updatevcf.tsv.gz -c CHROM,POS,AN HGDP.TGP.merged.haplo.chr22.bcf -Ob -o HGDP.TGP.merged.haplo.chr22.updated.bcf
#$bcftools index HGDP.TGP.merged.haplo.chr22.updated.bcf

#${bcftools} view -G HGDP.TGP.merged.haplo.${chrlist}.updated.bcf -Oz -o HGDP.TGP.merged.haplo.${chrlist}.updated.sites.vcf.gz
#${bcftools} index -f HGDP.TGP.merged.haplo.${chrlist}.updated.sites.vcf.gz

#${bcftools} query -f'%CHROM\t%POS\t%REF,%ALT\n' HGDP.TGP.merged.haplo.${chrlist}.updated.sites.vcf.gz | bgzip -c > HGDP.TGP.merged.haplo.${chrlist}.updated.sites.tsv.gz
#${tabix} -s1 -b2 -e2 HGDP.TGP.merged.haplo.${chrlist}.updated.sites.tsv.gz
