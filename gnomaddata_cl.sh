#!/bin/bash -l

# create allele frequency files at positions selected to use in calculation of genotype concordance

## reference_panel files - chr22
KGPchr22bcf=./reference_panel/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.bcf
HGDPTGPchr22bcf=HGDP_TGP/HGDP.TGP.merged.haplo.chr22.updated.bcf
bcftools=/usr/local/sw/bcftools-1.9/bcftools
#bcftools=/usr/local/sw/bcftools-1.11/bcftools

## main freq file: 
## gnomad/gnomad.genomes.r3.0.sites.chr22.vcf.gz

## get snp and freq 3.0
#$bcftools annotate -i 'TYPE="snp" && INFO/AF_nfe!="." ' -x^INFO/AC,^INFO/AN,^INFO/AF,^INFO/AC_nfe,^INFO/AN_nfe,^INFO/AF_nfe,^INFO/AN_nfe -Ob -o gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.bcf gnomad/gnomad.genomes.r3.0.sites.chr22.vcf.gz
#$bcftools index gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.bcf

## isec with kgp 
#$bcftools isec -n=2 -w1 -Ob -o gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf  gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.bcf $KGPchr22bcf
#$bcftools index gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf

# intersection of gnomad3.0 with both kgp and hgdptgp (isec1)
#$bcftools isec -n=2 -w1 gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.isec.bcf $HGDPTGPchr22bcf -Oz9 -o gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz
#$bcftools index gnomad/gnomad.genomes.r3.0.sites.chr22.tmp.kgp.hgdptgp.isec.vcf.gz
