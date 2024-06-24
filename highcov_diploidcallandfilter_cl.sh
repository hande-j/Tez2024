#!/bin/bash -l

#bcftools=/usr/local/sw/bcftools-1.11/bcftools
bcftools=/usr/local/sw/bcftools-1.9/bcftools
glimpsedir=$1
refdir=$2
ref=$refdir/hg38-v0-Homo_sapiens_assembly38.fasta 
reference_panel=$glimpsedir/reference_panel
datasetdir=$3
bamfileall=$4
name="$(awk '{print substr($0, 45)}' <<< ${bamfileall} | awk -F [_.-] '{print $1}')"
outdir=$glimpsedir/highcov_samples/${name}
maindir=$glimpsedir/highcov_samples
chrlist="chr22"
#mkdir -p samplename
#echo $name > samplename/$name.txt
samplename=$glimpsedir/highcov_samples/samplename/${name}.txt

##                   Diploid Genotyping

#while read i;
#do
#chrlist=$i
#BAM=${bamfileall}
#REFGEN=${ref}

## For KGP validaiton:
#VCF=${reference_panel}/1000GP.${chrlist}.sites.vcf.gz
#TSV=${reference_panel}/1000GP.${chrlist}.sites.tsv.gz
#OUT=${outdir}/GLIMPSE_validation/${name}.${chrlist}.validation.bcf

## for HGDPTGP validation:
#VCF=$glimpsedir/HGDP_TGP/HGDP.TGP.merged.haplo.chr22.sites.vcf.gz
#TSV=$glimpsedir/HGDP_TGP/HGDP.TGP.merged.haplo.chr22.sites.tsv.gz
#OUT=${outdir}/GLIMPSE_validation/${name}.${chrlist}.hgdptgp.validation.bcf

## Genotype calling with bcftools with quality filters 
#${bcftools} mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -q 30 -Q 20 -C 50 -T ${VCF} -r ${chrlist} ${BAM} -Ou | ${bcftools} call -Aim -C alleles -T ${TSV} -Ob -o ${OUT}
#${bcftools} index -f ${OUT}

## 	Reheader - files need sample names
#while read i;do
#chrlist=$i
#${bcftools} reheader -s ${samplename} ${outdir}/GLIMPSE_validation/${name}.${chrlist}.validation.bcf -o ${outdir}/GLIMPSE_validation/${name}.${chrlist}.validation.hdr.bcf
#${bcftools} index ${outdir}/GLIMPSE_validation/${name}.${chrlist}.validation.hdr.bcf
#done < ../chrlist

## Concat Each Sample
#list=$maindir/$name/$name.valid.bcflist
#ls -v $maindir/$name/GLIMPSE_validation/$name.*.hdr.bcf > $list 
#$bcftools concat -f $list -Ob -o $maindir/$name/GLIMPSE_validation/$name.all.validation.bcf
#$bcftools index $maindir/$name/GLIMPSE_validation/$name.all.validation.bcf



##        Filtering the Diploid Called Genotypes (Creating Validation Data) of High Coverage Samples

# 1&2:  DP and QUAL Filters:
# 1 DP Diltering min:8, max: 2*depth
# file:tenbamcovlist with columns sample,depth, minDP, maxDP:
#AKT16   8.58423 8       17
#BAR25   8.86314 8       17
#Nea2    10.8218 8       21
#Nea3    9.75588 8       19
#WC1     8.72987 8       17
#sf12    53.4582 8       106
#atp016  12.3625 8       24
#KK1     9.41379 8       19
#LBK     13.6523 8       27
#BOT2016 12.358  8       24

#name=$1
#depth=$2
#low=$3
#high=$4

#infile=$maindir/$name/GLIMPSE_validation/$name.all.validation.bcf
#outfile=$maindir/$name/GLIMPSE_validation/$name.all.validation.dpqual.bcf
#
#$bcftools view -i 'INFO/DP >= '$low' & INFO/DP <= '$high' & QUAL>30' -Ob -o $outfile $infile
#$bcftools index $outfile

# 3&4: get genome mask positions and remove repeat regions

# source for genome mask (GRCh38): https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed
# repeatmask source: UCSC Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1891779418_amgqN9j241alcLRk4CZMaNRtub4A
# all positions in the RepeatMasker track selected, head of file:
#chr1    67108753        67109046        L1P5    1892    +
#chr1    8388315 8388618 AluY    2582    -
#chr1    25165803        25166380        L1MB5   4085    +
#chr1    33554185        33554483        AluSc   2285    -
#chr1    41942894        41943205        AluY    2451    -
#chr1    50331336        50332274        HAL1    1587    +

secondir=$glimpsedir/genome_mask
genomemask=${secondir}/20160622.allChr.mask.bed
repeatmask=${secondir}/repeatmask.gz

#   get passed positions from the strict genome mask
#infile=$maindir/$name/GLIMPSE_validation/$name.all.validation.dpqual.bcf
#outfile=$maindir/$name/GLIMPSE_validation/$name.all.validation.dpqualmask38.bcf
#${bcftools} view -R ${genomemask} -Ob -o $outfile $infile
#${bcftools} index $outfile
#   remove repeatmask positions
#infile=$maindir/$name/GLIMPSE_validation/$name.all.validation.dpqualmask38.bcf
#outfile=$maindir/$name/GLIMPSE_validation/$name.all.validation.dpqualmask38norep.bcf
#${bcftools} view -T ^${repeatmask} -Ob -o $outfile $infile
#${bcftools} index $outfile

# merge the samples
#ls $maindir/*/GLIMPSE_validation/*.all.validation.dpqualmask38norep.bcf > $maindir/validationfilter2.bcflist
#$bcftools merge -l $maindir/validationfilter2.bcflist -Oz9 -o imputed_vcf/validation.allten.filtered2.vcf.gz
#$bcftools index imputed_vcf/validation.allten.filtered2.vcf.gz

# HumanOrigins_v2_hg38 Subset 
humanorigins_v2=$datasetdir/HumanOrigins_v2_hg38_pos.bed
#${bcftools} view ${maindir}/imputed_vcf/validation.allten.filtered2.vcf.gz -R ${humanorigins_v2} -Oz9 -o ${maindir}/imputed_vcf/validation.allten.filtered2.hov2.vcf.gz
#${bcftools} index ${maindir}/imputed_vcf/validation.allten.filtered2.hov2.vcf.gz
infile=${maindir}/imputed_vcf/validation.allten.filtered2.hov2.vcf.gz
outfile=${maindir}/imputed_vcf/validation.allten.filtered2.hov2.id.vcf.gz
#${bcftools} annotate --rename-chrs ../fiftysix/chromnames8to7.txt $infile -Ou | ${bcftools} annotate --set-id '%CHROM\:%POS' -Ou | ${bcftools} annotate --rename-chrs ../fiftysix/chromnames7to8.txt -Oz9 -o $outfile
#${bcftools} index $outfile
