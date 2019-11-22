#!/bin/bash

################## Merge two complete (LRR and BAF annotated) vcf files ###############################


if [ -v $1 ]
then
echo 'to run this script you need vcftools and bcftools installed globally' 
echo 'this script should be run like this:'
echo './Merge_vcf.sh <path folder> <path folder 1> <path folder 2>'
else

path=$1
path_p1=$2
path_p2=$3


# comm -12 $path_p1'SNPs_QC.txt' $path_p2'SNPs_QC.txt' > $path'comm_SNPs.txt'
awk 'NR==FNR{seen[$0]=1; next} seen[$0]' $path_p1'SNPs_QC.txt' $path_p2'SNPs_QC.txt' | cut -f3 > $path'comm_SNPs.txt'

vcftools --vcf $path_p1'data_BAF_LRR.vcf' --snps $path'comm_SNPs.txt' --recode --recode-INFO-all --out $path'vcf1_filt'
vcftools --vcf $path_p2'data_BAF_LRR.vcf' --snps $path'comm_SNPs.txt' --recode --recode-INFO-all --out $path'vcf2_filt'


bgzip $path'vcf1_filt.recode.vcf'
bgzip $path'vcf2_filt.recode.vcf'

tabix $path'vcf1_filt.recode.vcf.gz'
tabix $path'vcf2_filt.recode.vcf.gz'

bcftools merge $path'vcf1_filt.recode.vcf.gz' $path'vcf2_filt.recode.vcf.gz' > $path'data_BAF_LRR.vcf'

# check number of snps is the same as in comm_SNPs.txt
if [ $(grep -v "^#" $path'data_BAF_LRR.vcf' | wc -l) != $(wc -l < $path'comm_SNPs.txt') ]
then 
	echo "number of lines not correct, check SNPs intersection"
fi

fi

