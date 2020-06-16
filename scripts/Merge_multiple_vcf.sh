#!/bin/bash

################## Merge more than 2 complete (LRR and BAF annotated) vcf files ###############################


if [ -v $1 ]
then
echo 'to run this script you need vcftools and bcftools installed globally' 
echo 'this script should be run like this:'
echo './Merge_vcf.sh <path folder> <path folder 1> <path folder 2>'
else

path=$1
path_multiple=$2

# comm -12 $path_p1'SNPs_QC.txt' $path_p2'SNPs_QC.txt' > $path'comm_SNPs.txt'
readarray -t path_names < ${path_multiple}
cp ${path_names[0]}SNPs_QC.txt ${path}tot_SNPs.txt 

for i in ${path_names[@]:1}
do
echo ${i}
awk 'NR==FNR{seen[$0]=1; next} seen[$0]' ${path}tot_SNPs.txt ${i}SNPs_QC.txt > ${path}tmp
mv ${path}tmp ${path}tot_SNPs.txt
done

cut -f3 ${path}tot_SNPs.txt > ${path}comm_SNPs.txt
n=${#path_names[@]}
file=()

for i in $(seq ${n})
do
echo ${i}
p=$(eval echo "\${path_names[${i}-1]}")
vcftools --vcf ${p}data_BAF_LRR.vcf --snps ${path}comm_SNPs.txt --recode --recode-INFO-all --out ${path}vcf${i}_filt
bgzip ${path}vcf${i}_filt.recode.vcf
tabix ${path}vcf${i}_filt.recode.vcf.gz
file+=(${path}vcf${i}_filt.recode.vcf.gz)
done

bcftools merge ${file[@]} > ${path}data_BAF_LRR.vcf

# check number of snps is the same as in comm_SNPs.txt
if [ $(grep -v "^#" $path'data_BAF_LRR.vcf' | wc -l) != $(wc -l < $path'comm_SNPs.txt') ]
then 
	echo "number of lines not correct, check SNPs intersection"
fi

fi

