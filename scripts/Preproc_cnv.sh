#!/bin/bash

################## Filter SNPs and annotate vcf with LRR and BAF ###############################


if [ -v $1 ]
then
echo 'to run this script you need vcftools and bcftools installed globally' 
echo 'this script should be run like this:'
echo './Preproc_CNV.sh <path folder>' 
else


path=$1

cd $path

echo $path

# filter SNPs
cut  -f2 SNP_Table_filt.txt > temp.txt
tail -n +2 temp.txt > SNPs_QC.txt
rm temp.txt


cp plink.vcf_RefAltCorrected.vcf plink_RefAltCorrected.vcf

vcftools --vcf plink_RefAltCorrected.vcf --snps SNPs_QC.txt --recode --recode-INFO-all --out vcf_file_QC


LRR=LRR_table.txt # na saved as '.', position ordered
BAF=BAF_table.txt # na saved as '.', position ordered


#### annotate vcf file with BAF and LRR
# zip
bgzip $LRR
bgzip $BAF


# annotate the vcf file
tabix -s1 -S1 -b2 -e2 $LRR'.gz'
tabix -s1 -S1 -b2 -e2 $BAF'.gz'
echo '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="GS estimate of BAF">' > baf.hdr
echo '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="GS estimate of LRR">' > lrr.hdr
bcftools annotate -a BAF_table.txt.gz -h baf.hdr -c CHROM,POS,ID,FMT/BAF -Ov -o output.vcf vcf_file_QC.recode.vcf
bcftools annotate -a LRR_table.txt.gz -h lrr.hdr -c CHROM,POS,ID,FMT/LRR -Ov -o data_BAF_LRR.vcf output.vcf
rm output.vcf

# note:	duplicates arise when variants have the	same position in plink.vcf --> bcftools	use only first IDs
# delete duplicates (all, do not keep a copy)
mv data_BAF_LRR.vcf data_BAF_LRR_withdup.vcf
grep -v "^#" data_BAF_LRR_withdup.vcf | awk {'print $3'} > tmp
uniq -d tmp > dupl.list
vcftools --vcf data_BAF_LRR_withdup.vcf --exclude dupl.list  --recode-INFO-all --recode --out data_BAF_LRR
mv data_BAF_LRR.recode.vcf data_BAF_LRR.vcf
rm tmp

# update SNP_QC list 
mv SNPs_QC.txt SNPs_QC_withdup.txt 
grep -v "^#" data_BAF_LRR.vcf | cut -f1,2,3 > SNPs_QC.txt

fi

