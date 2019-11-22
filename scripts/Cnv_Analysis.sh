#!/bin/bash

##################################### Detect CNV in the pairwise comparison #####################

if [ -v $1 ]
then
echo 'to run this script you need bcftools globally installed' 
echo 'this script should be run like this:'
echo './Cnv_Analysis.sh <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered>'
else


path=$1
file=$2
sample_name=$3

cd $path

nrow=$(wc -l < $file)
thr=$(($nrow - 2))

pre=$(head -n 2 $file | tail -n 1)
post=$(tail -n $thr $file)

mkdir outdir_${sample_name}_${pre}
 
for post_sample in ${post} 
do
	mkdir outdir_${sample_name}_${pre}/output_${post_sample}
	
	bcftools cnv -c ${pre}  -s ${post_sample} -t ^0,25,26 -o outdir_${sample_name}_${pre}/output_${post_sample} data_BAF_LRR.vcf

done

fi

