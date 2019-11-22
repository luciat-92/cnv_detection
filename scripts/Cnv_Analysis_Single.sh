#!/bin/bash

##################################### Run CNV detection for a single line #####################

if [ -v $1 ]
then
echo 'to run this script you need bcftools (and R ?) globally installed' 
echo 'this script should be run like this:'
echo './Cnv_Analysis_Single.sh <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered>'
else

path=$1
file=$2
sample_name=$3

cd $path

id=$(tail -n 1 $file)
echo $id

mkdir outdir_${sample_name}_$id

bcftools cnv -s $id  -t ^0,25,26 -o outdir_${sample_name}_$id data_BAF_LRR.vcf

fi
