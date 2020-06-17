#!/bin/bash

##################################### Match sample base on genotype and create new annotation file for the merging #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: argparse, circlize, ComplexHeatmap, RColorBrewer, doParallel' 
echo 'this script should be run like this:'
echo './GT_match_merge.sh <path to Rscripts folder> <path output folder> <path project filelist> <path annotation filelist> <File with samples to consider> <ncores>'
echo 'order in <path project filelist> <path annotation filelist> must match'
else

path_Rscript=$1
path=$2
path_file=$3
ann_file=$4
samples_file=$5
ncores=$6

samples_names=$(cat $samples_file) 
readarray -t path_names < ${path_file}
readarray -t ann_names < ${ann_file}
n=${#path_names[@]}

GT_file=()
ann_file_tot=()

for i in $(seq ${n})
do
p=$(eval echo "\${path_names[${i}-1]}")
a=$(eval echo "\${ann_names[${i}-1]}")
GT_file+=(${p}GT_table.txt)
ann_file_tot+=(${p}${a})
done

Rscript $path_Rscript'GT_match_mergedvcf_run.R'  --fold_fun $path_Rscript  --ann_file ${ann_file_tot[@]}  --GT_file ${GT_file[@]} --samples_name $samples_names  --comm_SNPs_file $path'comm_SNPs.txt'  --ncores $ncores --outf  $path'Merged_'

fi
