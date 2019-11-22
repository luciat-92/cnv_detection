#!/bin/bash

##################################### Match sample base on genotype and create new annotation file for the merging #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: argparse, circlize, ComplexHeatmap, RColorBrewer, doParallel' 
echo 'this script should be run like this:'
echo './GT_match_merge.sh <path to Rscripts folder> <path output folder> <path project 1> <path project 2> <path annotation file 1> <path annotation file 2> <File with samples to consider> <ncores>'
else

path_Rscript=$1
path=$2
path_p1=$3
path_p2=$4
ann_file1=$5
ann_file2=$6
samples_file=$7
ncores=$8

samples_names=$(cat $samples_file) 

Rscript $path_Rscript'GT_match_mergedvcf_run.R'  --fold_fun $path_Rscript  --ann_file $ann_file1 $ann_file2  --GT_file $path_p1'GT_table.txt' $path_p2'GT_table.txt' --samples_name $samples_names  --comm_SNPs_file $path'comm_SNPs.txt'  --ncores $ncores --outf  $path'Merged_'

fi
