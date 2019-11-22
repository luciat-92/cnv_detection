#!/bin/bash

##################################### Plot LRR and BAF for single analysis #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: ggplot2, argparse, grid, gridExtra' 
echo 'this script should be run like this:'
echo './Cnv_Plot_LRR-BAF_single.sh <path Rscript folder> <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered> <path to annotation file> <chromosomes to be plotted>'
else

if [ -v $6 ]
then
chr=0
else
chr=$6
fi

# echo $chr

path_Rscript=$1
path=$2
file=$3
sample_name=$4
ann_file=$5


line=$(tail -1 $file)

Rscript $path_Rscript'plot_LRR-BAF_single_run.R' --inputfold ${path}outdir_${sample_name}_${line}/ --line_id ${line} --chr $chr --outf ${path}/outdir_${sample_name}_${line}/ --ann_file ${ann_file}


fi


