#!/bin/bash

##################################### Quality control for the CNV and plot difference in the CN detection  #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages: argparse, ggplot2, grid, gridExtra' 
echo 'this script should be run like this:'
echo './Cnv_Diff.sh <path R script> <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered> <path annotation file>' 
else


path_Rscript=$1
path=$2
file=$3
sample_name=$4
annotation_file=$5


nrow=$(wc -l < $file)
thr=$(($nrow - 2))

pre=$(head -n 2 $file | tail -n 1)
post=$(tail -n $thr $file)

sex=$(head -n 1 $file)
# echo $sex

[[ -d $path'summary_res/' ]] || mkdir $path'summary_res/'

Rscript $path_Rscript'summary_diffCN_1comp_run.R' --fold_fun $path_Rscript --inputdir $path --ann_file $annotation_file --post ${post} --pre ${pre} --sample_name ${sample_name} --sex $sex --outf $path'summary_res/forPlot_'

fi
