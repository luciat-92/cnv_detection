#!/bin/bash

##################################### Plot detected CNV for a Sample in different comparisons (only QC CNV)  #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages: argparse, ggplot2, grid, gridExtra' 
echo 'this script should be run like this:'
echo './Cnv_Plot.sh <path R script> <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered> ' 
else


path_Rscript=$1
path=$2
file=$3
sample_name=$4



nrow=$(wc -l < $file)
thr=$(($nrow - 2))

pre=$(head -n 2 $file | tail -n 1)
post=$(tail -n $thr $file)

sex=$(head -n 1 $file)
# echo $sex



Rscript $path_Rscript'plot_summary_cnv_QC_run.R'  --post ${post} --pre ${pre} --sample_name ${sample_name} --sex $sex --outf forPlot_cnv_${sample_name} --width_plot 20 --fold_fun $path_Rscript --inputdir $path

fi
