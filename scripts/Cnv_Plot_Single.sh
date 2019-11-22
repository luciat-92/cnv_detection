#!/bin/bash

##################################### Plot detected CNV for a Sample in a single analysis (only QC CNV)  #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages: argparse, ggplot2, grid, gridExtra' 
echo 'this script should be run like this:'
echo './Cnv_Plot_Single.sh <path R script> <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered> ' 
else


path_Rscript=$1
path=$2
file=$3
sample_name=$4


id=$(tail -n 1 $file)
sex=$(head -n 1 $file)
# echo $sex

Rscript $path_Rscript'plot_single_cnv_QC_run.R' --inputdir $path --fold_fun $path_Rscript  --line $id --sample_name ${sample_name} --sex $sex --outf $path

fi
