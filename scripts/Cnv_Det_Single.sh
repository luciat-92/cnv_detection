#!/bin/bash

##################################### Quality control for the CNV and plot the delition and duplication detected  #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages: argparse, ggplot2, grid, gridExtra' 
echo 'this script should be run like this:'
echo './Cnv_Det_Single.sh <path R script> <path folder> <path .txt file with SentrixBarcode_SentrixPosition id> <name sample considered>' 
else


path_Rscript=$1
path=$2
file=$3
sample_name=$4

id=$(tail -n 1 $file)
sex=$(head -n 1 $file)
# echo $sex

[[ -d $path'summary_res/' ]] || mkdir $path'summary_res/'

Rscript $path_Rscript'summary_CNSingleAnalysis_run.R' --inputdir $path  --fold_fun $path_Rscript  --line $id --sample_name ${sample_name} --sex $sex --outf $path'summary_res/'

fi
