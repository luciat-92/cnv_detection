#!/bin/bash

##################################### Match sample base on genotype and update annotation file #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: argparse, circlize, ComplexHeatmap, RColorBrewer' 
echo 'this script should be run like this:'
echo './GT_match.sh <path to Rscripts folder> <path output folder> <project_name> <path annotation file>'
else

path_Rscript=$1
path=$2
project_name=$3
annotation_file=$4

Rscript $path_Rscript'GT_match_run.R'  --comp_GT_file $path'GT_comp.txt' --Samples_table_file $path'Samples_Table_filt.txt' --ann_file $annotation_file  --outf $path$project_name'_'

fi

