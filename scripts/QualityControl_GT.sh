#!/bin/bash

################## Quality control for genotype data ###############################


if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: argparse, ggplot2, ggrepel, RColorBrewer, doParallel' 
echo 'this script should be run like this:'
echo './QualityControl_GT.sh <path to Rscripts folder> <manifest file> <file with Pseudo Autosomal Region> <path output folder> <n cores for the parallelization> <n header line manifest>'
echo '<n header line manifest> is the number of lines to be skipped in the manifest file, the colnames line has to be included. Default 7'
else

path_Rscript=$1
manifest_file=$2
PAR_file=$3
path=$4
ncores=$5

if [ -v $6 ]
then
nhead=7
else
nhead=$6
fi



manifest_name=$(echo $manifest_file | cut -d"/" -f4)
manifest_name=$(echo $manifest_name | cut -d"." -f1)
manifest_name=$(echo "$manifest_name" | tr - .)

### NOTE: if the manifest_name gives an error in QC_extract_LRR_BAF_run.R, check columnnames in SNP_Table.txt 

Rscript $path_Rscript'PAR_SNP_build_run.R' --outf $path --manifest_file $manifest_file --PAR_file $PAR_file --n_header_line $nhead

Rscript $path_Rscript'QC_extract_LRR_BAF_run.R' --outf $path --manifest_name $manifest_name --SNP_table_file $path'SNP_Table.txt' --Full_table_file $path'Full_Data_Table.txt' --Sample_table_file $path'Samples_Table.txt' --PAR_file $path'PAR_SNPs.txt' --ncores $ncores --fold_fun $path_Rscript

Rscript $path_Rscript'plot_sample_QC_run.R' --outf $path --sample_table_file $path'Samples_Table_filt.txt' --SNP_table_info_file $path'info_QC.txt'  

fi

