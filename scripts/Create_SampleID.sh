#!/bin/bash

##################################### For each considered sample create file with SentrixBarcode_SentrixPosition ID  #####################

if [ -v $1 ]
then
echo 'to run this script you need R with the following packages installed: argparse' 
echo 'this script should be run like this:'
echo './Create_SampleID.sh <path to Rscripts folder> <path output folder> <path annotation file> <File with samples to consider>'
else

path_Rscript=$1
path=$2
annotation_file=$3
samples_file=$4

samples=$(cat $samples_file) 
# echo ${samples[@]}

for sample_name_id in ${samples[@]}
do
	echo $sample_name_id
	Rscript $path_Rscript'create_annSample_run.R' --sample_name ${sample_name_id} --ann_file $annotation_file --outf $path'ID_'

done

fi
