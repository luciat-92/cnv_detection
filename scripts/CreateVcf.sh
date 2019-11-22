#!/bin/bash

if [ -v $1 ]
then
echo 'to run this script you need to have plink installed'
echo 'this script should be run like this:'
echo './CreateVcf.sh <path folder> <path installation plink> <path genome studio .map .ped files>'
else


path=$1
path_plink=$2
genomestudio_output=$3

cd $path

$path_plink'plink' --file $genomestudio_output  --recode vcf
fi

