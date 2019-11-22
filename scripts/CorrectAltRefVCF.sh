#!/bin/bash
######################################################################################
###############clean the vcf file to get rid of the XY and MT snps ###################
######################################################################################

if [ -v $1 ]
then
echo 'this script should be run like this:'
echo './CorrectAltRefVCF.sh <path to vcf> <path to reference alleles list file>'
else
vcf=$1
reference_alleles=$2
grep '#' $1 > $1'_vcfHeader'
######### correct the samples in the header ######
perl -pe 's/\t\d+_/\t/g' $1'_vcfHeader' | grep -v '##contig=<ID=25' | grep -v '##contig=<ID=26'| grep -v '##contig=<ID=0' > $1'_correctedvcfHeader'
# now correct the vcf body
grep -v '#' $1 | grep -v '^25'| grep -v '^26'| grep -v '^0'| sort -k1,1 -k2,2n  > $1'_vcfbody'
paste $2 $1'_vcfbody' > tmp
#control that the number of lines are the same
if [ "$(wc -l < tmp)" -eq "$(wc -l < $2)" ]
then
awk 'BEGIN {OFS="\t"}{if (($3==$7) && ($1==$4) && ($2==$5)) {print $0} else if (($3!=$7) && ($1==$4) && ($2==$5)) {$8=$7;$7=$3;print $0} else if (($1!=$4) || ($2!=$5)) {print "WARNING: " $1,$2 " in ref file does not match the vcf"}}' tmp | cut -f4- > tmp2
echo $(grep 'WARNING' tmp2)
cat $1'_correctedvcfHeader' tmp2 > $1'_RefAltCorrected.vcf'
else
echo 'Warning: No Matching number of lines, are you sure you are using the correct manifest file?'
fi
rm $1'_correctedvcfHeader'  $1'_vcfHeader' tmp $1'_vcfbody' tmp2
fi
