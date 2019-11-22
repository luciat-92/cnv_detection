#!/bin/bash
################## GETTING CORRECT ALT/REF ###############################

##################################################################
# prepare the manifest file and extract the reference allele     #
#             from the reference genome fasta file               #
##################################################################

if [ -v $1 ]
then
echo 'to run this script you need to have bedtools installed and available globally'
echo 'this script should be run like this:'
echo './GetRefAlt.sh <path to manifestFile.csv> <path to fasta> <number_of_header_lines> <numer_of_tail_lines>'
echo 'the header and tail lines are the ones at the begining and end that do not have a snp record but other info'
else
a=$(wc -l $1 | awk '{print $1}')
b=`expr $a - $3`
c=`expr $b - $4`
echo $c
echo $b
tail -$b $1 | head -$c | perl -pe 's/,/\t/g' | cut -f 10-11 | awk 'BEGIN {OFS="\t"}{print "chr"$1,$2-1,$2}' | grep -v "chrXY" | grep -v "chrMT" | grep -v "chr0" > $1'.bed'
bedtools getfasta -fi $2 -bed $1'.bed' > $1'.ref'
perl -pe 's/\n|-|:/\t/g' $1'.ref' | perl -pe 's/>chr/\n/g' | awk 'BEGIN {OFS="\t"}{if(NR>=2) {print $1,$3,$4}}' |perl -pe 's/^X/23/g' | perl -pe 's/^Y/24/g'|  sort -k1,1 -k2,2n > $1'.tmp'
awk '{if($3=="a") {$3="A";print $0} else if($3=="t") {$3="T";print $0} else if($3=="c") {$3="C";print $0} else if($3=="g") {$3="G";print $0} else {print $0}}' $1'.tmp' > $1'.tmp2'
perl -pe 's/ +/\t/g' $1'.tmp2' > $1'_reference_alleles'
rm  $1'.bed' $1'.ref' $1'.tmp' $1'.tmp2'
fi
