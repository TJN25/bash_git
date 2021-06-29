#!/bin/bash

##makes alignments and running alifoldz and r-scape
##GCA accession number.

usage(){
    echo "sraAlignAndFold.sh is a script for making a multiple sequence alignment and
    getting secondary structure information.  
Usage:
 fetchGenomeGCA.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-r	Folder location

"
}

while getopts "i:o:h" arg; do
case $arg in
	i) 
	FOLDER=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${FOLDER} ]; then
	FOLDER="./"
fi




mkdir -p "$FOLDER/alifold/post_script"
mkdir -p "$FOLDER/RNAAlifold"
mkdir -p alignments_rnaalifold


let "fileNum = 0"
for file in alignments/*.stk

do
lines=`wc -l < $file`
if (( $lines < 1));then
continue
fi

outname=`basename $file`


if [ -f "$FOLDER/RNAAlifold/$outname.rnaalifold" ]; then
	echo "Already exists: $file"
	continue
fi



nseqs=`esl-alistat $file | grep "Number of sequences" | cut -d ":" -f2`



length=`esl-alistat $file | grep "Alignment length:" | cut -d ":" -f2`
largest_length=`esl-alistat $file | grep "Largest" | cut -d ":" -f2`


if (( $length < 500 )); then

diffLength=`expr $largest_length - $length`

if (( $diffLength > $length ));then
echo "Alignment is poor: $file"
continue
fi

if (( $nseqs > 5000 )); then
echo "Skipping $file (length: $length, nseqs: $nseqs)"
echo "$file (length: $length, nseqs: $nseqs)" >> skipped_alignments.txt
continue
fi


echo "Running RNAAlifold on $file (length: $length, nseqs: $nseqs)"

> tmp.stk

esl-alimanip   --lnfract 0.8 --lxfract 1.2 --lmax 500 --detrunc 30 ${file} | esl-alimask -g --gapthresh 0.8 -p --pfract 0.5 --pthresh 0.5 --keepins - | grep -v "SS_cons" > tmp.stk


RNAalifold -p -r -d2 --SS_cons --noLP --color --aln-stk=${file} tmp.stk >> ./RNAAlifold/$outname.rnaalifold
cat alirna.ps > ./alifold/post_script/$outname.ps      
else
	echo "Skipping: $file. Length: $length"
fi

done


for file in alignments_R*; #change to _R for positive control
do
outname=`basename $file .stk.stk`

mv $file ./alignments_rnaalifold/$outname.stk

done
