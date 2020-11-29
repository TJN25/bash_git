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




echo "Running on: $file"


nseqs=`esl-alistat $file | grep "Number of sequences" | cut -d ":" -f2`



length=`esl-alistat $file | grep "Alignment length:" | cut -d ":" -f2`
largest_length=`esl-alistat $file | grep "Largest" | cut -d ":" -f2`


if (( $length < 500 )); then

diffLength=`expr $largest_length - $length`

if (( $diffLength > $length ));then
echo "Alignment is poor: $file"
continue
fi


echo "Running RNAAlifold on $file (length: $length, nseqs: $nseqs)"


esl-reformat  clustal $file  | RNAalifold --aln-stk=${file} >> ./RNAAlifold/$outname.rnaalifold
cat alirna.ps > ./alifold/post_script/$outname.ps      
else
	echo "Skipping: $file. Length: $length"
fi

done


for file in alignments_G*; #change to _R for positive control
do
outname=`basename $file .stk.stk`

mv $file ./alignments_rnaalifold/$outname.stk

done
