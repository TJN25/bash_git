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




mkdir -p "$FOLDER/rscape_out"



let "fileNum = 0"
for file in alignments_G*;
do

checkname=`basename $file .stk`
if [ -f "../rscape_out/${checkname}_1.R2R.sto" ]; then
echo "Already exists: $file"
continue
else

echo "Running R-scape on: $file"
	
fi

time R-scape --r2rall --outdir ../rscape_out/ $file > rsacpe.out

done

