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
mkdir -p "$FOLDER/rscape_out"
mkdir -p "$FOLDER/RNAAlifold"


let "fileNum = 0"
for file in alignments/*.stk

do

outname=`basename $file`

#echo "${FOLDER}alifold/$outname.alifold"

if [ -f "$FOLDER/alifold/$outname.alifold" ]; then
	echo "Already exists: $file"
	continue
fi




echo "Running on: $file"


start=`grep "GCA" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
end=`grep "GCA" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

if [[ $start == "" ]]; then
	start=`grep "NC_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NC_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	start=`grep "NZ_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NZ_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	head $file
fi




length=`expr $end - $start`

if (( $length < 0 )); then
	length=$(( -1 * $length ))
fi

if (( $length < 750 )); then
esl-reformat  clustal $file  | RNAalifold --aln-stk=${file} >> ./RNAalifold/$outname.rnaalifold
cat alirna.ps > ./alifold/post_script/$outname.ps      
esl-reformat  clustal $file  | alifoldz.pl > ./alifold/$outname.alifold         

R-scape --r2rall --outdir ./rscape_out/ $file

else
	echo "Skipping: $file"
fi

done
