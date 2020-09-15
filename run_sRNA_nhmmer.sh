#!/bin/bash

usage(){
    echo "run_sRNA_nhmmer.sh is a script for running nhmmer and sorting the results.  
Usage:
 run_sRNA_nhmmer.sh [opts] [input]

Required:	       
	-d	<database> The nucleotide database to be searched against as the nhmmer target
	-f	<folder> The folder containing the query file for nhmmer
	-e	<extension> The file extension of the query file in <folder>

Options:
	-h	Display this help
	-c	Check if the sequences are >50 nt (only works on fasta files)
"
}

check_seq_lengths=""
database=""
folder=""
extension=""

missing=""
exitTrue="F"
evalue="1e-5"
while getopts "d:f:E:e:ch" arg; do
case $arg in
	d) 
	database=${OPTARG};;
	f)
	folder=${OPTARG};;
	E)
	evalue=${OPTARG};;
	e)
	extension=${OPTARG};;	
	c)
	check_seq_lengths="T";;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [[ $database == "" ]]; then
	missing="$missing -d <database> "
	exitTrue="T"
fi


if [[ $folder == "" ]]; then
	missing="$missing -f <folder> "
	exitTrue="T"	
fi

if [[ $extension == "" ]]; then
	missing="$missing -e <extension> "
	exitTrue="T"
fi

if [[ $exitTrue == "T" ]];then
	echo "$missing not found. Use -h for more help."
	exit
fi


if [[ $extension == "stk" ]]; then
	check_seq_lengths=""
fi

echo "extension value: $extension"

if [[ $check_seq_lengths == "T" ]]; then

	echo "checking for short seqs"
	mkdir -p $folder/short_seqs
	
	let "fileNum = 0"
	for file in $folder/*.$extension
	do
	
	seqlength=`grep -v ^">" $file | wc -c`
	if (( $seqlength < 50 )); then
	mv $file  single_seqs/short_seqs/
	fi
	
	done

fi

echo "making directories in `pwd`"
mkdir -p alignments
mkdir -p hmm
mkdir -p output

echo "running nhmmer"
let "fileNum = 0"
for file in $folder/*.$extension

do

echo $file
outname=`basename $file`


if [ -f "output/$outname.res" ]; then
	echo "Already exists: $file"
	continue
fi




nhmmer -E $evalue --tblout output/$outname.tbl -A alignments/tmp.stk --tformat FASTA  $file $database > output/$outname.res

lines=`wc -l < alignments/tmp.stk`


if (( $lines > 0 )); then

esl-alimanip   --detrunc  60  alignments/tmp.stk > alignments/$outname.stk

hmmbuild hmm/$outname.hmm alignments/$outname.stk

else

echo "No hits found for $outname"

fi

done



