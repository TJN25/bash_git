#!/bin/bash

usage(){
    echo "run_sRNA_nhmmer.sh is a script for running nhmmer and sorting the results.  
Usage:
 run_sRNA_nhmmer.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-r	Folder location

"
}

while getopts "d:h" arg; do
case $arg in
	d) 
	database=${OPTARG};;

    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

echo "checking for short seqs"
mkdir -p single_seqs/short_seqs
let "fileNum = 0"
for file in single_seqs/*.fna

do

seqlength=`grep -v ^">" $file | wc -c`
if (( $seqlength < 50 )); then
mv $file  single_seqs/short_seqs/
fi

done

echo "making directories in `pwd`"
mkdir -p alignments
mkdir -p hmm
mkdir -p output


echo "running nhmmer"
let "fileNum = 0"
for file in single_seqs/*.fna

do


outname=`basename $file .fna`

nhmmer -E 1e-5 --tblout output/$outname.tbl -A alignments/tmp.stk --tformat FASTA --qformat FASTA $file $database > output/$outname.res

esl-alimanip   --detrunc  60  alignments/tmp.stk > alignments/$outname.stk

hmmbuild hmm/$outname.hmm alignments/$outname.stk > output/hmm.out



done




