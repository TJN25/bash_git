#!/bin/bash

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

while getopts "d:o:f:h" arg; do
case $arg in
	d) 
	database=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
	f) 
	out_folder=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done


mkdir -p short_seqs
let "fileNum = 0"
for file in *.fasta

do

seqlength=`grep -v ^">" $file | wc -c`
if (( $seqlength < 50 )); then
mv $file  short_seqs/
fi

done

for file in *.fasta

do

cat $file >> ${OUTPUT}.fna

done


mkdir -p $out_folder
let "fileNum = 0"
for file in *.fasta

do


outname=`basename $file .fasta`

nhmmer -E 1e-5 --tblout $out_folder/$outname.tbl -A $out_folder/tmp.stk --tformat FASTA --qformat FASTA $file $database > $out_folder/$outname.res

esl-alimanip   --detrunc  60  $out_folder/tmp.stk > $out_folder/$outname.stk

hmmbuild $out_folder/$outname.hmm $out_folder/$outname.stk

cat $out_folder/$outname.hmm >> $OUTPUT.hmm 

done

rm $out_folder/$OUTPUT.tbl.out
let "fileNum = 0"
for file in $out_folder/*.tbl

do

cat $file >> $out_folder/$OUTPUT.tbl.out

done


