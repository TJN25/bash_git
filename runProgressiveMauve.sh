#!/bin/bash

usage(){
    echo "runProgressiveMauve.sh is a script for aligning a set of genomes. 
Usage:
 runProgressiveMauve.sh [opts] [input] [output]

Options:
	-h	Display this help

Input	       
	-i	file containing GCA Accession list
	-o output name

"
}

method="shifted"

while getopts "i:o:m:h" arg; do
case $arg in
	i) 
	in_file=${OPTARG}
	;;
	o) 
	out_file=${OPTARG}
	;;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done


gca_list=""

while read p; do
  gca_list="$gca_list ~/phd/RNASeq/fasta_files/${p}.fna"
done < $in_file





code_to_run=`echo "progressiveMauve  --output=$out_file --output-guide-tree=$out_file.guide_tree --backbone-output=$out_file.backbone  $gca_list"`
eval $code_to_run