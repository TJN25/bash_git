#!/bin/bash

##Run rmfam over all fasta files

usage(){
    echo "run_rmfam_scan.sh is a script for running rmfam_scan over a set of fasta
    files.  
Usage:
 run_rmfam_scan.sh [opts] [extension]

Options:
	-h	Display this help

Input
	-i extension	       

"
}

while getopts "i:" arg; do
case $arg in
	i)
		input=${OPTARG}
		;;

	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done


selectRandomSRA.R -i $input

 mkdir random_sequences

mv ${input}_random* random_sequences/
