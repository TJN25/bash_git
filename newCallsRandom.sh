#!/bin/bash

while getopts "g:" arg; do
case $arg in
	g) 
	gff=${OPTARG};;
      
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

cp ../$gff.gff ./
randomDataReformat.R -g $gff


mkdir -p "gff_files"
cp ../$gff.gff ./gff_files
cp ${gff}_random_sra.gff ./gff_files

combine_gff_files.R -f gff_files/ -o ${gff}_random -r -g $gff.gff 