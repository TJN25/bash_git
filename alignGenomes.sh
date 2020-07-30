#!/bin/bash

##Aligns genomes according to a .phy or .nwk file
##GCA accession number.

usage(){
    echo "alignGenomes.sh is a script for aligning a set of genomes given an order to do this in. 
Usage:
 alignGenomes.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-i	input file

"
}

while getopts "i:o:h" arg; do
case $arg in
	f) 
	INFILE=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

