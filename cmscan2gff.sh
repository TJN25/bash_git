#!/bin/bash

##Downloads the genome sequence (and gff file if the -g flag is used) when given a 
##GCA accession number.

usage(){
    echo "cmscan2gff.sh is a script for converting the output of a cmscan (--fmt 2) into a gff file.  
Usage:
 cmscan2gff.sh [input]

Options:
	-h	Display this help

Input	       
	-c cmscan file (required)
	-o Output name
	-g GCA accession

"
}

while getopts "c:o:g:h" arg; do
case $arg in
	c) 
	CMSCAN=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
	g) 
	GCA=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${CMSCAN} ]; then
    echo "Error: No input specified." >&2
    usage
    exit 1
fi

if [ -z ${OUTPUT} ]; then

OUTPUT=$(basename $CMSCAN .tblout).gff

fi


echo "##gff-version 3" > $OUTPUT
echo "#!gff-spec-version 1.21" >> $OUTPUT
echo "#!processor R script (local) with manual add of top section" >> $OUTPUT
sleep 1
echo '#!genome-build' "$(esearch -db assembly -query $(basename $OUTPUT .gff) | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyName)" >> $OUTPUT
sleep 1
echo '#!genome-build-accession NCBI_Assembly:' "$(basename $OUTPUT .gff)" >> $OUTPUT
echo "#!annotation-source rfam (local version)" >> $OUTPUT
echo '##sequence-region' "$(cat $CMSCAN | tr -s ' ' | cut -d ' ' -f4 | sort | grep "_" | grep -v ".f" | uniq) 1 4941439" >> $OUTPUT
sleep 1
echo '##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' "esearch -db assembly -query $(basename $OUTPUT .gff) | efetch -format docsum | xtract -pattern DocumentSummary -element Taxid" >> $OUTPUT
sleep 1
