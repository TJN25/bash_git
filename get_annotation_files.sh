#!/bin/bash

usage(){
    echo "get_annotation_files.sh is a script for downloading a  GFF file from a GCA or GCF accession.  
Usage:
 get_annotation_files.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-g	Genome accession (required)
"
}

while getopts "g:h" arg; do
case $arg in
	g) 
	input_name=${OPTARG};;

    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${input_name} ]; then
    echo "Error: No input specified." >&2
    usage
    exit 1
fi
if [ ! -f ${input_name}.gff ];then
main_link=$(esearch -db assembly -query ${input_name} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)


file_name=`echo ${main_link} | rev | cut -d '/' -f1 | rev`

gffLink="${main_link}/${file_name}_genomic.gff.gz"


    
curl $gffLink > ${input_name}.gff.gz 
sleep 1
gunzip ${input_name}.gff.gz 

else

echo "${input_name}.gff aready exists"

fi