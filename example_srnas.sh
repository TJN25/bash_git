#!/bin/bash

usage(){
    echo "example_srnas.sh gets info read depths and neighbouring genes from example sRNA
Usage:
 example_srnas.sh	[opts] [extension]

Options:
    -h  Display this help

Input
    -s srna_id
    -i write info           
    -o output

"
}

writeinfo="F"

while getopts "s:o:ih" arg; do
case $arg in
    s)
        srna=${OPTARG}
        ;;
    i) 
        writeinfo='T'
        ;;
    o) 
        outname=${OPTARG}
        ;;
    h)
        usage
        exit
      ;;    
    \?) 
    echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [[ -z $srna ]]; then
echo 'Error: sRNA ID needed. Specify with -s <srna_id>'
echo ' '
echo 'Use -h for more help.'
echo ' '
exit
fi


if [ -z ${outname} ]; then
    outname=${srna}
fi

echo "Making folders"

mkdir -p ${outname}_example_files
mkdir -p tmp


current_folder=`pwd`

echo "Plotting genes"
get_example_srnas.R ${srna} ${outname} ${current_folder}

> ${current_folder}/${outname}_example_files/${outname}_read_values.txt


echo "Getting read depths"
read_depths_examples.sh ${srna} ${outname} ${current_folder}

echo "Plotting read depths"
examples_heatmaps.R ${srna} ${outname} ${current_folder}