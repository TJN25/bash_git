#!/bin/bash

##makes alignments and running alifoldz and r-scape
##GCA accession number.

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

while getopts "i:o:h" arg; do
case $arg in
	i) 
	FOLDER=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${FOLDER} ]; then
	FOLDER="./"
fi




mkdir -p "$FOLDER/alifold/post_script"



let "fileNum = 0"
for file in alignments_rnaalifold/*.stk


do

outname=`basename $file`
if [ -f "$FOLDER/alifold/$outname.alifold" ]; then
	#echo "$FOLDER/alifold/$file.alifold"
	echo "Already exists: $file"
	continue
else
	echo "Checking size: $file"
fi


lines=`wc -l < $file`
if (( $lines < 1 ));then
continue
fi


nseqs=`grep "#=" $file | cut -d ' ' -f2 | sort | uniq | wc -l`


start=`grep "GCA" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f1`
end=`grep "GCA" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f2`

if [[ $start == "" ]]; then
	start=`grep "NC_" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NC_" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	start=`grep "NZ_" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NZ_" $file | head -n 1 | cut -d " " -f1 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	head $file
fi




length=`expr $end - $start`

if (( $length < 0 )); then
	length=$(( -1 * $length ))
fi

if (( $length < 500 )); then

ID=`grep "NC_" $file | head -n 1 | cut -d " " -f2`

if [[ $ID == "" ]]; then
	ID=`grep "NZ_" $file | head -n 1 | cut -d " " -f2`
fi

if [[ $ID == "" ]]; then
	ID=`grep GCA_" $file | head -n 1 | cut -d " " -f2`
fi

if [[ $ID == "" ]]; then
	head $file
else
alignmentLength=`grep $ID $file | grep -v "#" | tr -s ' ' | cut -d ' ' -f2 | wc -c`

diffLength=`expr $alignmentLength - $length`

if (( $diffLength > $length ));then
echo "Alignment is poor: $file"
continue
fi
fi




echo "Running alifoldz.pl on $file (length: $length, nseqs: $nseqs)"

time esl-reformat  clustal $file  | alifoldz.pl > ./alifold/$outname.alifold         
 

else

	echo "Skipping: $file"


fi

done







