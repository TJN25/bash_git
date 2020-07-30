#!/bin/bash

##Aligns genomes according to a .phy or .nwk file
##GCA accession number.

usage(){
    echo "generateRandomSRACalls.sh is a script for aligning a set of genomes given an order to do this in. 
Usage:
 generateRandomSRACalls.sh [opts] [input] [output]

Options:
	-h	Display this help

Input	       
	-g	GCA Accession
	-m Method (shuffled, shifted or both)
	-o output name

"
}

method="shifted"

while getopts "g:o:m:h" arg; do
case $arg in
	g) 
	in_file=${OPTARG}
	;;
	o) 
	out_file=${OPTARG}
	;;
	m)
	method=${OPTARG}
	;;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [[ -z $in_file ]]; then
echo 'Error: GCA needed. Specify with -i <gca>'
echo ' '
echo 'Use -h for more help.'
echo ' '
exit
fi

if [[ -z $out_file ]]; then
	out_file=$in_file
fi

if [[ $method == "both" ]]; then

echo "Shuffling Data"

selectRandomSRA.R -i $in_file -m shifted -o ${in_file}_shifted -q
selectRandomSRA.R -i $in_file -m shuffled -o ${in_file}_shuffled -q
randomDataReformat.R -i ${in_file}_shifted -g $in_file -q
randomDataReformat.R -i ${in_file}_shuffled -g $in_file -q

echo "Combining shifted calls"

cd ./random_sequences/shifted/
mkdir -p gff_files
mv ../../${in_file}_shifted_random_sra.gff gff_files
rm ./gff_files/${in_file}.gff
#cp ../../${in_file}.gff gff_files
combine_gff_files.R -f gff_files/ -o ${in_file}_shifted_random -r ${in_file}.gff -g ${in_file}.gff
mv ../../${in_file}_shifted_random* ./


echo "Combining shuffled calls"

cd ../../random_sequences/shuffled/
mkdir -p gff_files
mv ../../${in_file}_shuffled_random_sra.gff gff_files
rm ./gff_files/${in_file}.gff
#cp ../../${in_file}.gff gff_files
combine_gff_files.R -f gff_files/ -o ${in_file}_shuffled_random -r ${in_file}.gff -g ${in_file}.gff
mv ../../${in_file}_shuffled_random* ./

fi

if [[ $method == "shifted" ]]; then

echo "Shuffling Data using 'shifted' method"
selectRandomSRA.R -i $in_file -m shifted -o ${in_file}_shifted -q

randomDataReformat.R -i ${in_file}_shifted -g $in_file -q

echo "Combining shifted calls"

cd ./random_sequences/shifted/
mkdir -p gff_files
mv ../../${in_file}_shifted_random_sra.gff gff_files
rm ./gff_files/${in_file}.gff
#cp ../../${in_file}.gff gff_files

combine_gff_files.R -f gff_files/ -o ${in_file}_shifted_random -r ${in_file}.gff -g ${in_file}.gff
mv ../../${in_file}_shifted_random* ./

fi

if [[ $method == "shuffled" ]]; then

echo "Shuffling Data"

selectRandomSRA.R -i $in_file -m shuffled -o ${in_file}_shuffled -q
randomDataReformat.R -i ${in_file}_shuffled -g $in_file -q

echo "Combining shuffled calls"

cd ./random_sequences/shuffled/
mkdir -p gff_files
mv ../../${in_file}_shuffled_random_sra.gff gff_files
rm ./gff_files/${in_file}.gff
#cp ../../${in_file}.gff gff_files
combine_gff_files.R -f gff_files/ -o ${in_file}_shuffled_random -r ${in_file}.gff -g ${in_file}.gff
mv ../../${in_file}_shuffled_random* ./

fi






