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
	-e extension	       
	-a	align
	-o output
	-k keep tmp files

"
}

align="F"
keep="F"
while getopts "e:ao:hk" arg; do
case $arg in
	e)
		extension=${OPTARG}
		;;
	a) 
		align="T"
		;;
    h)
		usage
		exit
      ;;    
    k)
		keep="T"
      	;;   
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done


let "fileNum = 0"

for file in *.$extension

do

if [[ $align == "T" ]]; then

echo "Aligning"
esl-reformat  clustal $file > $file.clustal

echo "Running rmfam_scan"

rmfam_scan.pl -g -f ~/phd/RNASeq/RMfam/scripts/RMfam.cm $file.clustal -o $file.out
#mv *.cmscan ${file}-aligned.cmscan
#mv *.tblout ${file}-aligned.tblout

else

echo "Running rmfam_scan"

rmfam_scan.pl -g -f ~/phd/RNASeq/RMfam/scripts/RMfam.cm $file -o $file

#mv *.cmscan $file.cmscan
#mv *.tblout $file.tblout

fi 

if [[ $keep == "T" ]]; then

mkdir -p "rmfam_scan_tmp"
mv *.mcl rmfam_scan_tmp
mv *.clustal rmfam_scan_tmp


else

rm *.mcl
rm *.clustal


fi

mkdir -p "rmfam_cmscan"
mkdir -p "rmfam_gff_files"
mkdir -p "rmfam_tblout"
mv *.tblout rmfam_tblout
mv *.cmscan rmfam_cmscan
mv *.gff rmfam_gff_files

#rmfam_scan.pl -f ~/phd/RNASeq/RMfam/scripts/RMfam.cm sra_enterics-serratia-1007_1.fasta -o sra_enterics-serratia-1007_1.tblout

done

cd rmfam_tblout

let "fileNum = 0"

for file in *.tblout

do

sed 's/  /__/g' $file | sed 's/ /--/g' | sed 's/__/ /g' | tr -s ' ' | sed 's/ --/ /g' | sed 's/-- / /g' | sed 's/--!/ !/g' | tr ' ' '\t' | sed 's/--/_/g' > $file.formatted.tblout

done

