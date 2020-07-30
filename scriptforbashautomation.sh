#!/bin/bash
cd 
cd Desktop/forming-arCOG-pssms/seqs/muscle/mafft
let "fileNum = 0"
for file in *

do

command /Users/thomasnicholson25/Downloads/ncbi-blast-2.2.31+/bin/psiblast -db /Users/thomasnicholson25/Desktop/python/ar14.fa -in_msa $file -out_pssm $file.pssm -num_iterations 4 &

done

mkdir pssm
mv *.pssm pssm


#hmmpred comparing hmms.
#take a look hmmer
