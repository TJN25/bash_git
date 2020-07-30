#!/bin/bash
echo "Script is starting"
let "fileNum=0"
for file in ../../../../DB/RefSeq83/bacteria/*/*/*_genomic.gbff
do 
echo $file
echo "Script is running"
parse_gbff.1.0.py $file
genome_protein_name_lookup.1.0.py $file
cut -d '"' -f1 tmp.protein.txt > tmp.protein.faa
hmmsearch --tblout tmp.hmm.out -E 1e-5 --cpu 3  --noali ../DB/CRISPR_models.hmm tmp.protein.faa >> output.log 
rm tmp.protein.faa
rm tmp.protein.txt
cat tmp.hmm.out >> ../Working/predicted_cas_genes_in_bacteria_of_refseq83.txt
#rm tmp.hmm.out
#rm output.log
done