esl-alimanip --lnfract 0.8 --lxfract 1.2 --lmin 50 --lmax 500 --detrunc 30 GCA_000017745.1_101.stk.stk



for file in *.stk
do
outname=`basename $file .stk`

esl-alimanip --lnfract 0.8 --lxfract 1.2 --lmin 50 --lmax 500 --detrunc 30 $file > ../alignments_filtered/$outname.stk

done

for file in *.stk

do 

outname=`basename $file .fna.stk`

grep ^"#=GS" $file | cut -d " " -f2 | cut -d "[" -f1 | sort | uniq | grep -v  $outname > list_to_move.txt

while read line; do mv $line.fna.stk duplicates; done < list_to_move.txt

done



for file in *.stk.stk;

do 

lines=`wc -l < $file`
if (( $lines < 1));then
echo "$file has no data."
continue
fi
echo $file
outname=`basename $file .stk`
esl-alimask -g --gapthresh 0.49 -p --pfract 0.5 --pthresh 0.5 --keepins $file | esl-alimanip   --lnfract 0.6 --lxfract 1.4 --lmin 50 --lmax 500 --detrunc 30 - > ../alignments_filtered/$outname
done



run_sRNA_nhmmer.sh -d /Users/thomasnicholson/phd/RNASeq/srna_seqs/version_1/positive_control.fna -f alignments_filtered/ -o overlapping_models -e stk -E 1e-5


nhmmer -E $evalue --tblout tmp.tbl -A tmp.stk --tformat FASTA  $file $database

for file in *.stk;
do 

ID=`echo $file | cut -d '.' -f1,2`

grep ^"#=GS" $file | sort | uniq | cut -d "[" -f1 | sed -e "s/$/   $ID/" >> ../predicted_genomic_sequence_matches.txt

done






for file in *.stk;
do
lines=`wc -l < $file`
if (( $lines < 1));then
echo "$file has no data."
continue
fi
done






for file in *.txt;
do
max_num="0"
outname=`basename $file _combined_list.txt`
echo $outname
> fasta/$outname.fa
while read line; 
do 

esl-reformat fasta ../alignments_filtered/$line.stk >> fasta/$outname.fa

current_num=`esl-alistat ../alignments_filtered/$line.stk | grep "Number of sequences" | cut -d ' ' -f 4`

if(( $current_num > $max_num ));
then

max_seq="$line.stk"
max_num="$current_num"
fi

done < $file

hmmbuild hmm/$outname.hmm ../alignments_filtered/$max_seq

hmmalign --informat fasta hmm/$outname.hmm fasta/$outname.fa | esl-alimask -g --gapthresh 0.8 -p --pfract 0.5 --pthresh 0.5 - | esl-alimanip   --lnfract 0.6 --lxfract 1.4 --lmin 50 --lmax 500 --detrunc 50 - > combined_alignments_2/$outname.stk

done



> ../../original_stats_3.txt
> ../../new_stats_3.txt
for file in *;
do 

esl-alistat ../../alignments_filtered/$file | sed 's/%//g' | sed 's/Alignment /Alignment_/g' | sed 's/Average /Average_/g' | sed 's/Format: /Format:_/g' | sed 's/Number of sequences/Number_of_sequences/g' | sed -e "s/$/ $file/" >> ../../original_stats_3.txt


lines=`wc -l < $file`
if (( $lines < 1));then
# echo $file
continue
fi

esl-alistat $file | sed 's/%//g' | sed 's/Alignment /Alignment_/g' | sed 's/Average /Average_/g' | sed 's/Format: /Format:_/g' | sed 's/Number of sequences/Number_of_sequences/g' | sed -e "s/$/ $file/" >> ../../new_stats_3.txt

done


while read line; 
do 

mv $line ids_worse_after_merge/

done < $file


while read line; 
do 

name=`basename $line .stk`
echo $name 

while read line_2;
do

echo $line_2

cp ../alignments_filtered/${line_2}.stk combined_alignments_2/
done < ${name}_combined_list.txt

done < $file



promptValue() {
  read -p "$1"": " val
   $val
}



for file in *;
do 

esl-alistat $file | sed 's/%//g' | sed 's/Alignment /Alignment_/g' | sed 's/Average /Average_/g' | sed 's/Format: /Format:_/g' | sed 's/Number of sequences/Number_of_sequences/g' | sed -e "s/$/ $file/" >> ../original_stats_all.txt

done


mkdir -p alignments_rnaalifold
for file in alignments_G*;
do
outname=`basename $file .stk.stk`

mv $file ./alignments_rnaalifold/$outname.stk

done



for file in *.stk;
do 

ID=`echo $file | cut -d '.' -f1,2 | cut -d "_" -f1,2`
ID_2=`echo $file | cut -d '.' -f1,2 `

grep ^"#=GS" $file | sort | uniq | cut -d "/" -f1 | cut -d ' ' -f2 | sed -e "s/$/   $ID   $ID_2/" >> ../query_target_pairs_2.txt

done


for file in *.rnaalifold;   do   if [ $file == *"\.stk\.rnaalifold" ]; then ID=`basename $file .stk.rnaalifold`; else ID=`basename $file .stk.rnaalifold`; fi; MFE=`grep "=" $file | rev | cut -d "(" -f1 | rev | cut -d "=" -f1`; echo "$ID $MFE" >> ../positive_control.rnaalifold; done
