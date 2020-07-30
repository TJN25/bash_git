#!/bin/bash
#This folder should be the directory and have the sequence of interest in it in the form of name.faa
#This calls other scripts, makes some simple adjustments to the files and sorts the final results into a folder.
#See each script for details on what they do or look at the ReadMe.txt file.
#./Virus-Host_Matcher_1.31.sh -i test.faa -d _ -s 0 -e 3 -t y -r n -c 1


keep_files='n'
sig_hmms_only='y'
output_file="n"


while getopts "i:t::" arg; do
  case $arg in
    i)
      file1_input=$OPTARG
      #echo $file1
      ;; 
    t)
      keep_files=$OPTARG
      #echo $keep_files
      ;;       
	      
  esac
done





fasta_file=fastafile.txt
echo $fasta_file



arVOG_hmm=arVOG-all.hmm



if [[ -z $file1 ]]; then
echo 'Error: Fasta File Needed'

elif [ $file_input = '-h' ]; then
echo 'Virus-Host_Matcher.sh: compares proteins from a number of contigs to hmm models and then scores'
echo 'each contig based on similarity to those models.'
echo 'Version 1.0 2016'
echo 'usage: <filename> <delimiter> <start_pos> <end_pos>'
echo '<filename> should be in a protein fasta format and include information about the contig or species'
echo '<delim> is the character that separates the contig name from the protein name/detail (if this is | then type __)'
echo '<start_pos> is the number of delimiters that you must pass to be left with just the contig name.'
echo '<end_pos> is the number of delimiters from the end needed to be left with just the contig'
echo 'For MG-RAST files the format would be <filename> _ 0 3 5 n y' 
echo 'which would keep all contigs with 5 or more proteins and remove the temp files once completed and would remove' 
echo 'some of the less significant models'

else


file1=$file_input
output_file=`basename $file1`
 cat $file1 > fastafile.txt #needed for running other scripts
chmod +x format_fastafile_to_remove_small_contigs.py
./format_fastafile_to_remove_small_contigs.py
 chmod +x remove_small_contigs.py 
 ./remove_small_contigs.py $delim $pos_start $pos_end $contig_cutoff
 echo 'Match '$file1' to HMM virus models'
 echo $arVOG_hmm
 echo 'Running HMMSearch on '$file1' and arVOG models' 
# each of the hmmsearch commands compare the input file '$file1' to one of the HMM model datasets
# the renaming is done to enable the files to be mre easily used with subsequent scripts
 hmmsearch --tblout arVOG_file.txt --noali $arVOG_hmm $fasta_file > $output_file.output_ar.txt 

echo 'Writing results into file'
chmod +x find_and_replace_space_with_tab.py
# the HMMsearch output includes a number of spaces and to work with the file more easily these are changed to one tab
./find_and_replace_space_with_tab.py
# only the protein name, e value and domain are kept from the output
grep 'arVOG_' arVOG_file.tab | cut -f1-3 > arVOG_list.txt
chmod +x generate_output_file.py
# a list of the protein, the lowest e value and the domain this e value is from are listed.
./generate_output_file.py
sort output_tmp.txt | uniq > output.txt
chmod +x significant_matches_output.py
# a list that includes only those proteins that a match was found for
./significant_matches_output.py



echo 'Calculating Scores'

chmod +x format_for_scoring.py
# This step assigns a minimun e value to the e values that are below the threshold in order to ensure
#that the score calculation is not weighhted to much based on the very small e values.
./format_for_scoring.py
# Calculates a score indicating which domain the genome likely belongs to
Rscript score_matches_for_contigs.R

#chmod +x format_output_scores.py
#./format_output_scores.py

keep_files_=$keep_files

fi

# these steps are needed to rename the files and put them in a folder

if [[ $keep_files_ = 'y' ]]; then

mkdir $output_file.temp.folder
mv output.txt output_tmp.txt output_sig_hits_only.txt temp_file.txt scores_calc.txt output_scores.txt proteins_and_contigs_names.txt fastafile.txt output_with_species.txt $output_file.temp.folder
mv arVOG_file.tab $output_file.temp.folder
mv arVOG_file.txt $output_file.temp.folder
mv arVOG_list.txt $output_file.temp.folder
mv proteinnames_list.txt fastafile_subset.txt $output_file.temp.folder
mv $output_file.output_ar.txt $output_file.temp.folder
command cat $output_file.temp.folder/output.txt > $output_file.txt
command cat $output_file.temp.folder/output_sig_hits_only.txt > $output_file.sig_hits_only.txt
command cat $output_file.temp.folder/output_with_species.txt > $output_file.output_with_species.txt
command cat $output_file.temp.folder/output_scores.txt > $output_file.scores.txt
command mkdir results.$output_file.folder
command mv $output_file.* results.$output_file.folder
echo 'Temp files being moved to temp folder'
echo 'Finished'

elif [[ $keep_files_ = 'n' ]]; then


mkdir $output_file.temp.folder
mv output.txt output_tmp.txt output_sig_hits_only.txt temp_file.txt scores_calc.txt output_scores.txt proteins_and_contigs_names.txt fastafile.txt output_with_species.txt $output_file.temp.folder
mv arVOG_file.tab $output_file.temp.folder
mv arVOG_file.txt $output_file.temp.folder
mv arVOG_list.txt $output_file.temp.folder
mv proteinnames_list.txt fastafile_subset.txt $output_file.temp.folder
mv $output_file.output_ar.txt $output_file.temp.folder
command cat $output_file.temp.folder/output.txt > $file1.txt
command cat $output_file.temp.folder/output_sig_hits_only.txt > $output_file.sig_hits_only.txt
command cat $output_file.temp.folder/output_with_species.txt > $output_file.output_with_species.txt
command cat $output_file.temp.folder/output_scores.txt > $output_file.scores.txt
command mkdir results.$output_file.folder
command mv $output_file.* results.$output_file.folder
rm results.$output_file.folder/$output_file.temp.folder/*
rmdir results.$output_file.folder/$output_file.temp.folder/

echo 'Temp files being moved to temp folder'

done

echo 'Finished'

else

echo "Did not complete the commands: If you wish to debug then include 'y' in the command."

fi