#1/bin/bash

srna=${1}
outname=${2}
current_folder=${3}


current_values=`echo "${current_folder}/${outname}_example_files/${outname}_read_values.txt"`
# echo $current_values
# > ${current_values}
while read line;
do
# echo $line
genome=`echo $line | cut -d ',' -f8`
srna=`echo $line | cut -d ',' -f10`
contig=`echo $line | cut -d ',' -f1`
genus=`echo $line | cut -d ',' -f9`
start=`echo $line | cut -d ',' -f3`
stop=`echo $line | cut -d ',' -f4`

if [[ $start  == "start" ]];then
continue
fi

tmpstart=$(($start - 5000))
tmpstop=$(($start + 5000))

start=$tmpstart
stop=$tmpstop

esl-seqstat -a ~/phd/RNASeq/sequences/${genome}.fna | grep "=" > ./tmp/tmp.stats
    current_count=0
while read stat_line;
do
chromosome=`echo "$stat_line" | grep $contig | wc -l | cut -d ' ' -f8`
if (( $chromosome > 0 ));
then
tmpstart=$(($start + $current_count))
tmpstop=$(($stop + $current_count))
if (( $tmpstart < $tmpstop  )); then
start=$tmpstart
stop=$tmpstop
else
stop=$tmpstart
start=$tmpstop
fi

if [[ ! -d ~/phd/RNASeq/genera/${genus}/${genome}.data/plot_files ]];
then
echo "Genome: ${genome} has no plot files"
continue
fi
#echo "~/phd/RNASeq/genera/${genus}/${genome}.data/plot_files $start $stop"
> ./tmp/test.plot
counter="1"
for plotfile in ~/phd/RNASeq/genera/${genus}/${genome}.data/plot_files/*.plot;
do
if [[ $plotfile == *"ncRNA"* ]]; then
continue
elif [[ $plotfile == *"rev"* ]]; then
continue
elif [[ $plotfile == *"fwd"* ]]; then
continue            
fi
# echo "$plotfile $genus $genome"
sed -n "${start},${stop}{p;${stop}q;}"  $plotfile > ./tmp/tmp.plot
if [[  $counter ==  "1" ]];
then
cat ./tmp/tmp.plot > ./tmp/test.plot
counter="2"
else
cat ./tmp/test.plot > ./tmp/tmp2.plot 
paste ./tmp/tmp2.plot ./tmp/tmp.plot > ./tmp/test.plot
fi
done
  #cp test.plot .
#summarise_sRNA_read_depths.R ./test.plot
#cat ./test.values | sed -e "s/$/ $ID $genus $genome_name $genome $contig $realstart $realstop $start $stop/" >> ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/${group}_read_depths_summary.txt
sRNA_examples_summary.R "${genome}_${start}" "${current_folder}"

cat ${current_values} > ./tmp/tmp2.values 
paste ./tmp/tmp2.values ./tmp/tmp.values > ${current_values}

else
length=`echo $stat_line | cut -d ' ' -f3`
current_count=$(($current_count + $length))
fi

done < ./tmp/tmp.stats

done < ${current_folder}/${outname}_example_files/${outname}_list.csv