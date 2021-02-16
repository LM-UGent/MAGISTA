#!/usr/bin/bash
# Note: the first argument should be the name of a file containing all the bins to evaluate 

MAGISTA_dir=~/Workspace/public/MAGISTA/

run_nucops2=${MAGISTA_dir}utils/nucops2/bin/nucops2
run_genosig=${MAGISTA_dir}utils/genosig/bin/genosig
run_distribution_summary=${MAGISTA_dir}distribution_summary.R
run_GenDisCal=GenDisCal

if [ -f tmp.list ]
then
 rm tmp.list
fi

echo "GC 1K"
echo "Bin.Id,bin.size,bin.part" > leads.csv
while read filename
do
   $run_nucops2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 1000 -q progress
   binpart=$($run_nucops2 fastasummary -q progress ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$($run_nucops2 fastasummary -q progress ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$filename,$bintot,$binpart" >> leads.csv

   $run_nucops2 winsplit -i $filename -o ${filename%.fa}_w1k_s500.fasta -w 1000 -s 500 -q progress
   $run_genosig generate contigs -i ${filename%.fa}_w1k_s500.fasta -t freq 1 --readable -m 500 | cut -d , -f 3 > ${filename%.fa}_w1k_s500_points.csv
   echo ${filename%.fa}_w1k_s500_points.csv
   rm ${filename%.fa}_w1k_s500.fasta
done < $1 > tmp.list
Rscript $run_distribution_summary tmp.list tmp.csv
while read f
do
   rm $f
done < tmp.list
rm tmp.list
cp tmp.csv tmp.GC.csv
cp leads.csv leads.GC.csv
paste -d , leads.csv tmp.csv > distribs_w1k_GC.csv
rm tmp.csv


echo "K4S 50K"
echo "Bin.Id,bin.size,bin.part" > leads.csv
while read filename
do
   $run_nucops2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 50000 -q progress
   binpart=$($run_nucops2 fastasummary -q progress ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$($run_nucops2 fastasummary -q progress ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$filename,$bintot,$binpart" >> leads.csv

   $run_nucops2 winsplit -i $filename -o ${filename%.fa}_w50k_s5k.fasta -w 50000 -s 5000 -q progress
   $run_genosig generate contigs -i ${filename%.fa}_w50k_s5k.fasta -t karl 4 -o ${filename%.fa}_w50k_s5k.sdb
   $run_GenDisCal ${filename%.fa}_w50k_s5k.sdb -b karl 4 -m SVC 0.05 -q all | cut -d , -f 4 > ${filename%.fa}_w50k_s5k_cmp.csv
   echo ${filename%.fa}_w50k_s5k_cmp.csv
   rm ${filename%.fa}_w50k_s5k.fasta ${filename%.fa}_w50k_s5k.sdb
done < $1 > tmp.list
Rscript $run_distribution_summary tmp.list tmp.csv
while read f
do
   rm $f
done < tmp.list
rm tmp.list
paste -d , leads.csv tmp.csv > distribs_w50k_k4s.csv 
rm tmp.csv

echo "F4C 5K"
echo "Bin.Id,bin.size,bin.part" > leads.csv
while read filename
do
   $run_nucops2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 5000 -q progress
   binpart=$($run_nucops2 fastasummary -q progress ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$($run_nucops2 fastasummary -q progress ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$filename,$bintot,$binpart" >> leads.csv

   $run_nucops2 winsplit -i $filename -o ${filename%.fa}_w5k_s1k.fasta -w 5000 -s 1000 -q progress
   $run_genosig generate contigs -i ${filename%.fa}_w5k_s1k.fasta -t freq 4 -o ${filename%.fa}_w5k_s1k.sdb
   $run_GenDisCal ${filename%.fa}_w5k_s1k.sdb -b freq 4 -m corr -q all| cut -d , -f 4 > ${filename%.fa}_w5k_s1k_cmp.csv
   echo ${filename%.fa}_w5k_s1k_cmp.csv
   rm ${filename%.fa}_w5k_s1k.fasta ${filename%.fa}_w5k_s1k.sdb
done < $1 > tmp.list
Rscript $run_distribution_summary tmp.list tmp.csv
while read f
do
   rm $f
done < tmp.list
rm tmp.list
paste -d , leads.csv tmp.csv > distribs_w5k_f4c.csv
rm tmp.csv

echo "M3C 10K"
echo "Bin.Id,bin.size,bin.part" > leads.csv
while read filename
do
   $run_nucops2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 10000 -q progress
   binpart=$($run_nucops2 fastasummary -q progress ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$($run_nucops2 fastasummary -q progress ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$filename,$bintot,$binpart" >> leads.csv

   $run_nucops2 winsplit -i $filename -o ${filename%.fa}_w10k_s5k.fasta -w 10000 -s 5000 -q progress
   $run_genosig generate contigs -i ${filename%.fa}_w10k_s5k.fasta -t mmz 3 -o ${filename%.fa}_w10k_s5k_m3.sdb
   $run_GenDisCal ${filename%.fa}_w10k_s5k_m3.sdb -b exp 3 -q all -m corr | cut -d , -f 4 > ${filename%.fa}_w10k_s5k_m3c_cmp.csv
   echo ${filename%.fa}_w10k_s5k_m3c_cmp.csv
   rm ${filename%.fa}_w10k_s5k.fasta ${filename%.fa}_w10k_s5k_m3.sdb
done < $1 > tmp.list
Rscript $run_distribution_summary tmp.list tmp.csv
while read f
do
   rm $f
done < tmp.list
rm tmp.list
paste -d , leads.csv tmp.csv > distribs_w10k_m3c.csv
rm tmp.csv

echo "M4C 10K"
echo "Bin.Id,bin.size,bin.part" > leads.csv
while read filename
do
   $run_nucops2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 10000 -q progress
   binpart=$($run_nucops2 fastasummary -q progress ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$($run_nucops2 fastasummary -q progress ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$filename,$bintot,$binpart" >> leads.csv

   $run_nucops2 winsplit -i $filename -o ${filename%.fa}_w10k_s5k.fasta -w 10000 -s 5000 -q progress
   $run_genosig generate contigs -i ${filename%.fa}_w10k_s5k.fasta -t mmz 4 -o ${filename%.fa}_w10k_s5k_m4.sdb
   $run_GenDisCal ${filename%.fa}_w10k_s5k_m4.sdb -b exp 4 -m corr -q all | cut -d , -f 4 > ${filename%.fa}_w10k_s5k_m4c_cmp.csv
   echo ${filename%.fa}_w10k_s5k_m4c_cmp.csv
   rm ${filename%.fa}_w10k_s5k.fasta ${filename%.fa}_w10k_s5k_m4.sdb
done < $1 > tmp.list
Rscript $run_distribution_summary tmp.list tmp.csv
while read f
do
   rm $f
done < tmp.list
rm tmp.list
paste -d , leads.csv tmp.csv > distribs_w10k_m4c.csv
rm tmp.csv


