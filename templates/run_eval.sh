#!/bin/bash
# Note: the first argument should be the name of a file containing all the bins to evaluate 
<<<ENV>>>
nucops_path=<<<NUCOPS_PATH>>>
genosig_path=<<<GENOSIG_PATH>>>
GDC_path=<<<GDC_PATH>>>

if [ -f tmp.list ]
then
 rm tmp.list
fi

if [ ! -f distribs_w1k_GC.csv ]
then
 echo "GC 1K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 1000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w1k.fasta -w 1000 -s 1000 -q all
   if [ $bintot -gt 30000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w1k.fasta -k 30M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w1k.fasta
   fi
   ${genosig_path} generate contigs -i ${filename%.fa}_w1k.fasta -t freq 1 --readable -m 500 | cut -d , -f 3 > ${filename%.fa}_w1k_points.csv
   echo ${filename%.fa}_w1k_points.csv
   rm ${filename%.fa}_w1k.fasta
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 cp tmp.csv tmp.GC.csv
 cp leads.csv leads.GC.csv
 paste -d , leads.csv tmp.csv > distribs_w1k_GC.csv
 rm tmp.csv
fi

if [ ! -f  distribs_w50k_k4s.csv ]
then
 echo "K4S 50K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 50000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w50k.fasta -w 50000 -s 50000 -q all
   if [ $bintot -gt 300000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w50k.fasta -k 300M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w50k.fasta
   fi
   ${genosig_path} generate contigs -i ${filename%.fa}_w50k.fasta -t karl 4 -o ${filename%.fa}_w50k.sdb
   ${GDC_path} ${filename%.fa}_w50k.sdb -b karl 4 -m SVC 0.05 -q all | cut -d , -f 4 > ${filename%.fa}_w50k_cmp.csv
   ${GDC_path} ${filename%.fa}_w50k.sdb -b karl 4 -m SVC 0.05 -q all -g 0.01 -o PNG/${filename%.fa}_w50k_g.csv
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w50k_s5k_g.csv
   echo ${filename%.fa}_w50k_cmp.csv
   rm ${filename%.fa}_w50k.fasta ${filename%.fa}_w50k.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w50k_k4s.csv 
 rm tmp.csv
fi

if [ ! -f distribs_w100k_k4s.csv ]
then
 echo "K4S 100K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 100000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w100k.fasta -w 100000 -s 100000 -q all
   if [ $bintot -gt 300000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w100k.fasta -k 300M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w100k.fasta
   fi

   ${genosig_path} generate contigs -i ${filename%.fa}_w100k.fasta -t karl 4 -o ${filename%.fa}_w100k.sdb
   ${GDC_path} ${filename%.fa}_w100k.sdb -b karl 4 -m SVC 0.05 -q all | cut -d , -f 4 > ${filename%.fa}_w100k_cmp.csv
   ${GDC_path} ${filename%.fa}_w100k.sdb -b karl 4 -m SVC 0.05 -g 0.01 -o PNG/${filename%.fa}_w100k_g.csv -q all
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w100k_s5k_g.csv
   echo ${filename%.fa}_w100k_cmp.csv
   rm ${filename%.fa}_w100k.fasta ${filename%.fa}_w100k.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w100k_k4s.csv
 rm tmp.csv
fi


if [ ! -f distribs_w5k_f4c.csv ]
then
 echo "F4C 5K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 5000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w5k.fasta -w 5000 -s 5000 -q all
   if [ $bintot -gt 150000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w5k.fasta -k -q all 150M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w5k.fasta
   fi

   ${genosig_path} generate contigs -i ${filename%.fa}_w5k.fasta -t freq 4 -o ${filename%.fa}_w5k.sdb
   ${GDC_path} ${filename%.fa}_w5k.sdb -b freq 4 -m corr -q all| cut -d , -f 4 > ${filename%.fa}_w5k_cmp.csv
   ${GDC_path} ${filename%.fa}_w5k.sdb -b freq 4 -m corr -q all -g 0.001 -o PNG/${filename%.fa}_w5k_g.csv
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w5k_s1k_g.csv
   echo ${filename%.fa}_w5k_cmp.csv
   rm ${filename%.fa}_w5k.fasta ${filename%.fa}_w5k.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w5k_f4c.csv
 rm tmp.csv
fi

if [ ! -f distribs_w10k_m3c.csv ]
then
 echo "M3C 10K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 10000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w10k.fasta -w 10000 -s 10000 -q all
   if [ $bintot -gt 300000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w10k.fasta -k 300M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w10k.fasta
   fi

   ${genosig_path} generate contigs -i ${filename%.fa}_w10k.fasta -t mmz 3 -o ${filename%.fa}_w10k_m3.sdb
   ${GDC_path} ${filename%.fa}_w10k_m3.sdb -b exp 3 -q all -m corr | cut -d , -f 4 > ${filename%.fa}_w10k_m3c_cmp.csv
   ${GDC_path} ${filename%.fa}_w10k_m3.sdb -b exp 3 -q all -m corr -g 0.001 -o PNG/${filename%.fa}_w10k_m3c_g.csv
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w10k_s5k_m3c_g.csv
   echo ${filename%.fa}_w10k_m3c_cmp.csv
   rm ${filename%.fa}_w10k.fasta ${filename%.fa}_w10k_m3.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w10k_m3c.csv
 rm tmp.csv
fi

if [ ! -f distribs_w10k_m4c.csv ]
then
 echo "M4C 10K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 10000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w10k.fasta -w 10000 -s 10000 -q all
   if [ $bintot -gt 300000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w10k.fasta -k 300M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w10k.fasta
   fi

   ${genosig_path} generate contigs -i ${filename%.fa}_w10k.fasta -t mmz 4 -o ${filename%.fa}_w10k_m4.sdb
   ${GDC_path} ${filename%.fa}_w10k_m4.sdb -b exp 4 -m corr -q all | cut -d , -f 4 > ${filename%.fa}_w10k_m4c_cmp.csv
   ${GDC_path} ${filename%.fa}_w10k_m4.sdb -b exp 4 -m corr -q all -g 0.001 -o PNG/${filename%.fa}_w10k_m4c_g.csv
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w10k_s5k_m4c_g.csv
   echo ${filename%.fa}_w10k_m4c_cmp.csv
   rm ${filename%.fa}_w10k.fasta ${filename%.fa}_w10k_m4.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w10k_m4c.csv
 rm tmp.csv
fi
 
 
if [ ! -f distribs_w10k_m5c.csv ]
then
 echo "M5C 10K"
 echo "bin.size,bin.part" > leads.csv
 while read filename
 do
   ${nucops_path}2 "select" -i $filename -o ${filename%.fa}_selected_contigs.fasta -l 10000 -q all
   binpart=$(${nucops_path} fastasummary -q all --validnullinput ${filename%.fa}_selected_contigs.fasta -t)
   bintot=$(${nucops_path} fastasummary -q all --validnullinput ${filename} -t)
   rm ${filename%.fa}_selected_contigs.fasta
   echo "$bintot,$binpart" >> leads.csv

   ${nucops_path} winsplit -i $filename -o ${filename%.fa}_w10k.fasta -w 10000 -s 10000 -q all
   if [ $bintot -gt 300000000 ]
   then
     ${nucops_path} select -l 999 -q all -i ${filename%.fa}_w10k.fasta -k 300M -o ${filename%.fa}.fasta.tmp
     mv ${filename%.fa}.fasta.tmp ${filename%.fa}_w10k.fasta
   fi

   ${genosig_path} generate contigs -i ${filename%.fa}_w10k.fasta -t mmz 5 -o ${filename%.fa}_w10k_m5.sdb
   ${GDC_path} ${filename%.fa}_w10k_m5.sdb -b exp 5 -m corr -q all | cut -d , -f 4 > ${filename%.fa}_w10k_m5c_cmp.csv
   ${GDC_path} ${filename%.fa}_w10k_m5.sdb -b exp 5 -m corr -q all -g 0.001 -o PNG/${filename%.fa}_w10k_m5c_g.csv
   #Rscript ~/Scripts/XY_view.R PNG/${filename%.fa}_w10k_s5k_m5c_g.csv
   echo ${filename%.fa}_w10k_m5c_cmp.csv
   rm ${filename%.fa}_w10k.fasta ${filename%.fa}_w10k_m5.sdb
 done < $1 > tmp.list
 Rscript ~/Scripts/distribution_summary.R tmp.list tmp.csv
 while read f
 do
   rm $f
 done < tmp.list
 rm tmp.list
 paste -d , leads.csv tmp.csv > distribs_w10k_m5c.csv
 rm tmp.csv
 rm leads.csv
fi



