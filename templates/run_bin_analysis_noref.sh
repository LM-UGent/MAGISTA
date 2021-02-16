#!/bin/bash
fa_list=$1
prefix=$2

#required scripts/software
MAGISTA_dir=~/Workspace/public/MAGISTA/
run_checkm=checkm
run_checkm2csv=${MAGISTA_dir}checkm2csv.sh
run_eval=${MAGISTA_dir}run_eval.sh

#rm  distribs_*.csv

example_bin=$(cat $fa_list | tail -1)
ext=$(echo "$example_bin" | sed "s/^.*\\.f/f/g")
bin_dir=$(realpath "$example_bin" | sed "s/[^\\/]*$//g")
if [[ -f binned_and_aligned_contigs.fa ]]; then
 rm binned_and_aligned_contigs.fa
fi

if [[ -x "$(command -v ${run_checkm})" ]]; then
 has_checkm=yes
fi
if [[ has_checkm = yes ]]; then
 if [[ -f ${prefix}_checkm.csv ]]; then
  echo CheckM has already been run
 else
  $run_checkm lineage_wf -t 16 -x $ext "$bin_dir" "./${prefix}_checkm" > ${prefix}_checkm.txt
  $run_checkm2csv ${prefix}_checkm.txt > ${prefix}_checkm.tmp
  firstline=yes
  while read line
  do
   if [ $firstline == yes ]
   then
    firstline=no
    echo "$line"
   else
    echo "${prefix}_$line"
   fi
  done < ${prefix}_checkm.tmp > ${prefix}_checkm.csv
  rm ${prefix}_checkm.tmp
  rm ${prefix}_checkm.txt
 fi
else
 echo CheckM is not in the path and has not been run
fi

if [[ $(echo distribs_*.csv) =~ \* ]]; then
 echo "Run distribution generation"
 $run_eval $fa_list
 for f in distribs_*.csv
 do
  tail -n +2 $f > $f.contents
  suffix=${f#*_}
  suffix=${suffix%.*}
  cat $f | head -1 | sed "s/,/,$suffix./g" | sed 's/\"//g' > "$f.header"
  cat $f.header $f.contents > "$f"
  rm $f.header $f.contents
 done
else
  echo Evaluation already done - please remove all files beginning with distribs_ to redo this part of the analysis
fi

for f in distribs_*.csv
do
 mv "$f" "$f.tmp"
 cat "$f.tmp" | awk 'NR<2{print;next}{print| "sort -k1 -t ','"}' > "$f"
 rm "$f.tmp"
done

if [[ has_checkm = yes ]]; then
 paste -d , distribs_*.csv ${prefix}_checkm.csv | sed "s/Distance/1.0/g" > "${prefix}_MAGISTIC_input.csv"
else
 paste -d , distribs_*.csv | sed "s/Distance/1.0/g" > "${prefix}_MAGISTA_input.csv"
fi

