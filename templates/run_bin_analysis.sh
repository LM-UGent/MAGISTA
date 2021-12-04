#!/bin/bash
fa_list=$1
ref_list=$2
prefix=$3

#required scripts/software
<<<ENV>>>
run_checkm=checkm
run_checkm2csv=<<<MAGISTA_PATH>>>/checm2csv.sh
run_eval=<<<MAGISTA_PATH>>>/run_eval.sh

echo "Please make sure that your contigs do not contain special characters before running this Script"
echo "You can use '${nucops_path} renameseqs -i <file> -f {strictname} ' to achive this (requires nucops > v0.2.4) "

#rm  distribs_*.csv

example_bin=$(cat $fa_list | tail -1)
ext=$(echo "$example_bin" | sed "s/^.*\\.f/f/g")
bin_dir=$(realpath "$example_bin" | sed "s/[^\\/]*$//g")
if [[ -f binned_and_aligned_contigs.fa ]]; then
 rm binned_and_aligned_contigs.fa
fi

echo "Extension: $ext"
echo "Bin directory: $bin_dir"

if [[ -f ${prefix}_checkm.csv ]]; then
 echo CheckM has already been run
else
 mkdir "${bin_dir}/${prefix}_checkm_bins"
 while read line
 do
   cp "$line" "${bin_dir}/${prefix}_checkm_bins/"
 done < ${fa_list}
 $run_checkm lineage_wf -t 16 -x $ext "${bin_dir}/${prefix}_checkm_bins" "./${prefix}_checkm" > ${prefix}_checkm.txt
 rm -rf "${bin_dir}/${prefix}_checkm_bins"
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

echo "Run distribution generation"
$run_eval $fa_list
for f in distribs_*.csv
do
 tail -n +2 $f > $f.contents
 suffix=${f#*_}
 suffix=${suffix%.*}
 cat $f | head -1 | sed "s/,/,$suffix./g" | sed 's/\"//g'  > "$f.header"
 cat $f.header $f.contents > "$f"
 rm $f.header $f.contents
done

for f in distribs_*.csv
do
 mv "$f" "$f.tmp"
 cat "$f.tmp" | awk 'NR<2{print;next}{print| "sort -k3 -t ','"}' > "$f"
 rm "$f.tmp"
done

for f in distribs_*.csv
do
   mv "$f" "$f.tmp"
   cat "$f.tmp" | awk 'NR<2{print;next}{print| "sort -k3 -t ','"}' > "$f"
   rm "$f.tmp"
done
paste -d , ${prefix}_checkm.csv  distribs_*.csv | sed "s/Distance/1.0/g" > "${prefix}_summaries.csv"



