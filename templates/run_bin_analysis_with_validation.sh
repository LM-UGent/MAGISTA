#!/bin/bash
fa_list=$1
ref_list=$2
prefix=$3

#required scripts/software

<<<ENV>>>

run_checkm=checkm
run_checkm2csv=<<<MAGISTA_PATH>>>/checm2csv.sh
run_AMBERevalsum=<<<MAGISTA_PATH>>>/AMBER_eval_summary.sh
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
 source activate checkm_env
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
 source deactivate
fi

tempfile=binned_and_aligned_contigs.fa
while read line
do
 cat $line
done < "$fa_list" > $tempfile
realpath "$tempfile"
if [[ -f "${prefix}_amber_summary.csv" ]]; then
 echo amber_summary already exists
else
# echo "Re-running AMBER analysis even if it exists - edit $0 to change"
 $run_AMBERevalsum $tempfile $fa_list $ref_list $prefix
fi

cat "${prefix}_amber_summary.csv" | awk 'NR<2{print;next}{print| "sort -k1 -t ','"}' > "${prefix}_amber_summary.tmp"
mv "${prefix}_amber_summary.tmp" "${prefix}_amber_summary.csv"
rm $tempfile
#Remove cocmments once this is all done
#if [[ $(echo distribs_*.csv) =~ \* ]]; then
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
#else
#  "Evaluation alreaady done" check moved to the file itself
#  echo Evaluation already done
#fi


for f in distribs_*.csv
do
 mv "$f" "$f.tmp"
 cat "$f.tmp" | awk 'NR<2{print;next}{print| "sort -k3 -t ','"}' > "$f"
 rm "$f.tmp"
done

#for var in 99 # 95 90 80
#do
# paste -d , "blasteval_summary_${var}.csv" ${prefix}_checkm.csv  distribs_*.csv | sed "s/Distance/1.0/g" > "${prefix}_blasteval_and_distrib_summaries_${var}.csv"
#done

paste -d , "${prefix}_amber_summary.csv" ${prefix}_checkm.csv  distribs_*.csv | sed "s/Distance/1.0/g" > "${prefix}_amber_and_distrib_summaries_.csv"
rm distribs_*.csv


