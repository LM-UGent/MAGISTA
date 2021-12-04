magista_dir=<<<MAGISTA_PATH>>>
cd full_test
isfirst=1
while read prefix
do
 $magista_dir/run_bin_analysis.sh ${prefix}.list ref.${prefix%.*}.list test_${prefix}
 if [ $isfirst -eq 1 ]
 then
   cat test_${prefix}_*.csv > test_data.csv
 else
   tail -n +2 test_${prefix}_*.csv >> test_data.csv
 fi
 rm distrib*
 isfirst=0
done < <(ls *.list | grep -v ref | sort | uniq | sed 's:.list::g')
cd ..
