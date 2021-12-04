magista_dir=<<<MAGISTA_PATH>>>
cd full_training
isfirst=1
while read prefix
do
 $magista_dir/run_bin_analysis_with_validation.sh ${prefix}.list ref.${prefix%.*}.list validated_training_${prefix}
 if [ $isfirst -eq 1 ]
 then
   cat validated_training_${prefix}_amber_and_distrib_summaries_.csv> validated_training_data.csv
 else
   tail -n +2 validated_training_${prefix}_amber_and_distrib_summaries_.csv >> validated_training_data.csv
 fi
 rm distrib*
 isfirst=0
done < <(ls *.list | grep -v ref | sort | uniq | sed 's:.list::g')
cd ..
