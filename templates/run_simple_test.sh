_pwd=/home/ggoussar/test_samples/templates/MAGISTA/templates
cd <<<MAGISTA_PATH>>>/test
<<<MAGISTA_PATH>>>/run_bin_analysis.sh fa.list ref.list unvalidated_Bmock12_mc
<<<MAGISTA_PATH>>>/run_bin_analysis_with_validation.sh fa.list ref.list validated_Bmock12_mc;
cd "$_pwd"

