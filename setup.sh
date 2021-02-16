#!/bin/bash
_pwd="$(pwd)/"
cd utils/genosig
make
cd ../nucops2
make
cd "$_pwd"
cat templates/MAGISTA.R | sed s:MAGISTA_dir\=.*:MAGISTA_dir=\"${_pwd}\":g > MAGISTA.R
cat templates/run_bin_analysis_noref.sh | sed s:MAGISTA_dir=.*:MAGISTA_dir\=\"${_pwd}\":g > run_bin_analysis_noref.sh
cat templates/run_eval.sh | sed s:MAGISTA_dir\=.*:MAGISTA_dir=\"${_pwd}\":g > run_eval.sh
cp templates/checkm2csv.sh checkm2csv.sh
cp templates/distribution_summary.R distribution_summary.R
chmod a+x run_bin_analysis_noref.sh
chmod a+x run_eval.sh
chmod a+x checkm2csv.sh
Rscript MAGISTA.R FIRSTRUN

if [ $? -eq 0 ]; then
    echo "---"
    echo "Setup complete"
    echo "Please make sure that CheckM and GenDisCal are installed on your machine and located in your PATH when you run the bin analysis."
    if ! [ -x "$(command -v GenDisCal)" ]; then
        echo 'Note: GenDisCal is not in the current path'
    fi
    if ! [ -x "$(command -v checkm)" ]; then
        echo 'Note: checkm is not in the current path'
    fi
else
    echo "Setup Failed - please check that R is installed on your system and that the following libraries are available: "
    echo " randomForest"
    echo " dplyr"
    echo " moments"
fi

