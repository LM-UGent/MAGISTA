#!/bin/sh
manager=$1

curpath=$(pwd)

configure_tool(){
  default=$1
  alternate=$2
  if $(command -v $default &> /dev/null); then
   echo $default
  else
   echo $alternate
  fi
}

requires_manual_path_installation(){
  tool=$1
  executable=$2
  website=$3
  comment=$4
  if ! command -v $executable &> /dev/null; then
   echo "[ERR] $tool must be installed AND be available in your PATH. It is available at $website (NOTE: $comment)"
   return 1
  else
   echo "      $tool seems to be installed - version was not checked though"
   return 0
  fi
}

echo ""
echo "***"
echo "CHECKING DESIRED MANAGER"
echo "***"

manager_ok=0
active_env=
activate_env=
if [ .$manager == ".pip" ]; then manager_ok=1; fi
if [ .$manager == ".conda" ]; then
 active_env=$(conda info | grep "active environment" | sed "s/[^:]*:.//g")
 activate_env="conda activate $active_env"
 manager_ok=1;
fi
if [ .$manager == ".conda.source" ]; then
 active_env=$(conda info | grep "active environment" | sed "s/[^:]*:.//g")
 manager=conda
 activate_env="source activate $active_env"
 manager_ok=1;
fi

if [ $manager_ok -eq 0 ]; then
 echo "Please specify the package manager you want to use: $0 [pip|conda|conda.source]"
 echo "  (This is a reminder to activate your conda environment if you use one)"
 exit 1
fi



dependencies_met=1

echo "***"
echo "CHECKING DEPENDENCIES"
echo "***"
if [ .$manager == ".conda" ]; then
 requires_manual_path_installation "Conda" conda https://docs.conda.io/en/latest/miniconda.html ""
 if [ $? -eq 1 ]; then
  dependencies_met=0
  manager_ok=0
 fi
fi
requires_manual_path_installation python python http://python.org/ "Python3"
if [ $? -eq 1 ]; then
 if [ .$manager == ".conda" ]; then
  if [manager_ok -eq 0 ]; then
   dependencies_met=0
   manager_ok=0
  else
   $manager install python=3
  fi
 fi
fi
requires_manual_path_installation R Rscript https://www.r-project.org  "tested with v4.0.5"
if [ $? -eq 1 ]; then
 if [ .$manager == ".conda" ]; then
  if [manager_ok -eq 0 ]; then 
   dependencies_met=0
  else
   $manager instal R=4
  fi
 fi
fi
requires_manual_path_installation hmmer hmmsearch http://hmmer.org/ ">=v3.1b1"
if [ $? -eq 1 ]; then
 dependencies_met=0
fi
requires_manual_path_installation prodigal prodigal https://github.com/hyattpd/Prodigal "2.60 or >=2.6.1"
if [ $? -eq 1 ]; then
 dependencies_met=0
fi
requires_manual_path_installation pplacer pplacer http://matsen.fhcrc.org/pplacer/ ">=1.1"
if [ $? -eq 1 ]; then
 dependencies_met=0
fi


if [ $dependencies_met -eq 0 ]; then
 echo "Some dependencies must be installed manually - we recommend that you set up a dedicated conda environment for this"
 exit 1
else
 echo "All required dependencies were met"
fi

echo ""
echo "***"
echo "INSTALLING ADDITIONAL DEPENDENCIES"
echo "***"


checkm_path=$(configure_tool checkm ${curpath}/dependencies/checkm)
if [ $checkm_path == "${curpath}/dependencies/checkm" ]; then
  echo "Configuring checkm and dependencies"
  ${curpath}/dependencies/install_checkm.sh $manager
fi
${curpath}/dependencies/install_checkm.sh $manager

printf "install.packages(c('tidyverse','patchwork','reshape2','extrafont','randomForest','randomForestExplainer','caret','moments'),repos='http://cran.us.r-project.org')" | R --no-save

nucops_path=$(configure_tool nucops "${curpath}/dependencies/nucops")
if [ "$nucops_path" == "${curpath}/dependencies/nucops" ]; then
 echo "   Compiling nucops"
 cd "${curpath}/dependencies/nucops_dir"
 make
 cp bin/nucops "${curpath}/dependencies/nucops"
 cd "${curpath}"
fi
GDC_path=$(configure_tool GenDisCal ${curpath}/dependencies/GenDisCal)
if [ "$nucops_path" == "${curpath}/dependencies/GenDisCal" ]; then
 echo "   Compiling GenDisCal"
 cd "${curpath}/dependencies/GenDisCal_dir"
 make
 cp bin/nucops "${curpath}/dependencies/GenDisCal"
 cd "${curpath}"
fi
genosig_path=$(configure_tool genosig ${curpath}/dependencies/genosig)
if [ "$nucops_path" == "${curpath}/dependencies/genosig" ]; then
 echo "   Compiling genosig"
 cd "${curpath}/dependencies/genosig_dir"
 make
 cp bin/genosig "${curpath}/dependencies/genosig"
 cd "${curpath}"
fi


echo ""
echo "***"
echo "UPDATING PATHS"
echo "***"


for f in templates/*.*
do
  cp $f ${curpath}/$(basename $f)
  echo $f
  cat $f | sed "s:<<<MAGISTA_PATH>>>:$curpath:g" | sed "s:<<<GDC_PATH>>>:$GDC_path:g" | sed "s:<<<NUCOPS_PATH>>>:$nucops_path:g" | sed "s:<<<GENOSIG_PATH>>>:$genosig_path:g" | sed "s:<<<ENV>>>:$activate_env:g" > ${curpath}/$(basename $f)
done

for d in full_test full_training
do
  for f in ${d}/*.list
  do
    echo $f
    while read line
    do
      realpath "${curpath}/${d}/${line}"
    done < $f > "${d}/tmp.tmp"
    mv "${d}/tmp.tmp" $f
  done
done

echo ""
echo "***"
echo "INSTALLATION COMPLETE"
echo "***"

