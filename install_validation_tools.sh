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
if [ .$manager == ".pip" ]; then manager_ok=1; fi
if [ .$manager == ".conda" ]; then manager_ok=1; fi

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

requires_manual_path_installation QUAST metaquast.py "https://cab.spbu.ru/software/quast/"
if [ $? -eq 1 ]; then
  dependencies_met=0
fi

if [ $dependencies_met -eq 0 ]; then
 echo "Some dependencies must be installed manually - we recommend that you set up a dedicated conda environment for this"
 exit 1
else
 echo "All required dependencies were met"
fi

if [ $manager_ok -eq 0 ]; then
 echo "Please specify the package manager you want to use: $0 [pip|conda]"
 echo "  (This is a reminder to activate your conda environment if you use one)"
 exit 1
fi

echo ""
echo "***"
echo "INSTALLING ADDITIONAL DEPENDENCIES"
echo "***"

$manager install cami-amber

echo ""
echo "***"
echo "INSTALLATION COMPLETE"
echo "***"

