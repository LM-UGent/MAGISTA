checkm_download_address=https://github.com/Ecogenomics/CheckM/archive/refs/tags/v1.1.2.zip
manager=$1

configure_tool(){
  default=$1
  alternate=$2
  if $(command -v $default &> /dev/null); then
   echo $default
  else
   echo $alternate
  fi
}
extraarg=""
if [ manager == conda ]; then
 export CONDA_ALWAYS_YES="true"
fi
cd dependencies
$manager install pysam
$manager install numpy
$manager install matplotlib
pip install checkm-genome
#wget -O checkm_dir.zip $checkm_download_address
#unzip checkm_dir.zip -d checkm_dir &> /dev/null
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzf checkm_data_2015_01_16.tar.gz
checkm data setRoot "$(pwd)"

if [ manager == conda ]; then
 unset CONDA_ALWAYS_YES
fi

cd ../

