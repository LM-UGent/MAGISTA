# MAGISTA
## Dependencies
Required:
 - python (you will likely have trouble if you do not use 3.6 to run both checkM and AMBER)
 - R
 - CheckM and its dependencies
Strongly recommended:
 - conda
## Installation
Extract the source code in the location of your choice.
If you have conda installed, we strongly recommend creating a virtual environment to use MAGISTA
```
conda create -n MAGISTA_env
conda activate MAGISTA_env
```
Next, run the following commands:
```
./install_dependencies.sh
Rscript ./MAGISTA.R FIRSTRUN
```


Installation is currently not completely tested - please report any trouble that you have with installation to g_goussarov@hotmail.com.