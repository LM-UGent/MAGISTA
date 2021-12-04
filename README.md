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

# MAGISTA

The MAGISTA software is intended to give de-novo evaluation of MAG (metagenome-assembled genome) quality. It estimates completeness and purity of MAGs using a random forest-based approach.

## Usage
MAGISTA consists of two separate subprograms: MAGISTA.R and run_bin_analysis.sh, both of which should be run from the command line.
### run_bin_analysis.sh
Usage: `run_bin_analysis.sh <file-list file> <prefix>`  
The `<file-list file>` should be a file with the absolute path of each bin that needs to be analysed, with one path per line  
`<prefix>` is a name that is prepended to the output.  
This script produces a csv file which can be parsed by MAGISTA.R  
If checkM is not in your PATH, the output will not be usable with the MAGISTIC argument
### MAGISTA.R
Usage: `Rscript MAGISTA.R <MAGISTA|MAGISTIC> <input file> [<output file>] [<model name>]`  
The first argument should be either `MAGISTA` or `MAGISTIC` depending on whether the checkM output should be included or not.  
`<input file>` should be the path to the file produced by run_bin_analysis_noref.sh  
`<output file>` should be the path to the output file.  
`<model name>` should be the name of the model, if you do not wish to use teh default models or create a new one.