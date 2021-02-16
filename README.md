# MAGISTA

The MAGISTA software is intended to give de-novo evaluation of MAG (metagenome-assembled genome) quality. It estimates completeness and purity of MAGs using a random forest-based approach.

## Installation
MAGISTA has only been tested on Linux. 

### Dependencies
Before running MAGISTA, please make sure the following dependencies are fullfilled:
 * Bash
 * GCC (should be openmp-compatible)
 * GenDisCal: https://github.com/LM-UGent/GenDisCal
 * CheckM (optional): https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm
 * R 4.0.0 or higher, and the following libraries:
   - randomForest
   - moments
   - tidyr
 
 In case you are unfamiliar with R, you can install these libraries by opening R from your terminal and typing the following commands:
 ```
 install.packages(c("randomForest","moments","tidyr"))
 ```
 You may need to select a mirror in the process
 
 If you run into trouble with the installation, please contact me at Gleb.Goussarov@sckcen.be
### Install
Once you have checked that all dependencies are met, simply go to the directory where this repositary is located nd run `setup.sh`
## Usage
MAGISTA consists of two separate subprograms: MAGISTA.R and run_bin_analysis_noref.sh, both of which should be run from the command line.
### run_bin_analysis_noref.sh
Usage: `run_bin_analysis_noref.sh <file-list file> <prefix>`  
The <file-list file> should be a file with the absolute path of each bin that needs to be analysed  
<prefix> is a name that is prepended to the output.  
This script produces a csv file which can be parsed by MAGISTA.R  
If checkM is not in your PATH, the output will not be usable with the MAGISTIC argument
### MAGISTA.R
Usage: `Rscript MAGISTA.R <MAGISTA|MAGISTIC> <input file> [<output file>] [<model name>]`  
The first argument should be either `MAGISTA` or `MAGISTIC` depending on whether the checkM output should be included or not.  
<input file> should be the path to the file produced by run_bin_analysis_noref.sh
<output file> should be the path to the output file.
<model name> should be the name of the model, if you do not wish to use teh default models or create a new one.
  
