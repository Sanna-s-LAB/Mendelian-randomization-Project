#!/bin/bash
#########################################################################
###### RUN CLUMPING
#########################################################################
Rscript -e "library(R.utils)"

script_r="/home/.../microbiome/Clumping.R"
folder="/home/.../microbiome/GWAS_final"
files=$(find "$folder" -type f -name "merged_*.txt")

for file in $files
do
	Rscript "$script_r" "$file"
done
