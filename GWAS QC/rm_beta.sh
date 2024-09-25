#!/bin/bash

####################################################################################
######## REMOVE -4<BETA<4 AND ADJUST BASE PAIR POSITION
####################################################################################

input_directory="/home/.../microbiome/GWAS_microbiome"
if [ ! -d "$input_directory" ]; then
  echo "Directory $input_directory non trovata."
  exit 1
fi

# Loop on all files
for input_file in "$input_directory"/*
do
  filename=$(basename "$input_file")
  output_file="/home/students/federica.grosso/nas/microbiome/GWAS_new/$filename"
  awk '{ if ($4 == "480000000" ) {$4 = "48000000";} print }' "$input_file" > "$output_file"
  Rscript RM_beta.R "$output_file"
done
