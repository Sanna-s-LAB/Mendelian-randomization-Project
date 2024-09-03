#!/bin/bash

# Definizione della directory di input
input_directory="/home/students/federica.grosso/nas/microbiome/GWAS_microbiome"

# Controllo se la directory esiste
if [ ! -d "$input_directory" ]; then
  echo "Directory $input_directory non trovata."
  exit 1
fi

# Loop su tutti i file nella directory
for input_file in "$input_directory"/*
do
  filename=$(basename "$input_file")
  output_file="/home/students/federica.grosso/nas/microbiome/GWAS_new/$filename"
  # Esegui il comando awk sul file corrente e reindirizza l'output in un nuovo file nella cartella di output
  awk '{ if ($4 == "480000000" ) {$4 = "48000000";} print }' "$input_file" > "$output_file"
  # Esecuzione dello script R passando solo il file di input come argomento
  Rscript RM_beta.R "$output_file"
done
