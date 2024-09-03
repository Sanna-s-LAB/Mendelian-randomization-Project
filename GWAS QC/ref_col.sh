#!/bin/bash
# Definizione della directory di input
input_directory="/home/students/federica.grosso/nas/microbiome/GWAS_new"
# Controllo se la directory esiste
if [ ! -d "$input_directory" ]; then
  echo "Directory $input_directory non trovata."
  exit 1
fi
# Loop su tutti i file nella directory
for input_file in "$input_directory"/new_*.tsv
do
  echo "Elaborazione del file: $input_file"
  # Crea un file temporaneo per ogni file di input
  awk '{print $3":"$4, $0}' "$input_file" > temp_file.txt
  awk '{ if (NR == 1) sub(/^chromosome:base_pair_location/, "SNP"); print }' temp_file.txt > "${input_file%.tsv}_SNP_file.txt"
  output_file="${input_file%.tsv}_SNP_file.txt"
  rm temp_file.txt
  Rscript merge.R "$output_file"
  echo "Operazioni completate. Il file finale è $output_file"
done
