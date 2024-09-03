#!/bin/bash
Rscript -e "library(R.utils)"
# Definisci il percorso dello script R
script_r="/home/students/federica.grosso/nas/microbiome/CLUMPING_micr_server.R"
# Cartella contenente i file da analizzare
#folder="/home/students/federica.grosso/nas/microbiome/GWAS_microbiome"
folder="/home/students/federica.grosso/nas/microbiome/GWAS_final"
# Trova tutti i file con estensione _buildGRCh37.tsv nella cartella e nelle sotto-cartelle
#files=$(find "$folder" -type f -name "*_buildGRCh37.tsv")
files=$(find "$folder" -type f -name "merged_*.txt")
# Itera attraverso i file trovati
for file in $files
do
	# Esegui lo script R sul file corrente
	Rscript "$script_r" "$file"
done
