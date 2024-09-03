#!/bin/bash
Rscript -e "library(R.utils)"
# Percorso del file R da eseguire
R_script="/home/students/federica.grosso/nas/microbiome/Results_CAD/Prova_CAD.R"
# Percorso della cartella contenente i file di esposizione
exposure_folder="/home/students/federica.grosso/nas/microbiome/Clumping_results_no_NA"
# Percorso del file di outcome
file_path_outcome="/home/students/federica.grosso/nas/microbiome/Outcomes/Coronary artery diseases/CAD_META.gz"
# Ciclo su tutti i file di esposizione nella cartella
for exposure_file in $exposure_folder/*.csv; do
    Rscript "$R_script" "$exposure_file" "$file_path_outcome"
done
