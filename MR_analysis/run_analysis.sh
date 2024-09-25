#!/bin/bash
Rscript -e "library(R.utils)"
# Percorso del file R da eseguire
R_script="/home/.../microbiome/Results_CAD/Prova_CAD.R"
# Percorso della cartella contenente i file di esposizione
exposure_folder="/home/.../microbiome/Clumping_results_no_NA"
# Percorso del file di outcome
file_path_outcome="/home/.../microbiome/Outcomes/Coronary artery diseases/CAD_META.gz"
# Ciclo su tutti i file di esposizione nella cartella
for exposure_file in $exposure_folder/*.csv; do
    Rscript "$R_script" "$exposure_file" "$file_path_outcome"
done

Rscript "/home/.../microbiome/Results_CAD/FDR_correction.R"
