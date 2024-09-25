#!/bin/bash
Rscript -e "library(R.utils)"

R_script="/home/.../microbiome/Results_CAD/Prova_CAD.R"
exposure_folder="/home/.../microbiome/Clumping_results_no_NA"
file_path_outcome="/home/.../microbiome/Outcomes/Coronary artery diseases/CAD_META.gz"

for exposure_file in $exposure_folder/*.csv; do
    Rscript "$R_script" "$exposure_file" "$file_path_outcome"
done

# FDR correction
Rscript "/home/.../microbiome/Results_CAD/FDR_correction.R"
