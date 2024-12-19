#!/bin/bash
Rscript -e "library(R.utils)"

R_script=".../MR.R"
exposure_folder=".../Clumping_results_no_NA"
file_path_outcome=".../Outcome_file.gz"

for exposure_file in $exposure_folder/*.csv; do
    Rscript "$R_script" "$exposure_file" "$file_path_outcome"
done

# FDR correction
Rscript ".../FDR_correction.R"
