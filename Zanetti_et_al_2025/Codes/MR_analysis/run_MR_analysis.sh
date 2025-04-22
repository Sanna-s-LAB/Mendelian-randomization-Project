#!/bin/bash

# Female files
outcomes_F=(~/HDL_female.txt ~/nonHDL_female.txt ~/LDL_female.txt ~/TC_female.txt ~/TG_female.txt)
exposures_F=(/mnt/dzanetti/UKB-PPP_SS/toshare_W_only/*F_only_regenie.gz)

for exp in "${exposures_F[@]}"; do
  for out in "${outcomes_F[@]}"; do
    exp_base=$(basename "$exp" .gz)
    out_base=$(basename "$out" .txt)
    log_file="${exp_base}_${out_base}_F.log"
    Rscript --vanilla run_analysis.R "$exp" "$out" "F" > "$log_file" 2>&1
  done
done

# Male files
outcomes_M=(~/HDL_male.txt ~/nonHDL_male.txt ~/LDL_male.txt ~/TC_male.txt ~/TG_male.txt)
exposures_M=(/mnt/dzanetti/UKB-PPP_SS/toshare_M_only/*M_only_regenie.gz)

for exp in "${exposures_M[@]}"; do
  for out in "${outcomes_M[@]}"; do
    exp_base=$(basename "$exp" .gz)
    out_base=$(basename "$out" .txt)
    log_file="${exp_base}_${out_base}_M.log"
    Rscript --vanilla run_analysis.R "$exp" "$out" "M" > "$log_file" 2>&1
  done
done
