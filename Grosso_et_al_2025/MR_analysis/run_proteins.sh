#!/bin/bash

if [ -z "$1" ]; then
  echo "Uso: $0 <server_id: 0|1|2>"
  exit 1
fi

SERVER_ID=$1
MAX_JOBS=30
TOTAL_SERVERS=3

exp_path="~/microbiome/Clumping_results_no_NA1"
file_path_outcome="~/microbiome/Outcomes/ProteinsUKBB2023_merged"
script_path="~/microbiome/Results_IP_UKBB2023/MR_UKBB2023.R"
results_dir="~/microbiome/Results_IP_UKBB2023"

mkdir -p "$results_dir"

pair_index=0

for outcome_file in "$file_path_outcome"/*.gz; do
  for exposure_file in "$exp_path"/*.csv; do

    if [ $((pair_index % TOTAL_SERVERS)) -ne $SERVER_ID ]; then
      ((pair_index++))
      continue
    fi

    exp_base=$(basename "$exposure_file" .csv)
    exp_base=${exp_base#clumping_}
    out_file=$(basename "$outcome_file" .gz)
	log_file="${results_dir}/log_file_${exp_base}__${out_file}.log"
	
    # Extract outcome OID name
    if [[ $out_file =~ OID[0-9]+ ]]; then
      OID="${BASH_REMATCH[0]}"
    else
      echo "Error: impossible to extract OID from $out_file"
      ((pair_index++))
      continue
    fi

    output_result="${results_dir}/MRallRES_${exp_base}_${OID}.csv"

    # If the MRallRES file exists, then it skip the analysis
    if [ -f "$output_result" ]; then
      echo "Skipped: $exp_base vs $OID (already run)"
      ((pair_index++))
      continue
    fi

    # Check the jobs running
    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
      sleep 1
    done

    (
      echo "Run: $exp_base vs $OID"
      Rscript "$script_path" "$exposure_file" "$outcome_file" \
        >& "$log_file"
    ) &

    ((pair_index++))

  done
done

wait
echo "Server $SERVER_ID: All jobs completed."
