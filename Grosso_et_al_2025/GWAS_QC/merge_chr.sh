#!/bin/bash

BASE_DIR="/home/res-fellows/federica.grosso/nas/microbiome/Outcomes/ProteinsUKBB2023_merged"

find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | while read subfolder; do
    echo "Unione file in: $subfolder"

    files=($(find "$subfolder" -maxdepth 1 -name "*.gz" | sort))
    [ ${#files[@]} -eq 0 ] && continue

    foldername=$(basename "$subfolder")
    output_file="$BASE_DIR/${foldername}.gz"

    {
        # Primo file completo (intestazione + dati)
        zcat "${files[0]}"

        # Altri file senza intestazione
        for f in "${files[@]:1}"; do
            zcat "$f" | tail -n +2
        done
    } | gzip > "$output_file"

    echo "Creato: $output_file"
done
