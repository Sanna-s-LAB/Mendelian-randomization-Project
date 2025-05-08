#!/bin/bash

# Directory dei file discovery
DISCOVERY_DIR="/home/res-fellows/federica.grosso/nas/microbiome/Outcomes/ProteinsUKBB2023_touse"
# Directory del file di mappatura
MAP_PREFIX="/home/res-fellows/federica.grosso/nas/metadata/olink_rsid_map_mac5_info03_b0_7_chr"

# Directory base per i file output
OUTPUT_BASE="/home/res-fellows/federica.grosso/nas/microbiome/Outcomes/ProteinsUKBB2023_merged"

# Ciclo su tutti i cromosomi (1–22)
for chr in {1..22}; do
    MAP_FILE="${MAP_PREFIX}${chr}_patched_v2.tsv.gz"

    # Controlla se il file di mappatura esiste
    if [[ ! -f "$MAP_FILE" ]]; then
        echo "Map file $MAP_FILE not found, skipping chromosome $chr"
        continue
    fi

    # Trova i file discovery del cromosoma corrente
    find "$DISCOVERY_DIR" -type f -name "*chr${chr}_*.gz" | while read -r DISCOVERY_FILE; do
        BASE_NAME=$(basename "$DISCOVERY_FILE" .gz)

        # Estrai la sottocartella relativa
        RELATIVE_SUBDIR=$(basename "$(dirname "$DISCOVERY_FILE")")
        OUTPUT_SUBDIR="${OUTPUT_BASE}/${RELATIVE_SUBDIR}"

        # Crea la sottocartella di output se non esiste
        mkdir -p "$OUTPUT_SUBDIR"

        OUTPUT_FILE="${OUTPUT_SUBDIR}/merged_${BASE_NAME}.gz"

        echo "Processing $DISCOVERY_FILE with map $MAP_FILE -> $OUTPUT_FILE"


        # Fusione + aggiunta colonna p-value
        zcat "$MAP_FILE" | awk -v id1=3 -v id2=1 -v rsid=4 '
        NR==FNR {
            key = $id2;
            value = $rsid;
            map[key] = value;
            next
        }
        {
            key = $id1;
            if (map[key]) {
                print $0, map[key];
            }
        }' - <(zcat "$DISCOVERY_FILE") | \
        awk -F' ' 'BEGIN {OFS=" "} {if (NR==1) print $0, "pvalue"; else print $0, 10^(-$13)}' | \
        gzip > "$OUTPUT_FILE"
    done
done
