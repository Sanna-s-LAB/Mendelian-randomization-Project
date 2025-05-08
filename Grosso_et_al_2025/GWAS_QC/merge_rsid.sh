#!/bin/bash

DISCOVERY_DIR="~/microbiome/Outcomes/ProteinsUKBB2023_touse"
MAP_PREFIX="~/metadata/olink_rsid_map_mac5_info03_b0_7_chr"

# Mother output directiory 
OUTPUT_BASE="~/microbiome/Outcomes/ProteinsUKBB2023_merged"

# For all chromosomes
for chr in {1..22}; do
    MAP_FILE="${MAP_PREFIX}${chr}_patched_v2.tsv.gz"

    if [[ ! -f "$MAP_FILE" ]]; then
        echo "Map file $MAP_FILE not found, skipping chromosome $chr"
        continue
    fi

    # Find chromosome file in the directory
    find "$DISCOVERY_DIR" -type f -name "*chr${chr}_*.gz" | while read -r DISCOVERY_FILE; do
        BASE_NAME=$(basename "$DISCOVERY_FILE" .gz)

        # Extract subdir
        RELATIVE_SUBDIR=$(basename "$(dirname "$DISCOVERY_FILE")")
        OUTPUT_SUBDIR="${OUTPUT_BASE}/${RELATIVE_SUBDIR}"

        # Create subdirectory if it doesn't exist
        mkdir -p "$OUTPUT_SUBDIR"

        OUTPUT_FILE="${OUTPUT_SUBDIR}/merged_${BASE_NAME}.gz"

        echo "Processing $DISCOVERY_FILE with map $MAP_FILE -> $OUTPUT_FILE"


        # Merge the rsid column and create the pvalue column
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
