#!/bin/bash

# Directory contenente i file CSV da unire
directory="/home/students/federica.grosso/nas/microbiome/Results_CAD"

# File di output
output_file="/home/students/federica.grosso/nas/microbiome/Results_CAD/merged_MR_CAD.csv"

# Crea un file temporaneo per memorizzare i dati uniti
temp_file=$(mktemp)

# Inizializza una variabile per memorizzare l'intestazione
header=""

# Trova tutti i file CSV nella directory specificata e uniscili in un unico file temporaneo
for file in "$directory"/*.csv; do
    # Se è il primo file, memorizza l'intestazione
    if [ -z "$header" ]; then
        header=$(head -n 1 "$file")
        echo "$header" > "$output_file"
    fi

    # Salta la prima riga (intestazione) e aggiungi il resto al file temporaneo
    tail -n +2 "$file" >> "$temp_file"
done

# Concatena il file temporaneo con l'intestazione nel file di output finale
cat "$temp_file" >> "$output_file"

# Rimuovi il file temporaneo
rm "$temp_file"

echo "Unione completata. File creato: $output_file"
