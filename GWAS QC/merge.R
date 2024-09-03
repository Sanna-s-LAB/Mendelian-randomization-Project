library(dplyr)
# Leggi il percorso del file dalla linea di comando
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

dat <- read.table(path, header=T)
dat1 <- read.table("/home/students/federica.grosso/nas/microbiome/Final_EUR_with_labels.txt", header=T)
# Unione dei dati
merged_data <- merge(dat, dat1, by.x = "SNP", all.x = T)

merged_data <- merged_data %>%
  mutate(diff_col = if_else(is.na(variant_id) | is.na(rsid), 1, if_else(variant_id == rsid, 0, 2))) %>%  # Aggiungi prima la colonna di differenza
  mutate(variant_id = if_else(is.na(variant_id), rsid, variant_id))  # Poi sostituisci gli NA

nome_file <- basename(path)
parti <- strsplit(nome_file, split = "_")[[1]]
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/GWAS_merged/merged_", parti[2], "_", parti[3], ".txt")

# Salvataggio del nuovo dataset
write.table(merged_data, file= output_file, sep="\t", row.names=FALSE, quote=F)
