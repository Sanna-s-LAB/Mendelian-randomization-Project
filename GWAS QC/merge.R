library(dplyr)
## MERGE GWAS and REF GWAS

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

dat <- read.table(path, header=T)
dat1 <- read.table("/home/.../microbiome/Final_EUR_with_labels.txt", header=T)

merged_data <- merge(dat, dat1, by.x = "SNP", all.x = T)

# Create diff_col to decide which SNP to keep
merged_data <- merged_data %>%
  mutate(diff_col = if_else(is.na(variant_id) | is.na(rsid), 1, if_else(variant_id == rsid, 0, 2))) %>%  
  mutate(variant_id = if_else(is.na(variant_id), rsid, variant_id))  

file_name <- basename(path)
parts <- strsplit(file_name , split = "_")[[1]]
output_file <- paste0("/home/.../microbiome/GWAS_merged/merged_", parts[2], "_", parts[3], ".txt")

write.table(merged_data, file= output_file, sep="\t", row.names=FALSE, quote=F)
