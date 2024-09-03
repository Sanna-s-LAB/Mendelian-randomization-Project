library(openxlsx)
library(readxl)
# CORRECTION

# Upload the dataset
dataset <- read.csv("/home/students/federica.grosso/nas/microbiome/Results_CAD/merged_MR_CAD.csv")

# Select rows with IVW
datasub <- subset(dataset, method == "Inverse variance weighted")

pval <- datasub$pval

# Adjust p-values with Benjamin-Hochberg
pval_adj <- p.adjust(pval, method = "BH", n = length(pval))
#pval_adj <0.05

# unire le colonnne

data2 <- data.frame(pval,pval_adj)

# Unione dei dati
merged_data <- merge(dataset, data2, by = "pval", all = TRUE)
# Ottieni l'indice della colonna pval nel dataset mergeato
pval_index <- grep("pval", colnames(merged_data))

# Sposta la colonna pval come penultima
merged_data <- merged_data[, c(setdiff(1:ncol(merged_data), pval_index), pval_index)]
write.xlsx(merged_data, "/home/students/federica.grosso/nas/microbiome/Results_CAD/merged_MR_CAD_corrected.xlsx")
