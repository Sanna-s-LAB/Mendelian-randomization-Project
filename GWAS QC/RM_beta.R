
# Leggi il percorso del file dalla linea di comando
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

# Lettura del file
data <- read.table(path, header=TRUE, sep="", fill=T)

# Filtraggio delle righe in base al valore della colonna beta
data_filtrato <- subset(data, beta > -4 & beta < 4)

nome_file <- basename(path)

output_file <- paste0("/home/students/federica.grosso/nas/microbiome/GWAS_new/new_", nome_file)

# Salvataggio del nuovo dataset
write.table(data_filtrato, file= output_file, sep="\t", row.names=FALSE, quote=F)
