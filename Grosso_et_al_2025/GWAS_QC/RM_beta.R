#######################################################################
##### REMOVE -4<beta<4
#######################################################################
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

data <- read.table(path, header=TRUE, sep="", fill=T)
data_filtered <- subset(data, beta > -4 & beta < 4)

nome_file <- basename(path)
output_file <- paste0(".../GWAS_new/new_", nome_file)
write.table(data_filtered, file= output_file, sep="\t", row.names=FALSE, quote=F)
