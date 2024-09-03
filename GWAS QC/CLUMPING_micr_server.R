require(R.utils)
require(TwoSampleMR)

# Leggi il percorso del file dalla linea di comando
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
suppressWarnings({
  exp_data <- read_exposure_data(
    filename = file_path,
    phenotype_col = "trait",
    sep = "\t",
    clump = F,
    snp_col = "rsid_to_use",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    #   units_col = "Units",
    pos_col = "base_pair_location",
    chr_col = "chromosome",
    #   samplesize_col = "TotalSampleSize"
  )
})

cd <- clump_data(
  exp_data,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 0.000005, #5e-6
  clump_p2 = 0.000005,
  bfile = "/home/students/federica.grosso/nas/pop_reference/EUR_phase3", #riferimento
  plink_bin = "/home/shared_tools/bin/plink"
)

nome_fileE <- basename(file_path)
nome_file_senza_ext <- tools::file_path_sans_ext(nome_fileE)
# Divide il nome del file basandosi sul carattere "_"
parti_nome <- strsplit(nome_file_senza_ext, "_")[[1]]
access <- paste(parti_nome[2], sep = "_")

#output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Clumping_results/clumping_", access, ".csv")
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Clumping_results_no_NA/clumping_", access, ".csv")
write.csv(cd, file = output_file, quote=F, row.names = F)
