################################################################
##### CLUMPING
################################################################
require(R.utils)
require(TwoSampleMR)

# Read command line
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]

# Read exposure
suppressWarnings({
  exp_data <- read_exposure_data(
    filename = file_path,
    phenotype_col = "trait",
    sep = "\t", # or "" if not tabulated
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
  bfile = "/home/.../pop_reference/EUR_phase3", #ref
  plink_bin = "/home/shared_tools/bin/plink"
)

nome_fileE <- basename(file_path)
nome_file_senza_ext <- tools::file_path_sans_ext(nome_fileE)
parti_nome <- strsplit(nome_file_senza_ext, "_")[[1]]
access <- paste(parti_nome[2], sep = "_")

output_file <- paste0(".../Clumping_results_no_NA/clumping_", access, ".csv")
write.csv(cd, file = output_file, quote=F, row.names = F)
