################################################################################
### Libraries
################################################################################
library(dplyr)
require(TwoSampleMR)
require(ieugwasr)
require(ggplot2)
require(MRPRESSO)
library(officer)
source("/home/students/federica.grosso/nas/MR_PRESSO_plots.R")

# Read path from command line
args <- commandArgs(trailingOnly = TRUE)
file_path_exposure <- args[1]
file_path_outcome <- args[2]

# File name exposure
file_nameE <- basename(file_path_exposure)
file_name_without_ext <- tools::file_path_sans_ext(file_nameE)
# Divide il nome del file basandosi sul carattere "_"
name_parts <- strsplit(file_name_without_ext, "_")[[1]]
exposure_name <- paste(name_parts[2])

# File name outcome
file_nameO <- basename(file_path_outcome)
outcome_parts <- strsplit(file_nameO, split = "_")[[1]]
outcome_name <- outcome_parts[1]

# Upload exposure file and outcome file
dataE<-read.csv(file_path_exposure, header=T)
suppressWarnings({
  dataO<-read_outcome_data(file_path_outcome, 
                           snps = dataE$SNP,
                           sep="\t", 
                           phenotype_col = "Outcome",
                           snp_col = "oldID",
                           beta_col = "Effect",
                           se_col="StdErr",
                           eaf_col = "Freq1",
                           effect_allele_col = "Allele1",
                           other_allele_col = "Allele2",
                           pval="P-value")
  dataO$outcome <- "CAD"
})

################################################################################
### MENDELIAN RANDOMIZATION ANALYSIS
################################################################################

# HARMONIZE
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/Harm_data_E_", exposure_name, "_O_", outcome_name, ".txt")
dat <- harmonise_data(dataE, dataO, action = 2)
write.table(dat, file = output_file, sep = "\t", quote=F)

# MR RESULTS
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/MR_E_", exposure_name, "_O_", outcome_name, ".csv")
mr_results <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
write.csv(mr_results, file = output_file, quote=F, row.names = F)
# wald with one SNP, significant p-value < 0.05

################################################################################
### SENSITIVITY ANALYSIS
################################################################################

# Heterogeneity test
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/Het_E_", exposure_name, "_O_", outcome_name, ".txt")
het<-mr_heterogeneity(dat)
write.table(het, file = output_file, sep = "\t", quote=F)

# Pleiotropy test
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/Ple_E_", exposure_name, "_O_", outcome_name, ".txt")
ple <- mr_pleiotropy_test(dat)
write.table(ple, file = output_file, sep = "\t", quote=F)

# Single SNP analysis
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/SINGLESNP_E_", exposure_name, "_O_", outcome_name, ".txt")
res_single <- mr_singlesnp(dat)
write.table(res_single, file = output_file, sep = "\t", quote=F)

# Leave-one-out analysis
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/LOO_E_", exposure_name, "_O_", outcome_name, ".txt")
res_loo <- mr_leaveoneout(dat)
write.table(res_loo, file = output_file, sep = "\t", quote=F)

# MR Steiger directionality test
# output_file <- paste0("/home/students/federica.grosso/nas/MRresults2/BMI/dir_E_", nome_proteina, "_O_", nome_outcome, ".txt")
#dir <- directionality_test(dat) # non posso farlo perchè non ho sample size
#write.table(dir, file = output_file, sep = "\t", quote=F, append = TRUE)

################################################################################
### PLOTS
################################################################################

## Scatter plot
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/plots_E_", exposure_name, "_O_", outcome_name, ".pdf")
pdf(output_file)

p1<-mr_scatter_plot(mr_results,dat)
print(p1)

## Forest plot
p2 <- mr_forest_plot(res_single)
print(p2)

## Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)
print(p3)

## Funnel plot
p4 <- mr_funnel_plot(res_single)
print(p4)

dev.off()

################################################################################
### MR EGGER
################################################################################

output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/MR_PRESSO/MR_PRESSO_E_", exposure_name, "_O_", outcome_name, ".txt")
mr_presso <- mr_presso(BetaOutcome="beta.outcome", 
                         BetaExposure="beta.exposure", 
                         SdOutcome="se.outcome", 
                         SdExposure="se.exposure", 
                         data=dat, 
                         OUTLIERtest = TRUE, 
                         DISTORTIONtest = TRUE, 
                         SignifThreshold = 0.05, 
                         NbDistribution = 1000, seed = 1)
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/MR_PRESSO/MR_PRESSO_E_", exposure_name, "_O_", outcome_name, ".docx")

# Save results in a file word
doc <- read_docx()
doc <- doc %>%
  body_add_par("MR Presso results:", style = "heading 1")
  
# 'Main MR results'
  doc <- doc %>%
    body_add_par("Main MR results:", style = "heading 2") %>%
    body_add(mr_presso$`Main MR results`)
  
# 'MR-PRESSO results'
  doc <- doc %>%
    body_add_par("MR-PRESSO results:", style = "heading 2") %>%
    body_add_par("Outlier test:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Outlier Test`)

# Global test RSS
  doc <- doc %>%
    body_add_par("Global test - RSSobs:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Global Test`$RSSobs)
  
# Global test p-value
  doc <- doc %>%
    body_add_par("Global test - Pvalue:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
  
#  Distortion test
  doc <- doc %>%
    body_add_par("Distortion Test - Outlier indices:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
  doc <- doc %>%
    body_add_par("Distortion Test - Distortion Coefficient:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`)
  
  doc <- doc %>%
    body_add_par("Distortion Test - Pval:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue)
print(doc, target = output_file)

# Scatter plot of MR-PRESSO (res dots for removed variants)
outliers_removed <- mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
dat$outlier_status <-  ifelse(seq_along(dat$SNP) %in% outliers_removed, "Removed", "Kept")

dat_kept <- subset(dat, outlier_status == "Kept")
dat_removed <- subset(dat, outlier_status == "Removed")
mr_results_kept <- mr(dat_kept, method_list = c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
output_file <- paste0("/home/students/federica.grosso/nas/microbiome/Results_CAD/MR_PRESSO/MR_PRESSO_E_", exposure_name, "_O_", outcome_name, ".pdf")

pdf(output_file)
p1 <- mr_scatter_plot_col(mr_results,dat,dat_kept,dat_removed)
print(p1)
p2 <- mr_scatter_plot_IVW(mr_results,mr_results_kept,dat,dat_kept,dat_removed)
print(p2)
dev.off()


