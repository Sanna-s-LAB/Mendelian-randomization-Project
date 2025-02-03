################################################################################
### Libraries
################################################################################
library(dplyr)
require(TwoSampleMR)
require(ieugwasr)
require(ggplot2)
require(MRPRESSO)
library(officer)
source(".../MR_PRESSO_plots.R")

# Read path from command line
args <- commandArgs(trailingOnly = TRUE)
file_path_exposure <- args[1]
file_path_outcome <- args[2]

# File name exposure
nome_fileE <- basename(file_path_exposure)
nome_file_senza_ext <- tools::file_path_sans_ext(nome_fileE)
# Divide name based on "_"
parti_nome <- strsplit(nome_file_senza_ext, "_")[[1]]
nome_exposure <- paste(parti_nome[2])

# File name outcome
nome_fileO <- basename(file_path_outcome)
parti_outcome <- strsplit(nome_fileO, split = "_")[[1]]
nome_outcome <- parti_outcome[1]

# Path where to save output files
base_path <- ".../Results/"

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
output_file <- paste0(base_path,"Harm_data_E_", nome_exposure, "_O_", nome_outcome, ".txt")
dat <- harmonise_data(dataE, dataO, action = 2)
write.table(dat, file = output_file, sep = "\t", quote=F)

# MR RESULTS
#output_file <- paste0(output_file <- paste0(base_path,"MR_E_", nome_exposure, "_O_", nome_outcome, ".csv")
mr_results <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
#write.csv(mr_results, file = output_file, quote=F, row.names = F)
# wald with one SNP, significant p-value < 0.05

################################################################################
### SENSITIVITY ANALYSIS
################################################################################

# Heterogeneity test
#output_file <-  paste0(base_path,"Het_E_", nome_exposure, "_O_", nome_outcome, ".txt")
het<-mr_heterogeneity(dat, method_list="mr_ivw")
#write.table(het, file = output_file, sep = "\t", quote=F)

# Pleiotropy test
#output_file <-  paste0(base_path,"Results_AMD/Ple_E_", nome_exposure, "_O_", nome_outcome, ".txt")
ple <- mr_pleiotropy_test(dat)
#write.table(ple, file = output_file, sep = "\t", quote=F)

# Single SNP analysis
output_file <- paste0(base_path,"SINGLESNP_E_", nome_exposure, "_O_", nome_outcome, ".txt")
res_single <- mr_singlesnp(dat)
write.table(res_single, file = output_file, sep = "\t", quote=F)

# Leave-one-out analysis
output_file <- paste0(base_path,"LOO_E_", nome_exposure, "_O_", nome_outcome, ".txt")
res_loo <- mr_leaveoneout(dat)
write.table(res_loo, file = output_file, sep = "\t", quote=F)

################################################################################
### PLOTS
################################################################################

## Scatter plot
output_file <- paste0(base_path,"plots_E_", nome_exposure, "_O_", nome_outcome, ".pdf")
pdf(output_file)

p1<-mr_scatter_plot(mr_results,dat)
#title("Scatter Plot")
print(p1)

## Forest plot
p2 <- mr_forest_plot(res_single)
#title("Forest Plot")
print(p2)

## Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)
#title("Leave-One-Out Plot")
print(p3)

## Funnel plot
## Asymmetry in a funnel plot is useful for gauging the reliability of a particular MR analysis
p4 <- mr_funnel_plot(res_single)
#title("Funnel Plot")
print(p4)

dev.off()

################################################################################
### MR-PRESSO
################################################################################
mr_presso <- tryCatch({
  result <- mr_presso(BetaOutcome="beta.outcome", 
                         BetaExposure="beta.exposure", 
                         SdOutcome="se.outcome", 
                         SdExposure="se.exposure", 
                         data=dat, 
                         OUTLIERtest = TRUE, 
                         DISTORTIONtest = TRUE, 
                         SignifThreshold = 0.05, 
                         NbDistribution = 1000, seed = 1)
  result  
  }, error = function(e) {
    # Manage error
    cat("Errore nell'esecuzione di mr_presso:", e$message, "\n")
    NULL  # If there is an error
  })

################################################################################
### Results output creation
################################################################################

  output_file <- paste0(base_path,"MRallRES_E_", nome_exposure, "_O_", nome_outcome, ".csv")

# Add MR-PRESSO row
  if (!is.null(mr_presso)) {
    mr_results[nrow(mr_results) + 1, ] <- mr_results[nrow(mr_results), ] 
    
    if (!is.null(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)){
      mr_results[nrow(mr_results), "nsnp"] <- mr_results$nsnp[1]- length(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
    } else {
      mr_results[nrow(mr_results), "nsnp"] <- mr_results$nsnp[1]
    }
    mr_results[nrow(mr_results), "method"] <- "MR-PRESSO"
    mr_results[nrow(mr_results), "b"] <- mr_presso$`Main MR results`$`Causal Estimate`[1]
    mr_results[nrow(mr_results), "se"] <- mr_presso$`Main MR results`$Sd[1]
    mr_results[nrow(mr_results), "pval"] <- mr_presso$`Main MR results`$`P-value`[1]
  } else { # If MR-PRESSO doesn't detect outliers
    mr_results[nrow(mr_results) + 1, ] <- mr_results[nrow(mr_results), ] 
    mr_results[nrow(mr_results), "method"] <- "MR-PRESSO"
    mr_results[nrow(mr_results), "b"] <- NA
    mr_results[nrow(mr_results), "se"] <- NA
    mr_results[nrow(mr_results), "pval"] <- NA
  }

# Add pleiotropy and heterogeneity to the file
 if (!is.null(ple) && nrow(ple) > 0 && ncol(ple) > 0){
  ple <- ple %>%
    rename_with(~paste0(., ".ple"))
	while(nrow(ple) < nrow(mr_results)) {
  ple[nrow(ple) + 1, ] <- NA
	}
  } else {
  ple <- data.frame(ple = rep(NA, nrow(mr_results)))
}
  if (!is.null(het) && nrow(het) > 0 && ncol(het) > 0) {
  het <- het %>%
    rename_with(~ paste0(.,".het"))
  while(nrow(het) < nrow(mr_results)) {
	het[nrow(het) + 1, ] <- NA
}
 } else  {
  het <- data.frame(het = rep(NA, nrow(mr_results)))
}

# Save scatter plot with Inverse Variance Weighted, Weighted Median and MR-PRESSO
mr_results <- mr_results[mr_results$method != "MR Egger", ]
p4 <- mr_scatter_plot(mr_results,dat)
print(p4)
ggsave("/.../SP_CAD.png", plot = p1[[1]], width= 6, height=6 , dpi = 300)

# Add variants column  
  final <- cbind(mr_results,ple,het)
  final <- final[, !colnames(final) %in% c("id.exposure.ple", "id.outcome.ple","outcome.ple","exposure.ple","method.het","id.exposure.het", "id.outcome.het","outcome.het","exposure.het")]
  filtered_variants <- dat$SNP[dat$mr_keep != FALSE]
  all_filtered_variants <- paste(filtered_variants, collapse = "; ")
  final$IV_list <- rep(all_filtered_variants, nrow(final))

if (!is.null(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) &&
    !all(is.na(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) &&
    !identical(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, "No significant outliers") &&
    length(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) > 0) 
{
    outliers_indices <- mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    filtered_variants_no_outliers <- filtered_variants[-outliers_indices]
    all_filtered_variants_no_outliers <- paste(filtered_variants_no_outliers, collapse = "; ")
    final$IV_list[nrow(final)] <- all_filtered_variants_no_outliers
} else {
    final$IV_list[nrow(final)] <- all_filtered_variants
}
	
write.csv(final, file = output_file, quote=F, row.names = F)

################################################################################################################################################
##### MR-PRESSO single results
######################################################################################################################
 output_file <- paste0(base_path,"MR_PRESSO/MR_PRESSO_E_", nome_exposure, "_O_", nome_outcome, ".docx") 
   doc <- read_docx()
   doc <- doc %>%
    body_add_par("MR Presso results:", style = "heading 1")
 
  if (!is.null(mr_presso$`Main MR results`)) {
  doc <- doc %>%
    body_add_par("Main MR results:", style = "heading 2") %>%
    body_add(mr_presso$`Main MR results`)}
  
  if (!is.null(mr_presso$`MR-PRESSO results`$`Outlier Test`)) {
  doc <- doc %>%
    body_add_par("MR-PRESSO results:", style = "heading 2") %>%
    body_add_par("Outlier test:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Outlier Test`)}
  
  if (!is.null(mr_presso$`MR-PRESSO results`$`Global Test`$RSSobs)) {
  doc <- doc %>%
    body_add_par("Global test - RSSobs:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Global Test`$RSSobs)}
  
  if (!is.null(mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue)) {
  doc <- doc %>%
    body_add_par("Global test - Pvalue:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue)}

if (!is.null(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) &&
    !all(is.na(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) &&
    !identical(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, "No significant outliers") &&
    length(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) > 0) {
  doc <- doc %>%
    body_add_par("Distortion Test - Outlier indices:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)}

if (!is.null(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`) &&  !any(is.na(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`))) {
  doc <- doc %>%
    body_add_par("Distortion Test - Distortion Coefficient:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`)}

if (!is.null(mr_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue) && !any(is.na(mr_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue))) {
  doc <- doc %>%
    body_add_par("Distortion Test - Pval:", style = "heading 3") %>%
    body_add(mr_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue)}
  # Save on word document
  print(doc, target = output_file)
  

# Scatter plot of MR-PRESSO (res dots for removed variants)
outliers_removed <- mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
dat$outlier_status <-  ifelse(seq_along(dat$SNP) %in% outliers_removed, "Removed", "Kept")
dat_kept <- subset(dat, outlier_status == "Kept")
dat_removed <- subset(dat, outlier_status == "Removed")
mr_results_kept <- mr(dat_kept, method_list = c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
output_file <-paste0(base_path,"MR_PRESSO/MR_PRESSO_E_", nome_exposure, "_O_", nome_outcome, ".pdf")
pdf(output_file)
p1 <- mr_scatter_plot_col(mr_results,dat,dat_kept,dat_removed)
print(p1)
p2 <- mr_scatter_plot_IVW(mr_results,mr_results_kept,dat,dat_kept,dat_removed)
print(p2)

dev.off()

rm(list = ls())
