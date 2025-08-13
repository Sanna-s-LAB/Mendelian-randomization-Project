#!/usr/bin/env Rscript
library("TwoSampleMR")
library("devtools")
library("plyr")
library("psych")
library("ggplot2")
library("pander")
library("knitr")
library("MRPRESSO")
library("data.table")
library("plyr")
library("dplyr")
library("ieugwasr")
library("LDlinkR")  
library("R.utils")
library("stringr")



# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Please provide exactly 3 arguments: exposure_file, outcome_file, sex (F or M)")
}

exp_data <- args[1]
out <- args[2]
sex <- args[3]

if( sex =="F") {
  base_path="~/Women/final_results_for_paper/"
} else if (sex=="M"){
  base_path="~/Men/final_results_for_paper/"
}

basename <- basename(out)
outcome_name <- sub("_.*", "", basename)
oid <- sub(".*(OID\\d+).*", "\\1", exp_data)



# === Find proxies using LDlinkR ===

exposure_snps <- trimws(as.character(exp_data$SNP))
outcome_snps <- trimws(as.character(out$rsid))

# Identify missing SNPs
missing_snps <- setdiff(exposure_snps, outcome_snps)

# Report results
cat("Missing SNPs in outcome:", length(missing_snps), "\n")
print(missing_snps)


n_proxies_found <- 0  # default value

# === Check if any missing SNPs have potential proxies ===
if (length(missing_snps) == 0) {
  cat("No missing SNPs, skipping proxy search.\n")
} else {
  # === Find proxies using LDlinkR ===
  proxy_results <- list()

  for (snp in missing_snps) {
    cat("Checking SNP:", snp, "\n")
    
    proxy_df <- NULL  

    tryCatch({
      proxy_df <- LDproxy(snp, pop = "EUR", r2d = "r2", token = "XXXX")
    }, error = function(e) {
      cat("Error during LDproxy() call for SNP", snp, ":", e$message, "\n")
    })

    # If proxy_df is still NULL or doesn't contain expected columns, skip
    if (is.null(proxy_df) || !"R2" %in% names(proxy_df) || !"RS_Number" %in% names(proxy_df)) {
      cat("No valid proxy data returned for", snp, "\n")
      proxy_results[[snp]] <- NA
      next
    }

    # Filter to R2 ≥ 0.8
    proxy_df_filtered <- proxy_df %>% filter(R2 >= 0.8)

    if (nrow(proxy_df_filtered) == 0) {
      cat("No proxies with R2 ≥ 0.8 for SNP:", snp, "\n")
      proxy_results[[snp]] <- NA
      next
    }

    # Match to outcome_snps 
    matching_proxies <- proxy_df_filtered$RS_Number[
      proxy_df_filtered$RS_Number %in% outcome_snps 
    ]

    if (length(matching_proxies) == 0) {
      cat("No valid matching proxies for", snp, "\n")
      proxy_results[[snp]] <- NA
      next
    }

    cat("Valid proxies for", snp, ":", paste(matching_proxies, collapse = ", "), "\n")
    proxy_results[[snp]] <- matching_proxies
    cat("Selected proxy for", snp, ":", proxy_results[[snp]], "\n")
  }

  # Count successful proxies
  n_proxies_found <- sum(!is.na(unlist(proxy_results)))
  cat("Proxies found and matched:", n_proxies_found, "\n")
}


# === Replace exposure data rows using valid proxies, if any were found ===

if (n_proxies_found > 0) {

  # Extract OID from exposure file path
  oid <- str_extract(basename(exp), "OID[0-9]+")

  # Build correct proxy file path using OID
  proxy_file_path <- file.path(
    paste0("combined_", oid, "_UKBB_proteomics_female_only.tsv.gz")
  )

  if (!file.exists(proxy_file_path)) {
    stop("Proxy file does not exist for OID: ", oid, "\nExpected file: ", proxy_file_path)
  }

  cat("Using proxy file:", proxy_file_path, "\n")

  # Read only necessary columns from the proxy file
  proxy_df <- fread(proxy_file_path, select = c(
    "rsid", "effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "p_value", "n"
  ))

  # For every proxy in proxy_results, update matching row in exp_data


  proxy_log <- data.frame(
  original_snp = character(),
  proxy_used = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(exp_data))) {
  original_snp <- as.character(exp_data$SNP[i]) 

  if (is.na(original_snp) || original_snp == "") {
    cat("Skipping invalid SNP at row", i, ": NA or empty\n")
    next
  }

  if (!(original_snp %in% names(proxy_results))) {
    cat("SNP not found in proxy_results at row", i, ":", original_snp, "\n")
    next
  }

 proxy_candidates <- proxy_results[[original_snp]]

if (!is.null(proxy_candidates) && length(proxy_candidates) > 0) {
  # Filter proxies that exist in the exposure proxy file
  proxy_rows <- proxy_df %>% filter(rsid %in% proxy_candidates)

  if (nrow(proxy_rows) == 0) {
    cat("No proxy SNPs found in exposure proxy file for", original_snp, "\n")
    next
  }

# Filter proxies with p-value ≤ 5e-8
significant_proxies <- proxy_rows %>% filter(p_value <= 5e-8)

if (nrow(significant_proxies) == 0) {
  cat("No significant proxies (p ≤ 5e-8) found for", original_snp, "\n")
  next
}

  # Select proxy with lowest p-value
  proxy_row <- significant_proxies %>% arrange(p_value) %>% slice(1)
  

  exp_data$SNP[i] <- proxy_row$rsid
  proxy_log <- rbind(proxy_log, data.frame(
  original_snp = original_snp,
  proxy_used = proxy_row$rsid,
  stringsAsFactors = FALSE
))

  exp_data$effect_allele.exposure[i]  <- proxy_row$effect_allele
  exp_data$other_allele.exposure[i]   <- proxy_row$other_allele
  exp_data$beta.exposure[i]           <- proxy_row$beta
  exp_data$se.exposure[i]             <- proxy_row$standard_error
  exp_data$eaf.exposure[i]            <- proxy_row$effect_allele_frequency
  exp_data$pval.exposure[i]           <- proxy_row$p_value
  exp_data$pval_origin.exposure[i]    <- proxy_row$p_value
  exp_dat$samplesize.exposure[i]     <- proxy_row$n

  cat("Replaced SNP", original_snp, "with best proxy", proxy_row$rsid, "at row", i, "\n")
    } else {
      cat("Proxy candidates not found for SNP", original_snp, "at row", i, "\n")
    }
  }


  # Save list of proxies used
  log_file <- paste0("proxy_used_", oid, "_TC", ".tsv")
  fwrite(proxy_log, log_file, sep = "\t")

  cat("Saved proxy usage log to:", log_file, "\n")

} else {
  cat("No valid proxies found — skipping proxy replacement step.\n")
}


#
# 1. Clumping of exposure data
#
data <- data.table::fread( 
  exp_data, 
  header=TRUE, sep="\t", stringsAsFactors=FALSE, 
)
data <- as.data.frame(data)
data <- plyr::rename(data, c("Pval" = "pval"))

exp <- ld_clump(
  dat = data,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.00000005,
  pop = "EUR",
  bfile = "~/POP_reference/EUR_phase3",
  plink_bin = "~/plink"
)
exp$Phenotype_new_2 <- oid

#
# 2. Format exposure and outcome data
#
exposure_dat <- format_data(
  exp,
  type = "exposure",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "pval",
  samplesize_col = "TotalSampleSize",
  phenotype_col = "Phenotype_new_2"
)
exposure_dat$data_source.exposure <-"textfile"


dim <- dim(exposure_dat)
print(dim)
names <- names(exposure_dat)
print(names) 
str<- str(exposure_dat)
print(str)
head<- head(exposure_dat)
print(head)
tail<- tail(exposure_dat)
print(tail)


outcome_dat<- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = out,
  sep = '\t',
  snp_col = 'rsid',
  beta_col = 'EFFECT_SIZE',
  se_col = 'SE',
  effect_allele_col = 'ALT',
  other_allele_col = 'REF',
  samplesize_col = 'N',
  eaf_col = "POOLED_ALT_AF",
  pval_col = 'pvalue',
  phenotype_col = "Phenotype"
)

dim <- dim(outcome_dat)
print(dim)

#
# 3. Harmonization
#
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dim<- dim(dat)
print(dim)

data.table::fwrite( 
  dat, 
  file= paste0(base_path,oid, ".dat_",outcome_name,"_",sex,".txt"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
  showProgress=FALSE
)

#
# 4. MR analysis
#
mr_results<- mr(dat,  method_list=c("mr_egger_regression",    
                                    "mr_weighted_median",   
                                    "mr_ivw", "mr_wald_ratio"))
print(mr_results)
data.table::fwrite( 
  mr_results, 
  file= paste0(base_path,oid, ".mr_res_",outcome_name,"_",sex,".txt"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
  showProgress=FALSE
)


#
# 5. Sensitivity analyses
#
het <- mr_heterogeneity(dat)

data.table::fwrite( 
  het, 
  file= paste0(base_path,oid, ".het_",outcome_name,"_",sex,".txt"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
  showProgress=FALSE
)



ple <- mr_pleiotropy_test(dat)

data.table::fwrite( 
  ple, 
  file= paste0(base_path,oid, ".ple_",outcome_name,"_",sex,".txt"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
  showProgress=FALSE
)



res_leaveoneout <- tryCatch({
  mr_leaveoneout(dat)
}, error = function(e) {
  cat("Errore in mr_leaveoneout:", e$message, "\n")
  data.frame()  # Restituisce un data.frame vuoto in caso di errore
})
data.table::fwrite( 
  res_leaveoneout, 
  file= paste0(base_path,oid, ".leaveoneout_",outcome_name,"_",sex,".txt"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
  showProgress=FALSE
)




################################################################################
### PLOTS
################################################################################

## Create PDF
output_file <- paste0(base_path, "plot_", oid, "_",outcome_name, ".pdf")
pdf(output_file)

## Scatter plot
tryCatch({
  p1 <- mr_scatter_plot(mr_results, dat)
  print(p1)
}, error = function(e) {
  cat("Errore nel grafico Scatter:", e$message, "\n")
  plot.new()
  text(0.5, 0.5, "Errore nel grafico Scatter", cex = 1.5)
})

## Leave-one-out plot
tryCatch({
  if (nrow(res_leaveoneout) > 0) {
    p2 <- mr_leaveoneout_plot(res_leaveoneout)
    print(p2)
  } else {
    cat("Leave-one-out plot non generato: dati assenti.\n")
    plot.new()
    text(0.5, 0.5, "Dati Leave-one-out non disponibili", cex = 1.5)
  }
}, error = function(e) {
  cat("Errore nel grafico Leave-one-out:", e$message, "\n")
  plot.new()
  text(0.5, 0.5, "Errore nel grafico Leave-one-out", cex = 1.5)
})

dev.off()





################################################################################
### MR-PRESSO
################################################################################
dat1 <- dat %>% filter(mr_keep)
mr_presso <- tryCatch({
  result <- mr_presso(BetaOutcome="beta.outcome", 
                      BetaExposure="beta.exposure", 
                      SdOutcome="se.outcome", 
                      SdExposure="se.exposure", 
                      data=dat1, 
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
### Final table of results
################################################################################

output_file <- paste0(base_path,"MRallRES_", oid, "_",outcome_name, ".csv")

# Add MR-PRESSO row
if (!is.null(mr_presso)) {
  mr_results[nrow(mr_results) + 1, ] <- mr_results[nrow(mr_results), ] 
  outliers_indices <- mr_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  
  if (!is.null(outliers_indices)) {
    if (!"All SNPs considered as outliers" %in% outliers_indices) {
      mr_results[nrow(mr_results), "nsnp"] <- mr_results$nsnp[1] - length(outliers_indices)
    } else {
      mr_results[nrow(mr_results), "nsnp"] <- 0
    }
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

# Add ple and heterogeneity to the file
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

mr_presso
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
  print(outliers_indices)
  # Ensure it's numeric
  if (is.numeric(outliers_indices) && length(outliers_indices) > 0) {
    filtered_variants_no_outliers <- filtered_variants[-outliers_indices]
    all_filtered_variants_no_outliers <- paste(filtered_variants_no_outliers, collapse = "; ")
    final$IV_list[nrow(final)] <- all_filtered_variants_no_outliers
  } else if (outliers_indices == "All SNPs considered as outliers") 
  {
    final$IV_list[nrow(final)]<- NA
  }
  else
  {
    final$IV_list[nrow(final)] <- all_filtered_variants
  }
} else {
  final$IV_list[nrow(final)] <- all_filtered_variants
}

write.csv(final, file = output_file, quote=F, row.names = F)
