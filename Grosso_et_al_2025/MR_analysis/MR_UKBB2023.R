#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)

exp =args[1]
out=args[2]

library("TwoSampleMR")
library("devtools")
library("plyr")
library("lattice")
library("psych")
library("ggplot2")
library("pander")
library("knitr")
library("MRPRESSO")
library("data.table")
library("plyr")
library("dplyr")



exposure_dat <- read.csv(exp, header=T)

setwd("~/microbiome/Results_IP_UKBB2023/")

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

outcome_name <- sub(".*(OID\\d+).*", "\\1", out)

file_nameE <- basename(exp)
nome_file_senza_ext <- tools::file_path_sans_ext(file_nameE)

name_pars <- strsplit(nome_file_senza_ext, "_")[[1]]
exposure_name <- paste(name_parts[2])

outcome_dat<- read_outcome_data(
  snps = exposure_dat$rsid,
  filename = out,
  sep = ' ',
  snp_col = 'rsid',
  beta_col = 'BETA',
  se_col = 'SE',
  effect_allele_col = 'ALLELE1',
  other_allele_col = 'ALLELE0',
  samplesize_col = 'N',
  eaf_col = "A1FREQ",
  pval_col = 'pvalue'
)
outcome_dat$outcome <- outcome_name

dim <- dim(outcome_dat)
print(dim)

names(exposure_dat)[names(exposure_dat) == "rsid"] <- "SNP"

base_path<- "~/microbiome/Results_IP_UKBB2023/"




dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dim<- dim(dat)
print(dim)

  data.table::fwrite( 
    dat, 
    file= paste0(base_path,exposure_name,"_",outcome_name,"_dat",".txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
    showProgress=FALSE
  )

mr_results<- mr(dat,  method_list=c("mr_egger_regression",    
                                                "mr_weighted_median",   
                                                "mr_ivw", "mr_wald_ratio"))
print(mr_results)


 data.table::fwrite( 
    mr_results, 
    file= paste0(base_path,exposure_name,"_",outcome_name, "_mr_res.txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
    showProgress=FALSE
  )





het <- mr_heterogeneity(dat)

data.table::fwrite( 
    het, 
    file= paste0(base_path,exposure_name,"_",outcome_name, "_het.txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
    showProgress=FALSE
  )



ple <- mr_pleiotropy_test(dat)

data.table::fwrite( 
    ple, 
    file= paste0(base_path,exposure_name,"_",outcome_name, "_ple.txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
    showProgress=FALSE
  )



res_leaveoneout <- tryCatch({
  mr_leaveoneout(dat)
}, error = function(e) {
  cat("Errore in mr_leaveoneout:", e$message, "\n")
  data.frame()  # Empty dataframe if there is an error
})
data.table::fwrite( 
    res_leaveoneout, 
    file= paste0(base_path,exposure_name,"_",outcome_name, "_leaveoneout.txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", 
    showProgress=FALSE
  )






################################################################################
### PLOTS
################################################################################



## Crea il file PDF
output_file <- paste0(base_path, "plot_",exposure_name,"_",outcome_name, ".pdf")
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
### Results output creation
################################################################################

  output_file <- paste0(base_path,"MRallRES_",exposure_name,"_",outcome_name,  ".csv")

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


 
