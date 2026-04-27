process_tables <- function(df) {
  df$Significant <- NA
  list1 <- list()
  v <- c("Inverse variance weighted", "Weighted median", "MR Egger")
  combinations <-  selected_pairs <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$fdr) & sub_df$method == "Inverse variance weighted" & sub_df$fdr < 0.05) &
             any(!is.na(sub_df$pval.ple) & sub_df$method == "MR Egger" & sub_df$pval.ple >= 0.05) &
             any(!is.na(sub_df$pval) & sub_df$method == "Weighted median" & sub_df$pval < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval < 0.05 | is.na(sub_df$pval)))| !("MR-PRESSO"%in% sub_df$method))
             & any("Inverse variance weighted" %in% sub_df$method & 
                 sub_df$method == "Inverse variance weighted" & 
                 sub_df$Q_pval.het >= 0.05))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr) & sub_df$fdr < 0.05 & sub_df$exp_var>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$fdr) & sub_df$fdr < 0.05 & sub_df$Q_pval.het >= 0.05)))) {
      
      sub_df$Significant <- "YES"
    } else {
      sub_df$Significant <- "NO"
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  final <- do.call(rbind, list1)
  final <- final[order(final$Significant, decreasing = T), ]
  rownames(final) <- NULL
  
  return(final)
}


explained_variance <- function(data, N)
{
  eaf = data$eaf.exposure
  MAF <- ifelse(eaf <= 0.5, eaf, 1-eaf)
  beta = data$beta.exposure
  se =data$se.exposure
  R2 = 2 * beta^2 * MAF * (1 - MAF) / (2 * beta^2 * MAF * (1 - MAF) + se^2 * 2 * N * MAF * (1 - MAF))
  R_2=sum(R2)
  return(R_2)
}

filter_leaveoneout <- function(not_sig, all, IV_list_column) {
  to_keep <- list()
  
  for (i in seq_len(nrow(not_sig))) {
    exposure_i <- not_sig$exposure[i]
    outcome_i <- not_sig$outcome[i]
    snp_i <- not_sig$SNP[i]
    
    rwo_sig <- subset(all, exposure == exposure_i & outcome == outcome_i & method == "MR-PRESSO")
    
    if (nrow(rwo_sig) > 0) {
      # Check IV_list
      present_snp <- any(sapply(rwo_sig[[IV_list_column]], function(s) snp_i %in% strsplit(s, ";\\s*")[[1]]))
      
      # If there isn't the snp in MR-PRESSO, we keep this couple
      if (!present_snp) {
        to_keep <- append(to_keep, list(c(exposure_i, outcome_i)))
      }
    }
  }
  
  if (length(to_keep) == 0) {
    message("No exposure-outcome pairs to keep: list_to_keep is empty.")
    return(NULL)
  }
  else {
  # Convert the list to a data frame
  list_to_keep <- do.call(rbind, to_keep)
  colnames(list_to_keep) <- c("exposure", "outcome")
  list_to_keep <- as.data.frame(list_to_keep)
  
  toremove_W <- anti_join(not_sig, list_to_keep, by = c("exposure", "outcome"))
  
  return(list(keep = list_to_keep, remove = toremove_W))
  }
}
