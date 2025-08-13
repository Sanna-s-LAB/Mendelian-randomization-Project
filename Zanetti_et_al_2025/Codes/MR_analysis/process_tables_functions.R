# function to compute explained variance
explained_variance <- function(data, N)
{
  eaf = data$eaf.exposure
  MAF <- ifelse(eaf <= 0.5, eaf, 1-eaf)
  beta = data$beta.exposure
  se =data$se.exposure
  R2 = 2 * beta^2 * MAF * (1 - MAF) / (2 * beta^2 * MAF * (1 - MAF) + se^2 * 2 * N * MAF * (1 - MAF))
  #R_2=sum(R2)
  return(R2)
}

# function to process the leave-one-out data
process_leave_one_out <- function(
    leave_one_out, 
    base_data, 
    type = "men",
    p_col = "p",
    snp_col = "SNP",
    exposure_col = "exposure",
    outcome_col = "outcome",
    method_col = "method",
    iv_list_col = "IV_list.men",  
    loo_col_suffix = "loo"      
) {
  # remove usefulness columns
  leave_one_out$id.exposure <- NULL
  leave_one_out$id.outcome <- NULL
  leave_one_out$samplesize <- NULL
  leave_one_out <- na.omit(leave_one_out)
  leave_one_out <- leave_one_out[!(leave_one_out[[snp_col]] == "All"),]
  
  # Step 1: create column loo.[type]
  loo_with_label <- leave_one_out
  loo_col_name <- paste0(loo_col_suffix, ".", type)
  
  loo_with_label[[loo_col_name]] <- ifelse(loo_with_label[[p_col]] < 0.05, "no", loo_with_label[[snp_col]])
  
  # Step 2: group by exposure and outcome
  grouped <- split(loo_with_label, list(loo_with_label[[exposure_col]], loo_with_label[[outcome_col]]), drop = TRUE)
  
  # Step 3: apply on each group
  result_list <- lapply(grouped, function(df) {
    exposure <- unique(df[[exposure_col]])
    outcome <- unique(df[[outcome_col]])
    
    if (all(df[[loo_col_name]] == "no")) {
      loo_value <- "no"
    } else {
      loo_value <- paste(df[[loo_col_name]][df[[loo_col_name]] != "no"], collapse = ",")
    }
    
    data.frame(
      exposure = exposure, 
      outcome = outcome, 
      tmp_loo = loo_value, 
      stringsAsFactors = FALSE
    )
  })
  
  loo_summary <- do.call(rbind, result_list)
  rownames(loo_summary) <- NULL
  names(loo_summary)[names(loo_summary) == "tmp_loo"] <- loo_col_name
  
  # merge by base_data
  merged <- merge(base_data, loo_summary, by = c(exposure_col, outcome_col), all.x = TRUE)
  
  # Get combinations of exposure-outcome
  combinations <- unique(merged[c(exposure_col, outcome_col)])
  
  list_out <- list()
  for (i in seq_len(nrow(combinations))) {
    exp <- combinations[[exposure_col]][i]
    out <- combinations[[outcome_col]][i]
    sub_df <- subset(merged, merged[[exposure_col]] == exp & merged[[outcome_col]] == out)
    
    if (any(!is.na(sub_df[[loo_col_name]]) & 
            sub_df[[method_col]] == "MR-PRESSO" & 
            sub_df[[loo_col_name]] != "no")) {
      
      sub_sub <- subset(sub_df, sub_df[[method_col]] == "MR-PRESSO")
      
      # Divide SNPs from loo and IV_list
      loo_snps <- unlist(strsplit(sub_sub[[loo_col_name]], ","))
      iv_snps <- unlist(strsplit(sub_sub[[iv_list_col]], "; "))
      
      # Keep only the SNPs in IV_list
      valid_snps <- loo_snps[loo_snps %in% iv_snps]
      
      # Assign the result
      if (length(valid_snps) == 0) {
        sub_df[[loo_col_name]] <- "no"
      } else {
        sub_df[[loo_col_name]] <- paste(valid_snps, collapse = ",")
      }
    }
    list_out[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final_df <- do.call(rbind, list_out)
  rownames(final_df) <- NULL
  return(final_df)
}

# function to add significance column in main results (FDR corrected)
process_tables <- function(df) {
  df$SignificantBysex_F <- NA
  list1 <- list()
  v <- c("Inverse variance weighted", "Weighted median", "MR Egger")
  combinations <-  df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$fdr.women) & sub_df$method == "Inverse variance weighted" & sub_df$fdr.women < 0.05) &
             any(!is.na(sub_df$pval.ple.women) & sub_df$method == "MR Egger" & sub_df$pval.ple.women > 0.05) &
             any(!is.na(sub_df$loo.women) & sub_df$loo.women=="no") &
             any(!is.na(sub_df$pval.women) & sub_df$method == "Weighted median" & sub_df$pval.women < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.women < 0.05 | is.na(sub_df$pval.women)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.women) & sub_df$fdr.women < 0.05 & sub_df$exp_var.women>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.women == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$fdr.women) & sub_df$fdr.women < 0.05)))) {
      
      sub_df$SignificantBysex_F <- 1
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.women) & sub_df$fdr.women < 0.05 & sub_df$exp_var.women<0.01)){
      sub_df$SignificantBysex_F <- 4
    } else {
      sub_df$SignificantBysex_F <- 5
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  df <- do.call(rbind, list1)
  df$SignificantBysex_M <- NA
  list2 <- list()
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$fdr.men) & sub_df$method == "Inverse variance weighted" & sub_df$fdr.men < 0.05) &
             any(!is.na(sub_df$pval.men) & sub_df$method == "MR Egger" & sub_df$pval.ple.men > 0.05) &
             any(!is.na(sub_df$pval.men) & sub_df$method == "Weighted median" & sub_df$pval.men < 0.05) &
             any(!is.na(sub_df$loo.men) & sub_df$loo.men=="no") &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.men < 0.05 | is.na(sub_df$pval.men)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05 & sub_df$exp_var.men>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.men == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05)))) {
      
      sub_df$SignificantBysex_M <- 2
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05 & sub_df$exp_var.men<0.01)){
      sub_df$SignificantBysex_M <- 4
    } else {
      sub_df$SignificantBysex_M <- 5
    }
    
    list2[[paste0(exp, "_", out)]] <- sub_df
  }
  
  
  final <- do.call(rbind, list2)
  final$SignificantBysex_F <- as.numeric(final$SignificantBysex_F)
  final$SignificantBysex <- NA
  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M == 2, 3, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M==5, 1, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_M == 2 & final$SignificantBysex_F==5, 2, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==5, 5, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==2, 4, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==1 & final$SignificantBysex_M==4, 4, final$SignificantBysex)  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==5, 5, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==4, 5, final$SignificantBysex)  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==4, 5, final$SignificantBysex) 
  
  final$SignificantBysex_F <- NULL
  final$SignificantBysex_M <- NULL
  final <- final[order(final$SignificantBysex), ]
  rownames(final) <- NULL
  
  return(final)
}

# function to add sex-specific column in main results (FDR corrected)
process_tables_truly_sex_specific <- function(df, col_name = "SexSpecific") {

  list1 <- list()
  v <- c("Inverse variance weighted", "Weighted median", "MR Egger")
  combinations <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  df[[col_name]] <- "no"
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if (((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$fdr.women) & sub_df$method == "Inverse variance weighted" & sub_df$fdr.women < 0.05) &
             any(!is.na(sub_df$pval.ple.women) & sub_df$method == "MR Egger" & sub_df$pval.ple.women > 0.05) &
             any(!is.na(sub_df$loo.women) & sub_df$loo.women == "no") &
             any(!is.na(sub_df$pval.women) & sub_df$method == "Weighted median" & sub_df$pval.women < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.women < 0.05 | is.na(sub_df$pval.women))) | !("MR-PRESSO" %in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.women) & sub_df$fdr.women < 0.05 & sub_df$exp_var.women >= 0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.women == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$fdr.women) & sub_df$fdr.women < 0.05)))) & (
                (any("Inverse variance weighted" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Inverse variance weighted" ,"fdr.men"]) >= 0.05 )|
                 (any("Wald ratio" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Wald ratio","fdr.men"]) >= 0.05))))) {
      sub_df[[col_name]] <-"women-only" 
  
    } else if (((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$fdr.men) & sub_df$method == "Inverse variance weighted" & sub_df$fdr.men < 0.05) &
             any(!is.na(sub_df$pval.ple.men) & sub_df$method == "MR Egger" & sub_df$pval.ple.men > 0.05) &
             any(!is.na(sub_df$loo.men) & sub_df$loo.men == "no") &
             any(!is.na(sub_df$pval.men) & sub_df$method == "Weighted median" & sub_df$pval.men < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.men < 0.05 | is.na(sub_df$pval.men))) | !("MR-PRESSO" %in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05 & sub_df$exp_var.men >= 0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.men == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05))) )
        & (any("Inverse variance weighted" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Inverse variance weighted" ,"fdr.women"]) >= 0.05 )|
           (any("Wald ratio" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Wald ratio","fdr.women"]) >= 0.05)))) {
      
      sub_df[[col_name]] <-"men-only" 
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
 

  rownames(final) <- NULL
  return(final)
}

# function to create significance column on bidirectional results (not FDR corrected)
process_tables_bidirectional <- function(df) {
  df$SignificantBysex_F <- NA
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
         all(any(!is.na(sub_df$pval.women) & sub_df$method == "Inverse variance weighted" & !is.na(sub_df$pval.women) & sub_df$pval.women < 0.05) &
             any(!is.na(sub_df$loo.women) & sub_df$loo.women == "no") &
             any(!is.na(sub_df$pval.ple.women) & sub_df$method == "MR Egger" & sub_df$pval.ple.women > 0.05) &
             any(!is.na(sub_df$pval.women) & sub_df$method == "Weighted median" & sub_df$pval.women < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.women < 0.05 | is.na(sub_df$pval.women)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women) & sub_df$pval.women < 0.05 & sub_df$exp_var.women>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.women == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.women) & sub_df$pval.women < 0.05)))) {
      
      sub_df$SignificantBysex_F <- 1
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women) & sub_df$pval.women < 0.05 & sub_df$exp_var.women<0.01)){
      sub_df$SignificantBysex_F <- 4
    } else {
      sub_df$SignificantBysex_F <- 5
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  df <- do.call(rbind, list1)
  df$SignificantBysex_M <- NA
  list2 <- list()
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$pval.men) & sub_df$method == "Inverse variance weighted" & !is.na(sub_df$pval.men) & sub_df$pval.men < 0.05) &
             any(!is.na(sub_df$pval.men) & sub_df$method == "MR Egger" & sub_df$pval.ple.men > 0.05) &
             any(!is.na(sub_df$loo.men) & sub_df$loo.men == "no") &
             any(!is.na(sub_df$pval.men) & sub_df$method == "Weighted median" & sub_df$pval.men < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.men < 0.05 | is.na(sub_df$pval.men)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.men) & sub_df$pval.men < 0.05 & sub_df$exp_var.men>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.men == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.men) & sub_df$pval.men < 0.05)))) {
      
      sub_df$SignificantBysex_M <- 2
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$fdr.men) & sub_df$fdr.men < 0.05 & sub_df$exp_var.men<0.01)){
      sub_df$SignificantBysex_M <- 4
    } else {
      sub_df$SignificantBysex_M <- 5
    }
    
    list2[[paste0(exp, "_", out)]] <- sub_df
  }
  
  
  final <- do.call(rbind, list2)
  final$SignificantBysex_F <- as.numeric(final$SignificantBysex_F)
  final$SignificantBysex <- NA
  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M == 2, 3, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M==5, 1, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_M == 2 & final$SignificantBysex_F==5, 2, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==5, 5, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==2, 4, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==1 & final$SignificantBysex_M==4, 4, final$SignificantBysex)  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==5, 5, final$SignificantBysex)
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==4, 5, final$SignificantBysex)  
  final$SignificantBysex <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==4, 5, final$SignificantBysex) 
  
  final$SignificantBysex_F <- NULL
  final$SignificantBysex_M <- NULL
  final <- final[order(final$SignificantBysex), ]
  rownames(final) <- NULL
  
  return(final)
}

# function to add significance column in bidirectional results
process_tables_after_bidirectional <- function(df) {
  
  list1 <- list()
  
  combinations <-  df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  df$SignificantBysex_new <-NA
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if (all(sub_df$SignificantBysex == 5)){
      sub_df$SignificantBysex_new <-6 
    } 
    
    if (all(any(sub_df$SignificantBysex == 1 | sub_df$SignificantBysex == 3 | sub_df$SignificantBysex == 4) & any(sub_df$SignificantBysex_bid == 1 | sub_df$SignificantBysex_bid == 3))){
      sub_df$SignificantBysex_new <-5
    } else if (all(any(sub_df$SignificantBysex == 2 | sub_df$SignificantBysex == 3 | sub_df$SignificantBysex == 4) & any(sub_df$SignificantBysex_bid == 2 | sub_df$SignificantBysex_bid == 3))){
      sub_df$SignificantBysex_new <-5
    } else {
      sub_df$SignificantBysex_new <-  sub_df$SignificantBysex
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  rownames(final) <- NULL
  return(final)
}

# function to modify SignificantBysex column numbers in the Bidirectional MR table
column_bidirectional <- function(df) {
  
  list1 <- list()
  
  combinations <-  selected_pairs <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  df$SignificantBysex_bid <-NA
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if (all(sub_df$SignificantBysex == 1 & sub_df$SignificantBysex_main==3)){
      sub_df$SignificantBysex_bid <- "3_1"
    } else if (all(sub_df$SignificantBysex == 2 & sub_df$SignificantBysex_main==3)){
      sub_df$SignificantBysex_bid <- "3_2"
    } else if (all(sub_df$SignificantBysex == 3 & sub_df$SignificantBysex_main==3)){
      sub_df$SignificantBysex_bid <- "3_3"
    } else if (all(sub_df$SignificantBysex == 1 & sub_df$SignificantBysex_main==4)){
      sub_df$SignificantBysex_bid <- "4_1"
    } else if (all(sub_df$SignificantBysex == 2 & sub_df$SignificantBysex_main==4)){
      sub_df$SignificantBysex_bid <- "4_2"
    } else if (all(sub_df$SignificantBysex == 3 & sub_df$SignificantBysex_main==4)){
      sub_df$SignificantBysex_bid <- "4_3"
    } else {
      sub_df$SignificantBysex_bid <- sub_df$SignificantBysex
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  rownames(final) <- NULL
  return(final)
}

# function to create SignificantBysex.all in Lifelines table
process_tables_lifelines_all <- function(df) {
  df$SignificantBysex_F <- NA
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
         all(any(!is.na(sub_df$pval.women.all) & sub_df$method == "Inverse variance weighted" & sub_df$pval.women.all < 0.05) &
             any(!is.na(sub_df$pval.ple.women.all) & sub_df$method == "MR Egger" & sub_df$pval.ple.women.all > 0.05) &
             any(!is.na(sub_df$loo.women.all) & sub_df$loo.women.all=="no") &
             any(!is.na(sub_df$pval.women.all) & sub_df$method == "Weighted median" & sub_df$pval.women.all < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.women.all < 0.05 | is.na(sub_df$pval.women.all)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women.all) & sub_df$pval.women.all < 0.05 & sub_df$exp_var.women.all>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.women.all == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.women.all) & sub_df$pval.women.all < 0.05)))) {
      
      sub_df$SignificantBysex_F <- 1
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women.all) & sub_df$pval.women.all < 0.05 & sub_df$exp_var.women.all<0.01)){
      sub_df$SignificantBysex_F <- 4
    } else {
      sub_df$SignificantBysex_F <- 5
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  df <- do.call(rbind, list1)
  df$SignificantBysex_M <- NA
  list2 <- list()
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$pval.men.all) & sub_df$method == "Inverse variance weighted" & sub_df$pval.men.all < 0.05) &
             any(!is.na(sub_df$pval.ple.men.all) & sub_df$method == "MR Egger" & sub_df$pval.ple.men.all > 0.05) &
             any(!is.na(sub_df$loo.men.all) & sub_df$loo.men.all=="no") &
             any(!is.na(sub_df$pval.men.all) & sub_df$method == "Weighted median" & sub_df$pval.men.all < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.men.all < 0.05 | is.na(sub_df$pval.men.all)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.men.all) & sub_df$pval.men.all < 0.05 & sub_df$exp_var.men.all>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.men.all == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.men.all) & sub_df$pval.men.all < 0.05)))) {
      
      sub_df$SignificantBysex_M <- 2
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.men.all) & sub_df$pval.men.all < 0.05 & sub_df$exp_var.men.all<0.01)){
      sub_df$SignificantBysex_M <- 4
    } else {
      sub_df$SignificantBysex_M <- 5
    }
    
    list2[[paste0(exp, "_", out)]] <- sub_df
  }
  
  
  final <- do.call(rbind, list2)
  final$SignificantBysex_F <- as.numeric(final$SignificantBysex_F)
  final$SignificantBysex.all <- NA
  
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M == 2, 3, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M==5, 1, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_M == 2 & final$SignificantBysex_F==5, 2, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==5, 5, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==2, 4, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==1 & final$SignificantBysex_M==4, 4, final$SignificantBysex.all)  
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==5, 5, final$SignificantBysex.all)
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==4, 5, final$SignificantBysex.all)  
  final$SignificantBysex.all <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==4, 5, final$SignificantBysex.all) 
  
  final$SignificantBysex_F <- NULL
  final$SignificantBysex_M <- NULL
  
  rownames(final) <- NULL
  
  return(final)
}

# function to create SignificantBysex.nostatins in Lifelines table
process_tables_lifelines_nostatins <- function(df) {
  df$SignificantBysex_F <- NA
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
         all(any(!is.na(sub_df$pval.women.nostatins) & sub_df$method == "Inverse variance weighted" & sub_df$pval.women.nostatins < 0.05) &
             any(!is.na(sub_df$pval.ple.women.nostatins) & sub_df$method == "MR Egger" & sub_df$pval.ple.women.nostatins > 0.05) &
             any(!is.na(sub_df$loo.women.nostatins) & sub_df$loo.women.nostatins=="no") &
             any(!is.na(sub_df$pval.women.nostatins) & sub_df$method == "Weighted median" & sub_df$pval.women.nostatins < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.women.nostatins < 0.05 | is.na(sub_df$pval.women.nostatins)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women.nostatins) & sub_df$pval.women.nostatins < 0.05 & sub_df$exp_var.women.nostatins>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.women.nostatins == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.women.nostatins) & sub_df$pval.women.nostatins < 0.05)))) {
      
      sub_df$SignificantBysex_F <- 1
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.women.nostatins) & sub_df$pval.women.nostatins < 0.05 & sub_df$exp_var.women.nostatins<0.01)){
      sub_df$SignificantBysex_F <- 4
    } else {
      sub_df$SignificantBysex_F <- 5
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  df <- do.call(rbind, list1)
  df$SignificantBysex_M <- NA
  list2 <- list()
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if ((all(v %in% sub_df$method) & 
         all(any(!is.na(sub_df$pval.men.nostatins) & sub_df$method == "Inverse variance weighted" & sub_df$pval.men.nostatins < 0.05) &
             any(!is.na(sub_df$pval.ple.men.nostatins) & sub_df$method == "MR Egger" & sub_df$pval.ple.men.nostatins > 0.05) &
             any(!is.na(sub_df$loo.men.nostatins) & sub_df$loo.men.nostatins=="no") &
             any(!is.na(sub_df$pval.men.nostatins) & sub_df$method == "Weighted median" & sub_df$pval.men.nostatins < 0.05) &
             (any(sub_df$method == "MR-PRESSO" & (sub_df$pval.men.nostatins < 0.05 | is.na(sub_df$pval.men.nostatins)))| !("MR-PRESSO"%in% sub_df$method)))) |
        ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.men.nostatins) & sub_df$pval.men.nostatins < 0.05 & sub_df$exp_var.men.nostatins>=0.01)) |
        (("Inverse variance weighted" %in% sub_df$method & 
          any(sub_df$nsnp.men.nostatins == 2 & sub_df$method == "Inverse variance weighted" & 
              !is.na(sub_df$pval.men.nostatins) & sub_df$pval.men.nostatins< 0.05)))) {
      
      sub_df$SignificantBysex_M <- 2
    } else if("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df$pval.men.nostatins) & sub_df$pval.men.nostatins < 0.05 & sub_df$exp_var.men.nostatins<0.01)){
      sub_df$SignificantBysex_M <- 4
    } else {
      sub_df$SignificantBysex_M <- 5
    }
    
    list2[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list2)
  final$SignificantBysex_F <- as.numeric(final$SignificantBysex_F)
  final$SignificantBysex.nostatins <- NA
  
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M == 2, 3, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F == 1 & final$SignificantBysex_M==5, 1, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_M == 2 & final$SignificantBysex_F==5, 2, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==5, 5, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==2, 4, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==1 & final$SignificantBysex_M==4, 4, final$SignificantBysex.nostatins)  
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==5, 5, final$SignificantBysex.nostatins)
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==5 & final$SignificantBysex_M==4, 5, final$SignificantBysex.nostatins)  
  final$SignificantBysex.nostatins <- ifelse(final$SignificantBysex_F==4 & final$SignificantBysex_M==4, 5, final$SignificantBysex.nostatins) 
  
  final$SignificantBysex_F <- NULL
  final$SignificantBysex_M <- NULL
  #final <- final[order(final$SignificantBysex.nostatins), ]
  rownames(final) <- NULL
  
  return(final)
}

# function to create SexSpecific all or nostatins in Lifelines table
process_tables_truly_sex_specific_lifelines <- function(df, col_name = "SexSpecific", suffix = "") {
  # Colonne con suffisso dinamico
  fdr.women <- paste0("pval.women", suffix)
  pval.ple.women <- paste0("pval.ple.women", suffix)
  pval.women <- paste0("pval.women", suffix)
  loo.women <- paste0("loo.women", suffix)
  exp_var.women <- paste0("exp_var.women", suffix)
  nsnp.women <- paste0("nsnp.women", suffix)
  
  fdr.men <- paste0("pval.men", suffix)
  pval.ple.men <- paste0("pval.ple.men", suffix)
  loo.men <- paste0("loo.men", suffix)
  pval.men <- paste0("pval.men", suffix)
  exp_var.men <- paste0("exp_var.men", suffix)
  nsnp.men <- paste0("nsnp.men", suffix)
  
  
  list1 <- list()
  v <- c("Inverse variance weighted", "Weighted median", "MR Egger")
  combinations <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  
  df[[col_name]] <- "no"
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if (((all(v %in% sub_df$method) & 
          all(any(!is.na(sub_df[[fdr.women]]) & sub_df$method == "Inverse variance weighted" & sub_df[[fdr.women]] < 0.05) &
              any(!is.na(sub_df[[pval.ple.women]]) & sub_df$method == "MR Egger" & sub_df[[pval.ple.women]] > 0.05) &
              any(!is.na(sub_df[[loo.women]]) & sub_df[[loo.women]] == "no") &
              any(!is.na(sub_df[[pval.women]]) & sub_df$method == "Weighted median" & sub_df[[pval.women]] < 0.05) &
              (any(sub_df$method == "MR-PRESSO" & (sub_df[[pval.women]] < 0.05 | is.na(sub_df[[pval.women]]))) | !("MR-PRESSO" %in% sub_df$method)))) |
         ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df[[fdr.women]]) & sub_df[[fdr.women]] < 0.05 & sub_df[[exp_var.women]] >= 0.01)) |
         (("Inverse variance weighted" %in% sub_df$method & 
           any(sub_df[[nsnp.women]] == 2 & sub_df$method == "Inverse variance weighted" & 
               !is.na(sub_df[[fdr.women]]) & sub_df[[fdr.women]] < 0.05)))) & (
                 (any("Inverse variance weighted" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Inverse variance weighted" ,fdr.men]) >= 0.05 )|
                  (any("Wald ratio" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Wald ratio",fdr.men]) >= 0.05))))) {
      sub_df[[col_name]] <-"women-only" 
      
    } else if (((all(v %in% sub_df$method) & 
                 all(any(!is.na(sub_df[[fdr.men]]) & sub_df$method == "Inverse variance weighted" & sub_df[[fdr.men]] < 0.05) &
                     any(!is.na(sub_df[[pval.ple.men]]) & sub_df$method == "MR Egger" & sub_df[[pval.ple.men]] > 0.05) &
                     any(!is.na(sub_df[[loo.men]]) & sub_df[[loo.men]] == "no") &
                     any(!is.na(sub_df[[pval.men]]) & sub_df$method == "Weighted median" & sub_df[[pval.men]] < 0.05) &
                     (any(sub_df$method == "MR-PRESSO" & (sub_df[[pval.men]] < 0.05 | is.na(sub_df[[pval.men]]))) | !("MR-PRESSO" %in% sub_df$method)))) |
                ("Wald ratio" %in% sub_df$method & any(sub_df$method == "Wald ratio" & !is.na(sub_df[[fdr.men]]) & sub_df[[fdr.men]] < 0.05 & sub_df[[exp_var.men]] >= 0.01)) |
                (("Inverse variance weighted" %in% sub_df$method & 
                  any(sub_df[[nsnp.men]] == 2 & sub_df$method == "Inverse variance weighted" & 
                      !is.na(sub_df[[fdr.men]]) & sub_df[[fdr.men]] < 0.05))) )
               & (any("Inverse variance weighted" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Inverse variance weighted" ,fdr.women]) >= 0.05 )|
                  (any("Wald ratio" %in% sub_df$method & na.omit(sub_df[sub_df$method == "Wald ratio",fdr.women]) >= 0.05)))) {
      
      sub_df[[col_name]] <-"men-only" 
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  
  rownames(final) <- NULL
  return(final)
}

# Function for Cochran's Q with suffix
calc_cochran_Q <- function(df, b_women_col, se_women_col, b_men_col, se_men_col, suffix) {
  b_women <- df[[b_women_col]]
  se_women <- df[[se_women_col]]
  b_men <- df[[b_men_col]]
  se_men <- df[[se_men_col]]
  
  w_women <- 1 / se_women^2
  w_men   <- 1 / se_men^2
  
  diff_beta <- abs(b_women) - abs(b_men)
  beta_weighted <- (b_women * w_women + b_men * w_men) / (w_women + w_men)
  
  Q_pvalue <- 1 - pchisq(
    w_women * (b_women - beta_weighted)^2 +
      w_men * (b_men - beta_weighted)^2,
    df = 1
  )
  
  df[[paste0("Diff Beta WOMEN-MEN.", suffix)]] <- diff_beta
  df[[paste0("Q COCHRAN_WEIGHTED BETA.", suffix)]] <- beta_weighted
  df[[paste0("Q_Cochran_pvalue.", suffix)]] <- Q_pvalue
  
  return(df)
}

# Function to add the columns "Consistency_GLGC" which says if GLGC and Lifelines results are consistent
column_GLGC_sexspecific <- function(df) {
  
  list1 <- list()
  
  combinations <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  
  df$Consistency_GLGC_all_sexspecific <- NA
  df$Consistency_GLGC_nostatins_sexspecific <- NA
  df$Consistency_GLGC_Lifelines_sexspecific <- NA
  df$Consistency_QCochran_all_sexspecific <- NA
  df$Consistency_QCochran_nostatins_sexspecific <- NA
  df$Consistency_QCochran_Lifelines_sexspecific <- NA
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    ivw_row <- sub_df[sub_df$method %in% c("Inverse variance weighted", "Wald ratio"), ]
    
    if (
      all(sub_df$SexSpecific.all == sub_df$SexSpecific.nostatins) &
      (
        (
          all(sub_df$SexSpecific.nostatins == "women-only") &
          all(sub_df$SexSpecific %in% c("women-only")) &
          any(sign(na.omit(ivw_row$b.women.nostatins) )== sign(na.omit(ivw_row$b.women))) &
          any(sign(na.omit(ivw_row$b.women.all)) == sign(na.omit(ivw_row$b.women)))
        ) |
        (
          all(sub_df$SexSpecific.nostatins == "men-only") &
          all(sub_df$SexSpecific %in% c("men-only")) &
          any(sign(na.omit(ivw_row$b.men.nostatins)) == sign(na.omit(ivw_row$b.men))) &
          any(sign(na.omit(ivw_row$b.men.all)) == sign(na.omit(ivw_row$b.men)))
        )
      )
    ) {
      sub_df$Consistency_GLGC_Lifelines_sexspecific <- "yes"
    } else {
      sub_df$Consistency_GLGC_Lifelines_sexspecific <- "no"
    }
    
    if (
      all(sub_df$SexSpecific.all == sub_df$SexSpecific) &
      (
        (
          all(sub_df$SexSpecific.all ==  "women-only") &
          any(sign(na.omit(ivw_row$b.women.all)) == sign(na.omit(ivw_row$b.women)))
        ) |
        (
          all(sub_df$SexSpecific.all == "men-only") &
          any(sign(na.omit(ivw_row$b.men.all)) == sign(na.omit(ivw_row$b.men)))
        )
      )
    ) {
      sub_df$Consistency_GLGC_all_sexspecific <- "yes"
    } else {
      sub_df$Consistency_GLGC_all_sexspecific <- "no"
    }
    
    if (
      all(sub_df$SexSpecific.nostatins == sub_df$SexSpecific) &
      (
        (
          all(sub_df$SexSpecific.nostatins ==  "women-only") &
          any(sign(na.omit(ivw_row$b.women.nostatins)) == sign(na.omit(ivw_row$b.women)))
        ) |
        (
          all(sub_df$SexSpecific.nostatins ==  "men-only") &
          any(sign(na.omit(ivw_row$b.men.nostatins)) == sign(na.omit(ivw_row$b.men)))
        )
      )
    ) {
      sub_df$Consistency_GLGC_nostatins_sexspecific <- "yes"
    } else {
      sub_df$Consistency_GLGC_nostatins_sexspecific <- "no"
    }
    
    if (all(!is.na(ivw_row$Q_Cochran_pvalue.all) & !is.na(ivw_row$Q_Cochran_pvalue))) {
      if (ivw_row$Q_Cochran_pvalue.all < 0.05 & ivw_row$Q_Cochran_pvalue < 0.05) {
        sub_df$Consistency_QCochran_all_sexspecific <- "yes"
      } else {
        sub_df$Consistency_QCochran_all_sexspecific <- "no"
      }
    }
    
    #### Q Cochran NOSTATINS
    if (all(!is.na(ivw_row$Q_Cochran_pvalue.nostatins) & !is.na(ivw_row$Q_Cochran_pvalue))) {
      if (ivw_row$Q_Cochran_pvalue.nostatins < 0.05 & ivw_row$Q_Cochran_pvalue < 0.05) {
        sub_df$Consistency_QCochran_nostatins_sexspecific <- "yes"
      } else {
        sub_df$Consistency_QCochran_nostatins_sexspecific <- "no"
      }
    }
    
    #### Q Cochran LIFELINES
    if (all(
      !is.na(ivw_row$Q_Cochran_pvalue.all) &
      !is.na(ivw_row$Q_Cochran_pvalue.nostatins) &
      !is.na(ivw_row$Q_Cochran_pvalue)
    )) {
      if (
        ivw_row$Q_Cochran_pvalue.all < 0.05 &
        ivw_row$Q_Cochran_pvalue.nostatins < 0.05 &
        ivw_row$Q_Cochran_pvalue < 0.05
      ) {
        sub_df$Consistency_QCochran_Lifelines_sexspecific <- "yes"
      } else {
        sub_df$Consistency_QCochran_Lifelines_sexspecific <- "no"
      }
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  rownames(final) <- NULL
  return(final)
}





