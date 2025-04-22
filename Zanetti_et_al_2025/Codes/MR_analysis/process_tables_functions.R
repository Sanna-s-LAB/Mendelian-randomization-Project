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

# function to add significance column in main results (FDR corrected)
process_tables <- function(df) {
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
         all(any(!is.na(sub_df$fdr.women) & sub_df$method == "Inverse variance weighted" & sub_df$fdr.women < 0.05) &
             any(!is.na(sub_df$pval.ple.women) & sub_df$method == "MR Egger" & sub_df$pval.ple.women > 0.05) &
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
  
  combinations <-  selected_pairs <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  df$SignificantBysex_CISTRANS_new <-NA
  df$SignificantBysex_CIS_new <-NA
  df$SignificantBysex_TRANS_new <-NA
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    if (all(sub_df$SignificantBysex_CISTRANS == 5)){
      sub_df$SignificantBysex_CISTRANS_new <-6 
    } 
    if(all(!is.na(sub_df$SignificantBysex_CIS) & all(sub_df$SignificantBysex_CIS == 5))){
      sub_df$SignificantBysex_CIS_new <-6
    } 
    if (all(!is.na(sub_df$SignificantBysex_TRANS) & all(sub_df$SignificantBysex_TRANS == 5))){
      sub_df$SignificantBysex_TRANS_new <- 6
    }
    
    if (all(any(sub_df$SignificantBysex_CISTRANS == 1 | sub_df$SignificantBysex_CISTRANS == 3 | sub_df$SignificantBysex_CISTRANS == 4) & any(sub_df$SignificantBysex_bid == 1 | sub_df$SignificantBysex_bid == 3))){
      sub_df$SignificantBysex_CISTRANS_new <-5
      sub_df$SignificantBysex_CIS_new <-5
      sub_df$SignificantBysex_TRANS_new <-5
    } else if (all(any(sub_df$SignificantBysex_CISTRANS == 2 | sub_df$SignificantBysex_CISTRANS == 3 | sub_df$SignificantBysex_CISTRANS == 4) & any(sub_df$SignificantBysex_bid == 2 | sub_df$SignificantBysex_bid == 3))){
      sub_df$SignificantBysex_CISTRANS_new <-5
      sub_df$SignificantBysex_CIS_new <-5
      sub_df$SignificantBysex_TRANS_new <-5
    } else {
      sub_df$SignificantBysex_CISTRANS_new <-  sub_df$SignificantBysex_CISTRANS
      sub_df$SignificantBysex_CIS_new <-  sub_df$SignificantBysex_CIS
      sub_df$SignificantBysex_TRANS_new <-  sub_df$SignificantBysex_TRANS
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  
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
  
  return(final)
}

# function to select couples to remove after leave-one-out
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
  
  # Convert the list to a data frame
  list_to_keep <- do.call(rbind, to_keep)
  colnames(list_to_keep) <- c("exposure", "outcome")
  list_to_keep <- as.data.frame(list_to_keep)
  
  toremove_W<- anti_join(not_sig,list_to_keep, by=c("exposure", "outcome"))
  
  return(list(keep = list_to_keep, remove = toremove_W))
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

# Function to add the columns "Consistency_GLGC" which says if GLGC and Lifelines results are consistent
column_GLGC <- function(df) {
  
  list1 <- list()
  
  combinations <- df %>%
    dplyr::select("exposure", "outcome") %>%
    distinct()
  
  df$Consistency_GLGC_all <- NA
  df$Consistency_GLGC_nostatins <- NA
  df$Consistency_GLGC_Lifelines <- NA
  
  for (i in 1:nrow(combinations)) {
    exp <- combinations$exposure[i]
    out <- combinations$outcome[i]
    sub_df <- subset(df, exposure == exp & outcome == out)
    
    ivw_row <- sub_df[sub_df$method %in% c("Inverse variance weighted", "Wald ratio"), ]
    
     if (
        all(sub_df$SignificantBysex.all == sub_df$SignificantBysex.nostatins) &
        (
          (
            all(sub_df$SignificantBysex.nostatins == 1) &
            all(sub_df$SignificantBysex_CISTRANS %in% c(1)) &
            any(sign(ivw_row$b.women.nostatins) == sign(ivw_row$b.women)) &
            any(sign(ivw_row$b.women.all) == sign(ivw_row$b.women))
          ) |
          (
            all(sub_df$SignificantBysex.nostatins == 2) &
            all(sub_df$SignificantBysex_CISTRANS %in% c(2)) &
            any(sign(ivw_row$b.men.nostatins) == sign(ivw_row$b.men)) &
            any(sign(ivw_row$b.men.all) == sign(ivw_row$b.men))
          )
        )
      ) {
        sub_df$Consistency_GLGC_Lifelines <- "yes"
      } else {
        sub_df$Consistency_GLGC_Lifelines <- "no"
      }
    
    if (
      all(sub_df$SignificantBysex.all == sub_df$SignificantBysex_CISTRANS) &
      (
        (
          all(sub_df$SignificantBysex.all == 1) &
          any(sign(ivw_row$b.women.all) == sign(ivw_row$b.women))
        ) |
        (
          all(sub_df$SignificantBysex.all == 2) &
          any(sign(ivw_row$b.men.all) == sign(ivw_row$b.men))
        )
      )
    ) {
      sub_df$Consistency_GLGC_all <- "yes"
    } else {
      sub_df$Consistency_GLGC_all <- "no"
    }
    
    if (
      all(sub_df$SignificantBysex.nostatins == sub_df$SignificantBysex_CISTRANS) &
      (
        (
          all(sub_df$SignificantBysex.nostatins == 1) &
          any(sign(ivw_row$b.women.nostatins) == sign(ivw_row$b.women))
        ) |
        (
          all(sub_df$SignificantBysex.nostatins == 2) &
          any(sign(ivw_row$b.men.nostatins) == sign(ivw_row$b.men))
        )
      )
    ) {
      sub_df$Consistency_GLGC_nostatins <- "yes"
    } else {
      sub_df$Consistency_GLGC_nostatins <- "no"
    }
    
    list1[[paste0(exp, "_", out)]] <- sub_df
  }
  
  final <- do.call(rbind, list1)
  rownames(final) <- NULL
  return(final)
}




