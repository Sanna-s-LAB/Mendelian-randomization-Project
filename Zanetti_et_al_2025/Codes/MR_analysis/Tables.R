#!/usr/bin/env Rscript

library(plyr) 
library(dplyr)
library(openxlsx)
library(tidyverse)
library(readxl)
library(tidyr)
library(readr)
source("~/Scripts_for_github/process_tables_functions.R")

##### CIS+TRANS - FINAL RESULT FOR PAPER - TABLES #####

###### MEN ######

# Upload file
directory <- "~/Men/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men",
                         "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men","IV_list.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" &
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) &
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()

# FDR correction 

result_IVW <- new_males %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_males %>%
  filter(method == "Wald ratio")

result_Egger <- new_males %>%
  filter(method == "MR Egger")

result_Egger$fdr.men <- c("NA")

result_WM <- new_males %>%
  filter(method == "Weighted median")

result_WM$fdr.men <-  c("NA")

result_PRESSO <- new_males %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.men <- c("NA")

names(result_IVW)

result_IVW$fdr.men<-p.adjust(result_IVW$pval.men, method="BH")

result_WR$fdr.men<-p.adjust(result_WR$pval.men, method="BH")

men_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(men_all)
men_all$fdr.men<- as.numeric(men_all$fdr.men)
men_all$pval.men[men_all$nsnp.men == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

###### WOMEN ######

# Upload file
directory <- "~/Women/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                         "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" &
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) &
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()

# FDR correction
result_IVW <- new_women %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_women %>%
  filter(method == "Wald ratio")

result_Egger <- new_women %>%
  filter(method == "MR Egger")

result_Egger$fdr.women <- c("NA")

result_WM <- new_women %>%
  filter(method == "Weighted median")

result_WM$fdr.women <-  c("NA")

result_PRESSO <- new_women %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.women <- c("NA")

names(result_IVW)

result_IVW$fdr.women<-p.adjust(result_IVW$pval.women, method="BH")

result_WR$fdr.women<-p.adjust(result_WR$pval.women, method="BH")

women_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(women_all)
women_all$fdr.women<- as.numeric(women_all$fdr.women)
women_all$pval.women[women_all$nsnp.women == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  ALL ######

# Upload file
directory <- "~/All/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_a <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_a <- merged_df_a %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_a <- merged_df_a[,-c(1,2)] # remove id
indices_na_exposure_a <- which(is.na(merged_df_a$exposure)) # check for failed analyses
new_all <- merged_df_a
colnames(new_all) <- c("outcome", "exposure", "method", "nsnp.all", "b.all","se.all","pval.all",
                       "egger_intercept.ple.all","se.ple.all","pval.ple.all","Q.het.all","Q_df.het.all","Q_pval.het.all","IV_list.all")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" &
        any(nsnp.all %in% c(1,2,3) & method != "MR-PRESSO")
    ) &
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                       Q.het.all[method == "Weighted median"], Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                          Q_df.het.all[method == "Weighted median"], Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_pval.het.all[method == "Weighted median"], Q_pval.het.all)
  ) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.all)
  ) %>%
  ungroup()


# FDR correction
result_IVW <- new_all %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_all %>%
  filter(method == "Wald ratio")

result_Egger <- new_all %>%
  filter(method == "MR Egger")

result_Egger$fdr.all <- c("NA")

result_WM <- new_all %>%
  filter(method == "Weighted median")

result_WM$fdr.all <-  c("NA")

result_PRESSO <- new_all %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.all <- c("NA")

names(result_IVW)

result_IVW$fdr.all<-p.adjust(result_IVW$pval.all, method="BH")

result_WR$fdr.all<-p.adjust(result_WR$pval.all, method="BH")

all_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(all_all)
all_all$fdr.all<- as.numeric(all_all$fdr.all)
all_all$pval.all[all_all$nsnp.all == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs


######  ST1. TABLE of ALL RESULTS (significant or not) #####

setwd("~/")
ALL_RESULTS <- merge(women_all,men_all, by=c("exposure","outcome", "method"),all = TRUE)
ALL_RESULTS <- merge(ALL_RESULTS, all_all, by=c("exposure","outcome", "method"), all=TRUE)
ALL_RESULTS <- ALL_RESULTS[order(ALL_RESULTS$exposure, ALL_RESULTS$outcome, ALL_RESULTS$method, decreasing = FALSE), ]
ALL_RESULTS <- ALL_RESULTS[!is.na(ALL_RESULTS$exposure), ]

to_keep_w <- women_all %>%
  select(exposure, outcome) %>%
  distinct()

to_keep_m <- men_all %>%
  select(exposure, outcome) %>%
  distinct()

t<- inner_join(to_keep_w,to_keep_m) # to keep only common associations

ALL_RESULTS1 <- inner_join(ALL_RESULTS,t)

write.xlsx(ALL_RESULTS1, "~/Results tables/All_results.xlsx")

# To count common associations
n<- ALL_RESULTS1 %>%
  select(exposure, outcome) %>%
  distinct()


##### CIS - FINAL RESULT FOR PAPER - TABLES #####

######  MEN #####

# Upload file
directory <- "~/Men/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_cis_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men","IV_list.men",
                         "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()

# FDR correction 
result_IVW <- new_males %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_males %>%
  filter(method == "Wald ratio")

result_Egger <- new_males %>%
  filter(method == "MR Egger")

result_Egger$fdr.men <- c("NA")

result_WM <- new_males %>%
  filter(method == "Weighted median")

result_WM$fdr.men <-  c("NA")

result_PRESSO <- new_males %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.men <- c("NA")

names(result_IVW)

result_IVW$fdr.men<-p.adjust(result_IVW$pval.men, method="BH")
result_WR$fdr.men<-p.adjust(result_WR$pval.men, method="BH")

men_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(men_all)
men_all$fdr.men<- as.numeric(men_all$fdr.men)
men_all$pval.men[men_all$nsnp.men == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  WOMEN #####

# Upload file
directory <- "~/Women/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_cis_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                         "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()

# FDR correction
result_IVW <- new_women %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_women %>%
  filter(method == "Wald ratio")

result_Egger <- new_women %>%
  filter(method == "MR Egger")

result_Egger$fdr.women <- c("NA")

result_WM <- new_women %>%
  filter(method == "Weighted median")

result_WM$fdr.women <-  c("NA")

result_PRESSO <- new_women %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.women <- c("NA")

names(result_IVW)

result_IVW$fdr.women<-p.adjust(result_IVW$pval.women, method="BH")

result_WR$fdr.women<-p.adjust(result_WR$pval.women, method="BH")

women_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(women_all)
women_all$fdr.women<- as.numeric(women_all$fdr.women)
women_all$pval.women[women_all$nsnp.women == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  ALL #####

# Upload file
directory <- "~/All/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_cis_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_a <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_a <- merged_df_a %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_a <- merged_df_a[,-c(1,2)] # remove id
indices_na_exposure_a <- which(is.na(merged_df_a$exposure)) # check for failed analyses

new_all <- merged_df_a
colnames(new_all) <- c("outcome", "exposure", "method", "nsnp.all", "b.all","se.all","pval.all",
                       "egger_intercept.ple.all","se.ple.all","pval.ple.all","Q.het.all","Q_df.het.all","Q_pval.het.all","IV_list.all")

new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" &
        any(nsnp.all %in% c(1,2,3) & method != "MR-PRESSO")
    ) &
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                       Q.het.all[method == "Weighted median"], Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                          Q_df.het.all[method == "Weighted median"], Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_pval.het.all[method == "Weighted median"], Q_pval.het.all)
  ) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.all)
  ) %>%
  ungroup()


# FDR correction 
result_IVW <- new_all %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_all %>%
  filter(method == "Wald ratio")

result_Egger <- new_all %>%
  filter(method == "MR Egger")

result_Egger$fdr.all <- c("NA")

result_WM <- new_all %>%
  filter(method == "Weighted median")

result_WM$fdr.all <-  c("NA")

result_PRESSO <- new_all %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.all <- c("NA")

names(result_IVW)

result_IVW$fdr.all<-p.adjust(result_IVW$pval.all, method="BH")

result_WR$fdr.all<-p.adjust(result_WR$pval.all, method="BH")

all_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(all_all)
all_all$fdr.all<- as.numeric(all_all$fdr.all)

all_all$pval.all[all_all$nsnp.all == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  ST2. TABLE of ALL RESULTS (significant or not) #####

setwd("~/")
ALL_RESULTS <- merge(women_all,men_all, by=c("exposure","outcome", "method"),all = TRUE)
ALL_RESULTS <- merge(ALL_RESULTS, all_all, by=c("exposure","outcome", "method"), all=TRUE)
ALL_RESULTS <- ALL_RESULTS[order(ALL_RESULTS$exposure, ALL_RESULTS$outcome, ALL_RESULTS$method, decreasing = FALSE), ]
ALL_RESULTS <- ALL_RESULTS[!is.na(ALL_RESULTS$exposure), ]

to_keep_w <- women_all %>%
  select(exposure, outcome) %>%
  distinct()

to_keep_m <- men_all %>%
  select(exposure, outcome) %>%
  distinct()

t<- inner_join(to_keep_w,to_keep_m) # to keep only common associations

ALL_RESULTS1 <- inner_join(ALL_RESULTS,t)

write.xlsx(ALL_RESULTS1, "~/Results tables/All_results_CIS.xlsx")

# To count common associations
n<- ALL_RESULTS1 %>%
  select(exposure, outcome) %>%
  distinct()


### TRANS - FINAL RESULT FOR PAPER - TABLES ####

###### MEN ######

# Upload file
directory <- "~/Men/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_trans_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men","IV_list.men",
                         "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()

# FDR correction 
result_IVW <- new_males %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_males %>%
  filter(method == "Wald ratio")

result_Egger <- new_males %>%
  filter(method == "MR Egger")

result_Egger$fdr.men <- c("NA")

result_WM <- new_males %>%
  filter(method == "Weighted median")

result_WM$fdr.men <-  c("NA")

result_PRESSO <- new_males %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.men <- c("NA")

names(result_IVW)

result_IVW$fdr.men<-p.adjust(result_IVW$pval.men, method="BH")

result_WR$fdr.men<-p.adjust(result_WR$pval.men, method="BH")

men_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(men_all)
men_all$fdr.men<- as.numeric(men_all$fdr.men)
men_all$pval.men[men_all$nsnp.men == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  WOMEN #####

# Upload file
directory <- "~/Women/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_trans_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                         "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()

# FDR correction 
result_IVW <- new_women %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_women %>%
  filter(method == "Wald ratio")

result_Egger <- new_women %>%
  filter(method == "MR Egger")

result_Egger$fdr.women <- c("NA")

result_WM <- new_women %>%
  filter(method == "Weighted median")

result_WM$fdr.women <-  c("NA")

result_PRESSO <- new_women %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.women <- c("NA")

names(result_IVW)

result_IVW$fdr.women<-p.adjust(result_IVW$pval.women, method="BH")

result_WR$fdr.women<-p.adjust(result_WR$pval.women, method="BH")

women_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(women_all)
women_all$fdr.women<- as.numeric(women_all$fdr.women)
women_all$pval.women[women_all$nsnp.women == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  ALL #####

# Upload file
directory <- "~/All/cis_trans/final_results_for_paper"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_trans_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_a <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_a <- merged_df_a %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_a <- merged_df_a[,-c(1,2)] # remove id
indices_na_exposure_a <- which(is.na(merged_df_a$exposure)) # check for failed analyses
new_all <- merged_df_a
colnames(new_all) <- c("outcome", "exposure", "method", "nsnp.all", "b.all","se.all","pval.all",
                       "egger_intercept.ple.all","se.ple.all","pval.ple.all","Q.het.all","Q_df.het.all","Q_pval.het.all","IV_list.all")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" &
        any(nsnp.all %in% c(1,2,3) & method != "MR-PRESSO")
    ) &
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_all <- new_all %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                       Q.het.all[method == "Weighted median"], Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                          Q_df.het.all[method == "Weighted median"], Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_pval.het.all[method == "Weighted median"], Q_pval.het.all)
  ) %>%
  dplyr::mutate(
    Q.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.all),
    Q_df.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.all),
    Q_pval.het.all = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.all)
  ) %>%
  ungroup()


# FDR correction
result_IVW <- new_all %>%
  filter(method == "Inverse variance weighted")

result_WR <- new_all %>%
  filter(method == "Wald ratio")

result_Egger <- new_all %>%
  filter(method == "MR Egger")

result_Egger$fdr.all <- c("NA")

result_WM <- new_all %>%
  filter(method == "Weighted median")

result_WM$fdr.all <-  c("NA")

result_PRESSO <- new_all %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.all <- c("NA")

names(result_IVW)

result_IVW$fdr.all<-p.adjust(result_IVW$pval.all, method="BH")

result_WR$fdr.all<-p.adjust(result_WR$pval.all, method="BH")

all_all <- rbind(result_IVW, result_WR,result_Egger,result_WM, result_PRESSO)

str(all_all)
all_all$fdr.all<- as.numeric(all_all$fdr.all)

all_all$pval.all[all_all$nsnp.all == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

######  ST3. TABLE of ALL RESULTS (significant or not) #####

setwd("~/")
ALL_RESULTS <- merge(women_all,men_all, by=c("exposure","outcome", "method"),all = TRUE)
ALL_RESULTS <- merge(ALL_RESULTS, all_all, by=c("exposure","outcome", "method"), all=TRUE)
ALL_RESULTS <- ALL_RESULTS[order(ALL_RESULTS$exposure, ALL_RESULTS$outcome, ALL_RESULTS$method, decreasing = FALSE), ]
ALL_RESULTS <- ALL_RESULTS[!is.na(ALL_RESULTS$exposure), ]

to_keep_w <- women_all %>%
  select(exposure, outcome) %>%
  distinct()

to_keep_m <- men_all %>%
  select(exposure, outcome) %>%
  distinct()

t<- inner_join(to_keep_w,to_keep_m) # to keep only common associations

ALL_RESULTS1 <- inner_join(ALL_RESULTS,t)

write.xlsx(ALL_RESULTS1, "~/Results tables/All_results_TRANS.xlsx")

# To count common associations
n<- ALL_RESULTS1 %>%
  select(exposure, outcome) %>%
  distinct()



##### TABLES with columns ####

setwd("~/Results tables")
proteins <- read.table("~/ProteinsNames_new.csv",header=T,as.is=T,sep=",")
proteins<-subset(proteins,select =c("Assay","UniProt","exposure","Name"))

###### CIS+TRANS #######
ALL_RESULTS <- read_excel("All_results.xlsx")

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.men))%>%
  ungroup()

men <-  read.table("~/Men/final_results_for_paper/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)

Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.women))%>%
  ungroup()

women <- read.table("~/Women/final_results_for_paper/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )


colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# All results with variance explained columns for Wald ratio results
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Process tables with process_tables function to add the "SignificantBysex" column
final2 <- process_tables(ALL_RESULTS1)
final3 <-  merge(final2, proteins, by ="exposure")
write.xlsx(final3,"~/Results tables/Final_results_with_column.xlsx")

# Select only significant results
Sig_womenORmen <- final3 %>%
  dplyr::group_by(exposure, outcome) %>%
  filter(SignificantBysex %in% c(1,2,3,4)) %>%
  ungroup()

# Add Cochran's Q values
Sig_womenORmen$`Diff Beta WOMEN-MEN` <- abs(Sig_womenORmen$b.women)-abs(Sig_womenORmen$b.men)
Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA` <- (Sig_womenORmen$b.women*(1/Sig_womenORmen$se.women^2) 
                                             + Sig_womenORmen$b.men*(1/Sig_womenORmen$se.men^2)) / 
  ((1/Sig_womenORmen$se.women^2)+(1/Sig_womenORmen$se.men^2))
Sig_womenORmen$Q_Cochran_pvalue <- 1 - pchisq((1/Sig_womenORmen$se.women^2) * 
                                                (Sig_womenORmen$b.women- Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2 + 
                                                (1/Sig_womenORmen$se.men^2) * (Sig_womenORmen$b.men - Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2, df = 1)


Sig_womenORmen <- Sig_womenORmen[order(Sig_womenORmen$exposure, Sig_womenORmen$outcome, Sig_womenORmen$method, decreasing = FALSE), ]
write.xlsx(Sig_womenORmen, "~/Results tables/Significant_women_or_men_final.xlsx")

#Save significant proteins
Significant_proteins <- unique(Sig_womenORmen$exposure)
write.csv(Significant_proteins,"~/Results tables/Significant_proteins.csv")

# Count significant results
significant_count <- Sig_womenORmen  %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==1 & method == "Wald ratio" & !is.na(fdr.women) & fdr.women < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==1 & method == "Inverse variance weighted" & !is.na(fdr.women) & fdr.women < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==2 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==2 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==3 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==3 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==4 & method == "Wald ratio")
  )
significant_count

selected_pairs <- Sig_womenORmen %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

n <- Sig_womenORmen %>%
  dplyr::filter(SignificantBysex ==1) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)


###### CIS ######

all_CIS <- read_excel("~/Results tables/All_results_CIS.xlsx")
all_CIS <- all_CIS[!is.na(all_CIS$exposure), ]

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- all_CIS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.men))%>%
  ungroup()
men <-  read.table("~/Men/final_results_for_paper/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome",
                   "remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome",
                   "outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure",
                   "mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(all_CIS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_women_Wald <- all_CIS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.women))%>%
  ungroup()
women <- read.table("~/Women/final_results_for_paper/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# Final table with explained variance
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Process tables with process_tables function to add the "SignificantBysex" column
final2 <- process_tables(ALL_RESULTS1)
final3 <-  merge(final2, proteins, by ="exposure")
write.xlsx(final3,"~/Results tables/Final_results_CIS_with_column.xlsx")

# Here we take everything, since CIS analyses have been done only on CIS+TRANS significant results
Sig_womenORmen <-final3
Sig_womenORmen$`Diff Beta WOMEN-MEN` <- abs(Sig_womenORmen$b.women)-abs(Sig_womenORmen$b.men)
Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA` <- (Sig_womenORmen$b.women*(1/Sig_womenORmen$se.women^2) 
                                             + Sig_womenORmen$b.men*(1/Sig_womenORmen$se.men^2)) / 
  ((1/Sig_womenORmen$se.women^2)+(1/Sig_womenORmen$se.men^2))
Sig_womenORmen$Q_Cochran_pvalue <- 1 - pchisq((1/Sig_womenORmen$se.women^2) * 
                                                (Sig_womenORmen$b.women- Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2 + 
                                                (1/Sig_womenORmen$se.men^2) * (Sig_womenORmen$b.men - Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2, df = 1)


Sig_womenORmen <- Sig_womenORmen[order(Sig_womenORmen$exposure, Sig_womenORmen$outcome, Sig_womenORmen$method, decreasing = FALSE), ]
write.xlsx(Sig_womenORmen, "~/Results tables/Significant_CIS_women_or_men_final.xlsx")

# Count significant results
significant_count <- Sig_womenORmen  %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==1 & method == "Wald ratio" & !is.na(fdr.women) & fdr.women < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==1 & method == "Inverse variance weighted" & !is.na(fdr.women) & fdr.women < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==2 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==2 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==3 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==3 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==4 & method == "Wald ratio")
  )
significant_count

selected_pairs <- Sig_womenORmen %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

######TRANS ##########

all_TRANS <- read_excel("~/Results tables/All_results_TRANS.xlsx")
all_TRANS <- all_TRANS[!is.na(all_TRANS$exposure), ]

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- all_TRANS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.men))%>%
  ungroup()
men <-  read.table("~/Men/final_results_for_paper/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(all_TRANS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_women_Wald <- all_TRANS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.women))%>%
  ungroup()
women <- read.table("~/Women/final_results_for_paper/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# Final table with explained variance
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Process tables with process_tables function to add the "SignificantBysex" column
final2 <- process_tables(ALL_RESULTS1)
final3 <-  merge(final2, proteins, by ="exposure")
write.xlsx(final3,"Final_results_TRANS_with_column.xlsx")

# Here we take everything, since TRANS analyses have been done only on CIS+TRANS significant results
Sig_womenORmen <-final3
Sig_womenORmen$`Diff Beta WOMEN-MEN` <- abs(Sig_womenORmen$b.women)-abs(Sig_womenORmen$b.men)
Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA` <- (Sig_womenORmen$b.women*(1/Sig_womenORmen$se.women^2) 
                                             + Sig_womenORmen$b.men*(1/Sig_womenORmen$se.men^2)) / 
  ((1/Sig_womenORmen$se.women^2)+(1/Sig_womenORmen$se.men^2))
Sig_womenORmen$Q_Cochran_pvalue <- 1 - pchisq((1/Sig_womenORmen$se.women^2) * 
                                                (Sig_womenORmen$b.women- Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2 + 
                                                (1/Sig_womenORmen$se.men^2) * (Sig_womenORmen$b.men - Sig_womenORmen$`Q COCHRAN_WEIGHTED BETA`)^2, df = 1)

Sig_womenORmen <- Sig_womenORmen[order(Sig_womenORmen$exposure, Sig_womenORmen$outcome, Sig_womenORmen$method, decreasing = FALSE), ]
write.xlsx(Sig_womenORmen, "~/Results tables/Significant_TRANS_women_or_men_final.xlsx")

# Count significant results
significant_count <- Sig_womenORmen  %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==1 & method == "Wald ratio" & !is.na(fdr.women) & fdr.women < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==1 & method == "Inverse variance weighted" & !is.na(fdr.women) & fdr.women < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==2 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==2 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==3 & method == "Wald ratio" & !is.na(fdr.men) & fdr.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex==3 & method == "Inverse variance weighted" & !is.na(fdr.men) & fdr.men < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex==4 & method == "Wald ratio")
  )
significant_count

selected_pairs <- Sig_womenORmen %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()


##### BIDIRECTIONAL MR ######

# MEN 

# Upload file
directory <- "~/Men/BidirectionalMR"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men",
                         "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men","IV_list.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()

men_all <-new_males

men_all$pval.men[men_all$nsnp.men == 0] <- 1  # set p-value to 1 when MR-PRESSO removes all SNPs

# WOMEN 

# Upload file
directory <- "~/Women/BidirectionalMR"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                         "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()

women_all <- new_women

women_all$pval.women[women_all$nsnp.women == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

ALL_RESULTS <- merge(women_all,men_all, by=c("exposure","outcome", "method"),all = TRUE)
write.xlsx(ALL_RESULTS, "~/Results tables/All_BIDIRECTIONAL_results.xlsx")

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(pval.men))%>%
  ungroup()

men <-  read.table("~/Men/BidirectionalMR/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(pval.women))%>%
  ungroup()
women <- read.table("~/Women/BidirectionalMR/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# All results with explained variance
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Process tables with process_tables_bidirectional function to add the "SignificantBysex" column
final2 <- process_tables_bidirectional(ALL_RESULTS1)

###### ST4. Final column for bidirectional MR results ######
all_sign <- read_excel("~/Results tables/Significant_women_or_men_final.xlsx")
sub <- subset(all_sign, select=c("exposure","outcome","method","SignificantBysex"))
colnames(sub)<- c("outcome","exposure","method","SignificantBysex_main")
final3 <- merge(final2, sub, by=c("exposure","outcome","method"), all.x = T)

final4 <- final3 %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex_main,
       .direction = "downup") %>%
  ungroup()

# process data to write the column in the final way
final5 <- column_bidirectional(final4)

colnames(proteins)<- c("Assay" ,   "UniProt",  "outcome" ,"Name" )
final6 <-  merge(final5, proteins, by ="outcome")

write.xlsx(final6,"~/Results tables/Final_results_BIDIR_with_3columns.xlsx")
final6$SignificantBysex <- NULL
final6$SignificantBysex_main <-NULL
write.xlsx(final6,"~/Results tables/Final_results_BIDIR_with_column.xlsx") # ST5

# Count significant results (to review)
Sig_womenORmen <- final6 %>%
  dplyr::group_by(exposure, outcome) %>%
  filter(SignificantBysex_bid %in% c(1,2,"3_1","3_2","3_3","4_1")) %>%
  ungroup()

Sig_womenORmen <- Sig_womenORmen[order(Sig_womenORmen$exposure, Sig_womenORmen$outcome, Sig_womenORmen$method, decreasing = FALSE), ]

significant_count <- Sig_womenORmen  %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid==1 & method == "Wald ratio" & !is.na(pval.women) & pval.women < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid==1 & method == "Inverse variance weighted" & !is.na(pval.women) & pval.women < 0.05)
  )
significant_count

significant_count <- Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid==2 & method == "Wald ratio" & !is.na(pval.men) & pval.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid==2 & method == "Inverse variance weighted" & !is.na(pval.men) & pval.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid=="3_1" & method == "Wald ratio" & !is.na(pval.men) & pval.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid=="3_1" & method == "Inverse variance weighted" & !is.na(pval.men) & pval.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid=="3_2" & method == "Wald ratio" & !is.na(pval.men) & pval.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid=="3_2" & method == "Inverse variance weighted" & !is.na(pval.men) & pval.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid=="3_3" & method == "Wald ratio" & !is.na(pval.men) & pval.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid=="3_3" & method == "Inverse variance weighted" & !is.na(pval.men) & pval.men < 0.05)
  )
significant_count

significant_count <-Sig_womenORmen %>%
  dplyr::summarise(
    wald_ratio_significant_count = sum(SignificantBysex_bid=="4_1" & method == "Wald ratio" & !is.na(pval.men) & pval.men < 0.05),
    IVW_ratio_significant_count = sum(SignificantBysex_bid=="4_1" & method == "Inverse variance weighted" & !is.na(pval.men) & pval.men < 0.05)
  )
significant_count

##### FINAL SIGNIFICANT CIS+TRANS, CIS, TRANS ####

setwd("~/Results tables")
all_sign <- read_excel("Significant_women_or_men_final.xlsx")
all_CIS <- read_excel("Significant_CIS_women_or_men_final.xlsx") 
all_TRANS <- read_excel("Significant_TRANS_women_or_men_final.xlsx")  
proteins <- read.table("~/ProteinsNames_new.csv",header=T,as.is=T,sep=",")
proteins<-subset(proteins,select =c("Assay","UniProt","exposure","Name"))

all_sign <- all_sign[,-c(8,9,11,12,13,14,20,21,23,24,25,26,
                         28,29,30,31,32,33,34,35,36,37,38,39,43,44,45)]
colnames(all_sign) <- c("exposure","outcome" ,"method","nsnp.women", "b.women","se.women","pval.women","pval.ple.women","fdr.women" 
                        ,"nsnp.men","b.men", "se.men","pval.men","pval.ple.men","fdr.men","exp_var.men","exp_var.women","SignificantBysex_CISTRANS"
                        ,"Diff_Beta_WOMEN_MEN_CISTRANS","Q_COCHRAN_WEIGHTED_BETA_CISTRANS","Q_Cochran_pvalue_CISTRANS")
all_sign <- all_sign[!(is.na(all_sign$nsnp.men) & is.na(all_sign$nsnp.women)), ]

all_CIS <- all_CIS[,-c(8,9,11,12,13,14,20,21,22,24,25,26,
                       28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,44,45)]
all_CIS <- all_CIS[!(is.na(all_CIS$nsnp.men) & is.na(all_CIS$nsnp.women)), ]
colnames(all_CIS) <- c("exposure","outcome" ,"method","nsnp.women.CIS", "b.women.CIS","se.women.CIS","pval.women.CIS","pval.ple.women.CIS","fdr.women.CIS" 
                       ,"nsnp.men.CIS","b.men.CIS", "se.men.CIS","pval.men.CIS","pval.ple.men.CIS","fdr.men.CIS","SignificantBysex_CIS","Diff_Beta_WOMEN_MEN_CIS","Q_COCHRAN_WEIGHTED_BETA_CIS","Q_Cochran_pvalue_CIS")

all_TRANS <- all_TRANS[,-c(8,9,11,12,13,14,20,21,22,24,25,26,
                           28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,44,45)]
all_TRANS <- all_TRANS[!(is.na(all_TRANS$nsnp.men) & is.na(all_TRANS$nsnp.women)), ]
colnames(all_TRANS) <- c("exposure","outcome" ,"method","nsnp.women.TRANS", "b.women.TRANS","se.women.TRANS","pval.women.TRANS","pval.ple.women.TRANS","fdr.women.TRANS" 
                         ,"nsnp.men.TRANS","b.men.TRANS", "se.men.TRANS","pval.men.TRANS","pval.ple.men.TRANS","fdr.men.TRANS","SignificantBysex_TRANS","Diff_Beta_WOMEN_MEN_TRANS","Q_COCHRAN_WEIGHTED_BETA_TRANS","Q_Cochran_pvalue_TRANS")

Sig<- merge(all_sign, all_CIS, by=c("exposure","outcome", "method"),all = TRUE)
Sig1 <- merge(Sig, all_TRANS, by=c("exposure","outcome", "method"),all= TRUE)
Sig2 <- merge(Sig1, proteins, by="exposure")

to_remove <- Sig2 %>%
  group_by(exposure, outcome) %>%
  filter(all(is.na(nsnp.men) & is.na(nsnp.women))) %>%
  distinct(exposure, outcome)

Sig2 <- Sig2 %>% 
  anti_join(to_remove, by = c("exposure", "outcome"))

save <- Sig2 %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex_CISTRANS,SignificantBysex_CIS, SignificantBysex_TRANS,
       .direction = "downup") %>%
  ungroup()

test <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% 1 ) %>%
  dplyr::filter(SignificantBysex_CIS %in% 1 ) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

test2 <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% 2 ) %>%
  dplyr::filter(SignificantBysex_CIS %in% 2 ) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

test3 <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% 1 ) %>%
  dplyr::filter(SignificantBysex_TRANS %in% 1 ) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

test4 <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% 2 ) %>%
  dplyr::filter(SignificantBysex_TRANS %in% 2 ) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()


save$Q_Cochran_pvalue_CIS <-NULL
save$Q_Cochran_pvalue_TRANS <- NULL
save$Q_Cochran_pvalue_CISTRANS <-NULL

save$Diff_Beta_WOMEN_MEN_CISTRANS <- NULL
save$Diff_Beta_WOMEN_MEN_CIS <- NULL
save$Diff_Beta_WOMEN_MEN_TRANS <- NULL

save$Q_COCHRAN_WEIGHTED_BETA_CISTRANS <- NULL
save$Q_COCHRAN_WEIGHTED_BETA_CIS <- NULL
save$Q_COCHRAN_WEIGHTED_BETA_TRANS <- NULL

write.xlsx(save,"TABLE_SIGNIFICANT_WITH_ALL_CIS_TRANS.xlsx")

##### ST6. SIGNIFICANT AFTER BIDIRECTIONAL ####
bid <- read_excel("~/Results tables/Final_results_BIDIR_with_3columns.xlsx")

bid1<- subset(bid,select=c("exposure","outcome","method","SignificantBysex", "SignificantBysex_bid"))
colnames(bid1)<-c("outcome","exposure","method","SignificantBysex_bid","SignificantBysex_bid_new")

new<- merge(Sig2,bid1,by=c("exposure","outcome","method"), all=T)

new1 <- new %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex_CISTRANS,SignificantBysex_CIS, SignificantBysex_TRANS,SignificantBysex_bid, SignificantBysex_bid_new, 
       .direction = "downup") %>%
  ungroup()

New2 <- process_tables_after_bidirectional(new1)

N <- New2 %>%
  dplyr::select("exposure", "outcome","method") %>%
  distinct()

n <- Sig2 %>%
  dplyr::select("exposure", "outcome","method") %>%
  distinct()

s <- anti_join(N,n)

New2_clean <- anti_join(New2, s, by = c("exposure", "outcome", "method"))


##### ST5. LEAVE-ONE-OUT #####

# Upload leaveoneout results
women <- read.table("~/Women/final_results_for_paper/All_leaveoneout_women.txt", header=T)
women$id.exposure <- NULL
women$id.outcome <- NULL
women$samplesize <- NULL

men <- read.table("~/Men/final_results_for_paper/All_leaveoneout_men.txt", header=T)
men$id.exposure <- NULL
men$id.outcome <- NULL
men$samplesize <- NULL

all<- read_excel("~/Results tables/Final_results_with_column.xlsx")
all<- subset(all,select = c("exposure","outcome","method","IV_list.women","IV_list.men"))

# Select leave-one-out results only for significant after bidirectional for women and men
sig <- New2_clean
sig_pairs_W <- sig %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% c(1)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

loo_W <- merge(women, sig_pairs_W, all.y = TRUE)

not_sig_W <- loo_W %>%
  filter(p >= 0.05) 

# Apply filter_leaveoneout function to keep only significant results or results where MR-PRESSO remove the problematic SNP
list_to_keep_W <- filter_leaveoneout(not_sig_W, all, "IV_list.women")$keep
toremove_W <- filter_leaveoneout(not_sig_W, all, "IV_list.women")$remove

list_to_keep_W$driving_SNP <- "excluded_by_presso"
toremove_W$driving_SNP <-"yes"
toremove_W <- subset(toremove_W, select=c("exposure","outcome","driving_SNP"))
a <- rbind(list_to_keep_W,toremove_W)

loo_W1 <- merge(loo_W,a, by=c("exposure","outcome"), all.x=T)
loo_W1$sex<-"F"
loo_W1$driving_SNP[is.na(loo_W1$driving_SNP)] <- "no"

sig_pairs_M <- sig %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% c(2)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

loo_M <- merge(men, sig_pairs_M, all.y = TRUE)

not_sig_M <- loo_M %>%
  filter(p >= 0.05) 

# Apply filter_leaveoneout function to keep only significant results or results where MR-PRESSO remove the problematic SNP
list_to_keep_M <- filter_leaveoneout(not_sig_M, all, "IV_list.men")$keep
toremove_M <- filter_leaveoneout(not_sig_M, all, "IV_list.men")$remove

list_to_keep_M$driving_SNP <- "excluded_by_presso"
toremove_M$driving_SNP <-"yes"
toremove_M <- subset(toremove_M, select=c("exposure","outcome","driving_SNP"))
a <- rbind(list_to_keep_M,toremove_M)
a <-distinct(a)

loo_M1 <- merge(loo_M,a, by=c("exposure","outcome"), all.x=T)
loo_M1$sex<-"M"
loo_M1$driving_SNP[is.na(loo_M1$driving_SNP)] <- "no"

final_LOO <- rbind(loo_W1,loo_M1)
proteins <- read.table("~/ProteinsNames_new.csv",header=T,as.is=T,sep=",")
proteins<-subset(proteins,select =c("Assay","UniProt","exposure","Name"))
final_LOO <- merge(final_LOO,proteins,by="exposure")
write.xlsx(final_LOO ,"~/Results tables/Table_leaveoneout.xlsx")

# Remove pairs with not significant SNPs in leave-one-out (11 men and 5 in women)
sig_filtered_W <- anti_join(sig, toremove_W, by = c("exposure", "outcome"))
sig_after_loo <- anti_join(sig_filtered_W, toremove_M, by = c("exposure", "outcome"))

save <- sig_after_loo
save$Q_Cochran_pvalue_CIS <-NULL
save$Q_Cochran_pvalue_TRANS <- NULL
save$Q_Cochran_pvalue_CISTRANS <-NULL
save$Diff_Beta_WOMEN_MEN_CISTRANS <- NULL
save$Diff_Beta_WOMEN_MEN_CIS <- NULL
save$Diff_Beta_WOMEN_MEN_TRANS <- NULL

save$Q_COCHRAN_WEIGHTED_BETA_CISTRANS <- NULL
save$Q_COCHRAN_WEIGHTED_BETA_CIS <- NULL
save$Q_COCHRAN_WEIGHTED_BETA_TRANS <- NULL
save$SignificantBysex_bid <- NULL

names(save)[names(save) == "SignificantBysex_bid_new"] <- "SignificantBysex_bid"

save <- save %>%
  group_by(exposure, outcome)%>%
  filter(SignificantBysex_CISTRANS_new != 5) 

save$SignificantBysex_CISTRANS_new <- NULL
save$SignificantBysex_CIS_new <- NULL
save$SignificantBysex_TRANS_new <- NULL
save$SignificantBysex_bid <- NULL
write.xlsx(save,"~/Results tables/Table_significant_after_leaveoneout.xlsx")

sig_pairs_after_bidir <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% c(1, 2)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

w <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS ==1) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

m <- save %>%
  dplyr::filter(SignificantBysex_CISTRANS ==2) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

cis_only <- save %>% 
  dplyr::filter(SignificantBysex_CISTRANS %in% c(3,4)) %>%
  dplyr::filter(SignificantBysex_CIS %in% c(1, 2)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

##### ST7. TABLE WITH COCHRAN'S Q STATISTICS #####

sig_after_loo$SignificantBysex_bid <- NULL

Sig_womenORmen <- sig_after_loo %>%
  dplyr::group_by(exposure, outcome) %>%
  filter(SignificantBysex_CISTRANS_new %in% c(1,2,3,4)) %>%
  ungroup()
Sig4 <- Sig_womenORmen[Sig_womenORmen$method %in% c("Inverse variance weighted", "Wald ratio"), ]
Sig4$SignificantBysex_bid_new<- NULL
Sig4$SignificantBysex_CISTRANS_new<- NULL
Sig4$SignificantBysex_CIS_new<- NULL
Sig4$SignificantBysex_TRANS_new<- NULL
write.xlsx(Sig4 ,"~/Results tables/Table_cochran_after_leaveoneout.xlsx")


##### LIFELINES #####

pairs_for_lifelines <-  Sig4 %>%
  dplyr::filter(SignificantBysex_CISTRANS %in% c(1,2) | (SignificantBysex_CISTRANS %in% c(3,4) & Q_Cochran_pvalue_CISTRANS<0.05) ) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

###### ALL ######

# MEN

directory <- "~/Lifelines/Men"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_all_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men",
                         "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men","IV_list.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males <- new_males %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()

# WOMEN

# Upload file
directory <- "~/Lifelines/Women"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_all_F\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                         "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women <- new_women %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()

###### NOSTATINS USERS ####

# MEN

# Upload file
directory <- "~/Lifelines/Men"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_nostatins_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males_N <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males_N) <- c("outcome", "exposure", "method", "nsnp.men", "b.men","se.men","pval.men",
                           "egger_intercept.ple.men","se.ple.men","pval.ple.men","Q.het.men","Q_df.het.men","Q_pval.het.men","IV_list.men")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_males_N <- new_males_N %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.men %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()


# Put heterogeneity analyses in the same row of their methods
new_males_N <- new_males_N %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                        Q.het.men[method == "Weighted median"], Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                           Q_df.het.men[method == "Weighted median"], Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Inverse variance weighted",
                             Q_pval.het.men[method == "Weighted median"], Q_pval.het.men)
  ) %>%
  dplyr::mutate(
    Q.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q.het.men),
    Q_df.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.men),
    Q_pval.het.men = ifelse( dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.men)
  ) %>%
  ungroup()


# WOMEN
# Upload file
directory <- "~/Lifelines/Women"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_nostatins_F\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses

new_women_N <- merged_df_w 
colnames(new_women_N) <- c("outcome", "exposure", "method", "nsnp.women", "b.women","se.women","pval.women",
                           "egger_intercept.ple.women","se.ple.women","pval.ple.women","Q.het.women","Q_df.het.women","Q_pval.het.women","IV_list.women")

# Remove MR-PRESSO when we have 1, 2 or 3 SNPs because the row is empty (MR-PRESSO doesn't work)
new_women_N <- new_women_N %>%
  group_by(outcome, exposure) %>%
  filter(
    !(
      method == "MR-PRESSO" & 
        any(nsnp.women %in% c(1,2,3) & method != "MR-PRESSO")
    ) & 
      !(is.na(exposure) | exposure == "" )
  ) %>%
  ungroup()

# Put heterogeneity analyses in the same row of their methods
new_women_N <- new_women_N %>%
  group_by(outcome, exposure) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                         Q.het.women[method == "Weighted median"], Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                            Q_df.het.women[method == "Weighted median"], Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Inverse variance weighted",
                              Q_pval.het.women[method == "Weighted median"], Q_pval.het.women)
  ) %>%
  dplyr::mutate(
    Q.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q.het.women),
    Q_df.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_df.het.women),
    Q_pval.het.women = ifelse(dplyr::n() > 2 & method == "Weighted median", NA, Q_pval.het.women)
  ) %>%
  ungroup()


# Create final LIFELINES table
all <- merge(new_women,new_males,all=T)
colnames(all) <- c("outcome","exposure", "method", "nsnp.women.all" ,"b.women.all" ,"se.women.all" 
                   , "pval.women.all" ,"egger_intercept.ple.women.all" ,"se.ple.women.all" ,"pval.ple.women.all"   
                   , "Q.het.women.all","Q_df.het.women.all","Q_pval.het.women.all", "IV_list.women.all" , "nsnp.men.all"
                   , "b.men.all", "se.men.all", "pval.men.all","egger_intercept.ple.men.all", "se.ple.men.all", "pval.ple.men.all",
                   "Q.het.men.all", "Q_df.het.men.all", "Q_pval.het.men.all", "IV_list.men.all")
all$outcome <- sub("_.*| .*", "", all$outcome)
nostatins <- merge(new_women_N,new_males_N,all=T)
colnames(nostatins ) <- c("outcome","exposure", "method", "nsnp.women.nostatins" ,"b.women.nostatins" ,"se.women.nostatins" 
                          , "pval.women.nostatins" ,"egger_intercept.ple.women.nostatins" ,"se.ple.women.nostatins" ,"pval.ple.women.nostatins"   
                          , "Q.het.women.nostatins","Q_df.het.women.nostatins","Q_pval.het.women.nostatins", "IV_list.women.nostatins" , "nsnp.men.nostatins"
                          , "b.men.nostatins", "se.men.nostatins", "pval.men.nostatins","egger_intercept.ple.men.nostatins", "se.ple.men.nostatins", "pval.ple.men.nostatins",
                          "Q.het.men.nostatins", "Q_df.het.men.nostatins", "Q_pval.het.men.nostatins", "IV_list.men.nostatins")
nostatins$outcome <- sub("_.*| .*", "", nostatins$outcome)

setwd("~")
tot <- merge(all,nostatins, by=c("exposure","outcome","method"), all=T)

ALL_RESULTS <- tot

# Compute explained variance for Men
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.men.all)))%>%
  ungroup()
men <-  read.table("~/Lifelines/Men/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men.all" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men3$outcome <- sub("_.*| .*", "", men3$outcome)
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men.all"), all.x=TRUE)
men_WR.sig <- men_WR.sig %>% distinct()
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men.all")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.men.nostatins)))%>%
  ungroup()
men <-  read.table("~/Lifelines/Men/dat.males.txt",as.is=T,sep="\t")
colnames(men) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men.nostatins" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men3$outcome <- sub("_.*| .*", "", men3$outcome)
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men.nostatins"), all.x=TRUE)
men_WR.sig <- men_WR.sig %>% distinct()
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men.nostatins")
Sig_men_2 <- merge(Sig_men,men_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Compute explained variance for Women
Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.women.all)))%>%
  ungroup()
women <- read.table("~/Lifelines/Women/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(women3)<- c("IV_list.women.all" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women3$outcome <- sub("_.*| .*", "", women3$outcome)
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women.all"), all.x=TRUE)
women_WR.sig <- women_WR.sig %>% distinct()
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women.all")
Sig_women <- merge(Sig_men_2,women_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.women.nostatins)))%>%
  ungroup()
women <- read.table("~/Lifelines/Women/dat.females.txt",header=T,as.is=T,sep="\t")
colnames(women) <- c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","beta.exposure","beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic","ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome","outcome","mr_keep.outcome","pval_origin.outcome","data_source.outcome","samplesize.exposure","se.exposure","pval.exposure","id.exposure","exposure","mr_keep.exposure","pval_origin.exposure","data_source.exposure","action","SNP_index","mr_keep")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(women3)<- c("IV_list.women.nostatins" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
women3$outcome <- sub("_.*| .*", "", women3$outcome)
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","outcome","IV_list.women.nostatins"), all.x=TRUE)
women_WR.sig <- women_WR.sig %>% distinct()
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women.nostatins")


# Final table with explained variance
ALL_RESULTS1 <- merge(Sig_women,women_WR.sig1, by=c("exposure","outcome","method"),all = T)
#ALL_RESULTS1 <-subset(ALL_RESULTS1, method %in% c("Inverse variance weighted","Wald ratio"))


# Add SignificantBysex columns all and nostatins
lif <- process_tables_lifelines_all(ALL_RESULTS1)
lif2 <- process_tables_lifelines_nostatins(lif)

proteins <- read.table("~/ProteinsNames_new.csv",header=T,as.is=T,sep=",")
proteins<-subset(proteins,select =c("Assay","UniProt","exposure","Name"))

final3 <-  merge(lif2, proteins, by ="exposure")

table<-read_excel("~/Results tables/Table_significant_after_bidir.xlsx")

sub <- subset(table,select=c("exposure","outcome","method","b.women","b.men","SignificantBysex_CISTRANS"))

final5 <- merge(final3,sub,by=c("exposure","outcome","method"),all.x=T)

final6 <- final5 %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex_CISTRANS, SignificantBysex.all, SignificantBysex.nostatins,
       .direction = "downup") %>%
  ungroup()

###### ST8. LIFELINES FINAL TABLES ####
final7 <- column_GLGC(final6)
write.xlsx(final7,"~/Results tables/Final_LIFELINES_results_with_column_GLGC.xlsx")

# Count significant results
n <- final7 %>%
  dplyr::filter(SignificantBysex.all ==1) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)


n <- final6 %>%
  dplyr::filter(SignificantBysex.all ==2) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)

n <- final6 %>%
  dplyr::filter(SignificantBysex.all %in% c(3,4)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)


n <- final6 %>%
  dplyr::filter(SignificantBysex.nostatins ==1) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)

n <- final6 %>%
  dplyr::filter(SignificantBysex.nostatins ==2) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)

n <- final6 %>%
  dplyr::filter(SignificantBysex.nostatins %in% c(3,4)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n)

n <- final6 %>%
  dplyr::filter(SignificantBysex.nostatins %in% c(1,2) & SignificantBysex.all %in% c(1,2)) %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n) 

n <- final7 %>%
  dplyr::filter(Consistency_GLGC_Lifelines =="yes")  %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()
nrow(n) 
