#!/usr/bin/env Rscript

library(plyr) 
library(dplyr)
library(openxlsx)
library(tidyverse)
library(readxl)
library(tidyr)
library(readr)
library(data.table)
source("~/Tables_new/process_tables_functions.R")

##### ST1. FINAL RESULT FOR PAPER - TABLES #####

###### MEN ######

# Upload file
directory <- "~/MR_onlycis_UKBB/MR_res_proxies/male"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males)[4:ncol(new_males)] <- paste0(colnames(new_males)[4:ncol(new_males)], ".men")

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

result_main <- new_males %>%
     filter(method == "Inverse variance weighted" | method=="Wald ratio")

result_Egger <- new_males %>%
  filter(method == "MR Egger")

result_Egger$fdr.men <- c("NA")

result_WM <- new_males %>%
  filter(method == "Weighted median")

result_WM$fdr.men <-  c("NA")

result_PRESSO <- new_males %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.men <- c("NA")

result_main$fdr.men<-p.adjust(result_main$pval.men, method="BH")

men_all <- rbind(result_main,result_Egger,result_WM, result_PRESSO)

str(men_all)
men_all$fdr.men<- as.numeric(men_all$fdr.men)
men_all$pval.men[men_all$nsnp.men == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

leave_one_out <- read.table("~/MR_onlycis_UKBB/MR_res_proxies/male/leaveoneout_males.txt", header=T)
# Merge all the leaveoneout files
system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_males.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_males.txt'")
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
final_men <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = men_all,
  type = "men",
  iv_list_col = "IV_list.men"
)
###### WOMEN ######

# Upload file
directory <- "~/MR_onlycis_UKBB/MR_res_proxies/female"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women)[4:ncol(new_women)] <- paste0(colnames(new_women)[4:ncol(new_women)], ".women")

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

result_main <- new_women %>%
  filter(method == "Inverse variance weighted" | method=="Wald ratio")

result_Egger <- new_women %>%
  filter(method == "MR Egger")

result_Egger$fdr.women <- c("NA")

result_WM <- new_women %>%
  filter(method == "Weighted median")

result_WM$fdr.women <-  c("NA")

result_PRESSO <- new_women %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.women <- c("NA")

result_main$fdr.women<-p.adjust(result_main$pval.women, method="BH")
women_all <- rbind(result_main,result_Egger,result_WM, result_PRESSO)

str(women_all)
women_all$fdr.women<- as.numeric(women_all$fdr.women)
women_all$pval.women[women_all$nsnp.women == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs

leave_one_out <- read.table("~/MR_onlycis_UKBB/MR_res_proxies/female/leaveoneout_females.txt", header=T)
# Merge all the leaveoneout files
system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_females.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_females.txt'")
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
final_women <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = women_all,
  type = "women",
  iv_list_col = "IV_list.women"
)
######  ALL ######

# Upload file
directory <- "~/MR_onlycis_UKBB/MR_res_proxies/MA"

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

result_main <- new_all %>%
  filter(method == "Inverse variance weighted" | method=="Wald ratio")


result_Egger <- new_all %>%
  filter(method == "MR Egger")

result_Egger$fdr.all <- c("NA")

result_WM <- new_all %>%
  filter(method == "Weighted median")

result_WM$fdr.all <-  c("NA")

result_PRESSO <- new_all %>%
  filter(method == "MR-PRESSO")

result_PRESSO$fdr.all <- c("NA")

result_main$fdr.all<-p.adjust(result_main$pval.all, method="BH")
all_all <- rbind(result_main,result_Egger,result_WM, result_PRESSO)

str(all_all)
all_all$fdr.all<- as.numeric(all_all$fdr.all)
all_all$pval.all[all_all$nsnp.all == 0] <- 1 # set p-value to 1 when MR-PRESSO removes all SNPs


######  TABLE of ALL RESULTS (significant or not) #####

setwd("~/Tables_new")
ALL_RESULTS <- merge(final_women,final_men, by=c("exposure","outcome", "method"),all = TRUE)
ALL_RESULTS <- merge(ALL_RESULTS, all_all, by=c("exposure","outcome", "method"), all=TRUE)
ALL_RESULTS <- ALL_RESULTS[order(ALL_RESULTS$exposure, ALL_RESULTS$outcome, ALL_RESULTS$method, decreasing = FALSE), ]
ALL_RESULTS <- ALL_RESULTS[!is.na(ALL_RESULTS$exposure), ]

write.xlsx(ALL_RESULTS, "~/Tables_new/AllCIS_results.xlsx")

###### TABLES with columns ####

setwd("~/Tables_new")
proteins <- data.table::fread("~/olink_protein_map_3k_v1.tsv")
proteins<-subset(proteins,select =c("OlinkID","Assay","UniProt","olink_target_fullname"))
proteins <- proteins[!duplicated(proteins), ]

colnames(proteins) <- c("exposure","Assay","UniProt","olink_target_fullname")

ALL_RESULTS <- read_excel("AllCIS_results.xlsx")

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.men))%>%
  ungroup()

# Merge all the leaveoneout files
#system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_males.txt && tail -n +2 -q *dat*.txt >> dat_males.txt'")

men <-  read.table("~/MR_onlycis_UKBB/MR_res_proxies/male/dat_males.txt",as.is=T,sep="\t", header = T)
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure" )
men3 <- men3[!duplicated(men3), ]
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)

men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)

Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(fdr.women))%>%
  ungroup()

women <- read.table("~/MR_onlycis_UKBB/MR_res_proxies/female/dat_females.txt",header=T,as.is=T,sep="\t")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure" )
colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"   )
women3 <- women3[!duplicated(women3), ]
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# All results with variance explained columns for Wald ratio results
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)


# Apply function to add significance columns
final <- process_tables(ALL_RESULTS1)
final <- merge(final, proteins, by="exposure")
final1 <- process_tables_truly_sex_specific(final)
write.xlsx(final1,"~/Tables_new/ST1_all_results.xlsx")

##### ST2. LEAVE-ONE-OUT TABLE ######

leave_one_out_m <- read.table("~/MR_onlycis_UKBB/MR_res_proxies/male/leaveoneout_males.txt", header=F)
colnames(leave_one_out_m) <- c("exposure",	"outcome",	"id.exposure",	"id.outcome",	"samplesize",	"SNP",	"b",	"se",	"p")
leave_one_out_m$sex <- "M"

leave_one_out_f <- read.table("~/MR_onlycis_UKBB/MR_res_proxies/female/leaveoneout_females.txt", header=F)
colnames(leave_one_out_f) <- c("exposure",	"outcome",	"id.exposure",	"id.outcome",	"samplesize",	"SNP",	"b",	"se",	"p")
leave_one_out_f$sex <- "F"

LOO <- rbind(leave_one_out_m,leave_one_out_f)
write.xlsx(LOO,"~/Tables_new/ST2_leave_one_out.xlsx")


##### ST3. SIGNIFICANT RESULTS #####

sig <- final1 %>%
  dplyr::filter(SignificantBysex %in% c(1,2,3,4)) 

write.xlsx(sig,"~/Tables_new/ST3_significant_results.xlsx")

##### Save list of significant relationships per outcome ####

sig_women <- final1 %>%
  dplyr::filter(SignificantBysex %in% c(1,3,4)) 
exposure_list <- tapply(sig_women$exposure, sig_women$outcome, function(x) unique(x))

for (exposure_name in names(exposure_list)) {
  exposure_data <- exposure_list[[exposure_name]]
  
  filename <- paste0("~/Tables_new/list_women_", exposure_name)
  
  writeLines(as.character(exposure_data), con = filename)
}

sig_men <- final1 %>%
  dplyr::filter(SignificantBysex %in% c(2,3,4)) 
exposure_list <- tapply(sig_men$exposure, sig_men$outcome, function(x) unique(x))

for (exposure_name in names(exposure_list)) {

  exposure_data <- exposure_list[[exposure_name]]
  
  filename <- paste0("~/Tables_new/list_men_", exposure_name)
  
  writeLines(as.character(exposure_data), con = filename)
}

##### ST4. BIDIRECTIONAL MR ######

# MEN 

# Upload file
directory <- "~/MR_onlycis_UKBB/BidirectionalMR/male" # da cambiare

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 
colnames(new_males)[4:ncol(new_males)] <- paste0(colnames(new_males)[4:ncol(new_males)], ".men")


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


# merge all the leaveoneout files
#system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_males.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_males.txt'")

leave_one_out <- read.table("~/MR_onlycis_UKBB/BidirectionalMR/male/leaveoneout_males.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]

men_all <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = men_all,
  type = "men",
  iv_list_col = "IV_list.men")

# WOMEN 

# Upload file
directory <- "~/MR_onlycis_UKBB/BidirectionalMR/female" 

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 
colnames(new_women)[4:ncol(new_women)] <- paste0(colnames(new_women)[4:ncol(new_women)], ".women")


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

# merge all the leaveoneout files
#system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_females.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_females.txt'")

leave_one_out <- read.table("~/MR_onlycis_UKBB/BidirectionalMR/female/leaveoneout_females.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
women_all <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = women_all,
  type = "women",
  iv_list_col = "IV_list.women"
)

ALL_RESULTS <- merge(women_all,men_all, by=c("exposure","outcome", "method"),all = TRUE)
write.xlsx(ALL_RESULTS, "~/Tables_new/All_BIDIRECTIONAL_results.xlsx")

# Select only Wald ratio results to compute variance explained
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(pval.men))%>%
  ungroup()

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_males.txt && tail -n +2 -q *dat*.txt >> dat_males.txt'")
men <-  read.table("~/MR_onlycis_UKBB/BidirectionalMR/male/dat_males.txt",as.is=T,sep="\t", header = T)
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  )
colnames(men3)<- c("IV_list.men" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure" )
men3 <- men3[!duplicated(men3), ]
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","IV_list.men"), all.x=TRUE)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)

men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)

Sig_women_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & !is.na(pval.women))%>%
  ungroup()

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_females.txt && tail -n +2 -q *dat*.txt >> dat_females.txt'")
women <- read.table("~/MR_onlycis_UKBB/BidirectionalMR/female/dat_females.txt",header=T,as.is=T,sep="\t")
women2 <- women %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE")
women3<- women2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure" )
colnames(women3)<- c("IV_list.women" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"   )
women3 <- women3[!duplicated(women3), ]
women_WR.sig <- merge(Sig_women_Wald, women3, by=c("exposure","IV_list.women"), all.x=TRUE)
women_WR.sig$exp_var <- explained_variance(women_WR.sig,women_WR.sig$samplesize.exposure)
women_WR.sig1 <- data.frame(women_WR.sig$exposure,women_WR.sig$outcome,women_WR.sig$method,women_WR.sig$exp_var)
colnames(women_WR.sig1)<-c("exposure","outcome","method","exp_var.women")

# All results with variance explained columns for Wald ratio results
ALL_RESULTS1 <- merge(Sig_men,women_WR.sig1, by=c("exposure","outcome","method"),all = T)

# Process tables with process_tables_bidirectional function to add the "SignificantBysex" column
final2 <- process_tables_bidirectional(ALL_RESULTS1)

###### Final column for bidirectional MR results ######
all_sign <- read_excel("~/Tables_new/ST3_significant_results.xlsx")
sub <- subset(all_sign, select=c("exposure","outcome","method","SignificantBysex"))
colnames(sub)<- c("outcome","exposure","method","SignificantBysex_main")
final3 <- merge(final2, sub, by=c("exposure","outcome","method"), all=T)

final4 <- final3 %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex_main, SignificantBysex,
       .direction = "downup") %>%
  ungroup()

# process data to write the column in the final way
final5 <- column_bidirectional(final4)

setwd("~/Tables_new")
proteins <- data.table::fread("~/olink_protein_map_3k_v1.tsv")
proteins<-subset(proteins,select =c("OlinkID","Assay","UniProt","olink_target_fullname"))
proteins <- proteins[!duplicated(proteins), ]

colnames(proteins)<- c( "outcome", "Assay" ,   "UniProt"  ,"Name" )
final6 <-  merge(final5, proteins, by ="outcome")

write.xlsx(final6,"~/Tables_new/Final_results_BIDIR_with_3columns.xlsx")
final6$SignificantBysex <- NULL
final6$SignificantBysex_main <-NULL
write.xlsx(final6,"~/Tables_new/ST5_Final_results_BIDIR_with_column.xlsx") # ST5

##### ST5. SIGNIFICANT AFTER BIDIRECTIONAL AND COCHRAN'S ####
bid <- read_excel("~/Tables_new/Final_results_BIDIR_with_3columns.xlsx")

bid1<- subset(bid,select=c("exposure","outcome","method","SignificantBysex", "SignificantBysex_bid"))
colnames(bid1)<-c("outcome","exposure","method","SignificantBysex_bid","SignificantBysex_bid_new")

new<- merge(all_sign,bid1,by=c("exposure","outcome","method"), all=T)

new1 <- new %>%
  group_by(outcome, exposure) %>%
  fill(SignificantBysex, SignificantBysex_bid, SignificantBysex_bid_new, 
       .direction = "downup") %>%
  ungroup()

New2 <- process_tables_after_bidirectional(new1)

N <- New2 %>%
  dplyr::select("exposure", "outcome","method") %>%
  distinct()

n <- all_sign %>%
  dplyr::select("exposure", "outcome","method") %>%
  distinct()

s <- anti_join(N,n)

New2_clean <- anti_join(New2, s, by = c("exposure", "outcome", "method"))


Sig_womenORmen <- New2_clean %>%
  dplyr::group_by(exposure, outcome) %>%
  filter(SignificantBysex_new %in% c(1,2,3,4)) %>%
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

Sig4 <- Sig_womenORmen
Sig4$SignificantBysex_bid_new<- NULL
Sig4$SignificantBysex_new<- NULL
write.xlsx(Sig4 ,"~/Tables_new/ST5_Table_after_bidir_and_cochran_all_methods.xlsx")
Sig4 <- Sig4 %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio"))
write.xlsx(Sig4 ,"~/Tables_new/ST5_Table_after_bidir_and_cochran.xlsx")


##### ST6. LIFELINES #####

pairs_for_lifelines <-  Sig4 %>%
     dplyr::filter(SexSpecific %in% c("men-only", "women-only") | Q_Cochran_pvalue<=0.05)%>%
     dplyr::select("exposure", "outcome") %>%
     distinct()

exposure_list <- tapply(pairs_for_lifelines$exposure, pairs_for_lifelines$outcome, function(x) unique(x))

for (exposure_name in names(exposure_list)) {
  exposure_data <- exposure_list[[exposure_name]]
  
  filename <- paste0("~/Tables_new/list_lifelines_", exposure_name)
  
  writeLines(as.character(exposure_data), con = filename)
}


###### ALL ######

# MEN

directory <- "~/MR_onlycis_UKBB/Lifelines/male"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_all_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 

colnames(new_males)[4:ncol(new_males)] <- paste0(colnames(new_males)[4:ncol(new_males)], ".men")
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

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_males.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_males.txt'")
leave_one_out <- read.table("~/MR_onlycis_UKBB/Lifelines/male/leaveoneout_males.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
new_males <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = new_males,
  type = "men",
  iv_list_col = "IV_list.men"
)

# WOMEN

# Upload file
directory <- "~/MR_onlycis_UKBB/Lifelines/female"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_all_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses
new_women <- merged_df_w 

colnames(new_women)[4:ncol(new_women)] <- paste0(colnames(new_women)[4:ncol(new_women)], ".women")
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

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_females.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_females.txt'")
leave_one_out <- read.table("~/MR_onlycis_UKBB/Lifelines/female/leaveoneout_females.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
new_women <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = new_women,
  type = "women",
  iv_list_col = "IV_list.women"
)
###### NOSTATINS USERS ####

# MEN

# Upload file
directory <- "~/MR_onlycis_UKBB/Lifelines/male"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_nostatins_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df <- merged_df %>% dplyr::select(-any_of(c("ple", "het")))
merged_df <- merged_df[,-c(1,2)] # remove id
new_males_N <- merged_df
indices_na_exposure <- which(is.na(merged_df$exposure)) # check for failed analyses, 

colnames(new_males_N)[4:ncol(new_males_N)] <- paste0(colnames(new_males_N)[4:ncol(new_males_N)], ".men")
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

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_males.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_males.txt'")
leave_one_out <- read.table("~/MR_onlycis_UKBB/Lifelines/male/leaveoneout_males.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
new_males_N <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = new_males_N,
  type = "men",
  iv_list_col = "IV_list.men"
)
# WOMEN
# Upload file
directory <- "~/MR_onlycis_UKBB/Lifelines/female"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = ".*_nostatins_M\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)
merged_df_w <- rbind.fill(df_list)

# Remove columns "ple" or "het"
merged_df_w  <- merged_df_w  %>% dplyr::select(-any_of(c("ple", "het")))
merged_df_w  <- merged_df_w [,-c(1,2)] # remove id
indices_na_exposure_w <- which(is.na(merged_df_w$exposure)) # check for failed analyses

new_women_N <- merged_df_w 

colnames(new_women_N)[4:ncol(new_women_N)] <- paste0(colnames(new_women_N)[4:ncol(new_women_N)], ".women")

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

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *leaveoneout*.txt | head -n 1) > leaveoneout_females.txt && tail -n +2 -q *leaveoneout*.txt >> leaveoneout_females.txt'")
leave_one_out <- read.table("~/MR_onlycis_UKBB/Lifelines/female/leaveoneout_females.txt", header=T)
leave_one_out$id.exposure <- NULL
leave_one_out$id.outcome <- NULL
leave_one_out$samplesize <- NULL
leave_one_out <- na.omit(leave_one_out)
leave_one_out <- leave_one_out[!(leave_one_out$SNP == "All"),]
new_women_N <- process_leave_one_out(
  leave_one_out = leave_one_out,
  base_data = new_women_N,
  type = "women",
  iv_list_col = "IV_list.women"
)
# Create final LIFELINES table
all <- merge(new_women,new_males,all=T)
colnames(all)[4:ncol(all)] <- paste0(colnames(all)[4:ncol(all)], ".all")
all$outcome <- sub("_.*| .*", "", all$outcome)
nostatins <- merge(new_women_N,new_males_N,all=T)
colnames(nostatins)[4:ncol(nostatins)] <- paste0(colnames(nostatins)[4:ncol(nostatins)], ".nostatins")
nostatins$outcome <- sub("_.*| .*", "", nostatins$outcome)

setwd("~/Tables_new")
tot <- merge(all,nostatins, by=c("exposure","outcome","method"), all=T)

ALL_RESULTS <- tot

# Compute explained variance for Men
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.men.all)))%>%
  ungroup()

# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_males.txt && tail -n +2 -q *dat*.txt >> dat_males.txt'")
men <-  read.table("~/MR_onlycis_UKBB/Lifelines/male/dat_males.txt",as.is=T,sep="\t", header = T)
men2 <- men %>%
  filter(remove == "FALSE" &  ambiguous == "FALSE") 
men3<- men2 %>% dplyr::select("SNP" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
colnames(men3)<- c("IV_list.men.all" , "beta.exposure", "eaf.exposure" ,"samplesize.exposure"  ,"se.exposure" , "exposure"   ,"pval.exposure"  ,"outcome"   )
men3$outcome <- sub("_.*| .*", "", men3$outcome)
men_WR.sig <- merge(Sig_men_Wald, men3, by=c("exposure","outcome","IV_list.men.all"), all.x=TRUE)
men_WR.sig <- men_WR.sig %>% distinct()
men_WR.sig$samplesize.exposure <- as.numeric(men_WR.sig$samplesize.exposure)
men_WR.sig$eaf.exposure <- as.numeric(men_WR.sig$eaf.exposure)
men_WR.sig$beta.exposure <- as.numeric(men_WR.sig$beta.exposure)
men_WR.sig$se.exposure <- as.numeric(men_WR.sig$se.exposure)
men_WR.sig$exp_var <- explained_variance(men_WR.sig,men_WR.sig$samplesize.exposure)
men_WR.sig1 <- data.frame(men_WR.sig$exposure,men_WR.sig$outcome,men_WR.sig$method,men_WR.sig$exp_var)
colnames(men_WR.sig1)<-c("exposure","outcome","method","exp_var.men.all")
Sig_men <- merge(ALL_RESULTS,men_WR.sig1, by=c("exposure","outcome","method"),all = T)
Sig_men_Wald <- ALL_RESULTS %>%
  group_by(exposure, outcome) %>%
  filter(
    method == "Wald ratio" & (!is.na(pval.men.nostatins)))%>%
  ungroup()
# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_males.txt && tail -n +2 -q *dat*.txt >> dat_males.txt'")
men <-  read.table("~/MR_onlycis_UKBB/Lifelines/male/dat_males.txt",as.is=T,sep="\t", header=T)
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
# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_females.txt && tail -n +2 -q *dat*.txt >> dat_females.txt'")
women <- read.table("~/MR_onlycis_UKBB/Lifelines/female/dat_females.txt",header=T,as.is=T,sep="\t")
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
# merge all the dat files
# system("bash -c 'head -n 1 $(ls *dat*.txt | head -n 1) > dat_females.txt && tail -n +2 -q *dat*.txt >> dat_females.txt'")
women <- read.table("~/MR_onlycis_UKBB/Lifelines/female/dat_females.txt",header=T,as.is=T,sep="\t")
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

# Add SignificantBysex columns all and nostatins
lif3 <- process_tables_truly_sex_specific_lifelines(ALL_RESULTS1, col_name = "SexSpecific.all", suffix = ".all")
lif4 <- process_tables_truly_sex_specific_lifelines(lif3, col_name = "SexSpecific.nostatins", suffix = ".nostatins")

proteins <- data.table::fread("~/olink_protein_map_3k_v1.tsv")
proteins<-subset(proteins,select =c("OlinkID","Assay","UniProt","olink_target_fullname"))
proteins <- proteins[!duplicated(proteins), ]
colnames(proteins)<-c("exposure","Assay","UniProt","Name")

final3 <-  merge(lif4, proteins, by ="exposure")

table<-read_excel("~/Tables_new/ST5_Table_after_bidir_and_cochran.xlsx")

sub <-  table %>%
   dplyr::filter(SexSpecific %in% c("men-only", "women-only") | Q_Cochran_pvalue<=0.05)
n <- sub %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

table1<-read_excel("~/Tables_new/ST5_Table_after_bidir_and_cochran_all_methods.xlsx")

combinations_sub <- unique(sub[, c("exposure", "outcome")])

table_filtered <- merge(table1, combinations_sub, by = c("exposure", "outcome"))

df2_new <- anti_join(table_filtered, sub, by = c("exposure", "outcome", "method"))

df_final <- bind_rows(sub, df2_new)

sub <- subset(df_final,select=c("exposure","outcome","method","b.women","b.men","SexSpecific","Q_Cochran_pvalue"))

final5 <- merge(final3,sub,by=c("exposure","outcome","method"),all=T)

final6 <- final5 %>%
  group_by(outcome, exposure) %>%
  fill(SexSpecific.all,SexSpecific.nostatins, SexSpecific,
       .direction = "downup") %>%
  ungroup()
final6<-final6[!is.na(final6$SexSpecific),]

# Cochran's Q
final6 <- calc_cochran_Q(
  final6,
  b_women_col = "b.women.all", se_women_col = "se.women.all",
  b_men_col = "b.men.all", se_men_col = "se.men.all",
  suffix = "all"
)

final6 <- calc_cochran_Q(
  final6,
  b_women_col = "b.women.nostatins", se_women_col = "se.women.nostatins",
  b_men_col = "b.men.nostatins", se_men_col = "se.men.nostatins",
  suffix = "nostatins"
)


###### LIFELINES FINAL TABLES ####
final6 <- column_GLGC_sexspecific(final6)
write.xlsx(final6,"~/Tables_new/ST6_Final_LIFELINES_results_with_column_GLGC.xlsx")


##### ST5 - add column for opposite estimate #####


df <- read_excel("~/Tables_new/SuppTables.xlsx",sheet="ST5")

list1 <- list()

combinations <-  df %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()


for (i in 1:nrow(combinations)) {
  exp <- combinations$exposure[i]
  out <- combinations$outcome[i]
  sub_df <- subset(df, exposure == exp & outcome == out)
  
  if(any(sign(na.omit(sub_df$b.women))!=sign(na.omit(sub_df$b.men)))){
    sub_df$opposite_estimate <- "yes"
  } else {
    sub_df$opposite_estimate <- "no"
  }
  list1[[paste0(exp, "_", out)]] <- sub_df
}

df <- do.call(rbind, list1)

write.xlsx(df,"~/Tables_new/ST5.xlsx")




##### ST7  #####


df <- read_excel("~/Tables_new/SuppTables.xlsx",sheet="ST5")

sig <- df %>%
  dplyr::filter(SexSpecific %in% c("men-only", "women-only") | Q_Cochran_pvalue<=0.05) %>%
dplyr::select("Assay","UniProt","olink_target_fullname" ) %>%
  distinct()

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(biomaRt)
library(rentrez)
library(dplyr)

genes <- sig$Assay

# Step 1: Get full gene names using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("hgnc_symbol", "description"),
  filters = "hgnc_symbol",
  values = genes,
  mart = ensembl
)

# Step 2: Get function summaries using rentrez
get_gene_summary <- function(symbol) {
  query <- paste0(symbol, "[Gene Name] AND Homo sapiens[Organism]")
  search <- entrez_search(db = "gene", term = query, retmax = 1)
  if (length(search$ids) == 0) return(NA)
  summary <- entrez_summary(db = "gene", id = search$ids[[1]])
  return(summary$summary)
}

# Step 3: Create data frame of summaries
gene_summaries <- sapply(genes, get_gene_summary)
summary_df <- data.frame(hgnc_symbol = names(gene_summaries), Function = gene_summaries)

# Step 4: Merge full name + function into one table
final_df <- merge(gene_info, summary_df, by = "hgnc_symbol", all = TRUE)

colnames(final_df) <- c("Gene", "Full_Name", "Function")

print(final_df)

write.csv(final_df, "ST7.csv", row.names = FALSE)









