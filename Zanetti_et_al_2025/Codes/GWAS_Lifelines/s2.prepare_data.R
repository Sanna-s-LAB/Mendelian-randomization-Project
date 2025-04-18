library(dplyr)

setwd("/groups/umcg-lifelines/tmp01/projects/ov23_0790/data/phenotypes")

#
# Get phenotypes
#

d <- read.delim("1a_v_2_results.csv", as.is = T, check.names = F,  sep = ",")
d<- d[,c("project_pseudo_id", "age", "gender", "cholesterol_result_all_m_1", "hdlchol_result_all_m_1", "ldlchol_result_all_m_1", "triglyceride_result_all_m_1")]
colnames(d) <- c("id",  "age", "gender","TC", "HDC", "LDC", "TG")

# Set missing data to NA
d[d == ""] = NA
d[d == "NA"] = NA
d[d=="$5"] <- NA
d[d=="$6"] <- NA


# remove columns with NA in all lipids 
d_non_na <- d %>% 
  filter_at(vars(TC, HDC, LDC, TG), all_vars(!is.na(.))) %>%
  mutate_at(vars(TC, HDC, LDC, TG), as.numeric)

# Get non-HDL cholesterol
d_non_na$nonHDC <- d_non_na$TC - d_non_na$HDC

d_non_na %>% select(gender, TC, HDC, LDC, TG) %>% group_by(gender) %>% summarize_all(mean)

rm(d)

#
# add statin info
#
med <- read.delim("1a_v_1_results.csv", as.is = T, check.names = F,  sep = ",")

atc_codes_cols <- colnames(med)[grepl("atc_code_adu_c.*",colnames(med))]
#Statins: C10AA0
atc <- med[,c("project_pseudo_id", atc_codes_cols)]
atc[atc == ""] = NA

atc_non_na <- atc %>% 
  filter_at(atc_codes_cols, all_vars(!is.na(.)))

rm(med)

statin_users <- atc_non_na %>%
  rowwise() %>%
  filter(any(startsWith(as.character(c_across(everything())), "C10AA0")))

#
# add genotyping batch:
#
geno_ids <- read.delim("/groups/umcg-lifelines/tmp01/projects/ov23_0790/data/genotypes/linkage_1+2.actual_genotyped.txt", as.is = T, check.names = F, header = F, sep = "\t")
colnames(geno_ids) <- c("genotype_id", "phenotype_id", "genotyping_batch")

#
# merge all 3 tables
#
d_non_na$statin_user <- 0
d_non_na[d_non_na$id %in% statin_users$project_pseudo_id, "statin_user"] <- 1

merged <- left_join(d_non_na, geno_ids, by = c("id" = "phenotype_id"))

write.table(merged, file = "pheno_data_combined.txt", sep = "\t", quote = F, row.names = F)


#
# Split per group (sex , statins). Adjust for covariates and normalize
#
merged$TG <- log(merged$TG)
merged <- merged[!is.na(merged$genotyping_batch),]
merged$id <- merged$genotype_id
merged$genotype_id <- NULL

merged$group <- paste0(merged$gender, "_", merged$statin_user)


regress_covariates <- function(d, lipid){
  lm_formula = as.formula(paste(lipid, "~ age + I(age^2)"))
  res <- residuals(lm(lm_formula, data = d))
  
  return(res)
}

run_INT <- function(d) {
  d_norm <- apply(d, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  return(d_norm)
}

tmp_merged_adj <- data.frame()
lipids <- c("TC", "HDC", "LDC", "TG", "nonHDC")
for (gr in unique(merged$group)){
  data_subs <- merged[merged$group == gr,]
  pheno_adj <- data.frame(matrix(nrow = nrow(data_subs), ncol = length(lipids) + 1))
  colnames(pheno_adj) <- c("IID", lipids)
  pheno_adj$IID <- data_subs$id
  
  # Regress out covairates
  for (lipid in lipids){
    pheno_adj[,lipid] <- regress_covariates(data_subs, lipid)
  }
  
  # tmp df for plotting
  tmp_merged_adj <- rbind(tmp_merged_adj, cbind(pheno_adj, gr))
  
  # Inverse-rank transform
  pheno_adj[2:ncol(pheno_adj)] <- run_INT(pheno_adj[2:ncol(pheno_adj)])
  
  # Format for fastGWA
  pheno_adj <- cbind(a=0, pheno_adj)
  colnames(pheno_adj)[1] <- "FID"
  
  write.table(pheno_adj, file = paste0("pheno_adj_", gr,".txt"), sep = "\t", quote = F, row.names = F)
}


#
# Plot the phenotypes
#
library(ggplot2)
library(patchwork)

plot_lst <- list()
boxplots <- list()
for (pheno in lipids){
  pheno2 <- ensym(pheno) 
  plot_lst[[pheno]] <- ggplot(tmp_merged_adj, aes (x = !!pheno2, fill =  gr, color = gr)) + 
    geom_histogram(alpha = 0.3, position = "identity") + 
    scale_color_brewer(palette = "Set2") + 
    scale_fill_brewer(palette = "Set2") +
    theme_bw() + 
    ggtitle(pheno)
  
  boxplots[[pheno]] <- ggplot(tmp_merged_adj, aes (y = !!pheno2, fill =  gr, color = gr)) + 
    geom_boxplot(alpha = 0.3) + 
    scale_color_brewer(palette = "Set2") + 
    scale_fill_brewer(palette = "Set2") +
    theme_bw() + 
    ggtitle(pheno)
}

pdf("all_lipids_adj_hist.pdf", width = 10, height = 15)
(plot_lst[[1]] + plot_lst[[2]] ) / (plot_lst[[3]] + plot_lst[[4]] ) / (plot_lst[[5]] + plot_spacer() )
dev.off()

pdf("all_lipids_adj_boxplots.pdf", width = 10, height = 15)
(boxplots[[1]] + boxplots[[2]] ) / (boxplots[[3]] + boxplots[[4]] ) / (boxplots[[5]] + plot_spacer())
dev.off()


#
# Split only by sex, for any statin usage status
#

tmp_merged_adj <- data.frame()
lipids <- c("TC", "HDC", "LDC", "TG", "nonHDC")
for (gender in unique(merged$gender)){
  data_subs <- merged[merged$gender == gender,]
  pheno_adj <- data.frame(matrix(nrow = nrow(data_subs), ncol = length(lipids) + 1))
  colnames(pheno_adj) <- c("IID", lipids)
  pheno_adj$IID <- data_subs$id
  
  # Regress out covariates
  for (lipid in lipids){
    pheno_adj[,lipid] <- regress_covariates(data_subs, lipid)
  }
  
  # Inverse-rank transform
  pheno_adj[2:ncol(pheno_adj)] <- run_INT(pheno_adj[2:ncol(pheno_adj)])
  
  # Format for fastGWA
  pheno_adj <- cbind(a=0, pheno_adj)
  colnames(pheno_adj)[1] <- "FID"
  
  write.table(pheno_adj, file = paste0("pheno_adj_", gender,"_all.txt"), sep = "\t", quote = F, row.names = F)
}





