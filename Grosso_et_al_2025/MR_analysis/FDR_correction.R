##########################################################################
# FDR correction and creation of a single xlsx file for the single outcome
##########################################################################

# Load necessary libraries
library(dplyr)
library(openxlsx)
library(readxl)

# Directory containing CSV files
directory <- ".../Results/"
name_exp <- "exp_name"

# Get the list of CSV files in the directory
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty list to store the data frames
df_list <- list()

# Read each CSV file and store it in the list
for (file in files) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  df_list[[file]] <- df
}

# Merge all data frames by row-binding and retaining all columns
merged_df <- bind_rows(df_list, .id = "source_file")

# Remove columns that contain "ple" or "het" in their names
merged_df <- merged_df %>% select(-any_of(c("ple", "het")))

# Save the merged file temporarily
output_file_xlsx <- paste0(directory, "merged_MR_",name_exp,".xlsx")
write.xlsx(merged_df, output_file_xlsx)

# Now apply the p-value correction
# Subset the dataset for the required methods 
datasub <- subset(merged_df, method %in% c("Inverse variance weighted", "Wald ratio"))

# Extract the p-values
pval <- datasub$pval

# Adjust p-values using the Benjamin-Hochberg method
pval_adj <- p.adjust(pval, method = "BH", n = length(pval))

# Create a new data frame for the p-values and adjusted p-values
data2 <- data.frame(pval = pval, pval_adj = pval_adj)

# Merge the original dataset with the adjusted p-values, matching on 'pval'
merged_data <- merge(merged_df, data2, by = "pval", all = TRUE)

# Write the final dataset to an Excel file
final_output <- paste0(directory, "merged_MR_",name_exp,"_corrected.xlsx")
write.xlsx(merged_data, final_output)

cat("Merging, column removal, and p-value adjustment completed. File created:", final_output, "\n")
