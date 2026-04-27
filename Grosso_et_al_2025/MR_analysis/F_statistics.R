source("../functions.R")

directory <- "../Clumping_results_no_NA1"
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.csv, header = TRUE, stringsAsFactors = FALSE)

N <- 7738

# global F-statistics
calc_F_global <- function(R2, N, k) {
  if (R2 >= 1) return(Inf)
  return((R2 * (N - k - 1)) / ((1 - R2) * k))
}

# Extract file names
GCS_names <- tools::file_path_sans_ext(basename(files))

# Compute total R2 and F global for every intrument
results <- mapply(function(df, name) {
  R2_vector <- explained_variance(df, N)
  R2_total <- sum(R2_vector, na.rm = TRUE)
  k <- nrow(df)-1
  F_global <- calc_F_global(R2_total, N, k)
  return(data.frame(GCS = name, R2 = R2_total, F_stat = F_global, snps=k))
}, df_list, GCS_names, SIMPLIFY = FALSE)

result_df <- do.call(rbind, results)

write.csv(result_df, file = file.path(directory, "explained_variance_summary.csv"), row.names = FALSE)

print(result_df)
