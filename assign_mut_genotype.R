# threshold_mutation_analysis.R
# Author: Amit Sud
# Date: 1st May 2025
# Description: Assigns mutated status to single cells based on thresholds, and compares mutation prevalence in tumor vs normal cells.
# Input: 
#   - obs_data: A dataframe with cell-level metadata (sample_id, *_MUT.calls, *_difference_proportion, etc.)
#   - expression_data: A gene expression matrix/dataframe (barcodes as rownames, gene names as columns)
# Output: 
#   - A dataframe summarizing optimal thresholds per gene and sample based on mutation and expression status

library(dplyr)
library(tidyr)
library(purrr)

# ---- Helper Function: Determine mutated status based on threshold ----
assign_mutated_status <- function(value, threshold) {
  value > threshold
}

# ---- Main Function ----
run_threshold_analysis <- function(obs_data, expression_data) {
  results_list <- list()

  for (sample_id in unique(obs_data$sample_id)) {
    message("Processing sample: ", sample_id)
    sample_data <- obs_data %>% filter(sample_id == !!sample_id)

    mut_cols <- grep("_MUT\\.calls$", colnames(sample_data), value = TRUE)

    for (mut_col in mut_cols) {
      gene <- sub("_MUT\\.calls$", "", mut_col)
      diff_col <- paste0(gene, "_difference_proportion")

      if (!diff_col %in% colnames(sample_data)) next

      valid_cells <- sample_data %>% filter(.data[[diff_col]] != 1)

      if (!gene %in% colnames(expression_data)) next

      gene_expr <- expression_data[rownames(valid_cells), gene, drop = FALSE]
      valid_cells$gene_expressed <- gene_expr[, 1] > 0

      for (threshold in 0:3) {
        mut_flag <- assign_mutated_status(valid_cells[[mut_col]], threshold)
        mut_colname <- paste0(gene, "_mutated_", threshold)
        valid_cells[[mut_colname]] <- mut_flag

        tumor_cells <- valid_cells %>%
          filter(Merged_Cluster_Name_scevan_broad_tumor_normal == "Tumor", gene_expressed)
        normal_cells <- valid_cells %>%
          filter(Merged_Cluster_Name_scevan_broad_tumor_normal == "Normal", gene_expressed)

        tumor_pct <- if (nrow(tumor_cells) > 0) {
          mean(tumor_cells[[mut_colname]]) * 100
        } else { 0 }

        normal_pct <- if (nrow(normal_cells) > 0) {
          mean(normal_cells[[mut_colname]]) * 100
        } else { 0 }

        results_list[[length(results_list) + 1]] <- data.frame(
          sample_id = sample_id,
          gene = gene,
          threshold = threshold,
          tumor_mutated_percentage = tumor_pct,
          normal_mutated_percentage = normal_pct
        )
      }
    }
  }

  # Combine and filter
  results_df <- bind_rows(results_list)

  nonzero_tumor <- results_df %>%
    group_by(sample_id, gene) %>%
    filter(sum(tumor_mutated_percentage) > 0) %>%
    ungroup()

  optimal_thresholds_df <- nonzero_tumor %>%
    filter(tumor_mutated_percentage > 0, normal_mutated_percentage <= 5) %>%
    group_by(sample_id, gene) %>%
    arrange(threshold) %>%
    slice_head(n = 1) %>%
    ungroup()

  return(optimal_thresholds_df)
}
