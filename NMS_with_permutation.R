# Title: Calculate Mixing Score Between Cell Types Using Spatial Coordinates
# Author: [Your Name]
# Date: [Date]
#
# Description:
# This script computes a "mixing score" to quantify spatial proximity between two cell types
# in spatial transcriptomics data. It uses fixed-radius nearest neighbor (frNN) analysis and 
# permutation testing to assess statistical significance.
#
# Inputs:
# - A directory of composite CSV files, each containing cell-level spatial data and annotations
#   (default: input_directory, pattern: *_TRA_*_composite.csv or *_TRB_*_composite.csv)
#   Each CSV must include: 
#   - Cell.X.Position
#   - Cell.Y.Position
#   - Cell.Type (annotated cell type)
#
# Outputs:
# - A CSV file for each input, storing the computed mixing score across radii (10 to 400 pixels)
#   for multiple reference-target cell type combinations.
#   Saved in output_directory with filenames like: `<sample>_<gene>_<clonotype>_TRA_mixing_scores.csv`
#
# Dependencies:
# dbscan, ggplot2, dplyr, vroom

library(dbscan)
library(ggplot2)
library(dplyr)
library(vroom)

# --------- Configuration ---------
input_directory <- "path/to/input/"
output_directory <- "path/to/output/"
dir.create(output_directory, showWarnings = FALSE)
file_pattern <- "_TRA_.*_composite.csv$"  # or "_TRB_.*_composite.csv$"
# ---------------------------------

# --------- Mixing Score Function ---------
mixing_score_manual <- function(data, reference_celltype, target_celltype, radius, feature_colname, n_permutations = 1000) {
  reference_cells <- data %>% filter(!!sym(feature_colname) == reference_celltype)
  target_cells <- data %>% filter(!!sym(feature_colname) == target_celltype)
  
  if (nrow(reference_cells) == 0 || nrow(target_cells) == 0) {
    cat("[WARNING] No cells found for reference or target.\n")
    return(NULL)
  }

  reference_coords <- reference_cells %>% select(Cell.X.Position, Cell.Y.Position)
  target_coords <- target_cells %>% select(Cell.X.Position, Cell.Y.Position)
  
  reference_interactions <- frNN(as.matrix(reference_coords), eps = radius, sort = FALSE)
  target_interactions <- frNN(as.matrix(target_coords), eps = radius, query = as.matrix(reference_coords), sort = FALSE)

  ref_ref_count <- sum(lengths(reference_interactions$id)) / 2
  ref_target_count <- sum(lengths(target_interactions$id))
  
  if (ref_ref_count == 0) {
    cat("[WARNING] No reference-reference interactions found.\n")
    return(NULL)
  }

  mixing_score <- ref_target_count / ref_ref_count
  normalized_mixing_score <- 0.5 * mixing_score * (nrow(reference_cells) - 1) / nrow(target_cells)

  # Permutation test: randomly shuffle the cell type labels to generate null distribution
  set.seed(8)
  permuted_scores <- replicate(n_permutations, {
    permuted_data <- data %>%
      mutate(permuted_celltype = sample(!!sym(feature_colname)))
    
    perm_ref <- permuted_data %>% filter(permuted_celltype == reference_celltype)
    perm_tar <- permuted_data %>% filter(permuted_celltype == target_celltype)
    
    if (nrow(perm_ref) > 0 && nrow(perm_tar) > 0) {
      perm_ref_coords <- perm_ref %>% select(Cell.X.Position, Cell.Y.Position)
      perm_tar_coords <- perm_tar %>% select(Cell.X.Position, Cell.Y.Position)
      
      perm_ref_ref <- frNN(as.matrix(perm_ref_coords), eps = radius, sort = FALSE)
      perm_ref_tar <- frNN(as.matrix(perm_tar_coords), eps = radius, query = as.matrix(perm_ref_coords), sort = FALSE)

      ref_ref_count <- sum(lengths(perm_ref_ref$id)) / 2
      ref_tar_count <- sum(lengths(perm_ref_tar$id))

      if (ref_ref_count > 0) {
        return(0.5 * (ref_tar_count / ref_ref_count) * (nrow(perm_ref) - 1) / nrow(perm_tar))
      }
    }
    return(0)
  })

  p_value <- mean(permuted_scores >= normalized_mixing_score, na.rm = TRUE)

  return(data.frame(
    Reference_Cell_Type = reference_celltype,
    Target_Cell_Type = target_celltype,
    Radius = radius,
    Mixing_Score = mixing_score,
    Normalized_Mixing_Score = normalized_mixing_score,
    P_Value_One_Sided = p_value,
    Reference_Reference_Interactions = ref_ref_count,
    Reference_Target_Interactions = ref_target_count,
    Number_of_Reference_Cells = nrow(reference_cells),
    Number_of_Target_Cells = nrow(target_cells)
  ))
}
# ---------------------------------

# --------- Main Processing Loop ---------
composite_files <- list.files(input_directory, pattern = file_pattern, full.names = TRUE)

if (length(composite_files) == 0) {
  stop("[ERROR] No composite files found in the input directory.")
}

for (file_path in composite_files) {
  file_name <- basename(file_path)
  cat("[INFO] Processing:", file_name, "\n")

  # Extract sample/gene/clonotype
  parts <- strsplit(file_name, "_TRA_")[[1]]
  if (length(parts) != 2) next
  sample_gene <- strsplit(parts[1], "_")[[1]]
  gene <- tail(sample_gene, 1)
  sample_id <- paste(head(sample_gene, -1), collapse = "_")
  clonotype <- gsub("_composite.csv", "", parts[2])
  
  # Output file
  output_file <- file.path(output_directory, paste0(sample_id, "_", gene, "_", clonotype, "_TRA_mixing_scores.csv"))
  if (file.exists(output_file)) {
    cat("[INFO] Already exists, skipping:", output_file, "\n")
    next
  }

  # Read input
  data <- read.csv(file_path)

  # Define cell type pairs
  pairs <- list(
    c(paste0(gene, "_MUT"), clonotype),
    c(paste0(gene, "_MUT"), "other_clonotype"),
    c(paste0(gene, "_WT"), clonotype),
    c(paste0(gene, "_WT"), "other_clonotype")
  )

  results <- do.call(rbind, lapply(pairs, function(pair) {
    do.call(rbind, lapply(seq(10, 400, by = 10), function(radius) {
      mixing_score_manual(data, pair[1], pair[2], radius, feature_colname = "Cell.Type", n_permutations = 1000)
    }))
  }))

  # Save
  if (!is.null(results) && nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    cat("[INFO] Saved results to:", output_file, "\n")
  }
}

cat("[INFO] All files processed.\n")
# ---------------------------------
