# TCR UMI Error Correction and Processing Pipeline
# Author: Amit Sud
# Date: 1st May 2025
# Description: This script processes output from NanoRanger scTCR-seq to correct UMI errors,
#              annotate and filter TCR clonotypes, and export a cleaned result table.

# ---- Load required libraries ----
library(dplyr)
library(stringdist)
library(ggplot2)
library(gridExtra)

# ---- Load input data ----
# bcumi: nanoranger output file ending in "bcumi.csv.gz"
# clones: nanoranger file ending in "clones.txt.gz"
# cell_barcodes_cellbender: cell barcode file from cellbender output (tab-separated)
# clones_productive: nanoranger file ending in "clones.txt.gz" with only productive clones

bcumi <- vroom("<path_to_bcumi.csv.gz>")
clones <- vroom("<path_to_clones.txt.gz>")
cell_barcodes_cellbender <- vroom("<path_to_barcodes.tsv>", delim = "\t", col_names = FALSE)
clones_productive <- vroom("<path_to_productive_clones.txt.gz>")

# ---- Aggregate read counts by UMI ----
aggregated_bcumi <- bcumi %>%
  group_by(bc, umi, chains, cloneId) %>%
  summarise(reads = n(), .groups = 'drop')

# ---- Define UMI network-based error correction function ----
network_umi_error_correction <- function(aggregated_bcumi) {
  results <- tibble(
    bc = character(),
    founder_umi = character(),
    chains = character(),
    cloneId = character(),
    grouped_umis = character(),
    reads = integer()
  )
  
  for (bc in unique(aggregated_bcumi$bc)) {
    bc_subset <- aggregated_bcumi %>% filter(bc == !!bc)
    for (chain in unique(bc_subset$chains)) {
      chain_subset <- bc_subset %>% filter(chains == !!chain)
      for (cloneId in unique(chain_subset$cloneId)) {
        clone_subset <- chain_subset %>% filter(cloneId == !!cloneId)
        
        while (nrow(clone_subset) > 0) {
          distances <- stringdistmatrix(clone_subset$umi, clone_subset$umi, method = "lv")
          matches <- distances <= 2
          match_counts <- rowSums(matches) - 1

          if (max(match_counts) == 0) {
            results <- bind_rows(results, clone_subset %>%
              mutate(founder_umi = umi, grouped_umis = "") %>%
              select(bc, founder_umi, chains, cloneId, grouped_umis, reads))
            break
          }

          founder_index <- which(match_counts == max(match_counts))[1]
          group_indices <- which(matches[founder_index,])
          total_reads <- sum(clone_subset$reads[group_indices])
          grouped <- unique(clone_subset$umi[group_indices])

          results <- bind_rows(results, tibble(
            bc = bc,
            founder_umi = clone_subset$umi[founder_index],
            chains = chain,
            cloneId = as.character(cloneId),
            grouped_umis = paste(grouped, collapse = ", "),
            reads = total_reads
          ))

          clone_subset <- clone_subset[-group_indices, ]
        }
      }
    }
  }
  return(results)
}

# ---- Run correction ----
corrected_bcumi <- network_umi_error_correction(aggregated_bcumi)

# ---- Merge with clone info ----
clones <- clones %>% filter(chains %in% c("TRA", "TRB"))
corrected_bcumi_2 <- corrected_bcumi %>%
  group_by(bc, cloneId, chains) %>%
  summarise(umis = n_distinct(founder_umi), reads = sum(reads), .groups = 'drop') %>%
  filter(chains %in% c("TRA", "TRB"))

merged_df <- corrected_bcumi_2 %>%
  inner_join(clones, by = c("cloneId", "chains")) %>%
  mutate(barcode = paste0(bc, "-1")) %>%
  group_by(barcode) %>%
  mutate(contig_id = paste0(barcode, "_contig_", dense_rank(cloneId))) %>%
  ungroup() %>%
  mutate(is_cell = barcode %in% cell_barcodes_cellbender$X1)

# ---- Plot UMIs per barcode for each chain ----
line_colors <- c("red", "blue", "green", "pink")

process_and_plot <- function(chain_value) {
  chain_data <- corrected_bcumi_2 %>% filter(chains == chain_value) %>%
    mutate(rank = dense_rank(desc(umis)))

  ggplot(chain_data, aes(x = log(rank), y = log(umis))) +
    geom_point() +
    labs(
      x = "Log of bc Rank (ranked by number of reads)",
      y = "Log of reads",
      title = paste("Log Scale of reads per bc-cloneId-chains Pair for", chain_value)
    ) +
    theme_minimal() +
    geom_hline(yintercept = log(c(20,10,5,2)), linetype = "dashed", color = line_colors)
}

plots <- lapply(unique(corrected_bcumi_2$chains), process_and_plot)
grid.arrange(grobs = plots, ncol = 2)

# ---- Annotate high confidence cells ----
# change this threshold to reflect the distribution of number of UMIs per cell barcode
merged_df <- merged_df %>%
  mutate(high_confidence = umis > 2 & chains %in% c("TRA", "TRB"))

# ---- Rename and structure to match downstream expectations ----
merged_df <- merged_df %>%
  rename_with(~ c(
    "barcode", "is_cell", "contig_id", "high_confidence", "length", "chain",
    "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
    "fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
    "fwr3", "fwr3_nt", "cdr3", "cdr3_nt", "fwr4", "fwr4_nt",
    "reads", "umis", "raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id"
  ), everything()) %>%
  mutate(
    productive = raw_clonotype_id %in% clones_productive$cloneId,
    raw_clonotype_id = paste0("clonotype", raw_clonotype_id)
  ) %>%
  mutate(across(
    c(fwr1, fwr1_nt, cdr1, cdr1_nt, fwr2, fwr2_nt, cdr2, cdr2_nt, fwr3, fwr3_nt, cdr3, cdr3_nt, fwr4, fwr4_nt),
    toupper
  )) %>%
  select(
    "barcode", "is_cell", "contig_id", "high_confidence", "length", "chain",
    "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
    "fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
    "fwr3", "fwr3_nt", "cdr3", "cdr3_nt", "fwr4", "fwr4_nt",
    "reads", "umis", "raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id"
  )

# ---- Export final output (customize this path for your use case) ----
write.table(
  merged_df,
  file = "TCRseq_processed_output.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ---- End of Script ----
