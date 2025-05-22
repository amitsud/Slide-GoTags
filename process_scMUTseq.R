# NanoRanger Genotyping UMI Correction and GEX Matching Pipeline
# Author: Amit Sud
# Date: 1st May 2025
# Description: This script processes mutation-level data from NanoRanger, applies UMI error correction,
#              and matches results to gene expression UMIs from CellRanger output to assign genotypes.

# ---- Required Libraries ----
library(dplyr)
library(stringdist)
library(parallel)
library(rhdf5)
library(vroom)

# ---- Input files (placeholders to be customized) ----
cellbender_bc <- "<path_to_cellbender_barcodes.tsv>"
h5file <- "<path_to_molecule_info.h5>"
gene <- "<target_gene_name>"
got_df <- vroom("<path_to_pileup_output.tsv>")
cellbender_filtered_h5 <- "<path_to_cellbender_filtered.h5>"

# ---- Load barcode list from CellBender ----
seurat_bcs <- scan(file = cellbender_bc, what = "character", quiet = TRUE)
seurat_bcs <- gsub("-.*", "", seurat_bcs)

# ---- Function: Decode numeric UMI values to nucleotide sequences ----
decode_umi <- function(umi_numeric) {
  int_to_nucleotide <- c("A", "C", "G", "T")
  numeric_to_binary <- function(x) as.character(intToBits(x))[1:20]
  binary_to_nucleotide <- function(bits) {
    decoded <- ""
    for (pos in seq(1, length(bits), by = 2)) {
      bit_pair <- sum(as.integer(bits[pos:(pos + 1)]) * c(1, 2))
      decoded <- paste0(int_to_nucleotide[bit_pair + 1], decoded)
    }
    decoded
  }
  sapply(umi_numeric, function(x) binary_to_nucleotide(numeric_to_binary(x)), USE.NAMES = FALSE)
}

# ---- Function: Error correct UMIs within each barcode group ----
network_umi_error_correction <- function(bc_data) {
  results <- tibble(bc = character(), founder_umi = character(), REF_count = integer(), ALT_count = integer(), founder_number = integer(), umi = character())
  founder_number <- 0
  if (nrow(bc_data) == 0) return(results)
  
  while (nrow(bc_data) > 0) {
    if (nrow(bc_data) < 2) {
      bc_data$founder_umi <- bc_data$umi
      bc_data$founder_number <- NA
      results <- bind_rows(results, bc_data %>% select(-total))
      break
    }
    distances <- stringdist::stringdistmatrix(bc_data$umi, bc_data$umi, method = "lv")
    matches <- distances <= 2
    match_counts <- rowSums(matches) - 1
    
    if (length(match_counts) == 0 || all(is.na(match_counts)) || max(match_counts) == 0) {
      bc_data$founder_umi <- bc_data$umi
      bc_data$founder_number <- NA
      results <- bind_rows(results, bc_data %>% select(-total))
      break
    }
    
    max_matches <- max(match_counts)
    founder_indices <- which(match_counts == max_matches)
    
    for (founder_index in founder_indices) {
      group_indices <- which(matches[founder_index, ])
      results <- bind_rows(results, tibble(
        bc = bc_data$bc[founder_index],
        founder_umi = bc_data$umi[founder_index],
        REF_count = sum(bc_data$REF_count[group_indices]),
        ALT_count = sum(bc_data$ALT_count[group_indices]),
        founder_number = founder_number + 1,
        umi = paste(bc_data$umi[group_indices], collapse = ", ")
      ))
    }
    founder_number <- founder_number + 1
    bc_data <- bc_data[-unique(unlist(lapply(founder_indices, function(fi) which(matches[fi, ])))), ]
  }
  results
}

# ---- Load GEX molecules from CellRanger HDF5 ----
gex_molecules <- list(
  barcode_idx = h5read(h5file, "barcode_idx") + 1,
  barcodes = h5read(h5file, "barcodes"),
  umi_counts = h5read(h5file, "count"),
  feature_idx = h5read(h5file, "feature_idx") + 1,
  feature_name = h5read(h5file, "features/name"),
  feature_id = h5read(h5file, "features/id"),
  umi = h5read(h5file, "umi")
)
h5closeAll()

# ---- Filter for target gene ----
target_gene_idx <- which(gex_molecules$feature_name == gene)
target_gene_id <- gex_molecules$feature_id[target_gene_idx]
target_gene_entries <- which(gex_molecules$feature_idx == target_gene_idx)

df_10x <- tibble(
  BC = gex_molecules$barcodes[gex_molecules$barcode_idx[target_gene_entries]],
  UMI = unlist(mclapply(gex_molecules$umi[target_gene_entries], mc.cores = detectCores(), FUN = decode_umi)),
  counts = gex_molecules$umi_counts[target_gene_entries]
) %>% mutate(BC_UMI = paste0(BC, "_", UMI))

# ---- Prepare and filter mutation pileup data ----
got_df <- got_df %>%
  group_by(umi, bc) %>%
  summarise(total = n(),
            REF_count = sum(base == "A", na.rm = TRUE),
            ALT_count = sum(base == "G", na.rm = TRUE),
            perc_REFALT = (REF_count + ALT_count) / total,
            .groups = 'drop') %>%
  filter(perc_REFALT > 0.95) %>%
  select(bc, umi, REF_count, ALT_count, total)

# ---- Apply UMI correction across barcodes ----
bc_groups <- split(got_df, got_df$bc)
founder_results <- bind_rows(lapply(bc_groups, network_umi_error_correction))

# ---- Filter low-support UMIs ----
filtered_founder_results <- founder_results %>% filter(REF_count > 1 | ALT_count > 1)
filtered_founder_results <- filtered_founder_results %>%
  mutate(BC_UMI = paste0(bc, "_", founder_umi),
         Exact_BC_UMI_Gene_Match = BC_UMI %in% df_10x$BC_UMI,
         Approx_BC_UMI_Gene_Match = unlist(mclapply(BC_UMI, mc.cores = detectCores(), FUN = function(x){
           ain(x, df_10x$BC_UMI, method = "lv", maxDist = 2)
         })))

# ---- Collapsing multiple gene assignments ----
target_bc_idx <- which(gex_molecules$barcodes %in% got_df$bc)
molecules <- gex_molecules$barcode_idx %in% target_bc_idx

df_all_gene <- tibble(
  BC_IDX = gex_molecules$barcode_idx[molecules],
  BC = gex_molecules$barcodes[gex_molecules$barcode_idx[molecules]],
  UMI_bin = gex_molecules$umi[molecules]
)
df_all_gene <- df_all_gene %>%
  mutate(UMI = unlist(lapply(UMI_bin, decode_umi)),
         BC_UMI = paste0(BC, "_", UMI),
         gene_idx = gex_molecules$feature_idx[molecules],
         gene = gex_molecules$feature_name[gene_idx],
         count = gex_molecules$umi_counts[molecules])

to_keep <- df_all_gene$BC_UMI %in% filtered_founder_results$BC_UMI
df_all_gene_collapse <- df_all_gene[to_keep,]

# Resolve conflicts in BC_UMI to gene mapping
for (k in unique(df_all_gene_collapse$BC_UMI)) {
  target_rows <- which(df_all_gene_collapse$BC_UMI == k)
  if (length(target_rows) > 1) {
    sub_df <- df_all_gene_collapse[target_rows,]
    if (gene %in% sub_df$gene) {
      sub_df$gene <- ifelse(any(grepl("_CITE", sub_df$gene)),
                            paste0("Multiple_", gene, "_CITE"),
                            paste0("Multiple_", gene))
    } else {
      sub_df$gene <- "Multiple"
    }
    df_all_gene_collapse[target_rows[1],] <- sub_df[1,]
    df_all_gene_collapse <- df_all_gene_collapse[-target_rows[-1],]
  }
}

# ---- Join with corrected UMI table ----
filtered_founder_results_gene <- merge(filtered_founder_results, df_all_gene_collapse[,c("BC_UMI", "gene")], by = "BC_UMI", all.x = TRUE)
filtered_founder_results_gene <- filtered_founder_results_gene %>%
  mutate(TOTAL_count = REF_count + ALT_count,
         In_GEX = !is.na(gene),
         Gene_Group = case_when(
           Exact_BC_UMI_Gene_Match ~ "Exact",
           Approx_BC_UMI_Gene_Match ~ "Approx",
           In_GEX ~ "Other Gene",
           TRUE ~ "No Gene"
         ))

# ---- Determine minimum threshold from density peak ----
no_gene_counts <- filtered_founder_results_gene %>% filter(Gene_Group == "No Gene") %>% pull(TOTAL_count)
if (length(no_gene_counts) >= 2) {
  d <- density(log10(no_gene_counts))
  threshold <- 10^(optimize(approxfun(d$x,d$y),interval=c(0,3))$minimum)
} else {
  threshold <- NA
  print("Not enough data points to calculate density.")
}

filtered_founder_results_gene <- filtered_founder_results_gene %>%
  mutate(Keep = case_when(
    Gene_Group == "Exact" ~ TRUE,
    Gene_Group == "Other Gene" ~ FALSE,
    TOTAL_count > threshold ~ TRUE,
    TRUE ~ FALSE
  ))

# ---- Genotype Assignment ----
filtered_founder_results_gene_keep <- filtered_founder_results_gene %>% filter(Keep)
filtered_founder_results_gene_keep_umi_calls <- filtered_founder_results_gene_keep %>%
  mutate(UMI_genotype = ifelse(ALT_count > REF_count, "MUT", "WT"))

filtered_founder_results_gene_keep_bc_calls <- filtered_founder_results_gene_keep_umi_calls %>%
  group_by(bc) %>%
  summarise(
    WT.calls = sum(UMI_genotype == "WT"),
    MUT.calls = sum(UMI_genotype == "MUT")
  ) %>%
  mutate(bc = paste0(bc, "-1"))



##denoising part 1###
#read 10x input
umi_counts_input_df <- df_10x %>% group_by(BC) %>% summarise(counts = n_distinct(UMI))
colnames(umi_counts_input_df)[colnames(umi_counts_input_df) == "BC"] <- "barcode"
colnames(umi_counts_input_df)[colnames(umi_counts_input_df) == "counts"] <- "umi_count_input"
all_barcodes <- data.frame(barcode = seurat_bcs)
umi_counts_input_df <- merge(all_barcodes, umi_counts_input_df, by = "barcode", all.x = TRUE)
umi_counts_input_df$umi_count_input[is.na(umi_counts_input_df$umi_count_input)] <- 0
umi_counts_input_df$barcode <- paste(umi_counts_input_df$barcode, "-1", sep="")

# Read cellbender output
cellbender_output <- H5Fopen(cellbender_filtered_h5)
barcodes_output <- h5read(cellbender_output, "matrix/barcodes")
counts_output <- h5read(cellbender_output, "matrix/data")
indptr_output <- h5read(cellbender_output, "matrix/indptr")
indices_output <- h5read(cellbender_output, "matrix/indices")
gene_names_output <- h5read(cellbender_output, "matrix/features/name")

counts_spars_output <- sparseMatrix(i = indices_output + 1, 
                                    p = indptr_output, 
                                    x = counts_output, 
                                    dims = c(length(gene_names_output), length(barcodes_output)))

features_df_output <- data.frame(gene = gene_names_output)
str(features_df_output)

gene_index_output <- which(features_df_output$gene == gene)
umi_counts_output <- counts_spars_output[gene_index_output, ]
umi_counts_output_df <- data.frame(barcode = barcodes_output, umi_count_output = as.vector(umi_counts_output))

H5Fclose(cellbender_output)

# Merge the dataframes by 'barcode' and calculate difference
umi_counts_merged_df <- merge(umi_counts_input_df, umi_counts_output_df, by = "barcode")
umi_counts_merged_df$difference <- umi_counts_merged_df$umi_count_output - umi_counts_merged_df$umi_count_input

seurat_bcs <- paste0(seurat_bcs, "-1")
real_droplets_nanopore_df <- filtered_founder_results_gene_keep_bc_calls %>% filter(bc %in% seurat_bcs)

real_droplets_nanopore_df <- real_droplets_nanopore_df %>% left_join(umi_counts_merged_df, by = c("bc" = "barcode"))

# ---- Save Output ----
write.table(
  real_droplets_nanopore_df,
  file = "genotype_bc_calls_output.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# ---- End of Script ----
