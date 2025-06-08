# Author: Amit Sud  
# Date: 1st May 2025  
# Description: This script processes high-confidence TCR clonotype data (TRA and TRB),
#              filters for quality, reshapes the data to wide format, and merges it with an AnnData object.
# Input:
#   - TCR file (tab-delimited with clonotype information)
#   - AnnData object: adata (with cell barcodes as index)
# Output:
#   - Updates adata.obs with columns for TRA/TRB CDR3s, reads, and UMIs

import pandas as pd

# === Step 1: Load the TCR data ===
tcr_file = "path/to/tcr_clonotypes.txt"  # <-- Replace with actual path
tcr_data = pd.read_csv(tcr_file, sep="\t")

# === Step 2: Filter for high-confidence, productive TCRs for TRA and TRB chains ===
tcr_filtered = tcr_data[
    (tcr_data['high_confidence'] == True) &
    (tcr_data['productive'] == True) &
    (tcr_data['is_cell'] == True) &
    (tcr_data['chain'].isin(['TRA', 'TRB']))
][['barcode', 'chain', 'cdr3', 'reads', 'umis']]

# Keep the entry with the highest UMI count per barcode and chain
tcr_filtered = (
    tcr_filtered
    .sort_values('umis', ascending=False)
    .drop_duplicates(subset=['barcode', 'chain'])
)

# === Step 3: Pivot to wide format (1 row per barcode, TRA and TRB columns) ===
tcr_wide = tcr_filtered.pivot(index='barcode', columns='chain', values=['cdr3', 'reads', 'umis'])

# Flatten MultiIndex columns
tcr_wide.columns = [f"{chain}_{col}" for col, chain in tcr_wide.columns]

# === Step 4: Merge with AnnData object ===
merged_df = adata.obs.merge(tcr_wide, how='left', left_index=True, right_on='barcode')
merged_df.set_index(adata.obs.index, inplace=True)
adata.obs = merged_df
