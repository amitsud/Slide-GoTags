# Slide-GoTags
Code for [SlideGoTags](doi.org/)

The processing of sequencing data consists of four major steps:
1. Processing of the illumina sequencing data from the spatial library
2. Processing of the illumina sequencing data from the 10x single-nuclei RNA-Seq library to generate gene expression data
3. Processing of nanopore sequencing data from the single-nuclei T-cell receptor library
4. Processing of nanopore sequencing data from the single-nuclei genotyping library

# 1. Processing of the illumina sequencing reads from the spatial library

Please see the [slide-tags github repo](https://github.com/broadchenf/Slide-tags) for scripts on how to process the spatial library following illumina sequencing.

# 2. Processing of the illumina sequencing reads from the 10x single-nuclei RNA-Seq library to generate gene expression data

Please see [cellranger](https://github.com/10XGenomics/cellranger) for generation of the barcode and feature count matrix for single-nuclei.
Please see [cellbender](https://github.com/broadinstitute/CellBender) for denoising of the single-nuclei gene expression data.

# 3. Processing of ONT sequencing reads from the single-cell RNA-Seq T-cell receptor library

This consists of four steps:
1. Download and install [nanoranger](https://github.com/mehdiborji/nanoranger) and the T-cell receptor (TCR) reference library from [VDJ.online](https://vdj.online/library).
2. Run nanoranger to align 5' single-nuclei TCR-seq ONT sequencing reads and extract relevant barcode and UMI data (for an example script see TCR_nanoranger.sh).
3. UMI correction and quality control of single-nuclei TCR sequencing using this [R script](https://github.com/amitsud/Slide-GoTags/blob/main/process_scTCRseq.R).
4. Assignment of TCRA and TCRB to single-nuclei - user defined.

# 4. Processing of ONT sequencing reads from the single-cell RNA-Seq genotyping library (details included in this repository)

This consists of five steps:
1. Download and install [nanoranger](https://github.com/mehdiborji/nanoranger).
2. Run nanoranger to align 5' single-nuclei ONT genotyping reads and extract relevant barcode and UMI data.
3. UMI correction, quality control and matching of genotyping barcode/UMI with gene expresson barcode/UMI using this [R script](https://github.com/amitsud/Slide-GoTags/blob/main/process_scMUTseq.R).
4. Denoising and assignment of mutations to single-nuclei - user defined.
