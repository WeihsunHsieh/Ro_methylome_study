# Epigenomic and Transcriptomic Analyses of Black Raspberry During Fruit Ripening

This repository contains custom scripts and workflow descriptions used for the analysis of RNA-seq and whole-genome bisulfite sequencing (BS-seq) data in our study. 
The scripts are provided to ensure transparency and reproducibility of the analyses described in the manuscript.

---

* Paper: *submitted, currently under review*
* Submitted version on bioRxiv: https://doi.org/10.1101/2025.11.24.690302
* Code: https://github.com/

---

# Overview of Analytical Workflow

## BS-seq Analysis

1. Quality check of raw reads by FastQC, and adapter and quality trimming with Trimmomatic  
2. Mapping of bisulfite sequencing reads using BS-Seeker2  
3. Generating methylation calling by custom perl scripts
4. Identification of differentially methylated regions (DMRs) using CGmapTools  
5. Genomic annotation of DMRs using ChIPseeker  
6. Definition of differentially methylated genes (DMGs) by BEDTools

## RNA-seq Analysis

1. Quality check of raw reads by FastQC, and adapter and quality trimming with Trimmomatic  
2. Read alignment to the reference genome with HISAT2
3. Read counts per gene by HTSeq 
4. Differential expression analysis using DESeq2 
5. Identification of differentially expressed genes (DEGs)
6. WGCNA for identification of regulatory modules

## Integrative Analysis

- Overlap analysis between DMGs and DEGs using BEDTools
- Correlation and clustering analyses
- Principal component analysis (PCA)
- Methylation profiling visualization using deepTools

---

# Software and Tools

The following software packages were used:

- BLAST (v2.2.26)
- Blast2GO (v5.2.5)
- R (v4.5.0)
- DESeq2 (v1.34.0)
- BS-Seeker2 (v2.1.8)
- CGmapTools (v0.1.2)
- BEDTools (v2.30.0)
- deepTools (v3.5.4)
- ChIPseeker (v1.36.0 within Bioconductor v3.17)
- WGCNA 

Please refer to individual scripts for detailed command usage.

---

# BS-seq Processing

Bisulfite sequencing reads were aligned using BS-Seeker2.

DMRs were identified using CGmapTools with study-specific thresholds defined in the manuscript.

Genomic annotation of DMRs was performed using the ChIPseeker package in R.

Promoter regions were defined as 1.5 kb upstream of the transcription start site (TSS).

DMGs were defined as genes overlapping DMRs located in promoter or gene body regions.

---

# RNA-seq Processing

Quality control of raw reads using FastQC.  

Adapter trimming using Trimmomatic.  

Mapping to the reference genome using HISAT2.  

Gene-level read counting using HTSeq.  

Differential expression analysis was performed using DESeq2 in R.

Criteria for defining DEGs:

- |log2 fold change| ≥ 2 
- Adjusted p-value (FDR) < 0.05  

---

# BLAST Analysis

Protein sequence similarity searches were performed using BLASTP with the following parameters:

- E-value cutoff: 1e-4  
- Output format: tabular (format 8)  
- Maximum alignments reported: 1  
- Threads: 16  

Example command:

```
blastall -p blastp \
-i query.fasta \
-d protein_database \
-e 1e-4 \
-m 8 \
-a 16 \
-b 1 \
-v 1 \
-K 1 \
-o output.txt
```

Functional annotation was performed using Blast2GO with default settings.

---

# PCA Analysis

Principal component analysis (PCA) was performed using the `prcomp` function in R.

- RNA expression values were transformed using log10(x + 1).  
- Average methylation levels were transformed using log10(x + 1). 
- The first three principal components were retained.  
- 3D visualization was generated using the 'plotly' package.  

---

# Correlation Heatmaps

Pairwise correlation heatmaps were generated in R using normalized expression count values.

---

# Gene Expression Heatmaps

Gene expression heatmaps were generated using dChip software with default clustering parameters.

---

# DMG–DEG Overlap Analysis

Overlap between DMGs and DEGs was determined using BEDTools intersect.

---

# Methylation Profiling

Genome-wide methylation profiling plots were generated using deepTools:

- computeMatrix  
- plotProfile  

---

# GO Enrichment Analysis

GO enrichment was performed using the 'goseq' R package
with hypergeometric testing and Benjamini–Hochberg FDR correction.

Significance threshold: FDR < 0.05

---

# Co-expression Network Analysis (WGCNA)

Weighted gene co-expression network analysis (WGCNA) was performed in R using normalized expression values derived from DESeq2.

Analysis workflow including: 

Filtering low-expression genes,

Soft-threshold power selection, based on scale-free topology fit index (R² > 0.8),

Construction of adjacency and TOM matrices,

Module detection using dynamic tree cutting,

Module–trait correlation analysis.

---

# Repository Structure

The repository is organized by analysis category:

- RNA-seq related scripts  
- BS-seq related scripts  
- Integration analyses  
- Visualization scripts  

---

# Data Availability

All custom scripts used in this study are publicly available in this repository.  
Raw sequencing data are available under the accession number provided in the manuscript.
