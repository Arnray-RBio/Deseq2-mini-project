# RNA-seq DESeq2 Mini Project

## Overview
This is a beginner-friendly RNA-seq differential expression analysis project using a simulated count matrix. The aim is to practise the downstream RNA-seq workflow from count matrix to differential expression results and visualisation.

## Workflow
1. Create a fake RNA-seq count matrix and sample metadata
2. Run DESeq2 differential expression analysis
3. Generate PCA plot, volcano plot, and heatmap

## Project Structure
- Deseq2 scripts/: R scripts used in the analysis
- Data/: simulated count matrix and metadata
- Results/: DESeq2 output table
- Figures(PCA, volcano, heatmap)/: generated plots

## Requirements
R packages(install them before starting with the scripts) :
- DESeq2
- tidyverse
- pheatmap
- ggrepel
