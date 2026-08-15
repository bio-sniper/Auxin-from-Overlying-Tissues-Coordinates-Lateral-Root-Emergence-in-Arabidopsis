# Project Title

> Reproducible bioinformatics analysis for Auxin from Overlying Tissues Coordinates Lateral Root Emergence in Arabidopsis project.

## Overview

This repository contains the bioinformatics analysis scripts used in:

**Auxin from Overlying Tissues Coordinates Lateral Root Emergence in Arabidopsis**

The project investigates transcriptional genes involved in lateral root emergence in the Arabidopsis overlying tissues.

The analysis includes:

* RNA-seq data preprocessing
* Differential expression analysis
* Functional enrichment analysis
* Data visualization
* Projection of targeted genes onto public single-cell atlas

## Data

### Input data

The analysis uses Bulk RNA-seq  and scRNA-seq from:

* Species: *Arabidopsis thaliana*
* Tissue: primary roots and curved root segments with 48 h gravitropism treatment
* Experimental groups: 
  Group1:
  6-day old seedlings treated with estradiol 2 days.
  Induced vs. Mock (pWER:XVE:YUC1-TAA1; pSCR:XVE:YUC1-TAA1; pSHR:XVE:YUC1-TAA1;   pAPL:XVE:YUC1-  TAA1)
 
 Group2:
 6-day old seedlings with 48 h gravitropism 
 Before vs. Gravitropism (Col-0 and yuc8 yuc9)

* Biological replicates: 3
* Data source: Generated in this study

Raw sequencing data are available at:

* NCBI BioProject: [PRJXXXXXX]
* NCBI GEO: [GSExxxxxx]
* SRA: [SRRxxxxxx]

*Reference:
  Single-cell data used in this study referenced to *Vascular transcription factors guide plant epidermal   responses to limiting phosphate conditions* with GEO number: GSE141730.

*Processed expression data are provided as:

* FPKM matrix

> Large raw sequencing files and intermediate files are not stored in this repository.

## Analysis Workflow

The main analysis workflow is:

Gene expression matrix
       │
       ▼
Differential expression analysis
       │
       ├── DEG
       │
       ├── GO enrichment
       │
       ├──  Projection targeted gene set to single-cell atlas
       │
       ▼
Visualization
       │
       ▼
Final figures
```

## Repository Structure

```text
.
├── README.md
├── scripts/
├── 01_SCR_WER_DEG_GO.R
├── 02_Col_yuc8 yuc9_DEG_GO.R
├── 03_single_cell_cellwall_projection.R
├── data/
│   └── README.md
└── environment/
    └── packages.txt
```

## Analysis Scripts

### `Script for Auxin from Overlying Tissues Coordinates Lateral Root Emergence in Arabidopsis.R`

Performs differential expression analysis using limma.

Main steps:

1. Import expression matrix
2. Define experimental groups
3. Construct design matrix
4. Perform differential expression analysis
5. GO_enrichemt of up and downregulated genes
6. Clustering public single-cell atlas
7. Projection targeted cell-wall genes onto single-cell atlas

## Software and Packages

A complete package list is provided in:

`environment/packages.txt`

## Reproducibility

All analyses were performed using the scripts provided in this repository.

The scripts are numbered according to the recommended execution order:

01 → 02 → 03

Input files and sample names should follow the format described in the corresponding scripts.

For reproducibility, the analysis environment and package versions are documented in `environment/`.

## Results

The scripts generate the main analysis results and figures reported in the manuscript.

## Citation

If you use these scripts or analyses, please cite:

> Author: Lei Li (李磊), Xinfang Zhu, Zhaojun Ding, Xiangpei Kong
   Auxin from Overlying Tissues Coordinates Lateral Root Emergence in Arabidopsis. 
   Journal of Integrative Plant Biology, 
   2026.

## Contact

For questions regarding this repository or the analysis, please contact:

**Lei Li (李磊)**

Shandong University 

Email: lileipower415@gmail.com
