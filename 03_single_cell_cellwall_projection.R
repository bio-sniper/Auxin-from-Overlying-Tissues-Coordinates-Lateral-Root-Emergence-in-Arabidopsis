#################################################
#Project:
#Auxin from Overlying Tissues Coordinates Lateral Root Emergence

#Script:
#GO enrichment and single cell analysis
#Author:Lei Li

#Date:2026/07/15

#load required packages
library (tidyverse)
library(ggupset)
library(enrichplot)
library(DOSE)
library(scales)
library(org.At.tair.db)
library(clusterProfiler)
library(readxl)
library(edgeR)
library(Seurat)
# Purpose 3. Projection of pWER-specific cell wall DEGs onto a single cell atlas. 
# 3.1. Read single-cell reference atlas: 
## Reference: Vascular transcription factors guide plant epidermal responses to limiting phosphate conditions
# Dataset: GSE141730
# Downloaded from GEO and stored in the local working directory

file1 <- Read10X_h5(
  "GSM4212550_BDR3_filtered_gene_bc_matrices.h5"
)

file2 <- Read10X_h5(
  "GSM4212551_BDR4_filtered_gene_bc_matrices.h5"
)

file3 <- Read10X_h5(
  "GSM4212552_BDR5_filtered_gene_bc_matrices.h5"
)

obj1 <- CreateSeuratObject(counts = file1, project = "GSE141730", min.cells = 3, min.features = 200)
obj2 <- CreateSeuratObject(counts = file2, project = "GSE141730", min.cells = 3, min.features = 200)
obj3 <- CreateSeuratObject(counts = file3, project = "GSE141730", min.cells = 3, min.features = 200)

obj1$orig.ident <- "sample1"
obj2$orig.ident <- "sample2"
obj3$orig.ident <- "sample3"

# 3.2. Merge samples and perform the Seurat analysis workflow
# ========================================
merged <- merge(obj1, y = list(obj2, obj3),
                add.cell.ids = c("S1","S2","S3"),
                project = "GSE141730")

set.seed(seed = 21)

merged <- SCTransform(merged, verbose = FALSE)
merged <- RunPCA(merged)
merged <- RunUMAP(merged, dims = 1:30, spread = 0.8, seed.use = 21)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.245, random.seed = 21)

# 3.3. Set cell wall gene set.
# =======================================================
cell_wall_gene_set <- c ("AT4G13390",
                         "AT5G06640",
                         "AT4G08410",
                         "AT5G06630",
                         "RIC1",
                         "AT2G24980",
                         "AT4G08400",
                         "AT5G35190",
                         "EXT2",
                         "ATEXPA12",
                         "PME46",
                         "EXPA18")

# 3.3. Calculate cell wall gene set activity score.
# ====================================================
merged <- AddModuleScore(
  merged,
  features = list(cell_wall_gene_set),
  name = "CellWall"
)

FeaturePlot(
  merged,
  features = "CellWall1",
  pt.size = 2
)

ggsave("cell wall gene projection from WER down.tif",
       height = 200,
       width = 240,
       units = "mm",
       dpi = 400,
       device = "tiff")

# Output: Figure S2C

