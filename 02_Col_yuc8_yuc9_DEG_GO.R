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
library(limma)
library(org.At.tair.db)
library(clusterProfiler)
library(ggvenn)
library(readxl)
library(edgeR)

## Purpose 2: Analysis of Col vs. yuc8 yuc9 after 48 h gravitropism treatment.
# 2.1. Read expression data
# ============================================================

Col_data <- read.csv(
  "Col_0.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

yuc89_data <- read.csv(
  "yuc8_yuc9.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# ============================================================
# 2.2. Define experimental groups
# ============================================================

group <- factor(
  c(
    "Timepoint_0h",
    "Timepoint_0h",
    "Timepoint_0h",
    "Timepoint_48h",
    "Timepoint_48h",
    "Timepoint_48h"
  )
)

# ============================================================
# 2.3. Build design matrix
# ============================================================

design <- model.matrix(~ group)

colnames(design) <- levels(group)

design

# 2.4. Fit linear model
# ===========================================================
fit_Col <- lmFit(Col_data, design)
fit_yuc89 <- lmFit(yuc89_data, design)

fit_Col <- eBayes(fit_Col)
fit_yuc89 <- eBayes(fit_yuc89)

# 2.5. Extract gene list
# =========================================================

res_Col <- topTable(fit_Col, 
                    coef = 2, #48 h vs. 0 h
                    number=Inf,
                    adjust.method = "BH",
                    sort.by = "P")

deg_Col <- res_Col[which(res_Col$adj.P.Val < 0.05), ]

up_Col <- deg_Col[deg_Col$logFC > 0, ]

down_Col <- deg_Col[deg_Col$logFC < -0, ]

res_yuc89 <- topTable(fit_yuc89, 
                      coef = 2, 
                      number=Inf,
                      adjust.method = "BH",
                      sort.by = "P")

deg_yuc89 <- res_yuc89[which(res_yuc89$adj.P.Val < 0.05), ]

up_yuc89 <- deg_yuc89[deg_yuc89$logFC > 0, ]

down_yuc89 <- deg_yuc89[deg_yuc89$logFC < -0, ]

# 2.6. Analysis of upregulated genes.
# =======================================================
venn_data_up <- list(
  Col_up = rownames(up_Col),
  yuc89_up = rownames(up_yuc89)
)

ggvenn(
  venn_data_up,
  fill_color = c ("#FF7256","#66CD00"),
  show_percentage = FALSE,
  stroke_alpha = 0.5,
  stroke_size = 0.3
)

ggsave("overlapping of up regulated genes.tif",
       width = 100,
       height = 80,
       units = "mm",
       dpi = 400,
       device = "jpeg")

# 2.7. Analysis of down-regulated transcripts.
# ================================================
venn_data_down <- list(
  Col_down = rownames(down_Col),
  yuc89_down = rownames(down_yuc89)
)

ggvenn(
  venn_data_down,
  fill_color = c ("#FF7256","#66CD00"),
  show_percentage = FALSE,
  stroke_alpha = 0.5,
  stroke_size = 0.3
)

ggsave("overlapping of down regulated genes.tif",
       width = 100,
       height = 80,
       units = "mm",
       dpi = 400,
       device = "jpeg")

# 2.8. Analyze unique DEGs of the yuc8 yuc9 mutant.
# ==============================================================
yuc89_unique_up <- setdiff((rownames(up_yuc89)), rownames(up_Col))
yuc89_unique_down <- setdiff((rownames(down_yuc89)), rownames(down_Col))

yuc89_unique_up_go_enrichment <- enrichGO(
  gene = yuc89_unique_up,
  OrgDb = org.At.tair.db,
  ont = 'BP',
  keyType = 'TAIR',
  readable = FALSE,
  pvalueCutoff = 0.05)

write.csv(yuc89_unique_up_go_enrichment,
          file = 'yuc89 unique upregulated GO results.csv',
          row.names = F)

dotplot(yuc89_unique_up_go_enrichment)

ggsave('yuc89 unique up GO.tif',
       plot = last_plot(),
       device = 'tiff',
       width= 160,
       height= 120,
       units = 'mm',
       dpi = 400)

yuc89_unique_down_go_enrichment <- enrichGO(
  gene = yuc89_unique_down,
  OrgDb = org.At.tair.db,
  ont = 'BP',
  keyType = 'TAIR',
  readable = FALSE,
  pvalueCutoff = 0.05)

write.csv(yuc89_unique_down_go_enrichment,
          file = 'yuc89 unique downregulated GO results.csv',
          row.names = F)

dotplot(yuc89_unique_down_go_enrichment)

ggsave('yuc89 unique down GO.tif',
       plot = last_plot(),
       device = 'tiff',
       width= 160,
       height= 120,
       units = 'mm',
       dpi = 400)

# Output: Figure S3