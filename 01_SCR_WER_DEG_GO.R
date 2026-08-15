#################################################
#Project:
#Auxin from Overlying Tissues Coordinates Lateral Root Emergence

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

# Purpose 1. Transcriptional analysis of the est-induced lines
# 1.1. Read data: est-induced lines.
# ====================================================

SCR_data <- read.csv(
  "SCR_est.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

WER_data <- read.csv(
  "WER_est.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# 1.2. Define experimental groups
# ============================================================

group <- factor(
  c(
    "Mock",
    "Mock",
    "Mock",
    "Induced",
    "Induced",
    "Induced"
  )
)

# 1.3. Build design matrix
# ============================================================

design <- model.matrix(~ group)

colnames(design) <- levels(group)

design

# 1.4. Fit linear model
# ===========================================================
fit_SCR <- lmFit(SCR_data, design)
fit_WER <- lmFit(WER_data, design)

fit_SCR <- eBayes(fit_SCR)
fit_WER <- eBayes(fit_WER)

# 1.5. Extract gene list
# =========================================================

res_SCR <- topTable(fit_SCR, 
                    coef = 2, #Mock vs. Induced
                    number=Inf,
                    adjust.method = "BH",
                    sort.by = "P")

deg_SCR <- res_SCR[which(res_SCR$adj.P.Val < 0.05), ]

up_SCR <- deg_SCR[deg_SCR$logFC > 0, ]

down_SCR <- deg_SCR[deg_SCR$logFC < -0, ]

res_WER <- topTable(fit_WER, 
                    coef = 2, 
                    number=Inf,
                    adjust.method = "BH",
                    sort.by = "P")

deg_WER <- res_WER[which(res_WER$adj.P.Val < 0.05), ]

up_WER <- deg_WER[deg_WER$logFC > 0, ]

down_WER <- deg_WER[deg_WER$logFC < -0, ]

# 1.6. Analysis of upregulated genes.
# =======================================================
venn_data_up <- list(
  WER_up = rownames(up_WER),
  SCR_up = rownames(up_SCR)
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

# 1.7. Analysis of down-regulated transcripts.
# ================================================
venn_data_down <- list(
  WER_down = rownames(down_WER),
  SCR_down = rownames(down_SCR)
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

# 1.8. Analyze unique DEGs of the pWER lines.
# ==============================================================
WER_unique_up <- setdiff((rownames(up_WER)), rownames(up_SCR))
WER_unique_down <- setdiff((rownames(down_WER)), rownames(down_SCR))

WER_unique_up_go_enrichment <- enrichGO(
  gene = WER_unique_up,
  OrgDb = org.At.tair.db,
  ont = 'BP',
  keyType = 'TAIR',
  readable = FALSE,
  pvalueCutoff = 0.05)

write.csv(WER_unique_up_go_enrichment,
          file = 'WER unique upregulated GO results.csv',
          row.names = F)

dotplot(WER_unique_up_go_enrichment)

ggsave('WER unique up GO.tif',
       plot = last_plot(),
       device = 'tiff',
       width= 160,
       height= 120,
       units = 'mm',
       dpi = 400)

WER_unique_down_go_enrichment <- enrichGO(
  gene = WER_unique_down,
  OrgDb = org.At.tair.db,
  ont = 'BP',
  keyType = 'TAIR',
  readable = FALSE,
  pvalueCutoff = 0.05)

write.csv(WER_unique_down_go_enrichment,
          file = 'WER unique downregulated GO results.csv',
          row.names = F)

dotplot(WER_unique_down_go_enrichment)

ggsave('WER unique down GO.tif',
       plot = last_plot(),
       device = 'tiff',
       width= 160,
       height= 120,
       units = 'mm',
       dpi = 400)

# Output:Figure 1F and G
# ============================================================