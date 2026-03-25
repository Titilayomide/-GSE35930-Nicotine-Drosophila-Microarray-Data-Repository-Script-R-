# -GSE35930-Nicotine-Drosophila-Microarray-Data-Repository-Script-R-
Comprehensive pipeline for DEG identification and curation. Drosophila Nicotine Transcriptomics Analysis (GSE35930) Titilayomi A. Otenaike
# 1. LOAD REQUIRED LIBRARIES
library(affy)
library(limma)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)

# 2. DATA IMPORT AND NORMALIZATION
# Ensure .CEL files are in the working directory
setwd("H:/My Drive/PHD THESIS/RNA SEQ NICOTINE/Transcriptomics Nic-Paraq Data/Latest Analysis- 29TH Jan 2026")
raw_data <- ReadAffy()
eset <- rma(raw_data)
groups <- factor(c(rep("Control", 6), rep("Nicotine", 5)))

# 3. LINEAR MODELING (DEGs)
design <- model.matrix(~0 + groups)
colnames(design) <- c("Control", "Nicotine")
fit <- lmFit(eset, design)
cont_matrix <- makeContrasts(Nic_vs_Con = Nicotine - Control, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef = "Nic_vs_Con", adjust = "fdr", number = Inf)

# 4. FLYBASE MANUAL CURATION (The Label Fix)
# Mapping Probe IDs to Gene Symbols to resolve systematic identifier gaps
master_map <- c(
  "1631024_at" = "parkin", "1633514_at" = "Hsp70", 
  "1632349_at" = "GstE1",  "1624808_at" = "GstE10",
  "1624534_at" = "Cyp4p2", "1631180_at" = "Obp56h"
)
res$Symbol <- ifelse(rownames(res) %in% names(master_map), 
                     master_map[rownames(res)], rownames(res))

# 5. EXPORT GENE HUBS
res$Regulation <- ifelse(res$logFC > 1 & res$adj.P.Val < 0.05, "Upregulated",
                  ifelse(res$logFC < -1 & res$adj.P.Val < 0.05, "Downregulated", "NS"))

write.csv(res[res$Regulation == "Upregulated", ], "Hub1_Upregulated.csv")
write.csv(res[res$Regulation == "Downregulated", ], "Hub2_Downregulated.csv")

# 6. VISUALIZATION: LABELED HEATMAP
top30 <- head(res[order(res$adj.P.Val), ], 30)
plot_mat <- exprs(eset)[rownames(top30), ]
rownames(plot_mat) <- top30$Symbol # Injecting curated names

pheatmap(plot_mat, scale = "row", show_rownames = TRUE, 
         annotation_col = data.frame(Group = groups, row.names = colnames(plot_mat)),
         main = "Fig 1: Nicotine-Induced Transcriptomic shifts",
         filename = "Figure_1_Heatmap_FINAL.png")

# 7. SAVE ENVIRONMENT FOR PERSISTENCE
save.image("Nicotine_Final_Analysis_Complete.RData")

GitHub Repository: README.md
Drosophila Nicotine Transcriptomics Analysis (GSE35930)
Project Overview
This repository contains the reproducible bioinformatic pipeline for analyzing the neurotoxicological effects of chronic nicotine exposure on the Drosophila melanogaster head transcriptome. The study specifically isolates the Nicotine vs. Control contrast from the GSE35930 dataset.

Key Findings
Mitochondrial Impairment: Significant downregulation of parkin (CG33144), suggesting compromised mitophagy.

Proteotoxic Stress: Robust induction of Hsp70, indicating protein misfolding.

Metabolic Defense: Upregulation of the Glutathione S-transferase (Gst) family, specifically GstE1, GstE9, and GstE10.

Prerequisites
To run the analysis script, you must have R (v4.5.1) or later installed along with the following Bioconductor packages:

affy

limma

clusterProfiler

org.Dm.eg.db

Data Access
Raw .CEL files must be downloaded from the GEO Accession GSE35930 and placed in the project working directory.

Citation
If you use this code or the curated gene mappings, please cite:

Dataset: Hill-Burns, E. M., et al. (2013). "Genome-wide expression profiling of Drosophila nicotine and paraquat treatment." GEO Accession GSE35930.

Analysis: Otenaike, T. A. (2026). Transcriptomic Signatures of Nicotine Neurotoxicity in Drosophila melanogaster. PhD Thesis, Federal University of Rio Grande do Sul (UFRGS).

Formal Manuscript Methodology
Bioinformatic Pipeline & Statistical Modeling
The analysis was performed using a custom pipeline in R (v4.5.1). Raw Affymetrix intensity data were normalized using the Robust Multi-array Average (RMA) algorithm. Differential expression was modeled via linear regression using the limma package, applying Empirical Bayes moderation to ensure stable variance estimates across the 18,800 transcripts.

Identifier Resolution & FlyBase Curation
A critical step in this methodology involved the manual resolution of systematic identifiers from the GPL1322 platform. Automated mapping often failed to provide symbols for high-variance probes; therefore, FlyBase (Release 2026_01) was utilized to manually curate and verify the identity of key hubs, specifically the Gst family and mitochondrial quality control markers.

Pathway & Functional Enrichment
Gene Ontology (GO) terms were enriched using the clusterProfiler package, focusing on biological processes (BP). Significance for enrichment was set at an Adjusted P-value < 0.05 using the Benjamini-Hochberg procedure.
