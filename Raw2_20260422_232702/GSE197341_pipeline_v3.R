#!/usr/bin/env Rscript
# GSE197341 Full Seurat Pipeline - v3 (fixed cluster naming)
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(Matrix)
})

cat("=== GSE197341 Pipeline v3 ===\n")
cat("Time:", format(Sys.time()), "\n\n")

data_dir <- "/workspace/data/GSE197341"
out_deg  <- "/mnt/results/05_scRNA_DEG"
out_cd   <- "/mnt/results/06_cell_death"
dir.create(out_deg, showWarnings=FALSE, recursive=TRUE)
dir.create(out_cd,  showWarnings=FALSE, recursive=TRUE)

# Load data
cat("Loading data...\n")
mtx   <- readMM(gzcon(file(file.path(data_dir, "GSE197341_expression_matrix.mtx.gz"), "rb")))
genes <- read.table(gzfile(file.path(data_dir, "GSE197341_genes.txt.gz")), header=FALSE)$V1
cells <- read.table(gzfile(file.path(data_dir, "GSE197341_CellIDS.txt.gz")), header=FALSE)$V1
meta  <- read.csv(gzfile(file.path(data_dir, "GSE197341_metadata.txt.gz")), row.names=1)
rownames(mtx) <- genes; colnames(mtx) <- cells
cat("Matrix:", nrow(mtx), "x", ncol(mtx), "\n")

# Create Seurat object
cat("Creating Seurat object...\n")
seu <- CreateSeuratObject(counts=mtx, min.cells=3, min.features=200)
cat("Cells:", ncol(seu), "Genes:", nrow(seu), "\n")

# Add metadata via @meta.data
samples_vec <- meta[colnames(seu), "Samples"]
tp_vec <- ifelse(grepl("^Sham", samples_vec), "Sham",
          ifelse(grepl("^D3 ",  samples_vec), "D3",
          ifelse(grepl("^D7 ",  samples_vec), "D7",
          ifelse(grepl("^D14 ", samples_vec), "D14",
          ifelse(grepl("^D28 ", samples_vec), "D28", "Unknown")))))
seu@meta.data$Samples   <- samples_vec
seu@meta.data$timepoint <- tp_vec
seu@meta.data$condition <- ifelse(tp_vec == "Sham", "Sham", "Stroke")
cat("Timepoints:\n"); print(table(tp_vec))

# QC
cat("\nQC...\n")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^mt-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
cat("After QC:", ncol(seu), "cells\n")

# Pipeline
cat("Normalizing...\n")
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
cat("Scaling...\n")
seu <- ScaleData(seu, features=VariableFeatures(seu), verbose=FALSE)
cat("PCA...\n")
seu <- RunPCA(seu, npcs=30, verbose=FALSE)
cat("UMAP...\n")
seu <- RunUMAP(seu, dims=1:20, verbose=FALSE)
cat("Clustering...\n")
seu <- FindNeighbors(seu, dims=1:20, verbose=FALSE)
seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
cat("Clusters:", length(unique(seu$seurat_clusters)), "\n")
print(table(seu$seurat_clusters))

# Cell type annotation
cat("\nAnnotating cell types...\n")
markers_list <- list(
  Neuron          = c("Snap25","Rbfox3","Syt1","Tubb3"),
  Astrocyte       = c("Gfap","Aldh1l1","Aqp4","Slc1a3"),
  Microglia       = c("Cx3cr1","P2ry12","Tmem119","Csf1r"),
  Oligodendrocyte = c("Mbp","Plp1","Mog","Cnp"),
  OPC             = c("Pdgfra","Cspg4","Sox10","Olig2"),
  Endothelial     = c("Cldn5","Pecam1","Cdh5","Tie1"),
  Pericyte        = c("Pdgfrb","Rgs5","Notch3","Kcnj8"),
  Fibroblast      = c("Col1a1","Col1a2","Dcn","Lum"),
  Macrophage      = c("Cd68","Mrc1","Adgre1","Ccr2")
)

avail_markers <- intersect(unique(unlist(markers_list)), rownames(seu))
cat("Available markers:", length(avail_markers), "\n")

Idents(seu) <- "seurat_clusters"
avg_expr_cl <- AverageExpression(seu, features=avail_markers, assays="RNA")$RNA

# Fix: strip "g" prefix from column names (Seurat v5 adds "g" to numeric cluster names)
colnames(avg_expr_cl) <- sub("^g", "", colnames(avg_expr_cl))
cat("Avg expr cluster names:", colnames(avg_expr_cl), "\n")

# Score clusters
cluster_scores <- matrix(0, nrow=ncol(avg_expr_cl), ncol=length(markers_list),
                         dimnames=list(colnames(avg_expr_cl), names(markers_list)))
for (ct in names(markers_list)) {
  ct_m <- intersect(markers_list[[ct]], rownames(avg_expr_cl))
  if (length(ct_m) > 0)
    cluster_scores[, ct] <- colMeans(avg_expr_cl[ct_m, , drop=FALSE])
}
cat("\nCluster scores:\n"); print(round(cluster_scores, 3))

# Assign cell types - use numeric cluster IDs matching seurat_clusters
cluster_celltype <- apply(cluster_scores, 1, function(x) names(which.max(x)))
cat("\nCluster -> cell type:\n"); print(cluster_celltype)

# Map: seurat_clusters is factor with levels "0","1",...
seu@meta.data$cell_type <- cluster_celltype[as.character(seu@meta.data$seurat_clusters)]
cat("\nCell type distribution:\n"); print(table(seu@meta.data$cell_type))
cat("\nCell type by timepoint:\n"); print(table(seu@meta.data$cell_type, seu@meta.data$timepoint))

# DEG analysis
cat("\n=== DEG Analysis ===\n")
Idents(seu) <- "cell_type"

run_degs <- function(seu_sub, cell_label) {
  cat(sprintf("\n%s: %d cells\n", cell_label, ncol(seu_sub)))
  print(table(seu_sub@meta.data$timepoint))
  degs_list <- list()
  for (tp in c("D3","D7","D14","D28")) {
    cells_tp   <- rownames(seu_sub@meta.data)[seu_sub@meta.data$timepoint == tp]
    cells_sham <- rownames(seu_sub@meta.data)[seu_sub@meta.data$timepoint == "Sham"]
    cat(sprintf("  %s: %d vs Sham: %d\n", tp, length(cells_tp), length(cells_sham)))
    if (length(cells_tp) >= 10 && length(cells_sham) >= 10) {
      tryCatch({
        degs <- FindMarkers(seu_sub, ident.1=cells_tp, ident.2=cells_sham,
                            min.pct=0.1, logfc.threshold=0, test.use="wilcox", verbose=FALSE)
        degs$gene <- rownames(degs); degs$timepoint <- tp; degs$cell_type <- cell_label
        degs_list[[tp]] <- degs
        cat(sprintf("  -> %d DEGs\n", nrow(degs)))
      }, error=function(e) cat(sprintf("  -> Error: %s\n", e$message)))
    } else cat("  -> Skipped\n")
  }
  if (length(degs_list) > 0) do.call(rbind, degs_list)
  else data.frame(note=paste("No DEGs for", cell_label))
}

ec_cells    <- subset(seu, idents="Endothelial")
ec_degs_all <- run_degs(ec_cells, "Endothelial")

pc_cells    <- subset(seu, idents="Pericyte")
pc_degs_all <- run_degs(pc_cells, "Pericyte")

# Cell death module scores
cat("\n=== Module Scores ===\n")
apop_genes <- necro_genes <- ferro_genes <- NULL
tryCatch({
  library(msigdbr)
  sets_h  <- msigdbr(species="Mus musculus", category="H")
  sets_c2 <- msigdbr(species="Mus musculus", category="C2")
  apop_genes  <- sets_h[sets_h$gs_name=="HALLMARK_APOPTOSIS", "gene_symbol"][[1]]
  necro_genes <- sets_c2[sets_c2$gs_name=="REACTOME_REGULATED_NECROSIS", "gene_symbol"][[1]]
  ferro_genes <- sets_c2[sets_c2$gs_name=="WP_FERROPTOSIS", "gene_symbol"][[1]]
  cat("MSigDB loaded: Apop", length(apop_genes), "Necro", length(necro_genes), "Ferro", length(ferro_genes), "\n")
}, error=function(e) cat("msigdbr failed:", e$message, "\n"))

if (is.null(apop_genes)) {
  apop_genes  <- c("Casp3","Casp7","Casp8","Casp9","Bax","Bcl2","Bcl2l1","Cycs","Apaf1",
                   "Bid","Bad","Pmaip1","Bbc3","Mcl1","Xiap","Diablo","Tnf","Fasl","Fas",
                   "Tradd","Fadd","Ripk1","Cflar","Birc2","Birc3","Tp53","Mdm2","Cdkn1a",
                   "Gadd45a","Aifm1","Endog","Dffb","Dffa","Lmnb1","Acin1")
  necro_genes <- c("Ripk1","Ripk3","Mlkl","Fadd","Tradd","Tnfrsf1a","Tnf","Zbp1",
                   "Dak","Pgam5","Drp1","Cyld","Ikbkg","Casp8","Cflar","Traf2","Traf5")
  ferro_genes <- c("Gpx4","Slc7a11","Acsl4","Lpcat3","Alox15","Alox5","Ptgs2","Tfrc",
                   "Slc40a1","Hamp","Ftl1","Fth1","Hmox1","Nqo1","Nfe2l2","Keap1",
                   "Sat1","Dpp4","Cd44","Gls2","Prom2","Fancd2","Cisd1","Cisd2","Vdac2","Vdac3")
}

apop_genes  <- intersect(apop_genes,  rownames(seu))
necro_genes <- intersect(necro_genes, rownames(seu))
ferro_genes <- intersect(ferro_genes, rownames(seu))
cat("In dataset: Apop", length(apop_genes), "Necro", length(necro_genes), "Ferro", length(ferro_genes), "\n")

seu <- AddModuleScore(seu, features=list(apop_genes),  name="Apoptosis_score",  verbose=FALSE)
seu <- AddModuleScore(seu, features=list(necro_genes), name="Necroptosis_score", verbose=FALSE)
seu <- AddModuleScore(seu, features=list(ferro_genes), name="Ferroptosis_score", verbose=FALSE)

meta_df <- seu@meta.data
score_summary <- meta_df %>%
  group_by(cell_type, condition) %>%
  summarise(n_cells=n(),
    mean_apoptosis=mean(Apoptosis_score1,na.rm=TRUE), sd_apoptosis=sd(Apoptosis_score1,na.rm=TRUE),
    mean_necroptosis=mean(Necroptosis_score1,na.rm=TRUE), sd_necroptosis=sd(Necroptosis_score1,na.rm=TRUE),
    mean_ferroptosis=mean(Ferroptosis_score1,na.rm=TRUE), sd_ferroptosis=sd(Ferroptosis_score1,na.rm=TRUE),
    .groups="drop")
score_summary$dataset <- "GSE197341"

score_by_tp <- meta_df %>%
  group_by(cell_type, timepoint) %>%
  summarise(n_cells=n(),
    mean_apoptosis=mean(Apoptosis_score1,na.rm=TRUE),
    mean_necroptosis=mean(Necroptosis_score1,na.rm=TRUE),
    mean_ferroptosis=mean(Ferroptosis_score1,na.rm=TRUE),
    .groups="drop")
score_by_tp$dataset <- "GSE197341"

cat("\nModule score summary:\n"); print(as.data.frame(score_summary))

seu@meta.data$celltype_condition <- paste(seu@meta.data$cell_type, seu@meta.data$condition, sep="_")
Idents(seu) <- "celltype_condition"
all_pw_genes <- unique(c(apop_genes, necro_genes, ferro_genes))
avg_expr <- AverageExpression(seu, features=intersect(all_pw_genes, rownames(seu)),
                              assays="RNA", verbose=FALSE)$RNA

# Save outputs
cat("\n=== Saving ===\n")
write.csv(ec_degs_all, file.path(out_deg, "GSE197341_EC_DEG_by_timepoint.csv"), row.names=FALSE)
write.csv(pc_degs_all, file.path(out_deg, "GSE197341_pericyte_DEG_by_timepoint.csv"), row.names=FALSE)
write.csv(score_summary, file.path(out_cd, "GSE197341_module_scores.csv"), row.names=FALSE)
write.csv(score_by_tp,   file.path(out_cd, "GSE197341_module_scores_by_timepoint.csv"), row.names=FALSE)
write.csv(as.data.frame(avg_expr), file.path(out_cd, "GSE197341_avg_expr_pathway_genes.csv"), row.names=TRUE)
cat("CSVs saved\n")

Idents(seu) <- "cell_type"
p1 <- DimPlot(seu, group.by="cell_type", label=TRUE, repel=TRUE) +
      ggtitle("GSE197341 - dMCAO cortex (Nakamura 2025)") + theme_classic()
ggsave(file.path(out_deg, "GSE197341_UMAP_celltype.svg"), p1, width=8, height=6)
ggsave(file.path(out_deg, "GSE197341_UMAP_celltype.png"), p1, width=8, height=6, dpi=150)
p2 <- DimPlot(seu, group.by="timepoint") + ggtitle("GSE197341 - by timepoint") + theme_classic()
ggsave(file.path(out_deg, "GSE197341_UMAP_timepoint.png"), p2, width=8, height=6, dpi=150)
cat("UMAPs saved\n")

saveRDS(seu, file.path(data_dir, "GSE197341_seurat.rds"))
cat("Seurat object saved\n")

cat("\n=== COMPLETE ===\n")
cat("Time:", format(Sys.time()), "\n")
cat("Total cells:", ncol(seu), "\n")
cat("Cell types:\n"); print(sort(table(seu@meta.data$cell_type), decreasing=TRUE))
cat("Timepoints:\n"); print(table(seu@meta.data$timepoint))
if ("timepoint" %in% colnames(ec_degs_all)) {
  cat("EC DEGs per timepoint:\n"); print(table(ec_degs_all$timepoint))
}
if ("timepoint" %in% colnames(pc_degs_all)) {
  cat("Pericyte DEGs per timepoint:\n"); print(table(pc_degs_all$timepoint))
}
