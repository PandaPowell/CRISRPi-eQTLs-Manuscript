rm(list=ls())
library(data.table)
library(tidyverse)
options(bitmapType="cairo")

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

gold_genes <- fread("processed_data/cis_gold_genes.txt")

gold_genes_long <- gold_genes[
  , .(grna = trimws(unlist(strsplit(grna_target, "\\s*,\\s*")))),
  by = gold_gene
]

# (optional) drop exact duplicates
gold_genes_long <- unique(gold_genes_long)
gold_genes_long$grna_gene = paste(gold_genes_long$grna, gold_genes_long$gold_gene, sep = "_")
head(gold_genes_long)

# Run bedtools intersect
system("bash scripts/util/gold_gene_atac_intersection.sh")

## ---- libraries ----
# ---- packages ----
library(pheatmap)

ok_cu <- requireNamespace("ComplexUpset", quietly = TRUE)
ok_gg <- requireNamespace("ggplot2", quietly = TRUE)

# ---- inputs ----
egene_file   <- "processed_data/cres_with_grna_eqtls_interval.txt"
cgene_file   <- "processed_data/cres_with_grnas.txt"
overlap_bed  <- "ATAC_overlap/gold_genes_ATAC_peaks.bed"  # last column must be cell name

# ---- gene sets ----
egene_ids <- fread(egene_file) %>% filter(target_gene %in% gold_genes_long$grna_gene)
egene_ids <- unique(egene_ids$grna_target)
cgene_ids = fread(cgene_file) %>% filter(target_gene %in% gold_genes_long$grna_gene, significant == 1)
cgene_ids <-  unique(cgene_ids$grna_target)
gene_universe <- union(egene_ids, cgene_ids)

# ---- ATAC overlaps ----
atac <- fread(overlap_bed) %>%
  filter(V9 %in% c("B","CD4","CD8","Ery","mDC","Mega","Mono", "NK", "pDC") )# Keep most relevant cell tpyes 

# Identify columns (you’ve used Ensembl in V4; cell name was appended last)
gene_col <- "V4"
cell_col <- names(atac)[ncol(atac)]

# Keep only genes in the universe; keep ALL rows (even V9==".") to build full grid
atac <- atac[get(gene_col) %in% gene_universe]

# Presence per (gene, cell): 1 if ANY real overlap (V9 != "."), else 0
tab_bin <- atac[, .(present = as.integer(any(V8 != "."))),
                by = .(gene = get(gene_col), cell = get(cell_col))]

# ---- complete grid for ALL (gene_universe × cells) ----
all_cells <- sort(unique(atac[[cell_col]]))
grid <- CJ(gene = gene_universe, cell = all_cells, unique = TRUE)

tab_bin_full <- grid[tab_bin, on = .(gene, cell)]
tab_bin_full[is.na(present), present := 0L]

# ---- wide binary matrix ----
mat_bin_df <- dcast(tab_bin_full, gene ~ cell, value.var = "present", fill = 0)
mat_bin <- as.matrix(mat_bin_df[, -1, drop = FALSE])
rownames(mat_bin) <- mat_bin_df[[1]]
mode(mat_bin) <- "numeric"

# ---- gene classes (eGene / cGene / both) ----
g <- rownames(mat_bin)
gene_classes <- data.frame(
  gene  = g,
  class = ifelse(g %in% egene_ids & g %in% cgene_ids, "both",
                 ifelse(g %in% egene_ids, "eGene", "cGene")),
  stringsAsFactors = FALSE
)
gene_classes$class <- factor(gene_classes$class, levels = c("eGene","cGene","both"))

# ---- order rows by class then specificity (breadth ascending) ----
breadth  <- rowSums(mat_bin)
ord_rows <- order(gene_classes$class, breadth)   # fewest cells first
mat_ord  <- mat_bin[ord_rows, , drop = FALSE]
ann_row  <- gene_classes[ord_rows, , drop = FALSE]
ann_row$breadth <- breadth[ord_rows]

gap_rows <- cumsum(table(ann_row$class))
gap_rows <- gap_rows[gap_rows < nrow(mat_ord)]

# ---- heatmap: binary (0=white, 1=dark) ----
cols   <- c("white", "#08306b")
breaks <- c(-0.5, 0.5, 1.5)
ann_colors <- list(
  class   = c(eGene="#1b9e77", cGene="#d95f02", both="#7570b3"),
  breadth = colorRampPalette(c("#f7fbff","#08306b"))(100)
)

pheatmap(
  mat_ord,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  border_color = NA,
  color  = cols,
  breaks = breaks,
  annotation_row = data.frame(class=ann_row$class, breadth=ann_row$breadth, row.names=ann_row$gene),
  annotation_colors = ann_colors,
  gaps_row = as.integer(gap_rows)
)

# ---- UpSet: intersections; highlight genes unique to each cell (if ComplexUpset) ----
cells_ord <- names(sort(colSums(mat_bin), decreasing = TRUE))
df_cu <- as.data.frame(mat_bin)
df_cu <- df_cu[, cells_ord, drop = FALSE]
df_cu$gene <- rownames(mat_bin)

library(UpSetR)

# 0) Start from your binary matrix
df_upset <- as.data.frame(mat_bin)                 # genes x cells, 0/1
stopifnot(all(as.matrix(df_upset) %in% c(0,1)))
# Drop genes with no cell present (don’t affect intersections, but do slow things)
#df_upset <- df_upset[rowSums(df_upset) > 0, , drop = FALSE]

# 1) Order sets by size (desc)
set_sizes <- colSums(df_upset)
sets_ord  <- names(sort(set_sizes, decreasing = TRUE))
df_upset  <- df_upset[ , sets_ord, drop = FALSE]

egene_upset = df_upset[rownames(df_upset) %in% egene_ids,]
cat(sum(rowSums(egene_upset) == 0)/nrow(egene_upset)*100,"% of genes overlap no open Chromatin")
cat(sum(rowSums(egene_upset) == 1)/nrow(egene_upset)*100,"% of genes overlap one cell type")

cgene_upset = df_upset[rownames(df_upset) %in% cgene_ids,]
cat(sum(rowSums(cgene_upset) == 0)/nrow(cgene_upset)*100,"% of genes overlap no open Chromatin")
cat(sum(rowSums(cgene_upset) == 1)/nrow(cgene_upset)*100,"% of genes overlap one cell type")

png(filename = "plots/interval/atac_gold_upset_egene_cres.png", width = 10, height = 8, units = "in", res = 300)
# plot
UpSetR::upset(
  egene_upset,
  sets                = colnames(egene_upset),
  nsets               = ncol(egene_upset),
  nintersects         = 20,          # reduce if heavy
  order.by            = "freq",
  keep.order          = TRUE,
  mb.ratio            = c(0.6, 0.4),
  empty.intersections = "off",
  sets.x.label        = "Genes per cell type",
  mainbar.y.label     = "Intersection size",
  text.scale = c(3, 3, 2, 2, 3.3, 3),
  sets.bar.color      = "#4c78a8",
  matrix.color        = "#2a3f5f"
)
dev.off()

UpSetR::upset(
  cgene_upset,
  sets                = colnames(cgene_upset),
  nsets               = ncol(cgene_upset),
  nintersects         = 14,
  order.by            = "freq",
  keep.order          = TRUE,
  mb.ratio            = c(0.6, 0.4),
  empty.intersections = "off",
  sets.x.label        = "Genes per cell type",
  mainbar.y.label     = "Intersection size",
  sets.bar.color      = "#4c78a8",
  matrix.color        = "#2a3f5f"
)
