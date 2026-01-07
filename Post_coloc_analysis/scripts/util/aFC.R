rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

options(bitmapType = "cairo")

# ---- I/O (edit paths if needed) ----
setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/")

# Load Interval gene list (protein-coding only)
gene_list <- read_excel("Post_coloc_analysis/data/interval_genes.xlsx", sheet = 1) |>
  dplyr::filter(gene_biotype == "protein_coding")

# Load GTEx eGenes (Whole Blood) and intersect with Interval genes
afc <- data.table::fread("data/GTEx_Analysis_v10_eQTL_updated/Whole_Blood.v10.eGenes.txt.gz")
afc <- afc[qval <= 0.05, ]
afc <- afc[gene_name %in% gene_list$gene_name | gene_name %in% gene_list$feature_id]

# Drop NAs/Infs for a clean correlation
df <- subset(afc, is.finite(afc) & is.finite(slope))

# Correlation (change method="spearman" if preferred)
ct  <- cor.test(df$afc, df$slope, method = "pearson")
lab <- sprintf("Pearson r = %.2f", ct$estimate)

png("Post_coloc_analysis/plots/interval/afc_beta_cor.png", res = 300, units = "in", width = 6, height = 6)

ggplot(df, aes(x = afc, y = slope)) +
  geom_point(alpha = 0.6, size = 1.8) +
  labs(x = "Allelic fold change (aFC)", y = "eQTL effect size") +
  annotate(
    "label", x = Inf, y = Inf, label = lab,
    hjust = 1.05, vjust = 1.3, size = 3.6, label.size = 0
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "grey30")
  )
dev.off()


# ---- Compute CIs and selection flags (idempotent if already present) ----
if (!("lci" %in% colnames(afc)) || !("uci" %in% colnames(afc))) {
  afc[, `:=`(lci = afc - 1.96 * afc_se, uci = afc + 1.96 * afc_se)]
}
if (!("crosses_zero" %in% colnames(afc))) {
  afc$crosses_zero <- afc$lci <= 0 & afc$uci >= 0
}

beta_dec <- log2(0.75)  # −25% line on log2 scale
if (!("ci_class" %in% colnames(afc))) {
  afc[, ci_class := dplyr::case_when(
    afc < 0 & beta_dec >= lci & beta_dec <= uci ~ TRUE,
    TRUE ~ FALSE
  )]
}

# ---- Final set used for simulation (your criteria) ----
sel <- afc |>
  dplyr::filter(ci_class == TRUE, crosses_zero == FALSE, afc < 0)

stopifnot(nrow(sel) > 0)  # ensure we have data

# ---- Fixed target % changes: 10%, 15%, 20%, 25%, 30%, 50% (down-regulation) ----
pct_targets <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.50)
# Convert to log2(AFC) targets in the negative direction: log2(1 - p)
grid_vals <- log2(1 - pct_targets)   # e.g., 10% -> log2(0.90) ≈ -0.152

# ---- CI-based counts + mean slope per target ----
bin_summ <- lapply(grid_vals, function(gv) {
  idx <- sel$lci <= gv & sel$uci >= gv
  n_ci <- sum(idx, na.rm = TRUE)
  mu_slope <- if (n_ci > 0) mean(sel$slope[idx], na.rm = TRUE) else NA_real_
  data.frame(target_log2 = gv, n_ci_cover = n_ci, mean_slope = mu_slope)
}) |> dplyr::bind_rows()

to_pct <- function(x_log2) (2^x_log2 - 1) * 100
bin_summ <- bin_summ |>
  dplyr::mutate(target_pct = to_pct(target_log2)) |>
  dplyr::arrange(target_log2)

# ---- Panel D table (now with 4th column = mean eQTL slope) ----
tab <- bin_summ |>
  dplyr::transmute(
    log2AFC            = sprintf("%.3f", target_log2),
    pct_change         = sprintf("%+.1f%%", target_pct),
    N_CI_covers_target = n_ci_cover,
    mean_eQTL_slope    = sprintf("%.3f", mean_slope)
  ) |>
  dplyr::mutate(
    mean_eQTL_slope = ifelse(is.na(mean_eQTL_slope), "—", mean_eQTL_slope),
    row = dplyr::row_number()
  )

tab_long <- tab |>
  dplyr::mutate(dplyr::across(c(log2AFC, pct_change, N_CI_covers_target, mean_eQTL_slope), as.character)) |>
  tidyr::pivot_longer(
    cols = c(log2AFC, pct_change, N_CI_covers_target, mean_eQTL_slope),
    names_to = "col", values_to = "val"
  ) |>
  dplyr::mutate(
    col = factor(
      col,
      levels = c("log2AFC", "pct_change", "N_CI_covers_target", "mean_eQTL_slope"),
      labels = c("Target log2(AFC)", "Target % change", "N with CI covering target", "Mean eQTL slope (β)")
    ),
    y = max(row) - row + 1
  )

headers <- data.frame(
  col = levels(tab_long$col),
  y   = max(tab_long$y) + 1
)

pD <- ggplot() +
  geom_label(data = headers, aes(x = col, y = y, label = col),
             label.size = 0, alpha = 0.9) +
  geom_label(data = tab_long, aes(x = col, y = y, label = val),
             label.size = 0, alpha = 0.7) +
  xlab(NULL) + ylab(NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text         = element_blank(),
    panel.grid        = element_blank(),
    plot.title        = element_blank(),
    plot.subtitle     = element_blank(),
    plot.title.position = "plot"
  )

# ---- Assemble panels A | D and save ----
layout_AD <- pA | pD
if (!requireNamespace("svglite", quietly = TRUE)) {
  stop("Please install the 'svglite' package to export SVG: install.packages('svglite')")
}
ggsave(
  "Post_coloc_analysis/plots/interval/effect_size_selection_AD.svg",
  pD,
  width = 6, height = 6,
  device = svglite::svglite
)

# (Optional) PNG export for quick preview
# ggsave("SuppFig_effect_size_selection_AD.png", layout_AD, width = 11, height = 5.5, dpi = 300)