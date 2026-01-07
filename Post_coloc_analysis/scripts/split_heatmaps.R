# One/two-column per gene heatmap:
# - Only CRISPRi or only eQTL -> single column labeled "GENE"
# - Both -> two columns: "GENE (CRISPRi)" and "GENE (eQTL)"
# Colors: method-specific; 0 = light grey; hard jump at 1; cap at 5
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(svglite)
library(scales)
library(stringr)

# ---- Inputs expected -----------------------------------------------
# df_cri: columns row_lab (gene), gwas_label (trait), count (numeric)
# df_eqtl: columns row_lab (gene), gwas_label (trait), count (numeric)
#---------------------------------------
# Helper: build COUNT heatmap df for one method
#---------------------------------------
make_heat_df_counts <- function(df_labeled, gold_genes) {
  # gold gene's (target, label) combos
  base <- df_labeled %>%
    filter(ensembl_id %in% gold_genes$gold_gene) %>%
    mutate(row_lab = gene_name) %>%
    distinct(grna_target, row_lab, gwas_label)
  
  # for each (target,label), how many distinct genes are hit
  per_target_counts <- df_labeled %>%
    distinct(grna_target, gene_name, gwas_label) %>%
    group_by(grna_target, gwas_label) %>%
    summarise(n_target = n_distinct(gene_name), .groups = "drop")
  
  # sum all genes at the gold gene's targets (including the gold gene)
  counts <- base %>%
    left_join(per_target_counts, by = c("grna_target", "gwas_label")) %>%
    mutate(n_target = replace_na(n_target, 0L)) %>%
    group_by(row_lab, gwas_label) %>%
    summarise(count = sum(n_target), .groups = "drop")
  
  # complete grid and totals for ordering
  rows <- counts %>% distinct(row_lab)
  cols <- counts %>% distinct(gwas_label)
  counts %>%
    tidyr::complete(row_lab = rows$row_lab, gwas_label = cols$gwas_label) %>%
    group_by(row_lab) %>% mutate(row_tot = sum(replace_na(count, 0))) %>% ungroup() %>%
    group_by(gwas_label) %>% mutate(col_tot = sum(replace_na(count, 0))) %>% ungroup()
}

# Build data for both methods
df_cri  <- make_heat_df_counts(sig_ccres_labeled, gold_genes)
df_eqtl <- make_heat_df_counts(sig_ecres_labeled, gold_genes)

#---------------------------------------
# Shared ordering using the UNION of labels (keep CRISPRi order first)
#---------------------------------------
row_levels_cri  <- df_cri  %>% arrange(desc(row_tot)) %>% pull(row_lab)    %>% unique()
row_levels_eqtl <- df_eqtl %>% arrange(desc(row_tot)) %>% pull(row_lab)    %>% unique()
col_levels_cri  <- df_cri  %>% arrange(desc(col_tot)) %>% pull(gwas_label) %>% unique()
col_levels_eqtl <- df_eqtl %>% arrange(desc(col_tot)) %>% pull(gwas_label) %>% unique()

row_levels_all <- c(row_levels_cri, setdiff(row_levels_eqtl, row_levels_cri))
col_levels_all <- c(col_levels_cri, setdiff(col_levels_eqtl, col_levels_cri))

df_cri <- df_cri %>%
  mutate(row_lab = factor(row_lab, levels = row_levels_all),
         gwas_label = factor(gwas_label, levels = col_levels_all))

df_eqtl <- df_eqtl %>%
  mutate(row_lab = factor(row_lab, levels = row_levels_all),
         gwas_label = factor(gwas_label, levels = col_levels_all))

stopifnot(exists("df_cri"), exists("df_eqtl"))
stopifnot(all(c("row_lab","gwas_label","count") %in% names(df_cri)))
stopifnot(all(c("row_lab","gwas_label","count") %in% names(df_eqtl)))

# Tag methods and bind
df_cri$method  <- "CRISPRi"
df_eqtl$method <- "eQTL"
df_all <- bind_rows(df_cri, df_eqtl)

# ---- Drop genes with no nonzero counts in either method ------------
genes_keep <- df_all %>%
  group_by(row_lab) %>%
  summarise(any_present = any(replace_na(count, 0) > 0), .groups = "drop") %>%
  filter(any_present) %>%
  pull(row_lab)

df_all <- df_all %>% filter(row_lab %in% genes_keep)

# Which methods actually have signal (>0) per gene?
present_methods <- df_all %>%
  group_by(row_lab, method) %>%
  summarise(any_present = any(replace_na(count, 0) > 0), .groups = "drop") %>%
  filter(any_present)

# ---- Build column layout (1 or 2 per gene) -------------------------
gene_order <- sort(unique(present_methods$row_lab))
both_genes <- present_methods %>% count(row_lab, name = "n_methods") %>%
  filter(n_methods == 2) %>% pull(row_lab)

col_layout <- present_methods %>%
  mutate(col_id = if_else(row_lab %in% both_genes,
                          paste0(row_lab, " (", method, ")"),
                          row_lab)) %>%
  arrange(match(row_lab, gene_order),
          factor(method, levels = c("CRISPRi","eQTL"))) %>%
  mutate(col_index = row_number()) %>%
  dplyr::select(row_lab, method, col_id, col_index)

# Build plotting frame with method preserved
df_plot_long <- df_all %>%
  inner_join(col_layout, by = c("row_lab","method")) %>%   # adds col_id + keeps method
  dplyr::select(col_id, method, gwas_label, count)

# Complete grid (and reattach method for each col_id)
all_traits <- sort(unique(df_plot_long$gwas_label))
df_plot_long <- df_plot_long %>%
  tidyr::complete(col_id, gwas_label = all_traits, fill = list(count = 0)) %>%
  left_join(dplyr::distinct(col_layout, col_id, method), by = "col_id") %>%  # <- restore `method`
  mutate(
    col_id     = factor(col_id, levels = col_layout$col_id),
    gwas_label = factor(gwas_label, levels = all_traits),
    count_cap  = pmin(count, 5)
  )

# ---- Colors: 0 light grey; hard jump at 1; cap at 5 ----------------
zero_grey <- "#eeeeee"  # lighter grey for 0
na_col    <- "#f6f6f6"

# ramps ONLY used above 1
ramp_y_to_o   <- scales::colour_ramp(c("#ffd73e", "#e29421"))   # CRISPRi 1..3
ramp_o_to_r   <- scales::colour_ramp(c("#e29421", "#ce472e"))   # CRISPRi 3..5
ramp_b1_to_b2 <- scales::colour_ramp(c("#0099dc", "#4f46e5"))  # eQTL    1..5

map_colors_cri <- function(v) {
  v <- pmin(pmax(v, 0), 5)
  out <- rep(NA_character_, length(v))
  out[is.na(v)] <- na_col
  out[v == 0]   <- zero_grey
  out[v > 0 & v <= 1] <- "#ffd73e"        # hard jump at 1
  idx <- v > 1 & v <= 3; out[idx] <- ramp_y_to_o((v[idx] - 1) / 2)
  idx <- v > 3 & v <= 5; out[idx] <- ramp_o_to_r((v[idx] - 3) / 2)
  out
}
map_colors_eqtl <- function(v) {
  v <- pmin(pmax(v, 0), 5)
  out <- rep(NA_character_, length(v))
  out[is.na(v)] <- na_col
  out[v == 0]   <- zero_grey
  out[v > 0 & v <= 1] <- "#0099dc"        # hard jump at 1
  idx <- v > 1 & v <= 5; out[idx] <- ramp_b1_to_b2((v[idx] - 1) / 4)
  out
}

df_plot_long <- df_plot_long %>%
  mutate(fill_col = ifelse(method.x == "CRISPRi",
                           map_colors_cri(count_cap),
                           map_colors_eqtl(count_cap)))

# ---- Plot -----------------------------------------------------------
p_one_or_two <- ggplot(df_plot_long,
                       aes(x = col_id, y = gwas_label, fill = fill_col)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.98, height = 0.98) +
  scale_fill_identity() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Gene (method if both)", y = "GWAS label") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y     = element_text(size = 12),
    panel.grid      = element_blank(),
    plot.margin     = grid::unit(c(0,0,0,0), "pt")
  )

# ---- Legends (reflect the non-linear mapping) ----------------------
legend_vals <- data.frame(val = seq(0, 5, by = 0.02))
legend_vals$col_cri  <- map_colors_cri(legend_vals$val)
legend_vals$col_eqtl <- map_colors_eqtl(legend_vals$val)

leg_cri <- ggplot(legend_vals, aes(x = val, y = 1)) +
  geom_tile(aes(fill = col_cri), height = 1) +
  scale_fill_identity() +
  labs(title = "CRISPRi", x = NULL, y = NULL) +
  scale_x_continuous(breaks = c(0,1,3,5), labels = c("0","1","3","5+")) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = margin(0,5,0,5)
  )

leg_eqtl <- ggplot(legend_vals, aes(x = val, y = 1)) +
  geom_tile(aes(fill = col_eqtl), height = 1) +
  scale_fill_identity() +
  labs(title = "eQTL", x = NULL, y = NULL) +
  scale_x_continuous(breaks = c(0,1,3,5), labels = c("0","1","3","5+")) +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = margin(0,5,0,5)
  )

# ---- Assemble & save -----------------------------------------------
p_final <- p_one_or_two / (leg_cri | leg_eqtl) + plot_layout(heights = c(1, 0.18))

out_svg <- "plots/interval/gold_gene_gwas_trait_heatmap_gene_cols.svg"
ggsave(out_svg, p_final, device = svglite, width = 24, height = 9, bg = "white")

# Optional: crop SVG page to drawing (removes outer whitespace in editors)
# system2("inkscape", c("--export-area-drawing", "--export-type=svg",
#                       paste0("--export-filename=", sub(".svg$", "_cropped.svg", out_svg)),
#                       out_svg))
