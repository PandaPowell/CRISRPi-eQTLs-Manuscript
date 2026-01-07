# Single-column diagonal-split heatmap:
# - One column per gene (x)
# - One row per GWAS label (y)
# - Each cell split on the diagonal:
#     • lower-left triangle = CRISPRi
#     • upper-right triangle = eQTL
# Colors: method-specific; 0 = light grey; hard jump at 1; cap at 5

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(svglite)
library(scales)
library(stringr)

# --------- Inputs expected -------------------------------------------
# df_cri: columns row_lab (gene), gwas_label (trait), count (numeric)
# df_eqtl: columns row_lab (gene), gwas_label (trait), count (numeric)
stopifnot(exists("df_cri"), exists("df_eqtl"))
stopifnot(all(c("row_lab","gwas_label","count") %in% names(df_cri)))
stopifnot(all(c("row_lab","gwas_label","count") %in% names(df_eqtl)))

# Harmonize ordering across methods (optional: keep your previous ordering logic)
rows_all <- sort(unique(c(df_cri$row_lab, df_eqtl$row_lab)))
cols_all <- sort(unique(c(as.character(df_cri$gwas_label), as.character(df_eqtl$gwas_label))))

df_all <- bind_rows(
  df_cri  %>% mutate(method = "CRISPRi"),
  df_eqtl %>% mutate(method = "eQTL")
)

# Keep genes that have any nonzero signal in either method
genes_keep <- df_all %>%
  group_by(row_lab) %>%
  summarise(any_present = any(replace_na(count, 0) > 0), .groups = "drop") %>%
  filter(any_present) %>% pull(row_lab)

df_all <- df_all %>% filter(row_lab %in% genes_keep)

# Wide matrix: one row per (gene, trait), two columns for methods
df_wide <- df_all %>%
  dplyr::select(row_lab, gwas_label, method, count) %>%
  mutate(gwas_label = as.character(gwas_label)) %>%
  group_by(row_lab, gwas_label, method) %>%
  summarise(count = sum(count), .groups = "drop") %>%  # in case of duplicates
  pivot_wider(names_from = method, values_from = count, values_fill = 0) %>%
  complete(row_lab = genes_keep, gwas_label = cols_all, fill = list(CRISPRi = 0, eQTL = 0)) %>%
  mutate(
    gene = factor(row_lab, levels = genes_keep),
    trait = factor(gwas_label, levels = cols_all),
    x = as.integer(gene),
    y = as.integer(trait),
    cri_cap  = pmin(CRISPRi, 5),
    eqtl_cap = pmin(eQTL,   5)
  )

# --------- Color mapping ---------------------------------------------
zero_grey <- "#eeeeee"  # 0
na_col    <- "#f6f6f6"  # NA (unused here, but kept for completeness)

ramp_y_to_o   <- scales::colour_ramp(c("#ffd73e", "#e29421"))  # CRISPRi 1..3
ramp_o_to_r   <- scales::colour_ramp(c("#e29421", "#ce472e"))  # CRISPRi 3..5
ramp_b1_to_b2 <- scales::colour_ramp(c("#0099dc", "#4f46e5"))  # eQTL    1..5

map_colors_cri <- function(v) {
  v <- pmin(pmax(v, 0), 5)
  out <- rep(NA_character_, length(v))
  out[v == 0]   <- zero_grey
  out[v > 0 & v <= 1] <- "#ffd73e"              # hard jump at 1 (no fade)
  idx <- v > 1 & v <= 3; out[idx] <- ramp_y_to_o((v[idx] - 1)/2)
  idx <- v > 3 & v <= 5; out[idx] <- ramp_o_to_r((v[idx] - 3)/2)
  out
}
map_colors_eqtl <- function(v) {
  v <- pmin(pmax(v, 0), 5)
  out <- rep(NA_character_, length(v))
  out[v == 0]   <- zero_grey
  out[v > 0 & v <= 1] <- "#0099dc"              # hard jump at 1
  idx <- v > 1 & v <= 5; out[idx] <- ramp_b1_to_b2((v[idx] - 1)/4)
  out
}

df_wide <- df_wide %>%
  mutate(
    fill_cri  = map_colors_cri(cri_cap),
    fill_eqtl = map_colors_eqtl(eqtl_cap)
  )

# --------- Build triangle polygons per cell --------------------------
# For each (x,y), rectangle bounds:
df_rect <- df_wide %>%
  mutate(
    xmin = x - 0.5, xmax = x + 0.5,
    ymin = y - 0.5, ymax = y + 0.5
  )

# Lower-left triangle (CRISPRi): (xmin,ymin) -> (xmin,ymax) -> (xmax,ymin)
tri_cri <- df_rect %>%
  transmute(
    gene, trait, x, y, fill = fill_cri,
    xv1 = xmin, yv1 = ymin,
    xv2 = xmin, yv2 = ymax,
    xv3 = xmax, yv3 = ymin
  ) %>%
  tidyr::pivot_longer(
    cols = matches("^(xv|yv)\\d$"),
    names_to = c(".value", "idx"),
    names_pattern = "(xv|yv)(\\d)"
  )
# result: columns gene, trait, x, y, fill, idx (1..3), xv, yv

# Upper-right triangle (eQTL): (xmax,ymax) -> (xmin,ymax) -> (xmax,ymin)
tri_eqtl <- df_rect %>%
  transmute(
    gene, trait, x, y, fill = fill_eqtl,
    xv1 = xmax, yv1 = ymax,
    xv2 = xmin, yv2 = ymax,
    xv3 = xmax, yv3 = ymin
  ) %>%
  tidyr::pivot_longer(
    cols = matches("^(xv|yv)\\d$"),
    names_to = c(".value", "idx"),
    names_pattern = "(xv|yv)(\\d)"
  )

# (Small helper to avoid factor level headaches in the plot)
tri_cri  <- tri_cri  %>% mutate(gene = factor(gene, levels = genes_keep),
                                trait = factor(trait, levels = cols_all))
tri_eqtl <- tri_eqtl %>% mutate(gene = factor(gene, levels = genes_keep),
                                trait = factor(trait, levels = cols_all))

# --------- Plot: diagonal split tiles --------------------------------
p_heat <- ggplot() +
  # CRISPRi triangles (lower-left)
  geom_polygon(
    data = tri_cri,
    aes(x = xv, y = yv, group = interaction(gene, trait, fill), fill = fill),
    color = NA
  ) +
  # eQTL triangles (upper-right)
  geom_polygon(
    data = tri_eqtl,
    aes(x = xv, y = yv, group = interaction(gene, trait, fill), fill = fill),
    color = NA
  ) +
  # White borders for the outer squares
  geom_rect(
    data = df_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, color = "white", linewidth = 0.5
  ) +
  # Diagonal line to separate halves
  geom_segment(
    data = df_rect,
    aes(x = xmin, y = ymax, xend = xmax, yend = ymin),
    color = "white", linewidth = 0.5
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = seq_along(genes_keep),
    labels = genes_keep,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    breaks = seq_along(cols_all),
    labels = cols_all,
    expand = c(0,0)
  ) +
  labs(x = "Gold-gene CRE", y = "GWAS trait") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y     = element_text(size = 12),
    panel.grid      = element_blank(),
    plot.margin     = grid::unit(c(0,0,0,0), "pt")
  )

# --------- Legends (show the non-linear scales) ----------------------
legend_vals <- data.frame(val = seq(0, 5, by = 0.02))
legend_vals$col_cri  <- map_colors_cri(legend_vals$val)
legend_vals$col_eqtl <- map_colors_eqtl(legend_vals$val)

leg_cri <- ggplot(legend_vals, aes(x = val, y = 1)) +
  geom_tile(aes(fill = col_cri), height = 1) +
  scale_fill_identity() +
  labs(title = "Number of CRISPRi target genes", x = NULL, y = NULL) +
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
  labs(title = "Number of eQTL target genes", x = NULL, y = NULL) +
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

p_final <- p_heat / (leg_cri | leg_eqtl) + plot_layout(heights = c(1, 0.18))

# --------- Save (vector SVG) ----------------------------------------
out_svg <- "plots/interval/gold_gene_gwas_trait_heatmap_diagsplit.svg"
ggsave(out_svg, p_final, device = svglite, width = 22, height = 10, bg = "white")
# Optional: crop page to drawing using Inkscape CLI
# system2("inkscape", c("--export-area-drawing", "--export-type=svg",
#                       paste0("--export-filename=", sub(".svg$", "_cropped.svg", out_svg)),
#                       out_svg))
