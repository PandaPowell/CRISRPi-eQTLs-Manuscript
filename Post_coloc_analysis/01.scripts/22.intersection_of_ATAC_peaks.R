rm(list=ls())
options(bitmapType="cairo")
library(ggplot2)
library(tidyverse)
library(data.table)
library(cowplot)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

atac_overlap = fread("ATAC_overlap/intersecting_ATAC_peaks.bed") %>%
  mutate(cell = str_split_fixed(cell_peak,"_",3)[,1])

# what is the number of overlaps of each cell-type
overlaps_by_cell = atac_overlap %>%
  group_by(cell) %>%
  summarise("overlapping peaks" = n()) %>%
  filter(cell != ".")

cell_name_lookup <- tibble::tibble(
  cell = c("Bcell", "CD4", "CD8", "CLP", "CMP", "Erythro",
           "GMP1low", "GMP2mid", "GMP3high", "HSC", "LMPP",
           "MEGA", "MEP", "MPP", "Mono", "NK", "mDC", "pDC"),
  full_name = c("B cell", "CD4⁺ T cell", "CD8⁺ T cell", "Common Lymphoid Progenitor",
                "Common Myeloid Progenitor", "Erythroblast",
                "Granulocyte-Monocyte Progenitor 1 (low)",
                "Granulocyte-Monocyte Progenitor 2 (mid)",
                "Granulocyte-Monocyte Progenitor 3 (high)",
                "Hematopoietic Stem Cell", "Lymphoid-Primed Multipotent Progenitor",
                "Megakaryocyte", "Megakaryocyte-Erythroid Progenitor",
                "Multipotent Progenitor", "Monocyte", "Natural Killer cell",
                "Myeloid Dendritic Cell", "Plasmacytoid Dendritic Cell"),
  total_peaks = c(129345, 99751, 94910, 100051, 190630, 69331, 98850, 118348, 120918,
                 150131, 131616, 123911, 194128, 160459, 78552, 91710, 115588, 79616)
)

overlaps_by_cell_named <- overlaps_by_cell %>%
  left_join(cell_name_lookup, by = "cell") %>%
  mutate(ratio = `overlapping peaks` / total_peaks) %>%
  arrange(desc(ratio)) %>%
  mutate(full_name = factor(full_name, levels = rev(full_name)))

png("plots/ATACseq_intersections.png", width = 6, height = 4,units = "in", res = 300)
# 2. Plot the ratio as bars
#    We'll keep the factor order (full_name) that you already established
ggplot(overlaps_by_cell_named, aes(x = full_name, y = ratio)) +
  geom_col() +
  labs(
    title = "Proportion of Overlapping Peaks vs. Total Peaks",
    subtitle = "overlapping_peaks / total_peaks",
    x = "Cell",
    y = "Ratio"
  ) +
  theme_cowplot() +
  coord_flip()

dev.off()

## Now look into peaks that intersect CREs
atac_overlap = fread("ATAC_overlap/intersecting_ATAC_peaks_w_cres.bed", fill = T) %>%
  mutate(cell = str_split_fixed(cell_peak,"_",3)[,1]) %>%
  distinct(k562_peak, gRNA_target, cell, .keep_all = T)

# Load crispr and eqtl targets
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

total_cres = unique(cres_w_grnas$grna_target)
grna_cres_w_cgenes = unique(cres_w_grnas$grna_target[cres_w_grnas$significant == 1])
grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
cres_no_target = length(total_cres) - length(unique(c(grna_cres_w_cgenes,grna_cres_w_egenes)))
overlapping_cres = grna_cres_w_cgenes[grna_cres_w_cgenes %in% grna_cres_w_egenes]
overlapping_cres_cgene = unique(cres_w_grnas$target_gene[cres_w_grnas$grna_target %in% overlapping_cres & cres_w_grnas$significant == 1])
overlapping_cres_egene = unique(cres_w_grnas_egene$target_gene[cres_w_grnas_egene$grna_target %in% overlapping_cres])
overlapping_target_genes = overlapping_cres_cgene[overlapping_cres_cgene %in% overlapping_cres_egene]

atac_counts = atac_overlap %>%
  filter(cell_peak != ".") %>%
  group_by(gRNA_target) %>%
  summarise(n = n())

zero_counts = data.frame(gRNA_target = atac_overlap$gRNA_target[!atac_overlap$gRNA_target %in% atac_counts$gRNA_target],
                         n = 0)

atac_counts = bind_rows(atac_counts, zero_counts) %>%
  mutate(cre_w_cgene = ifelse(gRNA_target %in% grna_cres_w_cgenes, 1, 0),
         cre_w_egene = ifelse(gRNA_target %in% grna_cres_w_egenes, 1, 0),
         intersecting_cres = ifelse(cre_w_cgene == 1 & cre_w_egene == 1, 1, 0),
         non_intersecting_cres = ifelse((cre_w_cgene == 1 & cre_w_egene == 0) | (cre_w_cgene == 0 & cre_w_egene == 1), 1, 0))

median(atac_counts$n[atac_counts$cre_w_cgene == 1])
median(atac_counts$n[atac_counts$cre_w_egene == 1])

test = wilcox.test(atac_counts$n[atac_counts$cre_w_cgene == 1], atac_counts$n[atac_counts$cre_w_egene == 1])

data <- data.frame(
  Gene_type = c(rep("cGene\nCREs", length(atac_counts$n[atac_counts$cre_w_cgene == 1])), rep("eGene\nCREs", length(atac_counts$n[atac_counts$cre_w_egene == 1]))),
  overlaps = c(atac_counts$n[atac_counts$cre_w_cgene == 1], atac_counts$n[atac_counts$cre_w_egene == 1]))

colors <- c(
  "cGene\nCREs" = "#d35e60",
  "eGene\nCREs" = "#76c0c1"
)

svg("plots/figure_plots/ATAC_plot.svg", width = 4, height = 6)

ggplot(data, aes(x = Gene_type, y = overlaps, fill = Gene_type)) +
  geom_boxplot(width = 0.35, outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  ylim(0, 20) +
  annotate("text", x = 1.5, y = 20, label = paste0("P = ", formatC(test$p.value, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  labs(x = "", y = "Number of cell types with\noverlapping ATAC-seq peak") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )

dev.off()

# calculate proportion of zeros:
# Number of zeros in group 1
x1 <- sum(atac_counts$cre_w_cgene == 1 & atac_counts$n == 0)
# Total in group 1
n1 <- sum(atac_counts$cre_w_cgene == 1)

# Number of zeros in group 2
x2 <- sum(atac_counts$cre_w_egene == 1 & atac_counts$n == 0)
# Total in group 2
n2 <- sum(atac_counts$cre_w_egene == 1)

# Compare the two proportions using a two-proportion test
test_result <- prop.test(c(x1, x2), c(n1, n2))

test_result

contingency_table <- matrix(
  c(x1, n1 - x1,
    x2, n2 - x2),
  nrow = 2, byrow = TRUE
)

fisher.test(contingency_table)

data <- data.frame(
  Gene_type = c(rep("Both", length(atac_counts$n[atac_counts$intersecting_cres == 1])), 
                rep("cGene", length(atac_counts$n[atac_counts$cre_w_cgene == 1 & atac_counts$cre_w_egene == 0])),
                rep("eGene", length(atac_counts$n[atac_counts$cre_w_egene == 1 & atac_counts$cre_w_cgene == 0]))),
  overlaps = c(atac_counts$n[atac_counts$intersecting_cres == 1], 
               atac_counts$n[atac_counts$cre_w_cgene == 1 & atac_counts$cre_w_egene == 0],
               atac_counts$n[atac_counts$cre_w_egene == 1 & atac_counts$cre_w_cgene == 0]))

test = kruskal.test(overlaps ~ Gene_type, data = data)

colors <- c(
  "Both" = "#8e6c8a",
  "cGene" = "#d35e60",
  "eGene" = "#76c0c1"
)

svg("plots/figure_plots/intersecting_ATAC_plot.svg", width = 4, height = 6)

ggplot(data, aes(x = Gene_type, y = overlaps, fill = Gene_type)) +
  geom_boxplot(width = 0.35, outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  ylim(0, 20) +
  annotate("text", x = 2, y = 20, label = paste0("P = ", formatC(test$p.value, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  labs(x = "", y = "Number of cell types with\noverlapping ATAC-seq peak") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

dev.off()
