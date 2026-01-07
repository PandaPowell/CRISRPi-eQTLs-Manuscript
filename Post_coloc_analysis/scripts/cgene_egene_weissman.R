rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(cowplot)
library(ggplot2)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

weissman = fread("processed_data/Weissman_cis_gene_de.csv")
sig_weissman = weissman %>% filter(significant_bonferroni == 0)
# Load CRISPR eQTL and GRNA target data
cres_w_grnas_egene <- fread("processed_data/cres_with_grna_eqtls.txt") %>%
  arrange(desc(PP.H4.abf)) %>% 
  distinct(ensembl_id,.keep_all = T) %>%
  filter(ensembl_id %in% weissman$ensembl_id)

cres_w_grnas <- fread("processed_data/cres_with_grnas.txt") %>%
  filter(significant ==1) %>%
  distinct(ensembl_id,.keep_all = T) %>%
  filter(ensembl_id %in% weissman$ensembl_id)

a = sum(cres_w_grnas_egene$ensembl_id %in% sig_weissman$ensembl_id)
b = nrow(cres_w_grnas_egene) - a
c = sum(cres_w_grnas$ensembl_id %in% sig_weissman$ensembl_id)
d = nrow(cres_w_grnas) - c
table <- matrix(c(a, b, c, d), nrow = 2)
print(table)
test = fisher.test(table)

colors <- c(
  "cGenes" = "#d35e60",
  "eGenes" = "#76c0c1",
  "GWAS genes" = "#014d64"
)

# Proportion of genes with pLI > 0.9
data <- data.frame(
  Gene_type = factor(c("cGenes","eGenes"), levels = c("cGenes","eGenes")),
  Proportion = c(a/(a+b),c/(c+d)))
data
#png("plots/figure_plots/pli_plot.png", width = 4, height = 6,units = "in", res = 300)
#svg("plots/figure_plots/pli_plot.svg", width = 4, height = 6)

weiss_plot = ggplot(data, aes(x = Gene_type, y = Proportion, fill = Gene_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  labs(title = "", x = "", y = "Proportion of genes \n with pLi > 0.9", color = "Gene Type") +
  annotate("text", x = 1.5, y = 0.35, label = paste0("P = ", formatC(test$p.value, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )

weiss_plot

dev.off()
