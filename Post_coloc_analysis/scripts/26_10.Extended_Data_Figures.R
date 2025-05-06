rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(UpSetR)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

total_cres = unique(cres_w_grnas$grna_target)
grna_cres_w_cgenes = unique(cres_w_grnas$grna_target[cres_w_grnas$significant == 1])
target_cgenes = unique(cres_w_grnas$target_gene[cres_w_grnas$significant == 1])
cgenes = unique(cres_w_grnas$ensembl_id[cres_w_grnas$significant == 1])
grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
target_egenes = unique(cres_w_grnas_egene$target_gene)
egenes = unique(cres_w_grnas_egene$ensembl_id)
cres_no_target = unique(cres_w_grnas$grna_target[!cres_w_grnas$grna_target %in% c(grna_cres_w_cgenes,grna_cres_w_egenes)])
overlapping_cres = grna_cres_w_cgenes[grna_cres_w_cgenes %in% grna_cres_w_egenes]
overlapping_cres_cgene = unique(cres_w_grnas$target_gene[cres_w_grnas$grna_target %in% overlapping_cres & cres_w_grnas$significant == 1])
overlapping_cres_egene = unique(cres_w_grnas_egene$target_gene[cres_w_grnas_egene$grna_target %in% overlapping_cres])
overlapping_target_genes = overlapping_cres_cgene[overlapping_cres_cgene %in% overlapping_cres_egene]

# Extended Data Fig. Number of cGenes by CRISPRi study type
svg("plots/supp_figs/Number_cGenes_by_CRISPRi_studytype.svg", width = 6, height = 4)

cres_w_grnas %>% 
  filter(significant == 1) %>%
  group_by(data) %>%
  summarise(n = n()) %>%
  ggplot( aes(x = data, y = n)) +
  geom_bar(stat = "identity", fill = "#adadad", width = 0.5, position = position_dodge(width = 0.7)) +
  theme_cowplot() +
  labs(title = "", x = "CRISPRi dataset", y = "N significant CRE-genes") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),  # Rotate the x-axis labels
    axis.text.y = element_text(size = 14)
  ) + coord_flip()

dev.off()

############################################################################################
# Extended Data Fig. Number of eGenes by eQTL study type
# 
# Define the sc_qtls
colors <- c(
  "Bulk_tissue_eQTL" = "#014d64",
  "Sc_eQTL" = "#7ad2f6",
  "Bulk_cell_type_eQTL" = "#ebebeb"
)

sc_qtls <- c("NK cells", "B cells", "CD4 T cells", "CD8 T cells", "DC mean", "Mono cells", "other cells", "other T cells")

# Create the ggplot with ordered bars and angled x-axis labels
p1 = cres_w_grnas_egene %>%
  mutate(
    uniq_id = paste(eqtl, ensembl_id, grna_target, sep = "_"),
    eqtl_study = paste(eqtl, study_id, sep = "_")
  ) %>% 
  distinct(uniq_id, .keep_all = TRUE) %>% 
  group_by(eqtl_study) %>% 
  summarise(N = n()) %>%
  mutate(eqtl_study = trimws(gsub("_", " ", eqtl_study))) %>%
  arrange(N) %>%  # Arrange by N to ensure order
  mutate(eqtl_study = factor(eqtl_study, levels = unique(eqtl_study)),  # Use unique values for factor levels
         group = ifelse(eqtl_study %in% c("blood QTS000019", "GTEx", "blood QTS000029"), "Bulk_tissue_eQTL", 
                        ifelse(eqtl_study %in% sc_qtls, "Sc_eQTL", "Bulk_cell_type_eQTL"))) %>%
  ggplot( aes(x = eqtl_study, y = N, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = colors) +
  theme_cowplot() +
  labs(title = "", x = "eQTL dataset", y = "N associated CRE-genes", fill = "eQTL type") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),  # Rotate the x-axis labels
    axis.text.y = element_text(size = 12),
    legend.position = "top",  # Move the legend to the left side of the plot
  ) + coord_flip()

svg("plots/supp_figs/eGenes_by_eQTL_study_type.svg", width = 12, height = 8)
p1
dev.off()

# Extended Data Fig. Intersecting eGenes per eQTL dataset
svg("plots/supp_figs/Intersecting_eGenes_per_eQTL_dataset.svg", width = 12, height = 8)

overlapping_eqtls <- cres_w_grnas_egene %>%
  filter(target_gene %in% overlapping_cres_cgene)

# Create the ggplot with ordered bars and angled x-axis labels
overlapping_eqtls %>%
  mutate(
    uniq_id = paste(eqtl, ensembl_id, grna_target, sep = "_"),
    eqtl_study = paste(eqtl, study_id, sep = "_")
  ) %>% 
  distinct(uniq_id, .keep_all = TRUE) %>% 
  group_by(eqtl_study) %>% 
  summarise(N = n()) %>%
  mutate(eqtl_study = trimws(gsub("_", " ", eqtl_study))) %>%
  arrange(N) %>%  # Arrange by N to ensure order
  mutate(eqtl_study = factor(eqtl_study, levels = unique(eqtl_study)),  # Use unique values for factor levels
         group = ifelse(eqtl_study %in% c("blood QTS000019", "GTEx", "blood QTS000029"), "Bulk_tissue_eQTL", 
                        ifelse(eqtl_study %in% sc_qtls, "Sc_eQTL", "Bulk_cell_type_eQTL"))) %>%
  ggplot( aes(x = eqtl_study, y = N, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = colors) +
  theme_cowplot() +
  labs(title = "", x = "eQTL dataset", y = "N intersecting CRE-eGenes and cGenes", fill = "eQTL type") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),  # Rotate the x-axis labels
    axis.text.y = element_text(size = 12),
    legend.position = "top",  # Move the legend to the left side of the plot
  ) + coord_flip()

dev.off()
# 
# Extended Data Fig. eQTL power increases with MAF
# "eqtl_maf_power_barplot.png" 

# Extended Data Fig. CRISPRi power increases with gene expression and number of cells
#  Jaspers plot

# Extended Data Fig. Upset plot of eGenes by eQTL study
cres_data <- cres_w_grnas_egene %>%
  mutate(
    uniq_id   = paste(gene_name, grna_target, sep = "_"),
    eqtl_study = paste(eqtl, study_id, sep = "_")
  ) %>%
  distinct(grna_target, eqtl_study)   # keep only one row per (uniq_id, eqtl_study)

# 2) Pivot to wide format: each row is a gene/uniq_id, 
#    and each column is an eQTL dataset indicating 1 = present, 0 = absent
cres_wide <- cres_data %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from  = eqtl_study,
    values_from = present,
    values_fill = 0
  )

cres_wide <- as.data.frame(cres_wide)
# 3) Make an UpSet plot of intersections across eQTL datasets
#    By default, UpSetR expects all "set" columns to contain 0/1 membership flags.
#    Here we exclude 'uniq_id' itself from the sets.
set_cols <- setdiff(colnames(cres_wide), "grna_target")

svg("plots/supp_figs/Upset_CREs_eGenes_by_study.svg", width = 4, height = 6)

upset(
  data = cres_wide,
  sets = set_cols,
  nsets = length(set_cols),    # or fewer if you only want the top sets
  number.angles = 0, 
  order.by = "freq",
  mb.ratio = c(0.60, 0.40),
  text.scale = 1.2
)

dev.off

## Number of studies per CRE
studies_per_cre = cres_w_grnas_egene %>%
  mutate(eqtl_study = paste(eqtl, study_id, sep = "_")) %>%
  distinct(grna_target, gene_name, eqtl_study) %>%     # remove redundant gene-study-gRNA combinations
  distinct(grna_target, gene_name, .keep_all = T) %>% # Remove studies that target the same gene at a CRE
  group_by(grna_target) %>%                 
  summarise(n_studies = n()) %>%                       # count how many studies per gRNA
  arrange(desc(n_studies))

median(studies_per_cre$n_studies)

p_hist <- ggplot(studies_per_cre, aes(x = n_studies)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#4C72B0") +
  geom_vline(aes(xintercept = median(n_studies)), linetype = "dashed", color = "red", size = 1) +
  labs(
    x = "Number of eQTL studies per CRE (unique genes only)",
    y = "Number of CREs",
    title = NULL
  ) +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )
svg("plots/supp_figs/Number_of_studies_per_CRE.svg", width = 6, height = 4)

p_hist

dev.off()
# Extended Data Fig. Hi-C bar plots
# Load Hi-c data
hic = fread("Hi_C/K562/K562.hg19.AllInteractions.SP4.FDR0.1.txt")

intersect = fread("Hi_C/K562/cres_w_grnas_HiC_interactions.bed") %>%
  left_join(hic, "InteractorID")

# Load intersected Hi-c and grnas
cgene_hic = cres_w_grnas %>%
  filter(significant ==1) %>%
  left_join(intersect,"grna_target", relationship = "many-to-many") %>%
  filter(is.na(InteractorID)==F)

# Load intersected Hi-c and grnas
egene_hic = cres_w_grnas_egene %>%
  left_join(intersect,"grna_target", relationship = "many-to-many") %>%
  filter(is.na(InteractorID)==F)

# Fisher exact test
#                       cgene                      egene
# HiCoverlap              a                          b
# HiC non-overlap         c                          d

overlaping_cgenes = length(unique(cgene_hic$target_gene[cgene_hic$gene_name == cgene_hic$RefSeqName]))
overlaping_egenes = length(unique(egene_hic$target_gene[egene_hic$gene_name == egene_hic$RefSeqName]))
total_cgene_targets = length(unique(cres_w_grnas$target_gene[cres_w_grnas$significant == 1]))
total_egene_targets = length(unique(cres_w_grnas_egene$target_gene))

a = overlaping_cgenes
b = overlaping_egenes
c = total_cgene_targets-overlaping_cgenes
d = total_egene_targets-overlaping_egenes

print(matrix(c(a,b,c,d),nrow=2,ncol=2))
test = fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))

# Proportion of genes with pLI > 0.9
data <- data.frame(
  Gene_type = factor(c("cGenes","eGenes"), levels = c("cGenes","eGenes")),
  Proportion = c(a/(a+c),b/(b+d)))

# Generate a bar plot for gene distances
svg("plots/supp_figs/HiC_barplots.svg", width = 4, height = 6)

ggplot(data, aes(x = Gene_type, y = Proportion, fill = Gene_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = colors) +
  labs(title = "", x = "", y = "Proportion of CRE-genes") +
  annotate("text", x = 1.5, y = 0.25, label = paste0("P = ", formatC(test$p.value, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  theme_cowplot() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )

dev.off()
# 
# Extended Data Fig. Logistic regression forest plot
# 

# Extended Data Fig. trans-cGenes barplots
# 
# Load Morris trans-gene results
trans_network_john <- readxl::read_xlsx("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/science.adh7699_table_s4.xlsx", skip = 2) %>%
  filter(`Q-value` < 0.1) %>%
  filter(`Network CRE` %in% c("GFI1B (CRE-1)", "HHEX","IKZF1","NFE2 (CRE-1)", "RUNX1")) # remove miRNAs and take most significant cre when multiple are present

# Barplot of number of trans-genes per cis-gene
n_trans_assoc <- trans_network_john %>%
  group_by(`Network CRE`) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(study = "Morris") %>%
  rename(cis_gene = `Network CRE`)

n_trans_assoc$cis_gene = gsub(" \\(CRE-1\\)", "", n_trans_assoc$cis_gene)

# Load Gasperini tran effects
trans_network_gas <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt") %>%
  filter(qvalue_trans < 0.10)

### Remove trans-results from genes with +/-5mbp of gRNAs
trans_network_gas$TSS = ifelse(trans_network_gas$strand == 1, trans_network_gas$start_position, trans_network_gas$end_position)
trans_network_gas$tss_distance = with(trans_network_gas, ifelse(gsub("chr","",gRNA_chr) == chromosome_name, as.integer(gRNA_end)-TSS, NA))
trans_network_gas = trans_network_gas[abs(tss_distance) > 5e+06 | is.na(tss_distance) == T,]

# Load gasperini data linking grnas to cis-target genes
gas_cis_data <- fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt")

# Make barplot of number of associated genes
# link cis-gene to gRNA target in trans data
gas_cis_data <- gas_cis_data[hgnc_symbol_closest %in% c(gas_cis_genes, morrs_cis_genes),]
trans_network_gas <- trans_network_gas[gRNA_target %in% gas_cis_data$gRNA_target,] %>%
  left_join(gas_cis_data[,c("gRNA_target","hgnc_symbol_closest")], by = "gRNA_target")

# Barplot of number of trans-genes per cis-gene
n_trans_assoc_gas <- trans_network_gas %>%
  group_by(hgnc_symbol_closest) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(study = "Gasperini") %>%
  rename(cis_gene = hgnc_symbol_closest)

svg(filename = "plots/supp_figs/trans_cgene_barplot.svg", width = 6, height = 4)

colors <- c(
  "Gasperini" = "#7c260b",
  "Morris" = "#ed111a"
)

# Create bar plot
ggplot(bind_rows(n_trans_assoc, n_trans_assoc_gas), aes(x = cis_gene, y = n, fill = study)) +
  geom_col(
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.8,
    alpha = 0.8
  ) +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "",
       x = "cis-gene",
       y = "n trans-genes")


dev.off()

# Extended Data Fig. GFI1B correlation 
overlapping_genes= morrs_cis_genes[morrs_cis_genes %in% gas_cis_genes]

gene = "GFI1B"

# Filter crispri associations to cis-gene
trans_j <- trans_network_john[grep(gene, trans_network_john$`Network CRE`),]
trans_g = trans_network_gas[hgnc_symbol_closest == gene,]

comb_networks = trans_j %>%
  left_join(trans_g, by = c("Gene" = "hgnc_symbol")) %>%
  filter(is.na(log2fc) == F)

r = cor(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)
p = cor.test(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)$p.value
p = ifelse(p<2e-16, 2e-16, p)

svg(filename = "plots/supp_figs/cor_plots_GFI1B.svg", width = 4, height = 4)

ggplot(comb_networks, aes(x=`Log2 Fold Change`, y=log2fc)) + 
  geom_point(color = "black", alpha = 0.6) +
  annotate("text", x = -1.5, y = 1.5, label = paste0("r = ",round(r,2),", P < ", formatC(p, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  ylab("Trans-gene LogFC - Gasperini dataset") +
  xlab("Trans-gene LogFC - Morris dataset") +
  theme_cowplot()

dev.off()

# Extended Data Fig. GFI1B correlation 
gene = "HHEX"

# Filter crispri associations to cis-gene
trans_j <- trans_network_john[grep(gene, trans_network_john$`Network CRE`),]
trans_g = trans_network_gas[hgnc_symbol_closest == gene,]

comb_networks = trans_j %>%
  left_join(trans_g, by = c("Gene" = "hgnc_symbol")) %>%
  filter(is.na(log2fc) == F)

r = cor(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)
p = cor.test(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)$p.value
p = ifelse(p<2e-16, 2e-16, p)

svg(filename = "plots/supp_figs/cor_plot_HHEX.svg", width = 4, height = 4)

ggplot(comb_networks, aes(x=`Log2 Fold Change`, y=log2fc)) + 
  geom_point(color = "black", alpha = 0.6) +
  annotate("text", x = -0.5, y = 1, label = paste0("r = ",round(r,2),", P = ", formatC(p, format = "e", digits = 1)), size = 5, hjust = 0.5, family = "sans") +
  ylab("Trans-gene LogFC - Gasperini dataset") +
  xlab("Trans-gene LogFC - Morris dataset") +
  theme_cowplot()

dev.off()

# Extended Data Fig. trans-eGenes barplots

# Extended Data Fig enrichment forest plot
# 
# Extended Data Fig. Locus heatmaps for all gold genes
