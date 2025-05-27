rm(list=ls())
library(data.table)
library(tidyverse)

# Load crispr and eqtl targets
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

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
  Gene_type = factor(c("cGene","eGene"), levels = c("cGene","eGene")),
  Proportion = c(a/(a+c),b/(b+d)))

# Generate a bar plot for gene distances
png("plots/HiC_barplots.png", width = 6, height = 8,units = "in", res = 300)

ggplot(data, aes(x = Gene_type, y = Proportion, fill = Gene_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  theme_minimal() +
  labs(title = "", x = "", y = "Proportion of CRE-genes") +
  annotate("text", x = 1.5, y = 0.31, label = paste0("P-value: ", round(test$p.value, 2)), size = 7, hjust = 0.5) +
  theme_cowplot() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 24)
  )

dev.off()
