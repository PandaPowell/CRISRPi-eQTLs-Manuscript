rm(list=ls())
options(bitmapType="cairo")
# Load libraries
library(ggplot2)
library(pROC)
library(dplyr)
library(data.table)
library(tidyr)
library(cowplot)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

# Load Gold Standard Genes
gold_genes <- fread("cis_gold_genes.txt")

# Split grna_target into separate rows for better filtering
gold_genes <- gold_genes[, .(grna = unlist(strsplit(grna_target, ","))), by = gold_gene]

# Load cgenes and filter to cres with a gold-gene within 1 mb and a signifcant cgene
cres_w_grnas <- fread("cres_with_grnas.txt") 

cgene_cres = cres_w_grnas %>%
  filter(grna_target %in% gold_genes$grna) %>%
  mutate(gold_overlap = ifelse(ensembl_id %in% gold_genes$gold_gene,1,0)) %>%
  filter(is.na(pvalue)==F) %>%
  dplyr::select(grna_target, ensembl_id, target_gene, significant)

# Load all cis-genes
eqtl_distance = fread("eQTL_gene_distances/eqtl_cat_gene_distances.txt") %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  distinct(target_gene, .keep_all = T) %>%
  filter(grna_target %in% gold_genes$grna) %>%
  select(grna_target, ensembl_id, target_gene) %>%
  mutate(significant =0)

gold_cgene_cres = bind_rows(cgene_cres, eqtl_distance) %>%
  distinct(target_gene,.keep_all = T) %>%
  mutate(gold_overlap = ifelse(ensembl_id %in% gold_genes$gold_gene,1,0))

conf_matrix = table(gold = gold_cgene_cres$gold_overlap, cgene = gold_cgene_cres$significant)

## Extract values from confusion matrix
TN <- conf_matrix[1,1]  # True Negatives
FP <- conf_matrix[1,2]  # False Positives
FN <- conf_matrix[2,1]  # False Negatives
TP <- conf_matrix[2,2]  # True Positives

## Compute Sensitivity and Specificity
sensitivity <- TP / (TP + FN)  # aka Recall
specificity <- TN / (TN + FP)

## Print results
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")

# Load egenes and filter to cres with a gold-gene within 1 mb and a signifcant cgene
cres_w_grnas_egene <- fread("cres_with_grna_eqtls.txt") %>%
  filter(grna_target %in% gold_genes$grna) %>%
  distinct(target_gene, .keep_all = T) %>%
  select(grna_target, ensembl_id, target_gene) %>%
  mutate(significant = 1)
  
# Load all cis-genes
eqtl_distance = fread("eQTL_gene_distances/eqtl_cat_gene_distances.txt") %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  distinct(target_gene, .keep_all = T) %>%
  filter(grna_target %in% gold_genes$grna) %>%
  select(grna_target, ensembl_id, target_gene) %>%
  mutate(significant =0)

# K562 genes
cgene_cres = cres_w_grnas %>%
  filter(grna_target %in% gold_genes$grna) %>%
  filter(is.na(pvalue)==F) %>%
  select(grna_target, ensembl_id, target_gene, significant) %>%
  mutate(significant = 0)

gold_egene_cres = bind_rows(cres_w_grnas_egene,eqtl_distance, cgene_cres) %>%
  distinct(target_gene, .keep_all = T) %>%
  mutate(gold_overlap = ifelse(ensembl_id %in% gold_genes$gold_gene,1,0))

conf_matrix = table(gold = gold_egene_cres$gold_overlap, cgene = gold_egene_cres$significant)

## Extract values from confusion matrix
TN <- conf_matrix[1,1]  # True Negatives
FP <- conf_matrix[1,2]  # False Positives
FN <- conf_matrix[2,1]  # False Negatives
TP <- conf_matrix[2,2]  # True Positives

## Compute Sensitivity and Specificity
sensitivity <- TP / (TP + FN)  # aka Recall
specificity <- TN / (TN + FP)

## Print results
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")

par(pty = "s")
roc(gold_cgene_cres$gold_overlap, glm.fit$fitted.values, plot=TRUE)

# Compare models
comparison_table = table(model1 = gold_egene_cres$significant %in% gold_egene_cres$gold_overlap, 
                         model2 = gold_cgene_cres$significant %in% gold_cgene_cres$gold_overlap)


par(pty = "s")
roc(gold_cgene_cres$gold_overlap, gold_cgene_cres$significant, plot=TRUE)



roc(gold_cgene_cres$gold_overlap, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)


conf_matrix = table(gold = gold_cgene_cres$gold_overlap, cgene = gold_cgene_cres$significant)

## Convert the table into a data frame
df_conf <- as.data.frame(conf_matrix)
## By default, columns will be "Actual", "Predicted", and "Freq"

## Make a tile plot
ggplot(df_conf, aes(x = gold, y = cgene, fill = Freq)) +
  geom_tile() +
  ## Add the frequency text on top of each tile
  geom_text(aes(label = Freq), color = "white", size = 5) +
  ## Move x-axis to the top (optional)
  scale_x_discrete(position = "top") +
  ## Use a simple gradient from white to blue
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Confusion Matrix Heatmap") +
  theme_minimal()

cres_w_grnas_egene <- fread("cres_with_grna_eqtls.txt") %>%
  filter(grna_target %in% gold_genes$grna)
egenes <- unique(cres_w_grnas_egene$ensembl_id)

closest.genes <- fread("closest_gene_to_GWAS_CRE.txt") %>%
  filter(grna_target %in% gold_genes$grna)
closest.genes = unlist(closest.genes$ensembl_id)

# Function to compute binary labels for a given predicted gene list
get_binary_labels <- function(pred_list, true_list) {
  as.integer(pred_list %in% true_list)  # 1 if gene is in gold standard, 0 otherwise
}

# Compute binary labels for each method
labels_cgene <- get_binary_labels(cgenes, unique(gold_genes$gold_gene))
labels_egene <- get_binary_labels(egenes, unique(gold_genes$gold_gene))
labels_closest <- get_binary_labels(closest.genes, unique(gold_genes$gold_gene))

# Assign scores (ranking-based, replace with actual scores if available)
scores_cgene <- rev(seq_along(cgenes))
scores_egene <- rev(seq_along(egenes))
scores_closest <- rev(seq_along(closest.genes))

# Compute ROC curves
roc_cgene <- roc(labels_cgene, scores_cgene)
roc_egene <- roc(labels_egene, scores_egene)
roc_closest <- roc(labels_closest, scores_closest)

# Calculate AUC values
auc_cgene <- auc(roc_cgene)
auc_egene <- auc(roc_egene)
auc_closest <- auc(roc_closest)

# Convert ROC data into a tidy dataframe for ggplot2
roc_data <- bind_rows(
  data.frame(Method = "CRISPRi gene", TPR = roc_cgene$sensitivities, FPR = 1 - roc_cgene$specificities),
  data.frame(Method = "eQTL gene", TPR = roc_egene$sensitivities, FPR = 1 - roc_egene$specificities),
  data.frame(Method = "Closest gene", TPR = roc_closest$sensitivities, FPR = 1 - roc_closest$specificities)
)

# Define custom colors matching previous plot
custom_colors <- c("CRISPRi gene" = "#F8766D", "eQTL gene" = "#00BFC4", "Closest gene" = "purple")

# Generate AUC Labels
auc_labels <- c(
  paste0("CRISPRi gene (AUC = ", round(auc_cgene, 2), ")"),
  paste0("eQTL gene (AUC = ", round(auc_egene, 2), ")"),
  paste0("Closest gene (AUC = ", round(auc_closest, 2), ")")
)

#svg("plots/figure_plots/ROC_plot.svg", width = 6, height = 6)
png("plots/figure_plots/ROC_plot.png", width = 6, height = 6,units = "in", res = 300)
# Plot ROC curves using ggplot2
ggplot(roc_data, aes(x = FPR, y = TPR, color = Method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = custom_colors) +
  labs(title = "", x = "False Positive Rate", y = "True Positive Rate") +
  theme_cowplot() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  annotate("text", x = 0.8, y = 0.5, label = auc_labels[1], size = 5, hjust = 0.5, family = "sans", color = "#F8766D") +
  annotate("text", x = 0.8, y = 0.45, label = auc_labels[3], size = 5, hjust = 0.5, family = "sans", color = "purple") +
  annotate("text", x = 0.8, y = 0.40, label = auc_labels[2], size = 5, hjust = 0.5, family = "sans", color = "#00BFC4") +
  geom_abline(linetype = "dashed", color = "gray")  # Diagonal reference line
dev.off()

# Function to compute confusion matrix components directly from binary labels
calculate_confusion_matrix <- function(true_labels, predicted_labels) {
  TP <- sum(true_labels == 1 & predicted_labels == 1)  # True Positives
  FN <- sum(true_labels == 1 & predicted_labels == 0)  # False Negatives
  FP <- sum(true_labels == 0 & predicted_labels == 1)  # False Positives
  TN <- sum(true_labels == 0 & predicted_labels == 0)  # True Negatives
  
  return(data.frame(TP = TP, FN = FN, FP = FP, TN = TN))
}

# Since the true labels (`labels_cgene`, `labels_egene`, `labels_closest`) are already binary,
# we assume predicted labels are the same as the given labels for now (since there's no ranking).

conf_matrix_cgene <- calculate_confusion_matrix(labels_cgene, labels_cgene)
conf_matrix_egene <- calculate_confusion_matrix(labels_egene, labels_egene)
conf_matrix_closest <- calculate_confusion_matrix(labels_closest, labels_closest)

# Print confusion matrices
print("Confusion Matrix for CRISPRi Genes")
print(conf_matrix_cgene)

print("Confusion Matrix for eQTL Genes")
print(conf_matrix_egene)

print("Confusion Matrix for Closest Genes")
print(conf_matrix_closest)
