rm(list=ls())
library(data.table)
library(tidyverse)
library(gridExtra)
options(bitmapType="cairo")

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

interval = fread("cres_with_grna_eqtls_interval.txt") %>%
  filter(eqtl == "Interval") %>% 
  distinct(ensembl_id,.keep_all=T)

gtex = fread("cres_with_grna_eqtls_interval.txt") %>%
  filter(eqtl == "GTEx") %>% 
  distinct(ensembl_id,.keep_all=T)

gene_metrics = fread("data/gene_metrics_vgh_202407.tsv", fill=T)
colnames(gene_metrics) = c("ensemble_id", colnames(gene_metrics)[1:15])

metrics = colnames(gene_metrics)
plot_list <- list()  # Create an empty list to store plots

for (i in 2:length(gene_metrics)){
  
  m = metrics[i]
  cat("Running wilcox test for", m, "\n")
  
  # Extract cgene_metric
  interval_metric = gene_metrics[ensemble_id %in% interval$ensembl_id, ..i]
  interval_metric = na.omit(interval_metric)
  interval_metric <- unlist(interval_metric)
  interval_metric <- as.numeric(interval_metric)
  
  # Extract egene_metric
  gtex_metric = gene_metrics[ensemble_id %in% gtex$ensembl_id, ..i]
  gtex_metric = na.omit(gtex_metric)
  gtex_metric <- unlist(gtex_metric)
  gtex_metric <- as.numeric(gtex_metric)
  
  # Run Wilcoxon test
  test = wilcox.test(interval_metric, gtex_metric)
  print(test)
  
  # Combine the data into a data frame for plotting
  combined_data <- data.frame(
    Gene_type = c(rep("Interval", length(interval_metric)), rep("GTEx", length(gtex_metric))),
    Value = c(interval_metric, gtex_metric)
  )
  
  # Calculate the y-axis limits (just above the last outlier)
  stats_interval <- boxplot.stats(interval_metric)
  stats_gtex <- boxplot.stats(gtex_metric)
  
  # Max of the outliers for cgene and egene
  # Check if there are outliers; if not, use the max values from the metrics
  if (length(stats_interval$out) > 0 | length(stats_gtex$out) > 0) {
    max_outlier <- max(c(stats_interval$out, stats_gtex$out), na.rm = TRUE)
  } else {
    max_outlier <- max(c(interval_metric, gtex_metric), na.rm = TRUE)
  }
  
  # Set the lower limit as the minimum value and upper limit slightly above the max outlier
  min_value <- min(combined_data$Value, na.rm = TRUE)
  max_value <- max_outlier * 1.1  # 10% above the highest outlier
  
  # Create boxplot
  plot <- ggplot(combined_data, aes(x = Gene_type, y = Value, fill = Gene_type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "blue") +  # Adding individual values as jittered dots
    theme_cowplot() +
    annotate("text", x = 1.5, y = max_value, label = paste0("P-value: ", round(test$p.value, 4)), size = 5, hjust = 0.5) +
    labs(title = paste(m), x = "", y = m) + coord_cartesian(ylim = c(min_value, max_value)) 
  
  # Store the plot in the plot list
  plot_list[[i - 1]] <- plot
}

png("plots/interval/gene_metrics_boxplot_all_cgene_egene.png", width = 32, height = 10,units = "in", res = 300)
# Use grid.arrange to print all plots together
do.call(grid.arrange, c(plot_list, nrow = 2))
dev.off()

# Now just compare all eGenes not just those colocalized 

