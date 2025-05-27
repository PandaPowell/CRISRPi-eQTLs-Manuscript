rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load stingseq data
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`)

cor(sting_seq$`Log2 fold-change`[sting_seq$`Q-value (1 Mb)` < 0.1], abs(sting_seq$`TSS Distance`[sting_seq$`Q-value (1 Mb)`<0.1]))

# Subset your data first
subset_data <- sting_seq[sting_seq$`Q-value (100 kb)` < 0.05, ]

# Calculate the correlation
cor_estimate <- cor(subset_data$`Log2 fold-change`, abs(subset_data$`TSS Distance`))

# Create the scatter plot
p <- ggplot(subset_data, aes(x = `Log2 fold-change`, y = abs(`TSS Distance`))) +
  geom_point(alpha = 0.6) +  # scatter points
  annotate("text", 
           x = Inf, y = Inf, 
           label = paste0("cor = ", round(cor_estimate, 2)), 
           hjust = 1.1, vjust = 1.5, 
           size = 5) + 
  theme_minimal() +
  labs(x = "Log2 fold-change", y = "Absolute TSS Distance")

# Show the plot
print(p)

# Load and filter data
resample_results <- fread("CRISPR_data/resampling_results.txt") %>%
  unique() %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two", !is.na(target_site.start))

# Compute distance
resample_results$distance <- resample_results$target_site.start - resample_results$TSS

# Filter to rejected sites
rejected <- resample_results %>% filter(rejected == TRUE, abs(distance) < 100000)

# Compute correlation
cor_estimate <- cor(rejected$xi, abs(rejected$distance))

# Make scatter plot
ggplot(rejected, aes(x = xi, y = abs(distance))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("cor = ", round(cor_estimate, 2)),
           hjust = 1.1, vjust = 1.5,
           size = 5) +
  theme_minimal() +
  labs(x = "CRISPR xi (effect size)", y = "Absolute TSS Distance")
