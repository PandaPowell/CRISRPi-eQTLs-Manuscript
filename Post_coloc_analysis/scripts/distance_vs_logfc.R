rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load stingseq data
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`) %>%
  filter(`Q-value (1 Mb)` < 0.1) %>%
  mutate(tss_distance = abs(`TSS Distance`/10000))

# Temporarily rename the column in both model and new data
sting_seq2 <- sting_seq %>%
  rename(cells = `Cells bearing gRNAs`)

fit <- glm(`Log2 fold-change` ~ tss_distance + cells,
           family = "gaussian", data = sting_seq2)

new_data <- data.frame(
  tss_distance = seq(min(sting_seq2$tss_distance, na.rm = TRUE),
                     max(sting_seq2$tss_distance, na.rm = TRUE),
                     length.out = 100),
  cells = median(sting_seq2$cells, na.rm = TRUE)
)

new_data$predicted_lfc <- predict(fit, newdata = new_data)

# Plot
ggplot(sting_seq2, aes(x = tss_distance * -10000, y = `Log2 fold-change`)) + 
  geom_point(alpha = 0.5) +
  geom_line(data = new_data, aes(x = tss_distance * -10000, y = predicted_lfc), color = "blue") +
  labs(x = "TSS distance (bp)", y = "Log2 fold-change") +
  theme_minimal()
  
# Gasperini
resample_results = fread("CRISPR_data/resampling_results.txt") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two", is.na(target_site.start) == F) %>%
  left_join(cells_per_grna_df, by = c("grna_group" = "target.site")) %>%
  filter(rejected == TRUE) %>%
  mutate(tss_distance = abs(target_site.stop-TSS)/10000) %>%
  filter(xi < 0)

# Fit the model
fit2 <- glm(xi ~ tss_distance + cells_per_grna, family = "gaussian", data = resample_results)
summary(fit2)

# Create prediction data — fix cells_per_grna at its median
new_data <- data.frame(
  tss_distance = seq(min(resample_results$tss_distance, na.rm = TRUE),
                     max(resample_results$tss_distance, na.rm = TRUE),
                     length.out = 100),
  cells_per_grna = median(resample_results$cells_per_grna, na.rm = TRUE)
)

# Add predictions
new_data$predicted_xi <- predict(fit2, newdata = new_data)

# Plot with predicted line
ggplot(resample_results, aes(x = tss_distance * -10000, y = xi)) + 
  geom_point(alpha = 0.5) +
  geom_line(data = new_data, aes(x = tss_distance * -10000, y = predicted_xi), color = "blue") +
  labs(x = "TSS distance (bp)", y = "Effect size (xi)") +
  theme_minimal()

# Combined
# Make columns consistent in sting_seq2
sting_subset <- sting_seq2 %>%
  rename(
    xi = `Log2 fold-change`,
    cells_per_grna = cells
  ) %>%
  select(tss_distance, xi, cells_per_grna)

# Combine
combined <- rbind(
  resample_results[, .(tss_distance, xi, cells_per_grna)],
  sting_subset,
  fill = TRUE
)

fit_combined <- glm(xi ~ tss_distance + cells_per_grna, 
                    family = "gaussian", data = combined)
summary(fit_combined)

new_data_combined <- data.frame(
  tss_distance = seq(min(combined$tss_distance, na.rm = TRUE),
                     max(combined$tss_distance, na.rm = TRUE),
                     length.out = 100),
  cells_per_grna = median(combined$cells_per_grna, na.rm = TRUE)
)

new_data_combined$predicted_xi <- predict(fit_combined, newdata = new_data_combined)

svg("plots/supp_figs/logfc_vs_distance.svg", width = 6, height = 4)

ggplot(combined, aes(x = tss_distance * -10000, y = xi)) + 
  geom_point(alpha = 0.5) +
  geom_line(data = new_data_combined, 
            aes(x = tss_distance * -10000, y = predicted_xi), 
            color = "blue", size = 1) +
  labs(x = "TSS distance (bp)", y = "Effect size (xi)",
       title = "Regression of LogFC on distance")

dev.off()
