cgenes = fread("Post_coloc_analysis/cres_with_grnas.txt") %>%
  filter(significant == 1) %>%
  distinct(ensembl_id)

egenes = fread("Post_coloc_analysis/cres_with_grna_eqtls.txt") %>%
  distinct(ensembl_id)


# obtain list of count files
count.files = system("ls eQTL_catalogue/count_means/*.csv", intern = T)
counts = fread(count.files[1])
name = str_split_fixed(str_split_fixed(count.files[1], "/",4)[,3], "\\.",2)[,1]
colnames(counts)[2] = name

# Load and merge into one file
for (x in 2:length(count.files)){
  print(x)
  counts2 = fread(count.files[x])
  name = str_split_fixed(str_split_fixed(count.files[x], "/",4)[,3], "\\.",2)[,1]
  colnames(counts2)[2] = name
  counts = counts %>% left_join(counts2, "phenotype_id")
}

counts[,ensembl_id := str_split_fixed(counts$phenotype_id, "\\.",2)[,1]]
counts.cgenes = counts[ensembl_id %in% cgenes$ensembl_id,]
counts.cgenes <- counts.cgenes[, 2:(ncol(counts.cgenes))]
max_exp = apply(counts.cgenes[,1:(ncol(counts.cgenes)-1)], 1, max, na.rm=T)

# Reshape data to tidy format
tidy_data <- counts.cgenes %>%
  pivot_longer(cols = 1:(ncol(counts.cgenes)-1), names_to = "eqtl", values_to = "Mean_counts")

overlap.genes = unique(tidy_data$ensembl_id[tidy_data$ensembl_id %in% egenes$ensembl_id])
nonoverlap.genes = unique(tidy_data$ensembl_id[!tidy_data$ensembl_id %in% egenes$ensembl_id])
tidy_data$ensembl_id <- factor(tidy_data$ensembl_id, levels = c(overlap.genes,nonoverlap.genes))

png("Post_coloc_analysis/plots/cgene_heatmap.png", width = 24, height = 8,units = "in", res = 300)

heatmap_plot <- tidy_data %>%
  ggplot(aes(y = eqtl, x = ensembl_id, fill = Mean_counts)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#3C94FE","#FE4E3C"),  # Define the color scale
    values = scales::rescale(c(0, 1, 10)),  # Define custom intervals (e.g., 0, 5, 10)
    limits = c(0, 10),  # Define the range of the color scale
    oob = scales::squish  # Cap values outside the limits
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = ifelse(levels(tidy_data$ensembl_id) %in% nonoverlap.genes, "red","black")), 
    axis.text.y = element_text(size = 10),
    legend.key.height = unit(1, "cm")  # Adjust legend key height
  ) +
  labs(
    title = "Heatmap of Mean Counts",
    x = NULL,  # Remove x-axis title
    y = NULL,  # Remove y-axis title
    fill = "Mean Counts"  # Legend title
  )
heatmap_plot

dev.off()

# Now vice versa
K562_bulk = fread("/gpfs/commons/groups/lappalainen_lab/jmorris/210205_STINGseq-v2/data/K562_scRNAseq_bulkRNAseq.txt") %>%
  mutate(ensembl_id = str_split_fixed(ENSG, "\\.",2)[,1])

tidy_data2 = egenes %>%
  left_join(K562_bulk[,c("ensembl_id","TPM")], "ensembl_id") %>%
  mutate(TPM = ifelse(is.na(TPM) == T, 0, TPM))

overlap.genes = unique(tidy_data2$ensembl_id[tidy_data2$ensembl_id %in% cgenes$ensembl_id])
nonoverlap.genes = unique(tidy_data2$ensembl_id[!tidy_data2$ensembl_id %in% cgenes$ensembl_id])
tidy_data2$ensembl_id <- factor(tidy_data2$ensembl_id, levels = c(overlap.genes,nonoverlap.genes))

png("Post_coloc_analysis/plots/egene_heatmap.png", width = 24, height = 4,units = "in", res = 300)

heatmap_plot <- tidy_data2 %>%
  ggplot(aes(y = "K562", x = ensembl_id, fill = TPM)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#3C94FE","#FE4E3C"),  # Define the color scale
    values = scales::rescale(c(0, 1, 10)),  # Define custom intervals (e.g., 0, 5, 10)
    limits = c(0, 10),  # Define the range of the color scale
    oob = scales::squish  # Cap values outside the limits
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = ifelse(levels(tidy_data2$ensembl_id) %in% nonoverlap.genes, "red","black")), 
    axis.text.y = element_text(size = 10),
    legend.key.height = unit(1, "cm")  # Adjust legend key height
  ) +
  labs(
    title = "Heatmap of Mean Counts",
    x = NULL,  # Remove x-axis title
    y = NULL,  # Remove y-axis title
    fill = "Mean Counts"  # Legend title
  )
heatmap_plot

dev.off()