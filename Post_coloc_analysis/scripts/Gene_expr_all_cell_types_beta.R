rm(list=ls())

library(data.table)
library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(viridis)
library(cowplot)
library(gridExtra)
library(R.utils)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

cgene = read.table("cgenes.txt")
egene = read.table("egenes.txt")

df = data.frame(Gene = unique(c(cgene$V1,egene$V1)))
df$cgene = ifelse(df$Gene %in% cgene$V1, 1,0)
df$egene = ifelse(df$Gene %in% egene$V1, 1,0)

cell_names = c("B","CD4_T","CD8_T","DC","Mono","NK")

plot_list1 <- list()
plot_list2 = list()

for (i in 1:6){
  
  cell_name = cell_names[i]
  
  print(cell_name)
  
  onek = readRDS(paste0("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/cell_specific_rds/lvl1/",cell_name,"_raw_counts.RDS"))
  onek = CreateSeuratObject(counts = onek, project = "onek", min.cells = 3, min.features = 200)
  onek[["percent.mt"]] <- PercentageFeatureSet(onek, pattern = "^MT-")
  onek <- subset(onek, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  onek = subset(onek, features = c(cgene$V1,egene$V1))
  Idents(onek) = paste0(cell_name,"_cells")
  
  pseudo_counts = rowSums(onek[['RNA']]$counts)
  
  cell_df = data.frame(Gene = names(pseudo_counts),
                       sum_counts = pseudo_counts)
  cell_df$cgene = as.factor(ifelse(cell_df$Gene %in% cgene$V1, 1,0))
  cell_df$egene = as.factor(ifelse(cell_df$Gene %in% egene$V1, 1,0))
  
  test_counts = c(cell_df$sum_counts[cell_df$cgene ==1], cell_df$sum_counts[cell_df$egene ==1])
  test_group = factor(c(rep("cgene", length(cell_df$sum_counts[cell_df$cgene ==1])),
                        rep("egene", length(cell_df$sum_counts[cell_df$egene ==1]))))
  
  wilcox_test_result <- wilcox.test(test_counts ~ test_group, exact = FALSE)
  
  data <- data.frame(
    Counts = test_counts,
    GeneType = test_group
  )
  data$Log_counts = log(data$Counts)
  
  # Create the boxplot using ggplot2
  plot1 <- ggplot(data,aes(x = GeneType, y = Counts, fill = GeneType)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    theme_cowplot() + ylim(0,30000) + labs(title =cell_name, y = "sum_norm") +
    annotate("text", x = 1.5, y = 30000, label = paste0("P-value: ",round(wilcox_test_result$p.value,3)), size = 5, hjust = 0.5)
  
  plot2 = ggplot(data,aes(x = Counts, group = GeneType, fill = GeneType)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_cowplot() + xlim(0,50000) + labs(y = "sum_norm")
  
  plot_list1[[i]] <- plot1
  plot_list2[[i]] <- plot2
}

# Save the combined plots
png("plots/onek_boxplot_hist.png", width = 50, height = 16, units = "in", res = 300)
grid.arrange(grobs = c(plot_list1, plot_list2), ncol = 6, nrow = 2)
dev.off()

K562_bulk = fread("/gpfs/commons/groups/lappalainen_lab/jmorris/210205_STINGseq-v2/data/K562_scRNAseq_bulkRNAseq.txt") %>%
  mutate(cgene = ifelse(GENE %in% cgene$V1, 1,0), egene = ifelse(GENE %in% egene$V1, 1,0)) %>% rename(Gene = GENE)
K562_bulk$ENSG = str_split_fixed(K562_bulk$ENSG, "\\.", 2)[,1]

cgene_tpm = with(K562_bulk,zTPM[cgene == 1])
egene_tpm = with(K562_bulk,zTPM[egene == 1])
test_tpm = c(cgene_tpm, egene_tpm)
test_group = factor(c(rep("cgene", length(cgene_tpm)), rep("egene", length(egene_tpm))))

wilcox_test_result <- wilcox.test(test_tpm ~ test_group, exact = FALSE)

data <- data.frame(
  ztpm = test_tpm,
  GeneType = test_group
)

# Create the boxplot using ggplot2
plot1 <- ggplot(data,aes(x = GeneType, y =ztpm, fill = GeneType)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_cowplot() + labs(y = "Sum counts") +
  annotate("text", x = 1.5, y = 7, label = paste0("P-value: ",round(wilcox_test_result$p.value,3)), size = 5, hjust = 0.5)

plot3 = ggplot(data,aes(x = ztpm, group = GeneType, fill = GeneType)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_cowplot() + labs(y = "Sum counts")

# Save the combined plots
png("plots/histogram_boxplot_k562bulk.png", width = 12, height = 6, units = "in", res = 300)
grid.arrange(plot1, plot3, ncol = 2)
dev.off()

summits <- fread("data/k562_gene_dis_summits.bed") %>%
  rename(ENSG = gene) %>%
  left_join(K562_bulk, by = "ENSG") %>%
  drop_na() %>%
  mutate(distance = abs(TSS_pos - summit_pos)) %>%
  group_by(summit) %>%
  mutate(dis_rank = rank(distance, ties.method = "first"),
         count_rank = rank(TotalMolecules, ties.method = "first"),
         num_genes = n()) %>%
  ungroup() %>% filter(num_genes > 24)

summits_rank_plot = summits %>% filter(count_rank ==1, dis_rank < 25) %>% group_by(dis_rank) %>%
  summarise(n = n())

summits_count_plot = summits %>% filter(dis_rank < 25) %>% group_by(dis_rank) %>%
  summarise(sum_counts = sum(TotalMolecules))

plot1 = ggplot(summits_rank_plot, aes(x=dis_rank, y=n)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() +
  labs(title = "", x = "Distance rank", y = "Number of top expressed genes per CRE")

plot2 = ggplot(summits_count_plot, aes(x=dis_rank, y=sum_counts)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() +
  labs(title = "", x = "Distance rank", y = "Sum of gene counts per CRE")

# Save the combined plots
png("plots/rank_barplots.png", width = 12, height = 6, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

# summits where closest gene is most highly expressed
one = summits %>% filter(dis_rank == 1 & count_rank == 1)
# summits where second closest gene is most highly expressed
two = summits %>% filter(dis_rank == 2 & count_rank == 1)
# 3 or more
more = summits %>% filter(dis_rank >2 & count_rank == 1)

data <- data.frame(
  Rank = factor(c("1st","1st","2nd","2nd",">2nd", ">2nd"), levels = c("1st", "2nd", ">2nd")),
  method = c(rep(c("CRISPRi", "eQTL"),3)),
  percent = c(per_rank_one,median(eqtl_gene_dis$one),per_rank_two,median(eqtl_gene_dis$two),
              1-(per_rank_one+per_rank_two),median(eqtl_gene_dis$three)))

# Generate a bar plot for gene distances
tiff("plots/crispr_closest_gene_barplot.tiff", width = 6, height = 10,units = "in", res = 300)

ggplot(data, aes(x = Rank, y = percent, fill = method)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  theme_minimal() +
  labs(title = "", x = "Closest to CREs", y = "Proportion of total CREs") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

dev.off()


# Gasperini data
raw_data_dir_gasp = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/CRISPR_data/gasperini_2019/"

# URL of data
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="

# Gasperini et al results
all_deg_results_filename <- "GSE120861_all_deg_results.at_scale.txt"

# names of genes
genes_filename <- "GSE120861_at_scale_screen.genes.txt"

# names of cells
cells_filename <- "GSE120861_at_scale_screen.cells.txt"

# "reference cells" used by Gasperini et al for computational purposes
reference_cells_filename <- "GSE120861_50k_reference_cells.rds"

# all (gRNA, gene) pairs
gRNAgroup_pair_table_filename <- "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"

# list of gRNA groups used
gRNA_groups_filename <- "GSE120861_grna_groups.at_scale.txt"

# Monocle Cell Data Set object with all gRNA data
cds_filename <- "GSE120861_at_scale_screen.cds.rds"

# Expression data
expression_filename <- "GSE120861_at_scale_screen.exprs.mtx"

# list of files to download
filenames <- c(gRNAgroup_pair_table_filename,
               gRNA_groups_filename,
               all_deg_results_filename,
               genes_filename,
               cells_filename,
               reference_cells_filename,
               cds_filename,
               expression_filename)

options(timeout = 600)

# download files if not already present
for (filename in filenames) {
  if (!file.exists(paste0(raw_data_dir_gasp, "/", filename))) {
    cat(paste0("Downloading ", filename, "\n"))
    source <- paste0(remote, filename, ".gz")
    dest <- paste0(raw_data_dir_gasp, "/", filename, ".gz")
    download.file(source, dest)
    gunzip(paste0(dest))
  }
}

# Determine which cells have no gRNAs or NT gRNAs

# 2. Read monocole object
gc()
library(monocle, quietly = TRUE) # monocole required to load monocole CellDataSet object; we will extract only the cell-specific metadata
library(fst)
m_object <- readRDS(paste0(raw_data_dir_gasp, "/GSE120861_at_scale_screen.cds.rds"))

# 2a. save the coefficients of the mean-dispersion relationship
saveRDS(attr(m_object@dispFitInfo$blind$disp_func, "coefficients"), paste0(raw_data_dir_gasp, "/disp_coefficients.rds"))

# 2b. save the raw dispersions
write.fst(m_object@dispFitInfo$blind$disp_table, paste0(raw_data_dir_gasp, "/disp_table.fst"))

cell_metadata <- pData(m_object)
rm(m_object); gc()
covariates_cols <- 1:18
gRNA_cols <- 19:ncol(cell_metadata)

# 2c. Save the gRNA indicators
gRNA_indicators <- cell_metadata[,gRNA_cols]
write.fst(x = gRNA_indicators, compress = 0, path = paste0(raw_data_dir_gasp, "/gRNA_indicators.fst"))

# We perform a basic quality control for the Gasperini data; in particular, we restrict the cells in our analysis to those with at least 1 gRNA.
gRNA_indicator_matrix <- read.fst(paste0(raw_data_dir_gasp, "/gRNA_indicators.fst"))
# Remove guides that intersect finemapped snps
gas_targetsites = read.table("CRISPR_data/gasperini_finemap_grna_targetsites.txt")

gRNA_indicator_matrix = gRNA_indicator_matrix[,grep(paste(gas_targetsites$V1,collapse = "|"),colnames(gRNA_indicator_matrix))]

# Determine which cells have no finemap gRNA.
gRNA_indicator_counts <- apply(gRNA_indicator_matrix, 1, sum)
cells_to_keep <- which(gRNA_indicator_counts < 1) %>% as.integer()

# Filter gene expression matrix to only those cells
k562_scrna = readMM("CRISPR_data/gasperini_2019/GSE120861_at_scale_screen.exprs.mtx.gz")
gene_names = read.table("CRISPR_data/gasperini_2019/GSE120861_at_scale_screen.genes.txt")
rownames(k562_scrna) = gene_names$V1
print(nrow(k562_scrna))
print(ncol(k562_scrna))
k562_scrna2 = k562_scrna[,cells_to_keep]

pseudo_counts = rowSums(k562_scrna2)

k562_scrna2 = CreateSeuratObject(counts = k562_scrna2, project = "gasperini", min.cells = 3, min.features = 200)
Idents(k562_scrna2) = "k562"

pseduo_average = PseudobulkExpression(
  k562_scrna2,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  layer = "counts",
  method = "average",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE
)

cell_df = data.frame(Gene = rownames(k562_scrna),
                     sum_counts = pseudo_counts,
                     sum_norm = pseduo_average$RNA)

write.table(cell_df, "CRISPR_data/gasperini_2019/pseudobulk_scRNA.txt", quote = F, row.names =F, col.names = T)

# Load genes from morris et al
k562_genes = read.table("data/k562_genes.txt", header = T) %>% unique()
# Filter to genes in both studies
pseudobulk_scRNA = read.table("CRISPR_data/gasperini_2019/pseudobulk_scRNA.txt", header = T) %>%
  rename(ensembl_gene_id = Gene) %>% left_join(k562_genes, "ensembl_gene_id") %>% filter(is.na(gene)==F) %>%
  mutate(log_sum_counts = log(sum_counts)) %>%
  mutate(cgene = ifelse(gene %in% cgene$V1, 1,0), egene = ifelse(gene %in% egene$V1, 1,0)) %>%
  filter(cgene == 1 | egene == 1,cgene == 0 | egene == 0) %>%
  mutate(gene_type = case_when(
    cgene == 1 & egene == 0 ~ "cgene",
    cgene == 0 & egene == 1 ~ "egene"
  ))

cgene_counts = with(pseudobulk_scRNA,sum_counts[cgene == 1])
egene_counts = with(pseudobulk_scRNA,sum_counts[egene == 1])
test_counts = c(cgene_counts, egene_counts)
test_group = factor(c(rep("cgene", length(cgene_counts)), rep("egene", length(egene_counts))))

wilcox_test_result <- wilcox.test(test_counts ~ test_group, exact = FALSE)

data <- data.frame(
  Counts = test_counts,
  GeneType = test_group
)

# Create the boxplot using ggplot2
ggplot(data,aes(x = GeneType, y = Counts, fill = GeneType)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_cowplot() + ylim(0,500000)
  labs(y = "sum_norm")

pseudobulk_scRNA %>% ggplot(aes(x=sum_counts)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
  xlim(0,1000000)

pseudobulk_scRNA %>% ggplot(aes(x=log(sum_counts))) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

t.test(pseudobulk_scRNA$log_sum_counts[pseudobulk_scRNA$gene_type=="cgene"],pseudobulk_scRNA$log_sum_counts[pseudobulk_scRNA$gene_type=="egene"])

# Create and save the boxplot as a TIFF file
# Create the first plot
plot1 <- pseudobulk_scRNA %>%
  ggplot(aes(x = gene_type, y = sum_norm, fill = gene_type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_cowplot() +
  labs(y = "sum_norm")

# Create the second plot
plot2 <- K562_bulk %>%
  ggplot(aes(x = gene_type, y = TotalMolecules, fill = gene_type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_cowplot() + ylim(0,20000) +
  labs(y = "Gene counts")

# Save the combined plots
png("plots/boxplot_k562.png", width = 10, height = 8, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()


plot_3 = K562_bulk %>%
  ggplot(aes(x = zTPM, group = gene_type, fill = gene_type)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_cowplot() + xlim(-3,3) +
  labs(x = "zTPM")

plot_4 = K562_bulk %>%
  ggplot(aes(x = TotalMolecules, group = gene_type, fill = gene_type)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_cowplot() + xlim(0,30000) +
  labs(x = "Gene counts")

# Save the combined plots
png("plots/histogram_tpm.png", width = 8, height = 10, units = "in", res = 300)
grid.arrange(plot_3, plot_4, ncol = 1)
dev.off()