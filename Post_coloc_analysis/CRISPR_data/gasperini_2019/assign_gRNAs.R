library(sceptre)
library(sceptredata)
library(Matrix)
library(data.table)
library(fst)
suppressWarnings(library(tidyverse))

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/CRISPR_data/gasperini_2019/")

deg = fread("GSE120861_all_deg_results.at_scale.txt")
#scep = readRDS("GSE120861_at_scale_screen.cds.rds")

expression_mtx <- Matrix::readMM("GSE120861_at_scale_screen.exprs.mtx.gz")
expression_mtx <- as(expression_mtx * 1, "dgCMatrix")  # Convert logical to numeric sparse matrix
grna_matrix = t(read_fst("gRNA_indicators.fst"))
grna_matrix <- as(grna_matrix * 1, "dgCMatrix")
gene_names <- fread("GSE120861_at_scale_screen.genes.txt", header = F)
gRNA_target_df = fread("GSE120861_grna_groups.at_scale.txt", header = F)
#colnames(gRNA_target_df) = 

all_deg_results <- suppressWarnings(read_tsv("GSE120861_all_deg_results.at_scale.txt", col_types = "cddddddccccciiciiccl"))
pairs_to_analyze <- all_deg_results %>% rename(gene_id = ENSG, gRNA_id = gRNA_group) %>% select(gene_id, gRNA_id) %>% mutate(gene_id = factor(gene_id), gRNA_id = factor(gRNA_id)) %>% arrange()
write.fst(pairs_to_analyze, paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))

gRNAgroup_pair_table <- fread("GSE120861_gene_gRNAgroup_pair_table.at_scale.txt") %>%
  mutate(grna_target = str_split_fixed(gRNAgroup, "_", 2)[,1]) %>%
  mutate(grna_target = ifelse(gRNAgroup.chr == "NTC","non-targeting", grna_target)) %>%
  select(grna_id = gRNAgroup, grna_target, chr = gRNAgroup.chr, start = gRNAgroup.start, end = gRNAgroup.stop)
gRNAgroup_pair_table[gRNAgroup_pair_table == "NTC"] = NA
gRNAgroup_pair_table = gRNAgroup_pair_table[!duplicated(gRNAgroup_pair_table), ]

#high-MOI CRISPRi setup
# 1. import data
data(highmoi_example_data); data(grna_target_data_frame_highmoi)
sceptre_object_highmoi <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)

sceptre_object_highmoi <- import_data(
  response_matrix = expression_mtx,
  grna_matrix = grna_matrix,
  grna_target_data_frame = gRNAgroup_pair_table,
  moi = "high",
  response_names = gene_names$V1
)


# 2. set analysis parameters
positive_control_pairs <- construct_positive_control_pairs(sceptre_object_highmoi)
discovery_pairs <- construct_cis_pairs(sceptre_object_highmoi,
                                       positive_control_pairs = positive_control_pairs,
                                       distance_threshold = 5e6
)
sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object_highmoi,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = "left"
)