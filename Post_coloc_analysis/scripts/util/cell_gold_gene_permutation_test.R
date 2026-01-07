rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(readxl)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# eGene enrichment by cell-type
# Load Gold Standard Genes
gold_genes <- fread("processed_data/cis_gold_genes.txt")

# load protein coding genes
protein_coding_genes = read_lines("processed_data/ensembl_protein_coding_genes.txt")

# Load CRISPR eQTL and GRNA target data
cres_w_grnas_egene <- fread("processed_data/cres_with_grna_eqtls_interval.txt")
cres_w_grnas <- fread("processed_data/cres_with_grnas.txt")

# Split grna_target into separate rows for better filtering
gold_genes <- gold_genes[, .(grna = unlist(strsplit(grna_target, ","))), by = gold_gene]

cgenes = cres_w_grnas %>%
  filter(grna_target %in% gold_genes$grna,
         significant == 1, 
         ensembl_id %in% protein_coding_genes) %>%
  distinct(ensembl_id)

cgenes = unlist(cgenes, use.names = F)

egenes = cres_w_grnas_egene %>%
  filter(grna_target %in% gold_genes$grna, 
         ensembl_id %in% protein_coding_genes) %>%
  arrange(desc(PP.H4.abf)) %>%
  distinct(unique_id, .keep_all = T) %>%# remove duplicate eQTL colocs at different GWAS traits
  mutate(gold = ifelse(ensembl_id %in% gold_genes$gold_gene, 1, 0)) %>%
  mutate(eQTL_group = case_when(
    eqtl %in% c("B-cell_naive", "B_cells") ~ "B cells",
    eqtl %in% c("CD4_T-cell_anti-CD3-CD28", "CD4_T-cell_naive", "CD4_T_cells") ~ "CD4 T cells",
    eqtl %in% c("CD8_T-cell_anti-CD3-CD28", "CD8_T-cell_naive", "CD8_T_cells") ~ "CD8 T cells",
    eqtl %in% c("NK-cell_naive", "NK_cells") ~ "NK cells",
    eqtl %in% c("DC_cells") ~ "Dendritic cells",
    eqtl %in% c("T-cell", "Tfh_memory", "Th1-17_memory", "Th17_memory",
                "Th1_memory", "Th2_memory", "Treg_memory", "Treg_naive", "other_T_cells") ~ "Other T cells",
    eqtl %in% c("monocyte_LPS", "monocyte_Pam3CSK4", "monocyte_R848", "monocyte_naive", 
                "monocyte_CD16_naive", "monocyte","monocyte_IAV","Mono_cells") ~ "Monocytes",
    eqtl %in% c("macrophage_IFNg+Salmonella", "macrophage_Listeria",
                "macrophage_Salmonella", "macrophage_naive") ~ "Macrophages",
    eqtl %in% c("neutrophil") ~ "Neutrophils",
    eqtl %in% c("LCL", "LCL_naive", "LCL_statin", "MAGE") ~ "LCLs",
    eqtl %in% c("GTEx", "blood", "Interval") ~ "Blood",
    eqtl %in% c("other_cells") ~ "other",
    TRUE ~ "Unmapped"
  ))

# Prepare a binary membership matrix: ensembl_id x eQTL_group
upset_data <- egenes %>%
  filter(gold == 1) %>%
  distinct(ensembl_id, eQTL_group) %>%     # one row per gene/group
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = eQTL_group, values_from = present, values_fill = 0)

# Universe-aware permutation test for "exactly-one" rows
permtest_unique1_universe <- function(upset_data, universes, B = 10000, seed = 1,
                                      alternative = c("two.sided","greater","less")) {
  alternative <- match.arg(alternative)
  stopifnot(ncol(upset_data) >= 2)
  ids <- upset_data[[1]]
  X <- as.matrix(upset_data[, -1, drop = FALSE])
  storage.mode(X) <- "integer"     # ensure 0/1 integers, not doubles
  K <- ncol(X); N <- nrow(X)
  colnames(X) <- colnames(upset_data)[-1]
  
  # Build eligibility mask E (N x K): TRUE if gene i is in the universe for column j
  E <- matrix(FALSE, nrow = N, ncol = K, dimnames = list(ids, colnames(X)))
  for (j in seq_len(K)) {
    col_j <- colnames(X)[j]
    if (is.null(universes[[col_j]]))
      stop(sprintf("universes[['%s']] is missing", col_j))
    E[, j] <- ids %in% universes[[col_j]]
  }
  
  # Sanity: if any observed 1 falls outside eligibility, force to 0 (conservative) and warn
  outside <- (X == 1) & (!E)
  if (any(outside, na.rm = TRUE)) {
    warning("Some observed 1s were outside the provided universe for that column; setting to 0 for testing.")
    X[outside] <- 0L
  }
  
  # Observed statistic: number of rows with exactly one 1 across columns
  U_obs <- sum(rowSums(X) == 1)
  
  set.seed(seed)
  U_null <- numeric(B)
  
  # Precompute eligible indices and per-column totals to preserve in permutations
  eligible_idx <- lapply(seq_len(K), function(j) which(E[, j]))
  n1 <- vapply(
    seq_len(K),
    function(j) as.integer(sum(X[eligible_idx[[j]], j], na.rm = TRUE)),
    integer(1)
  )
  
  for (b in seq_len(B)) {
    Xp <- matrix(0L, nrow = N, ncol = K)
    for (j in seq_len(K)) {
      ei <- eligible_idx[[j]]
      if (length(ei) == 0L || n1[j] == 0L) next
      if (n1[j] > length(ei)) {
        stop(sprintf("Column '%s': requested %d ones but only %d eligible genes.",
                     colnames(X)[j], n1[j], length(ei)))
      }
      put1 <- sample(ei, n1[j], replace = FALSE)
      Xp[put1, j] <- 1L
    }
    U_null[b] <- sum(rowSums(Xp) == 1)
  }
  
  mu <- mean(U_null); sdv <- stats::sd(U_null)
  
  p_val <- switch(alternative,
                  "greater"   = (1 + sum(U_null >= U_obs)) / (B + 1),
                  "less"      = (1 + sum(U_null <= U_obs)) / (B + 1),
                  "two.sided" = (1 + sum(abs(U_null - mu) >= abs(U_obs - mu))) / (B + 1)
  )
  
  list(
    U_obs       = U_obs,
    mean_null   = mu,
    sd_null     = sdv,
    ci95_null   = as.numeric(quantile(U_null, c(0.025, 0.975))),
    p_value     = p_val,
    alternative = alternative,
    B           = B
  )
}

# upset_data: data.frame with first column 'ensembl_id' and 0/1 columns after
# Build universes as a named list of eligible gene IDs per column:
universes <- list(
  `Blood` = unique(egenes$ensembl_id[egenes$eQTL_group == "Blood"]),   # character vector of Ensembl IDs
  `Monocytes`        = unique(egenes$ensembl_id[egenes$eQTL_group == "Monocytes"]),
  `Neutrophils`         = unique(egenes$ensembl_id[egenes$eQTL_group == "Neutrophils"]),
  `LCLs` = unique(egenes$ensembl_id[egenes$eQTL_group == "LCLs"]),
  `Macrophages` = unique(egenes$ensembl_id[egenes$eQTL_group == "Macrophages"]),
  `Other T cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "Other T cells"]),
  `B cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "B cells"]),
  `NK cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "NK cells"]),
  `CD4 T cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "CD4 T cells"]),
  `CD8 T cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "CD8 T cells"]),
  `Dendritic cells` = unique(egenes$ensembl_id[egenes$eQTL_group == "Dendritic cells"]),
  `other` = unique(egenes$ensembl_id[egenes$eQTL_group == "other"])
)

res <- permtest_unique1_universe(upset_data, universes, B = 10000, seed = 1, alternative = "two.sided")
str(res)
