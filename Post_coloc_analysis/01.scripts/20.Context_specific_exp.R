rm(list=ls())
library(data.table)
library(tidyverse)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load stingseq data
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`)

# Load gasperini data, filter to targeted enhancers
resample_results = fread("CRISPR_data/resampling_results.txt") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two", is.na(target_site.start) == F)

# Load all genes expressed in whole blood
gtex.power <- fread("power_results/power.results.beta0.5GTExV8_blood.csv") %>%
  mutate(snp_gene = paste0(snp, ":", genes)) %>%
  distinct(snp_gene, .keep_all = T)
gtex.power$eQTL.power[is.na(gtex.power$eQTL.power)] = 0
# Remove genes with less than 6 counts, 6/755 = 0.00794702
gtex.power = gtex.power[count.mean > 1,]

eqtl.cat = fread("power_results/eqtl.power.combined.txt") %>% filter(count.mean > 1)

# Load cres
cres_w_grnas = fread("cres_with_grnas.txt") %>%
  filter(significant == 1)

cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

total_cres = unique(cres_w_grnas$grna_target)
cat("Total cres tested =",length(total_cres),"\n")

grna_cres_w_cgenes = unique(cres_w_grnas$grna_target[cres_w_grnas$significant == 1])
cat("Total number CREs with cgenes =",length(grna_cres_w_cgenes),"\n")

grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
cat("Total number CREs with egenes =",length(grna_cres_w_egenes),"\n")

cat("Total number CREs with a target gene =", length(unique(c(grna_cres_w_cgenes,grna_cres_w_egenes))),"\n")

cres_no_target = length(total_cres) - length(unique(c(grna_cres_w_cgenes,grna_cres_w_egenes)))
cat("CREs without target genes =",cres_no_target,"\n")

overlapping_cres = grna_cres_w_cgenes[grna_cres_w_cgenes %in% grna_cres_w_egenes]
cat("Total number of overlapping CRES =",length(overlapping_cres),"\n")

cres_w_grnas = cres_w_grnas %>%
  filter(grna_target %in% overlapping_cres)

cres_w_grnas_egene = cres_w_grnas_egene %>%
  filter(grna_target %in% overlapping_cres)

# What proportion of significant egenes are not expressed in k562
total_egenes = length(unique(cres_w_grnas_egene$ensembl_id))
overlap_egenes = length(unique(cres_w_grnas_egene$ensembl_id[cres_w_grnas_egene$ensembl_id %in% c(sting_seq$`Ensembl ID`, resample_results$gene_id)]))

overlap_egenes/total_egenes

# What proportion of significant cgenes are not expressed in GTEx
total_cgenes = length(unique(cres_w_grnas$ensembl_id))
overlap_cgenes = length(unique(cres_w_grnas$ensembl_id[cres_w_grnas$ensembl_id %in% eqtl.cat$ensembl_id]))

overlap_cgenes/total_cgenes

unique(cres_w_grnas$ensembl_id[!cres_w_grnas$ensembl_id %in% eqtl.cat$ensembl_id])
unique(cres_w_grnas$gene_name[!cres_w_grnas$gene_name %in% eqtl.cat$genes])
