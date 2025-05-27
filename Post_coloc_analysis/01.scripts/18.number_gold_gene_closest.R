rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(scales)
library(tidyverse)
library(cowplot)
library(gridExtra)

setDTthreads(8)
setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load cres
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

# Load gold standard causal gene list
gold_genes = fread("cis_gold_genes.txt", header = T)
colnames(gold_genes)[2] = "rare_genes"

# filter eQTLs and crispri overlapping gold standard gene list
cgenes = cres_w_grnas %>% 
  filter(ensembl_id %in% gold_genes$rare_genes, significant == T) %>% 
  distinct(grna_target, .keep_all = T)

egenes = cres_w_grnas_egene %>% 
  filter(ensembl_id %in% gold_genes$rare_genes) %>% 
  distinct(grna_target, .keep_all = T)

# Load eqtl distances
eqtl.cat.distances = fread("eQTL_gene_distances/eqtl_cat_gene_distances.txt") %>%
  mutate(unique_id = paste0(grna_target,ensembl_id, dataset_id),
         grna_target_dataset = paste0(grna_target,"_",dataset_id), 
         tss_distance = position-start) %>%
  arrange(tss_distance) %>% 
  distinct(unique_id, .keep_all = T) %>% # remove duplicates 
  group_by(grna_target_dataset) %>% 
  distinct(unique_id, .keep_all =T) %>% 
  mutate(dis_rank = rank(abs(tss_distance))) %>%
  ungroup() %>%
  filter(ensembl_id %in% gold_genes$rare_genes) %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  filter(dis_rank == 1, target_gene %in% cres_w_grnas_egene$target_gene) %>%
  distinct(gene_name, .keep_all=T)

sc.eqtl.distances = fread("eQTL_gene_distances/Onek1k_gene_distance.txt") %>%
  filter(ensembl_id %in% gold_genes$rare_genes) %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  filter(dis_rank == 1, target_gene %in% cres_w_grnas_egene$target_gene) %>%
  distinct(ensembl_id, .keep_all=T)

gtex.distance = fread("eQTL_gene_distances/GTEx_gene_distance.txt") %>%
  filter(ensembl_id %in% gold_genes$rare_genes) %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  filter(dis_rank == 1, target_gene %in% cres_w_grnas_egene$target_gene) %>%
  distinct(ensembl_id, .keep_all=T)

mage.distance = fread("eQTL_gene_distances/MAGE_gene_distance.txt") %>%
  filter(ensembl_id %in% gold_genes$rare_genes) %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_id)) %>%
  filter(dis_rank == 1, target_gene %in% cres_w_grnas_egene$target_gene) %>%
  distinct(ensembl_id, .keep_all=T)

length(unique(c(eqtl.cat.distances$gene_name, sc.eqtl.distances$gene_name, gtex.distance$gene_name, mage.distance$geneSymbol)))

closest_cgene = cres_w_grnas %>% distinct(target_gene, .keep_all=T) %>% # Remove duplicate target-genes for different eqtl/gwas
  group_by(grna_target) %>% mutate(dis_rank = rank(abs(tss_distance))) %>% 
  ungroup() %>%
  filter(significant ==1, ensembl_id %in% gold_genes$rare_genes) %>%
  arrange(dis_rank) %>%
  distinct(ensembl_id, .keep_all=T)

sum(closest_cgene$dis_rank == 1)
