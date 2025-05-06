# Load necessary libraries
rm(list = ls())
library(biomaRt)
library(readxl)
library(data.table)
library(tidyverse)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

# Get Ensembl v99 protein coding genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "chromosome_name", "transcription_start_site"), 
                              filters = "biotype", 
                              values = "protein_coding", 
                              mart = ensembl)
protein_coding_genes <- gene_list$ensembl_gene_id

# Load gencode data
annot_file = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/Gencode/gencode.v33lift37.GRCh38.genes.gtf"
# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$ensembl_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[4]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot.GRC37 <- annot[order(annot$chr_number, annot$start),c("chr_number","start","ensembl_id", "gene_name")]
annot.GRC37 = annot.GRC37[annot.GRC37$ensembl_id %in% protein_coding_genes,]

# Load Mendelian genes
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867420309995-mmc7.xlsx", "mendelian_genes.xlsx")
curated_genes = read_xlsx("mendelian_genes.xlsx") %>%
  rename(gene_name = Gene_Symbol_HGNC) %>%
  left_join(annot.GRC37, "gene_name") %>% drop_na() %>%
  distinct(ensembl_id, .keep_all = T) %>%
  dplyr::select(gene_name,chr_number,start,ensembl_id)

# UKBB WES blood cell trait genes
ukbb_genes = fread("data/genebass/blood_trait_burden_test_genes.csv", sep = ",")
ukbb_genes = ukbb_genes[grep("pLoF", ukbb_genes$annotation),] %>%
  distinct(gene_id, .keep_all = T) %>%
  rename(ensembl_id = gene_id) %>%
  left_join(annot.GRC37, "ensembl_id") %>% drop_na() %>%
  dplyr::select(gene_name,chr_number,start,ensembl_id)

#ok = fread("rare_genes_list.txt")
gw.rare.genes = rbind(ukbb_genes,curated_genes) %>%
  distinct(ensembl_id, .keep_all=T)

cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

cres_w_grnas = fread("cres_with_grnas.txt") 

sig.grna.targets = cres_w_grnas %>%
  distinct(grna_target, .keep_all=T) %>% 
  filter(significant == 1 | grna_target %in% cres_w_grnas_egene$grna_target)

#1:150544640

find.cis = function(chr, pos){
  cre.rare.genes = gw.rare.genes[chr_number == chr & inrange(start, pos-1e+06, pos+1e+06)]
  return(cre.rare.genes)
}

results = lapply(1:nrow(sig.grna.targets), function(x) find.cis(chr = sig.grna.targets$chr[x],
                                                            pos = sig.grna.targets$grna_pos[x]))

results2 = do.call(rbind, results) %>% 
  distinct(ensembl_id, .keep_all=T)

# load cgenes and egenes
cgenes = unique(cres_w_grnas$ensembl_id[cres_w_grnas$significant == T])
cgenes = cgenes[cgenes %in% protein_coding_genes]
egenes = unique(cres_w_grnas_egene$ensembl_id)
egenes = egenes[egenes %in% protein_coding_genes]

calculate_enrichment = function(mendel_genes){
  # Set the number of iterations
  n_iterations <- 10000
  
  # Initialize a vector to store the overlaps
  cgene_overlap_counts <- numeric(n_iterations)
  egene_overlap_counts <- numeric(n_iterations)
  
  # Calculate the overlap in each iteration
  set.seed(123) # for reproducibility
  for (i in 1:n_iterations) {
    # Randomly sample a set of genes with the same size as your GWAS list
    crandom_genes <- sample(protein_coding_genes, size = length(cgenes), replace = FALSE)
    erandom_genes = sample(protein_coding_genes, size = length(egenes), replace = FALSE)
    # Calculate the number of overlapping genes
    cgene_overlap_counts[i] <- sum(crandom_genes %in% mendel_genes)
    egene_overlap_counts[i] <- sum(erandom_genes %in% mendel_genes)
  }
  
  # Take the median of the random gene overlaps
  cgene_random_counts = median(cgene_overlap_counts)
  egene_random_counts = median(egene_overlap_counts)
  
  # Observed overlap
  cgene_observed_overlap <- sum(cgenes %in% mendel_genes)
  egene_observed_overlap <- sum(egenes %in% mendel_genes)
  
  a = cgene_observed_overlap # cgene in mendel set
  b = length(cgenes) - cgene_observed_overlap # cgenes not in mendel set
  c = cgene_random_counts
  d = length(crandom_genes) - cgene_random_counts
  
  cgene_table = matrix(c(a,b,c,d), nrow = 2)
  print(cgene_table)
  
  print(fisher.test(cgene_table, alternative = "greater"))
  
  a1 = egene_observed_overlap
  b1 = length(egenes) - egene_observed_overlap
  c1 = egene_random_counts
  d1 = length(erandom_genes) - egene_random_counts
  
  egene_table = matrix(c(a1,b1,c1,d1), nrow = 2)
  print(egene_table)
  
  print(fisher.test(egene_table, alternative = 'greater'))
  
  compare_table = matrix(c(a,b,a1,b1), nrow = 2)
  print(compare_table)
  
  print(fisher.test(compare_table))
  
  print(cgenes[cgenes %in% mendel_genes])
  print(egenes[egenes %in% mendel_genes])
}

calculate_enrichment(results2$ensembl_id)