rm(list=ls())
library(data.table)
library(tidyverse)
options(bitmapType="cairo")

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/")

# Load GRC 38 genes positions
annot_file = "Post_coloc_analysis/data/gencode.v47.annotation.gtf.gz"
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$ensembl_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[3]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot.GRC38 <- annot[order(annot$chr_number, annot$start),c("chr_number","start","ensembl_id", "gene_name")]

gold_genes <- fread("Post_coloc_analysis/processed_data/cis_gold_genes.txt")

data.sets = fread("eQTL_catalogue/eQTL_catalogue_datasets.txt") %>%
  filter(tissue_label != "fibroblast")

# # Function should loop through dataset id
# Then another function should loop through all variants corresponding to this dataset_id
# id = "QTD000016"
extract_variants = function(did, sid, sample_group){
  
  temp_df = data.sets[data.sets$dataset_id == did,]
  
  base_dir = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/"
  file_name = paste0(base_dir,sid,"/",did,"/",did,".all.tsv.gz")
  eqtl_sumstats = fread(file_name)
  
  genes_tpm = eqtl_sumstats[molecular_trait_id %in% gold_genes$gold_gene, c("molecular_trait_id","median_tpm")] %>%
    distinct() %>% 
    mutate(cell_type = sample_group)
  
  return(genes_tpm)
}

ids = unique(data.sets$dataset_id)

results = lapply(1:length(ids), function(x) extract_variants(ids[x],
                                                             unique(data.sets$study_id[data.sets$dataset_id == ids[x]]),
                                                             unique(data.sets$sample_group[data.sets$dataset_id == ids[x]])))

results2 = do.call(rbind, results) %>% 
  left_join(annot.GRC38, by=c("molecular_trait_id"="ensembl_id"))

fwrite(results2,"Post_coloc_analysis/processed_data/gold_gene_TPM.txt", quote = F, row.names = F)
