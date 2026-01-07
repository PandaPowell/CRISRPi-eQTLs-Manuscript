rm(list=ls())

libraries <- c("data.table",
               "tidyverse")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load GRC 38 genes positions
annot_file = "data/gencode.v47.annotation.gtf.gz"
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

name = c("B", "Mono", "NK","DC","CD4_T","CD8_T", "other", "other_T")

cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt") %>%
  filter(!eqtl %in% c(paste0(name,"_cells"), "GTEx", "MAGE")) %>%
  mutate(unique_id = paste0(grna_target,dataset_id)) %>%
  distinct(unique_id, .keep_all=T) %>% 
  arrange(dataset_id)

variants_datasets = cres_w_grnas_egene[,c("dataset_id","study_id","eqtl_hit_hg38","eqtl","grna_target")] %>% unique()
 
# # Function should loop through dataset id
# Then another function should loop through all variants corresponding to this dataset_id
# id = "QTD000016"
extract_variants = function(id, id2){
  
  print(id2)
  print(id)
  
  temp_df = variants_datasets[variants_datasets$dataset_id == id,]
  
  base_dir = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/"
  file_name = paste0(base_dir,id2,"/",id,"/",id,".all.tsv.gz")
  eqtl_sumstats = fread(file_name)
  
  extract_variants2 = function(variant_id, grna_target, eqtl){

    locus = eqtl_sumstats[variant == variant_id,]
    locus$grna_target = grna_target
    locus$eqtl = eqtl
    
    return(locus)
  }
  
  all_loci = lapply(1:nrow(temp_df), function(x) extract_variants2(temp_df$eqtl_hit_hg38[x],
                                                                   temp_df$grna_target[x],
                                                                   temp_df$eqtl[x]))
  
  all_loci = do.call(rbind, all_loci)
  all_loci$dataset_id = id
  
  return(all_loci)
}

ids = unique(variants_datasets$dataset_id)

results = lapply(1:length(ids), function(x) extract_variants(ids[x],
          unique(variants_datasets$study_id[variants_datasets$dataset_id == ids[x]])))

results2 = do.call(rbind, results) %>% 
  rename(ensembl_id = gene_id) %>%
  left_join(annot.GRC38, "ensembl_id")

fwrite(results2,"eQTL_gene_distances/eqtl_cat_gene_distances.txt", quote = F, row.names = F)