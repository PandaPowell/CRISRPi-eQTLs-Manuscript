rm(list=ls())

libraries <- c("data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "ggplot2",
               "biomaRt",
               "cowplot",
               "fst",
               "arrow")

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

all_cells = data.frame()

for (cell_name in name){
  
  print(cell_name)
  
  # cres with eqtls
  cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt") %>% filter(eqtl %in% paste0(cell_name,"_cells")) %>%
    distinct(grna_target, .keep_all=T) %>% arrange(eQTL_variant)
  
  eqtl_base_dir = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/lvl1Res/"
  sum_stats = paste0(eqtl_base_dir,cell_name,"_mean_mx/Cis_GW_allpairs/OneK1K_",cell_name,".cis_qtl_pairs.")
  
  gene_distance = function(GWAS_variant, eQTL_variant, eqtl.dir, grna_target){
    
    parts <- strsplit(eQTL_variant, ":")[[1]]
    chr = paste0("chr",parts[1])
    GWAS_pos = as.numeric(parts[2])
    eQTL_variant_alt = paste(parts[1],parts[2],parts[4],parts[3],sep=":")
    
    # load eqtl sum stats
    eqtl_base_dir = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/lvl1Res/"
    eqtl.sum.stats = as.data.table(read_parquet(paste0(eqtl.dir,chr,".parquet")))
    
    # filter sum stats to coloc snps calculate distance from start site
    gene_range = eqtl.sum.stats[variant_id %in% eQTL_variant | variant_id %in% eQTL_variant_alt,] %>%
      mutate(distance = start_distance) %>%
      mutate(dis_rank = rank(abs(distance))) %>%
      mutate(grna_target = grna_target, GWAS_variant = variant_id)
    
    if (nrow(gene_range) > 0) {
      return(gene_range)
    } else {
      # Create a row of NAs for all columns in avg_exp
      na_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(sum_stats)))
      colnames(na_row) <- colnames(sum_stats)
      # Add additional columns from mutate
      na_row$distance <- NA
      na_row$dis_rank <- NA
      na_row$grna_target <- grna_target
      na_row$GWAS_variant <- GWAS_variant
      return(na_row)
    }
    
  }
  
  ss_res = lapply(1:nrow(cres_w_grnas_egene), 
                  function(x) gene_distance(cres_w_grnas_egene$finemap_snp_intersect_grna[x],
                                            cres_w_grnas_egene$eQTL_variant[x],
                                            sum_stats,
                                            cres_w_grnas_egene$grna_target[x]))
  ss_res2 = do.call(rbind, ss_res)
  ss_res2$cell_type = paste0(cell_name,"_cells")
  
  all_cells = bind_rows(all_cells, ss_res2)
  
}

all_cells = all_cells %>% rename(gene_name = phenotype_id) %>%
  mutate(unique_id = paste0(all_cells$grna_target, all_cells$gene_name , all_cells$cell_type)) %>%
  left_join(annot.GRC38, "gene_name")

fwrite(all_cells,"eQTL_gene_distances/Onek1k_gene_distance.txt",quote = F, row.names = F)
