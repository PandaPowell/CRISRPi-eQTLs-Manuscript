rm(list=ls())

libraries <- c("data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "ggplot2",
               "fst",
               "arrow")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# cres with eqtls
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt") %>% filter(eqtl == "MAGE") %>%
  distinct(grna_target, .keep_all=T) %>% arrange(eQTL_variant)

gene_distance = function(chr, eqtl_variant, gwas_variant, sum_stats){
  
  # Calculate distance of gene in comparison to GWAS variant as this was intersected with grna range
  # Split the string by ":" & # Rearrange the parts
  parts <- strsplit(eqtl_variant, ":")[[1]]
  eqtl_variant_alt <- paste(parts[1], parts[2], parts[4], parts[3], sep=":")
  gwas_pos = as.numeric(parts[2])
  
  # filter sum stats to coloc snps calculate distance from start site
  gene_range = sum_stats[snp_hg19 %in% c(eqtl_variant,eqtl_variant_alt),] %>%
    rename(eQTL_variant_hg38 = snp_hg38, eQTL_variant_hg19 = snp_hg19) %>%
    distinct(ensemblID, .keep_all =T) %>%
    mutate(dis_rank = rank(abs(tss_distance)), GWAS_variant = gwas_variant)
  
  if(nrow(gene_range) > 0){
    return(gene_range)
  } else{
    cat(eQTL_variant,"not found\n")
  }
  
}

# Preload sumstats only once for the given chromosome
sum_stats_list = list()
  
# Helper function to handle gene distance calculations across coloc data
process_coloc_data = function(coloc_data, sum_stats_list) {
  
    result_list = list()
    
    for (i in 1:nrow(coloc_data)) {
      
      
      eqtl_variant = gsub("chr","",gsub("_",":",coloc_data$eQTL_variant[i]))
      chr = strsplit(eqtl_variant, ":")[[1]][1]
      grna_target = coloc_data$grna_target[i]
      gwas_variant = gsub("_",":",coloc_data$finemap_snp_intersect_grna[i])
      
      # Load sumstats only if not already loaded
      if (!chr %in% names(sum_stats_list)) {
        eqtl_dir = paste0("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/MAGE/eQTL_FastQTL_results_chr",chr,".hg19.txt")
        sum_stats_list[[chr]] = fread(eqtl_dir) %>%
          mutate(snp_hg19 = paste(gsub("chr", "", chr), upper_hg19, variantRef, variantAlt, sep = ":"))
      }
      
      # Compute gene distance using pre-loaded sum_stats
      sum_stats = sum_stats_list[[chr]]
      gene_res = gene_distance(chr, eqtl_variant, gwas_variant, sum_stats)
      gene_res$grna_target = grna_target
      
      if (!is.null(gene_res)) {
        result_list[[length(result_list) + 1]] = gene_res
      }
    }
    
    return(do.call(rbind, result_list))
}
  
all_res = process_coloc_data(cres_w_grnas_egene, sum_stats_list)
  
# Combine all results
all_res2 = all_res %>% unique() %>%
  mutate(ensembl_id = str_split_fixed(ensemblID, "\\.",2)[,1])

fwrite(all_res2,"eQTL_gene_distances/MAGE_gene_distance.txt",quote = F, row.names = F)
