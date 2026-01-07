library(fst)
library(data.table)
library(tidyverse)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

resample_results = read.fst("resampling_results.fst") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two") %>%
  dplyr::select(chr, target_site.start, target_site.stop, pair_id) %>% drop_na()

dis = resample_results$target_site.stop - resample_results$target_site.start


fwrite(resample_results, "resample_results.bed", sep = "\t", quote = F, col.names = F, row.names = F)
