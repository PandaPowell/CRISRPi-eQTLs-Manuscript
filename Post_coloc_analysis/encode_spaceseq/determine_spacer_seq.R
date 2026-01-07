library(data.table)
library(fst)
library(tidyverse)
setDTthreads(8)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/encode_spaceseq/")

# Load blast output
blast_out = fread("blast_results.out") %>% select(V1,V2,V9,V10) %>%
  rename(Target_Site_gRNA = V1, acc_number = V2, start = V9, stop = V10) %>%
  mutate(start.spacer = ifelse(start < stop, start,stop), 
         stop.spacer = ifelse(stop > start, stop,start)) %>%
  separate(acc_number, "\\.", into=c("tmp1", "tmp2")) %>%
  mutate(chr.spacer = paste0("chr",gsub("NC_0+","", tmp1))) %>%
  select(-tmp1,-tmp2)

# Rename X and Y chrm
blast_out$chr.spacer[blast_out$chr.spacer == "chr23"] = "chrX"
blast_out$chr.spacer[blast_out$chr.spacer == "chr24"] = "chrY"

# Load encode data file
supp_table = fread("../CRISPR_data/ENCODE_harmonised_CRISPRi_data_GRC38.tsv") %>% rename(Target_Site_gRNA = name)

# Join the table and check if blast positions falls within enhancer positions
joint_table = supp_table %>%
  left_join(blast_out,"Target_Site_gRNA") %>%
  filter(chrom == chr.spacer, 
         (chromEnd-start.spacer) > 0,
         (stop.spacer-chromStart) > 0) %>%
  mutate(distance = stop.spacer-start.spacer)

# Since distance between grna targets is not too large we can just take the position of the first
target_site_positions = joint_table %>%
  select(Target_Site_gRNA, start.spacer, stop.spacer)

# Load crispri data and join on spacer positions
spacer_results = supp_table %>%
  left_join(target_site_positions, "Target_Site_gRNA") %>%
  mutate(start.spacer = ifelse(is.na(start.spacer) == T, chromStart, start.spacer),
         stop.spacer = ifelse(is.na(stop.spacer) == T, chromEnd, stop.spacer)) %>%
  dplyr::select(chrom, start.spacer, stop.spacer, PerturbationTargetID) %>%
  distinct(PerturbationTargetID, .keep_all = T)

fwrite(spacer_results, "encode_gRNA_target_sites_GRCh38.bed", sep = "\t", quote = F, col.names = F, row.names = F)

