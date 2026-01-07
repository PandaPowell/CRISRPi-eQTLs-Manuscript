library(data.table)
library(fst)
library(tidyverse)
setDTthreads(8)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/gasperini_spaceseq")

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

# Load gasperini supp table 2
supp_table = fread("gRNA_targetsite_spacers.txt")

# Join the table and check if blast positions falls within enhancer positions
joint_table = supp_table %>% left_join(blast_out,"Target_Site_gRNA") %>%
  filter(chr.candidate_enhancer == chr.spacer, 
         (stop.candidate_enhancer-start.spacer) > 0,
         (stop.spacer-start.candidate_enhancer) > 0) %>%
  mutate(distance = stop.spacer-start.spacer)

if(length(unique(joint_table$Target_Site)) < 5779){
  stop("Target sites missing")
}

# quantify difference in gRNA spacer positions that target the same site
space_pos <- joint_table %>%
  group_by(Target_Site) %>%
  arrange(Target_Site, gRNA) %>% # Ensure the rows are sorted within each group
  mutate(diff.start = start.spacer - lag(start.spacer, default = first(start.spacer)),
         diff.stop = stop.spacer - lag(stop.spacer, default = first(stop.spacer)))

median(abs(space_pos$diff.start))
max(abs(space_pos$diff.start))

# Since distance between grna targets is not too large we can just take the position of the first
target_site_positions = joint_table[gRNA ==1 ,] %>% distinct(Target_Site, .keep_all=T) %>%
  select(target_site = Target_Site, start.spacer, stop.spacer)

# Load SCEPTRE gasperini data, filter to enhancer sites
resample_results = read.fst("../resampling_results.fst") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two") %>% 
  left_join(target_site_positions, "target_site") %>%
  dplyr::select(chr, start.spacer, stop.spacer, target_site) %>% drop_na() %>% unique()

fwrite(resample_results, "resample_gRNA_target_sites.bed", sep = "\t", quote = F, col.names = F, row.names = F)

