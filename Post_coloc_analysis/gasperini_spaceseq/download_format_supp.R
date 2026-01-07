library(readxl)
library(dplyr)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/gasperini_spaceseq")

supp = read_xlsx("mmc2.xlsx", sheet = 2)

# Filter to # 5779 enhancer sites
enhancer_site = supp[!supp$Category %in% c("NTC","TSS","Positive_control_to_globin_locus"),] %>%
  group_by(Target_Site) %>% mutate(gRNA = seq(1,n())) %>%
  ungroup() %>% mutate(Target_Site_gRNA = paste0(Target_Site,".",gRNA))

# Write file
write.table(enhancer_site, "gRNA_targetsite_spacers.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
