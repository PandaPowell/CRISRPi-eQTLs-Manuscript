rm(list=ls())
library(data.table)
library(tidyverse)
options(bitmapType="cairo")

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

# We need to load the CRISPRi data intersecting the finemapped SNPs.
# Load stingseq data
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`)

# Load gasperini data, filter to targeted enhancers
resample_results = fread("CRISPR_data/resampling_results.txt") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two", is.na(target_site.start) == F)

# Load encode reanalysed crispri data, without gasperini data
encode = fread("CRISPR_data/NoGasperini_crispri_data.tsv") %>%
  mutate(target_site = paste0(chrom_GRC37,":",chromStart_GRC37,"-",chromEnd_GRC37)) %>%
  rename(gene_name = measuredGeneSymbol)

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

# Load GRC 38 genes positions
annot_file = "data/gencode.v33.annotation.gtf.gz"
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
annot.GRC38.v33 <- annot[order(annot$chr_number, annot$start),c("chr_number","start","ensembl_id", "gene_name")]

# Merge the gene names back into the data.table
encode <- left_join(encode, annot.GRC38[,c("ensembl_id","gene_name")], by = "gene_name")

# Load sting-seq target sites that intersect a finemapped GWAS variant
ss_finemap = fread("00.intersect_data/sting_seq_credible_sets/merged_stingseq_credset.txt") %>%
  rename(SS_coord = SNP_coord) %>% 
  left_join(sting_seq,"SS_coord") %>%
  mutate(gwas = as.character(gwas), data = "stingseq") %>%
  distinct(SS_coord,Gene,gwas, .keep_all=T) %>%
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`, pvalue = `Q-value (1 Mb)`,significant, ensembl_id = `Ensembl ID`, gene_name = Gene,
                tss_distance = `TSS Distance`, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+")))

ss_finemap_bcx = fread("00.intersect_data/sting_seq_credible_sets/BCX/merged_credset.txt") %>%
  rename(SS_coord = SNP_coord) %>%
  left_join(sting_seq,"SS_coord") %>%
  mutate(data = "stingseq") %>%
  distinct(SS_coord,Gene,gwas, .keep_all=T) %>%
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`,  pvalue = `Q-value (1 Mb)`, significant, ensembl_id = `Ensembl ID`, gene_name = Gene,
                tss_distance = `TSS Distance`, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+")))

gas_finemap = fread("00.intersect_data/gasperini_credible_sets/merged_gasperini_credset_closest.txt") %>%
  left_join(resample_results, "target_site") %>%
  mutate(gwas = as.character(gwas), data= "gasperini") %>%
  distinct(pair_id,gwas, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS,
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

gas_finemap_bcx = fread("00.intersect_data/gasperini_credible_sets/BCX/merged_credset.txt") %>%
  left_join(resample_results, "target_site") %>%
  distinct(pair_id,gwas, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS, data = "gasperini",
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

encode_finemap = fread("00.intersect_data/encode_credible_sets/merged_encode_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(target_site,gene_name,gwas, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS) %>%
  mutate(gwas = as.character(gwas), data = paste(`Reference#`,V32,V33,V34,V35, sep = "_"), grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target = target_site, logfc = Effectsize, pvalue = pValue, significant = Significant, ensembl_id, 
                gene_name,tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

encode_finemap_bcx = fread("00.intersect_data/encode_credible_sets/BCX/merged_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(target_site,gene_name,gwas, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS, gwas = as.character(gwas), 
         data = paste(`Reference#`,V32,V33,V34,V35, sep = "_"), grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target = target_site, logfc = Effectsize, pvalue = pValue, significant = Significant, ensembl_id, 
                gene_name,tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

cres_w_grnas = bind_rows(ss_finemap,ss_finemap_bcx,gas_finemap,gas_finemap_bcx,encode_finemap,encode_finemap_bcx) %>%
  mutate(chr = as.integer(str_extract(grna_target, "[0-9]+(?=\\:)"))) %>%
  filter(tss_distance < 1000000, tss_distance > -1000000)

## Remove intersecting gRNAs targets and check distances
sites = cres_w_grnas %>% distinct(grna_target, .keep_all = T)

merged_sites = data.frame()

# Check if gRNA positions intersect if they do then we combine all the targets under the same CRE
# We do not create new CRE IDs however some CRE IDs are lost when they are assigned a new CRE id from another study
for (i in 1:22){
  
  distance = 4001
  
  chr_sites = sites[chr == i,]
  chr_sites = chr_sites[order(grna_pos),]
  chr_sites[, difference := grna_pos - shift(grna_pos, 1, type = "lag"),]
  chr_sites$difference[is.na(chr_sites$difference) == T] = 100000
  # Group variants with overlaping regions
  chr_sites[, group := cumsum(difference > distance) + 1]
  merged_sites_chr = chr_sites %>% arrange(desc(data)) %>% group_by(group) %>% 
    mutate(grna_target_merge = first(grna_target), grna_pos_merge = first(grna_pos)) %>%
    ungroup() %>% dplyr::select(grna_target , grna_target_merge,grna_pos, grna_pos_merge)
  
  merged_sites = bind_rows(merged_sites,merged_sites_chr)
  
}

# remove overlapping CREs, 
cres_w_grnas = cres_w_grnas %>% 
  left_join(merged_sites, "grna_target") %>% 
  dplyr::select(-grna_target, -grna_pos.x, -grna_pos.y) %>% 
  rename(grna_target = grna_target_merge, grna_pos = grna_pos_merge) %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id))

sig_ccres = cres_w_grnas %>%
  filter(significant ==1) %>%
  distinct(ensembl_id,gwas,grna_target,.keep_all = T)

# your original dictionary (named character vector)
dict <- c(
  '30000' = 'White blood cell count',
  '30010' = 'Red blood cell count',
  '30020' = 'Hb conc',
  '30030' = 'Haematocrit %',
  '30040' = 'Mean corpuscular volume',
  '30050' = 'Mean corpuscular Hb',
  '30060' = 'Mean corpuscular Hb conc',
  '30070' = 'Red blood cell distr width',
  '30080' = 'Platelet count',
  '30090' = 'Platelet crit',
  '30100' = 'Mean platelet volume',
  '30110' = 'Platelet distr width',
  '30120' = 'Lymphocyte count',
  '30130' = 'Monocyte count',
  '30140' = 'Neutrophil count',
  '30150' = 'Eosinophil count',
  '30160' = 'Basophil count',
  '30180' = 'Lymphocyte %',
  '30190' = 'Monocyte %',
  '30200' = 'Neutrophil %',
  '30210' = 'Eosinophil %',
  '30220' = 'Basophil %',
  '30240' = 'Reticulocyte %',
  '30250' = 'Reticulocyte count',
  '30260' = 'Mean reticulocyte volume',
  '30270' = 'Mean sphered cell volume',
  '30280' = 'Immature reticulocyte frac',
  '30290' = 'HLS reticulocyte %',
  '30300' = 'HLS reticulocyte count'
)

# upadate for already translated names
update <- c(
  'White blood cell (leukocyte) count' = 'White blood cell count',
  'Red blood cell (erythrocyte) count'= 'Red blood cell count',
  'Haemoglobin concentration' = 'Hb conc',
  'Haematocrit percentage' = 'Haematocrit %',
  'Mean corpuscular haemoglobin' = 'Mean corpuscular Hb',
  'Mean corpuscular haemoglobin concentration' = 'Mean corpuscular Hb conc',
  'Red blood cell (erythrocyte) distribution width' = 'Red blood cell distr width',
  'Mean platelet (thrombocyte) volume' = 'Mean platelet volume',
  'Platelet distribution width' = 'Platelet distr width',
  'Lymphocyte percentage' = 'Lymphocyte %',
  'Monocyte percentage' = 'Monocyte %',
  'Neutrophil percentage'  = 'Neutrophil %',
  'Eosinophil percentage' = 'Eosinophil %',
  'Basophil percentage' = 'Basophil %',
  'Reticulocyte percentage' = 'Reticulocyte %',
  'Immature reticulocyte fraction' = 'Immature reticulocyte frac',
  'High light scatter reticulocyte percentage' = 'HLS reticulocyte %',
  'High light scatter reticulocyte count' = 'HLS reticulocyte count'
)

# add mappings for your abbreviations
abbr <- c(
  BAS  = 'Basophil count',
  EOS  = 'Eosinophil count',
  HCT  = 'Haematocrit %',
  HGB  = 'Hb conc',
  LYM  = 'Lymphocyte count',
  MCHC = 'Mean corpuscular Hb conc',
  MCH  = 'Mean corpuscular Hb',
  MCV  = 'Mean corpuscular volume',
  MON  = 'Monocyte count',
  MPV  = 'Mean platelet volume',
  NEU  = 'Neutrophil count',
  PLT  = 'Platelet count',
  RBC  = 'Red blood cell count',
  RDW  = 'Red blood cell distr width',
  WBC  = 'White blood cell count'
)

lookup <- c(dict, abbr, update)  # combined named vector

sig_ccres_labeled <- sig_ccres %>%
  mutate(
    gwas = as.character(gwas),
    gwas_label = dplyr::recode(gwas, !!!lookup, .default = gwas),
    gwas_group = case_when(
      str_detect(gwas_label, regex(
        "White blood cell|Lymphocyte|Monocyte|Neutrophil|Eosinophil|Basophil|leukocyte",
        ignore_case = TRUE
      )) ~ "White blood cells",
      str_detect(gwas_label, regex(
        "Red blood cell|erythrocyte|Haemoglobin|Hb|Haematocrit|Hematocrit|Mean corpuscular|RDW|Reticulocyte|sphered cell",
        ignore_case = TRUE
      )) ~ "Red blood cells",
      str_detect(gwas_label, regex(
        "Platelet|thrombocyte|MPV|Platelet crit|Platelet distribution width",
        ignore_case = TRUE
      )) ~ "Platelets",
      TRUE ~ "Other"
    ),
    gwas_group = factor(gwas_group, levels = c("White blood cells","Red blood cells","Platelets","Other"))
  ) %>%
  distinct(ensembl_id,gwas_label,grna_target,.keep_all = T)

label_groups <- sig_ccres_labeled %>%
  distinct(gwas_label, gwas_group)

ccre_barplot_df <- sig_ccres_labeled %>%
  count(gwas_label, name = "n") %>%
  left_join(label_groups, by = "gwas_label") %>%
  arrange(desc(n)) %>%
  mutate(gwas_label = factor(gwas_label, levels = gwas_label))

ccre_barplot = ggplot(ccre_barplot_df, aes(x = gwas_label, y = n, fill = gwas_group)) +
  geom_bar(stat = "identity", color = "black", width = 0.5) +
  labs(title = "", x = "", y = "Number of CRE-genes") +
  scale_fill_manual(
    values = c(
      "White blood cells" = "white",
      "Red blood cells"   = "red",
      "Platelets"         = "yellow",
      "Other"             = "grey70"
    )
  ) +
  theme_cowplot() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(l = 20)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16)
  ) +
  coord_flip()

# barplot for egenes
cres_w_grnas_eqtl = fread("processed_data/cres_with_grna_eqtls_interval.txt")

sig_ecres_labeled <- cres_w_grnas_eqtl %>%
  mutate(
    gwas = as.character(gwas),
    gwas_label = dplyr::recode(gwas, !!!lookup, .default = gwas),
    gwas_group = case_when(
      str_detect(gwas_label, regex(
        "White blood cell|Lymphocyte|Monocyte|Neutrophil|Eosinophil|Basophil|leukocyte",
        ignore_case = TRUE
      )) ~ "White blood cells",
      str_detect(gwas_label, regex(
        "Red blood cell|erythrocyte|Hb|Hemoglobin|Haematocrit|Hematocrit|Mean corpuscular|RDW|Reticulocyte|sphered cell",
        ignore_case = TRUE
      )) ~ "Red blood cells",
      str_detect(gwas_label, regex(
        "Platelet|thrombocyte|MPV|Platelet crit|Platelet distribution width",
        ignore_case = TRUE
      )) ~ "Platelets",
      TRUE ~ "Other"
    ),
    gwas_group = factor(gwas_group, levels = c("White blood cells","Red blood cells","Platelets","Other"))
  ) %>%
  distinct(ensembl_id,gwas_label,grna_target,.keep_all = T)

label_groups <- sig_ecres_labeled %>%
  distinct(gwas_label, gwas_group)

ecre_barplot_df <- sig_ecres_labeled %>%
  count(gwas_label, name = "n") %>%
  left_join(label_groups, by = "gwas_label") %>%
  arrange(desc(n)) %>%
  mutate(gwas_label = factor(gwas_label, levels = gwas_label))

ecre_barplot = ggplot(ecre_barplot_df, aes(x = gwas_label, y = n, fill = gwas_group)) +
  geom_bar(stat = "identity", color = "black", width = 0.5) +
  labs(title = "", x = "", y = "Number of CRE-genes") +
  scale_fill_manual(
    values = c(
      "White blood cells" = "white",
      "Red blood cells"   = "red",
      "Platelets"         = "yellow",
      "Other"             = "grey70"
    )
  ) +
  theme_cowplot() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(l = 20)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16)
  ) +
  coord_flip()

library(gridExtra)
grid.arrange(ecre_barplot, ccre_barplot, ncol = 2)

df_counts <- full_join(
  ccre_barplot_df %>% rename(n_crispri = n),
  ecre_barplot_df %>% rename(n_eqtl   = n),
  by = "gwas_label") %>%
  mutate(
    n_crispri = coalesce(n_crispri, 0L),
    n_eqtl    = coalesce(n_eqtl, 0L)
  )

df_rank <- df_counts %>%
  mutate(
    rank_CRISPRi = dplyr::dense_rank(dplyr::desc(n_crispri)),
    rank_eQTL    = dplyr::dense_rank(dplyr::desc(n_eqtl)),
    rank_diff    = rank_CRISPRi - rank_eQTL
  ) %>%
  dplyr::select(gwas_label, gwas_group.x, n_crispri, n_eqtl, rank_CRISPRi, rank_eQTL, rank_diff)

library(cowplot)

rho   <- cor(df_rank$rank_eQTL, df_rank$rank_CRISPRi, method = "spearman")
cor.test(df_rank$rank_eQTL, df_rank$rank_CRISPRi, method = "spearman")
tau   <- cor(df_rank$rank_eQTL, df_rank$rank_CRISPRi, method = "kendall")

library(ggrepel)

df_plot <- df_rank %>%
  dplyr::filter(!is.na(rank_eQTL), !is.na(rank_CRISPRi))

png("plots/interval/rank_concordance.png",
    width = 9, height = 7, units = "in", res = 300)

ggplot(df_plot, aes(x = rank_eQTL, y = rank_CRISPRi)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(fill = gwas_group.x), shape = 21, size = 3, stroke = 0.6, color = "black") +
  geom_text_repel(aes(label = as.character(gwas_label)),
                  size = 3, max.overlaps = Inf, box.padding = 0.25,
                  point.padding = 0.15, segment.color = "grey70") +
  scale_x_reverse() + scale_y_reverse() +   # rank 1 at top-right
  scale_fill_manual(
    values = c(
      "White blood cells" = "white",
      "Red blood cells"   = "red",
      "Platelets"         = "yellow",
      "Other"             = "grey70"
    ),
    name = "GWAS group"
  ) +
  labs(x = "Rank (eQTL)", y = "Rank (CRISPRi)",
       title = "",
       subtitle = paste0("Spearman \u03C1 = ", round(rho, 3),
                         "   Kendall \u03C4 = ", round(tau, 3))) +
  theme_cowplot() +
  theme(legend.position = "top")

dev.off()

# test if CRISPRi maps to RBC or WBC more with Paired signed-rank on rank differences
library(stats)
# Prepare rank diffs (you already have this)
df_rank2 <- df_rank %>%
  filter(!is.na(rank_CRISPRi), !is.na(rank_eQTL)) %>%
  mutate(group = gwas_group.x,
         drank = rank_CRISPRi - rank_eQTL)  # negative => CRISPRi ranks higher

# Helper: one-sample signed-rank vs 0 (returns HL pseudo-median + CI)
run_signed_rank <- function(x, side = c("less","greater")) {
  side <- match.arg(side)
  wilcox.test(x, mu = 0, alternative = side, exact = FALSE, conf.int = TRUE)
}

# RBC: CRISPRi ranks higher → drank < 0
rbc <- df_rank2 %>% filter(trimws(group) == "Red blood cells") %>% pull(drank)
res_rbc <- run_signed_rank(rbc, side = "less")

# WBC: eQTL ranks higher → drank > 0
wbc <- df_rank2 %>% filter(trimws(group) == "White blood cells") %>% pull(drank)
res_wbc <- run_signed_rank(wbc, side = "greater")

res_rbc
res_wbc

## Average number of traits a CRISPRi snp is associated to
trait_counts_c = sig_ccres_labeled %>%
  distinct(grna_target, gwas_label, .keep_all = T) %>%
  group_by(grna_target) %>%
  summarise(n = n()) %>%
  mutate(method = "crispr")

mean(trait_counts_c$n)

trait_counts_e = sig_ecres_labeled %>%
  distinct(grna_target, gwas_label, .keep_all = T) %>%
  group_by(grna_target) %>%
  summarise(n = n()) %>%
  mutate(method = "eqtl")

mean(trait_counts_e$n)

counts_combined = rbind(trait_counts_c, trait_counts_e)

fit_df <- counts_combined %>%
  dplyr::add_count(grna_target, name = "n_grna") %>%
  dplyr::filter(n_grna == 1) %>%
  dplyr::select(-n_grna)

fit = glm(n ~ method, family = "poisson", data = fit_df)
summary(fit)

# Wilcox test
df <- fit_df %>%
  filter(!is.na(n), !is.na(method)) %>%
  mutate(method = factor(method, levels = c("crispr","eqtl")))

# Quick sample sizes
table(df$method)

# Distribution summary (medians & IQR)
df %>%
  group_by(method) %>%
  summarise(n_targets = n(),
            median_n  = median(n),
            IQR_n     = IQR(n),
            .groups = "drop")

# Collect all fine-mapped credible sets that intersect each CRE
ccres_labeled <- cres_w_grnas %>%
  distinct(ensembl_id,gwas,grna_target,.keep_all = T) %>%
  mutate(
    gwas = as.character(gwas),
    gwas_label = dplyr::recode(gwas, !!!lookup, .default = gwas),
    gwas_group = case_when(
      str_detect(gwas_label, regex(
        "White blood cell|Lymphocyte|Monocyte|Neutrophil|Eosinophil|Basophil|leukocyte",
        ignore_case = TRUE
      )) ~ "White blood cells",
      str_detect(gwas_label, regex(
        "Red blood cell|erythrocyte|Haemoglobin|Hb|Haematocrit|Hematocrit|Mean corpuscular|RDW|Reticulocyte|sphered cell",
        ignore_case = TRUE
      )) ~ "Red blood cells",
      str_detect(gwas_label, regex(
        "Platelet|thrombocyte|MPV|Platelet crit|Platelet distribution width",
        ignore_case = TRUE
      )) ~ "Platelets",
      TRUE ~ "Other"
    ),
    gwas_group = factor(gwas_group, levels = c("White blood cells","Red blood cells","Platelets","Other"))
  ) %>%
  distinct(ensembl_id,gwas_label,grna_target,.keep_all = T)

trait_counts_c2 = ccres_labeled %>%
  distinct(grna_target, gwas_label, .keep_all = T) %>%
  group_by(grna_target) %>%
  summarise(n = n()) %>%
  mutate(method = "crispr")

cre_n_traits = bind_rows(trait_counts_c2, trait_counts_e) %>% 
  filter(n>1) %>%
  add_count(grna_target, name = "freq") %>%
  filter(freq > 1) %>%
  dplyr::select(-freq) %>%
  pivot_wider(names_from = method, values_from = n) %>%
  mutate(discordance = crispr-eqtl) %>%
  filter(discordance>=0)

wilcox.test(cre_n_traits$crispr, cre_n_traits$eqtl)

sum(cre_n_traits$discordance>0)/nrow(cre_n_traits)

# does association hold from CREs with targets for both?
cre_n_traits_shared = counts_combined %>% 
  filter(n>1) %>%
  add_count(grna_target, name = "freq") %>%
  filter(freq > 1) %>%
  dplyr::select(-freq) %>%
  pivot_wider(names_from = method, values_from = n) %>%
  mutate(discordance = crispr-eqtl) %>%
  filter(discordance>=0)

mean(cre_n_traits_shared$crispr)
cre_n_traits_shared$eqtl
wilcox.test(cre_n_traits_shared$crispr, cre_n_traits_shared$eqtl)
sum(cre_n_traits_shared$discordance>0)/nrow(cre_n_traits_shared)

# Determine snps that differ in pleiotropy
snp_trait_counts = counts_combined %>%
  tidyr::pivot_wider(
    id_cols      = grna_target,
    names_from   = method,   # e.g., "crispr", "eqtl"
    values_from  = n,
    values_fill  = NA
  ) %>% 
  mutate(diff  = crispr - eqtl) %>%
  arrange(desc(abs(diff))) %>%
  filter(is.na(diff) == F)

# Heatmap, snps vs 
cres_w_grnas = fread("processed_data/cres_with_grnas.txt")
cres_w_grnas_egene = fread("processed_data/cres_with_grna_eqtls_interval.txt")

total_cres = unique(cres_w_grnas$grna_target)
cat("Total cres tested =",length(total_cres),"\n")
grna_cres_w_cgenes = unique(cres_w_grnas$grna_target[cres_w_grnas$significant == 1])
cat("Total number CREs with cgenes =",length(grna_cres_w_cgenes),"\n")
grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
cat("Total number CREs with egenes =",length(grna_cres_w_egenes),"\n")
cat("Total number CREs with a target gene =", length(unique(c(grna_cres_w_cgenes,grna_cres_w_egenes))),"\n")
cres_no_target = length(total_cres) - length(unique(c(grna_cres_w_cgenes,grna_cres_w_egenes)))
cat("CREs without target genes =",cres_no_target,"\n")
overlapping_cres = grna_cres_w_cgenes[grna_cres_w_cgenes %in% grna_cres_w_egenes]
cat("Total number of overlapping CRES =",length(overlapping_cres),"\n")

gold_genes <- fread("processed_data/cis_gold_genes.txt")

gold_genes <- gold_genes[, .(grna = unlist(strsplit(grna_target, ","))), by = gold_gene] %>%
  mutate(grna_gene_pair = paste0(grna,"_",gold_gene))

cgenes = cres_w_grnas %>%
  filter(grna_target %in% gold_genes$grna,
         significant == 1) %>%
  distinct(ensembl_id)

# Calculate the number of GWAS traits at snp-gene pairs
trait_counts_c = sig_ccres_labeled %>%
  filter(target_gene %in% gold_genes$grna_gene_pair, ensembl_id %in% cgenes$ensembl_id) %>%
  distinct(grna_target, gwas_label, ensembl_id, .keep_all = T) %>%
  group_by(target_gene) %>%
  summarise(n = n()) %>%
  mutate(method = "crispr")

median(trait_counts_c$n)

trait_counts_e = sig_ecres_labeled %>%
  filter(target_gene %in% gold_genes$grna_gene_pair) %>%
  distinct(grna_target, gwas_label,ensembl_id, .keep_all = T) %>%
  group_by(target_gene) %>%
  summarise(n = n()) %>%
  mutate(method = "eqtl")

median(trait_counts_e$n)

counts_combined = rbind(trait_counts_c, trait_counts_e)

fit_df <- counts_combined %>%
  dplyr::add_count(target_gene, name = "n_grna") %>%
  dplyr::filter(n_grna == 1) %>%
  dplyr::select(-n_grna)

fit = glm(n ~ method, family = "poisson", data = fit_df)
summary(fit)

plot_matrix = sig_ccres_labeled %>%
  filter(ensembl_id %in% gold_genes$gold_gene) %>%
  mutate(gold = ifelse(ensembl_id %in% gold_genes$gold_gene,1,0)) %>%
  distinct(grna_target, gene_name, gwas_label, .keep_all=T) %>%
  dplyr::select(grna_target, gene_name, gwas_label, gold)

df <- plot_matrix %>%
  mutate(gold = as.integer(gold > 0),
         row_lab = gene_name) %>%
  group_by(row_lab, gwas_label) %>%
  summarise(present = as.integer(any(gold == 1)), .groups = "drop") %>%
  group_by(row_lab) %>% mutate(row_tot = sum(present)) %>% ungroup() %>%
  group_by(gwas_label) %>% mutate(col_tot = sum(present)) %>% ungroup() %>%
  mutate(
    row_lab    = fct_reorder(row_lab, row_tot, .desc = TRUE),
    gwas_label = fct_reorder(gwas_label, col_tot, .desc = TRUE),
    present_num = as.numeric(present)   # 0/1 for gradient
  )

# Venn diagram
uniq_eqtls = sig_ecres_labeled %>%
  distinct(eQTL_variant, ensembl_id, .keep_all = T)

uniq_eqtls_no_gold_gene = uniq_eqtls %>%
  filter(!target_gene %in% gold_genes$grna_gene_pair)

mean(abs(uniq_eqtls_no_gold_gene$beta), na.rm=T)

uniq_eqtls_gold_gene = uniq_eqtls %>%
  filter(target_gene %in% gold_genes$grna_gene_pair)

mean(abs(uniq_eqtls_gold_gene$beta), na.rm=T)

wilcox.test(abs(uniq_eqtls_gold_gene$beta), abs(uniq_eqtls_no_gold_gene$beta))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# Data (abs betas; drop NA/Inf)
a <- abs(uniq_eqtls_gold_gene$beta); a <- a[is.finite(a)]
b <- abs(uniq_eqtls_no_gold_gene$beta); b <- b[is.finite(b)]

# Wilcoxon test
wt <- wilcox.test(a, b)
p_txt <- paste0("Wilcoxon p = ", format(wt$p.value, scientific = FALSE, digits = 6, trim = TRUE))

# Long format for plotting
df <- bind_rows(
  data.frame(group = "Gold gene",    abs_beta = a),
  data.frame(group = "No gold gene", abs_beta = b)
)

# Position for annotation
ymax <- max(df$abs_beta, na.rm = TRUE)

p <- ggplot(df, aes(x = group, y = abs_beta, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.35) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
  annotate("text", x = 1.5, y = 0.98 * ymax, label = p_txt, size = 4.2) +
  labs(x = NULL, y = "Absolute effect size |β|") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Save as PNG
ggsave("plots/interval/abs_beta_violin.png", p, width = 6, height = 4, dpi = 300)
