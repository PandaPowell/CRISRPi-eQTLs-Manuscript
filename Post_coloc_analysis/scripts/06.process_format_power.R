rm(list=ls())
library(tidyverse)
library(data.table)
library(fst)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load stingseq data
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`)

# Load gasperini data, filter to targeted enhancers
resample_results = fread("CRISPR_data/resampling_results.txt") %>%
  unique() %>% 
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

# Load sting-seq power
stingseq.p1 = fread("power_results/crispri_power/morris_large_power_per_pair.csv")
stingseq.p2 = fread("power_results/crispri_power/morris_small_power_per_pair.csv")
stingseq.p = rbind(stingseq.p1,stingseq.p2) %>%
  mutate(cre_gene = paste0(cre_pert,"_",gene)) %>%
  arrange(desc(fraction_sig)) %>%
  distinct()

# Load sting-seq target sites that intersect a finemapped GWAS variant
ss_finemap = fread("00.intersect_data/sting_seq_credible_sets/merged_stingseq_credset.txt") %>%
  rename(SS_coord = SNP_coord) %>% 
  left_join(sting_seq,"SS_coord") %>%
  mutate(pair_id = paste0(SS_coord,Gene)) %>% 
  mutate(gwas = as.character(gwas), data = "stingseq") %>%
  distinct(pair_id, .keep_all=T) %>% 
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`, pvalue = `Q-value (1 Mb)`, significant, 
                ensembl_id = `Ensembl ID`, gene_name = Gene, tss_distance = `TSS Distance`, finemap_snp_intersect_grna,
                sentinel_snp = finemap_snp, gwas, data, gRNAs) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+"))) %>%
  mutate(cre_pert = str_split_fixed(gRNAs, "__", 3)[,1]) %>%
  mutate(cre_pert = gsub("-[0-9]$", "", cre_pert), 
         cre_gene = paste0(cre_pert,"_",ensembl_id)) %>%
  left_join(stingseq.p, "cre_gene") %>%
  dplyr::select(-gRNAs) %>%
  rename(power = fraction_sig)

ss_finemap_bcx = fread("00.intersect_data/sting_seq_credible_sets/BCX/merged_credset.txt") %>%
  rename(SS_coord = SNP_coord) %>%
  left_join(sting_seq,"SS_coord") %>%
  mutate(pair_id = paste0(SS_coord,Gene), data = "stingseq") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`,  pvalue = `Q-value (1 Mb)`, significant, ensembl_id = `Ensembl ID`, gene_name = Gene,
                tss_distance = `TSS Distance`, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, gRNAs) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+"))) %>%
  mutate(cre_pert = str_split_fixed(gRNAs, "__", 3)[,1]) %>%
  mutate(cre_pert = gsub("-[0-9]$", "", cre_pert), 
         cre_gene = paste0(cre_pert,"_",ensembl_id)) %>%
  left_join(stingseq.p, "cre_gene") %>% dplyr::select(-gRNAs) %>% rename(power = fraction_sig)

# Load gasperini power calculation
gasperini.p = fread("power_results/crispri_power/gasperini_power_per_pair.csv") %>%
  mutate(cre_gene = paste0(cre_pert,"_",gene))

# combine with association results
gas_finemap = fread("00.intersect_data/gasperini_credible_sets/merged_gasperini_credset_closest.txt") %>%
  left_join(resample_results, "target_site") %>%
  mutate(gwas = as.character(gwas), data= "gasperini") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS,
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos, cre_pert = grna_group) %>%
  mutate(cre_gene = paste0(cre_pert,"_",ensembl_id)) %>%
  left_join(gasperini.p, "cre_gene") %>% rename(power = fraction_sig)

gas_finemap_bcx = fread("00.intersect_data/gasperini_credible_sets/BCX/merged_credset.txt") %>%
  left_join(resample_results, "target_site") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS, data = "gasperini",
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos, cre_pert = grna_group) %>%
  mutate(cre_gene = paste0(cre_pert,"_",ensembl_id)) %>%
  left_join(gasperini.p, "cre_gene") %>% rename(power = fraction_sig)

encode_finemap = fread("00.intersect_data/encode_credible_sets/merged_encode_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(name, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS) %>%
  mutate(gwas = as.character(gwas), data = "encode", grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target = target_site, logfc = Effectsize, pvalue = pValue, significant = Significant, ensembl_id, 
                gene_name,tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos, power = PowerAtEffectSize25)

encode_finemap_bcx = fread("00.intersect_data/encode_credible_sets/BCX/merged_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(name, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS, gwas = as.character(gwas), 
         data = "encode", grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target = target_site, logfc = Effectsize, pvalue = pValue, significant = Significant, ensembl_id, 
                gene_name,tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos, power = PowerAtEffectSize25)

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
  mutate(target_gene = paste0(grna_target,"_",gene_name),
         unique_id = paste0(grna_target,"_",gene_name,"_", effect_size)) %>%
  distinct(unique_id, .keep_all=T)

fwrite(cres_w_grnas,"cres_with_grnas_power.txt",quote = F, row.names = F)

cres_w_grnas = fread("processed_data/cres_with_grnas.txt")

# Load CRISPRi power
crispr.power = fread("processed_data/cres_with_grnas_power.txt")

target_cgenes = unique(cres_w_grnas$target_gene[cres_w_grnas$significant == 1])

# calculate proportion of CRE-gene pairs with power > 80% at different effect sizes for significant genes
crispr.power.sig = crispr.power[significant == 1,]

all_crispr_prop_sig = data.frame()

for (i in unique(crispr.power$effect_size[is.na(crispr.power$effect_size)==F])){
  
  print(i)
  
  tmp_df = crispr.power.sig[effect_size == i,]
  
  crispr_genes <- data.frame(
    eqtl = "CRISPRi",
    total_sig_cre_genes = length(target_cgenes),
    cre_gene_80 = nrow(tmp_df[power>=0.80]),
    "prop_power>80" = nrow(tmp_df[power>=0.80])/nrow(tmp_df),
    "mean" = mean(tmp_df$power),
    effect_size = i)
  
  all_crispr_prop_sig = bind_rows(all_crispr_prop_sig,crispr_genes)
  
}

all_crispr_prop_sig

# Compute the proportion safely (no dplyr needed)
all_crispr_prop_sig$proportion_of_sig_cres <- 
  all_crispr_prop_sig$cre_gene_80/all_crispr_prop_sig$total_sig_cre_genes

df <- all_crispr_prop_sig[is.finite(all_crispr_prop_sig$proportion_of_sig_cres), ]

# If effect_size is numeric but you want discrete bars, order as a factor
if (is.numeric(df$effect_size)) {
  df$effect_size <- factor(df$effect_size, levels = sort(unique(df$effect_size)))
}

p  = ggplot(df, aes(x = effect_size, y = proportion_of_sig_cres)) +
  geom_col(width = 0.65, fill = "grey70", color = "grey25") +
  ## 80% reference line
  geom_hline(yintercept = 0.8, linetype = "dotted", linewidth = 0.8, color = "grey20") +
  ## Label the 80% line
  annotate("text", x = 1, y = 0.82, label = "", hjust = 0, vjust = -0.2, size = 3.5) +
  ## Optional: add % labels on bars (no 'scales' pkg)
  geom_text(aes(label = paste0(round(proportion_of_sig_cres * 100), "%")),
            vjust = -0.4, size = 3.6) +
  ## Keep within [0,1] and leave headroom for labels
  coord_cartesian(ylim = c(0, 1.05), clip = "off") +
  ## Ticks at 0, .2, ..., 1.0 without using scales pkg
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(x = "Expression decrease", y = "Proportion of significant CREs") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.9),
    plot.margin = margin(10, 20, 10, 10)
  )

# save as SVG (vector)
ggplot2::ggsave(
  filename = "plots/interval/proportion_by_effect_size.svg",
  plot = p,
  device = "svg",
  width = 6, height = 4, units = "in", dpi = 300
)

# calculate proportion of CRE-gene pairs with power > 80% at different effect sizes
all_crispr_prop = data.frame()

for (i in unique(crispr.power$effect_size[is.na(crispr.power$effect_size)==F])){
  
  print(i)
  
  tmp_df = crispr.power[effect_size == i,]
  
  crispr_genes <- data.frame(
    eqtl = "CRISPRi",
    sig_cre_genes = length(target_cgenes),
    n_sig_genes = 160,
    "prop_power>80" = nrow(tmp_df[power>=0.80])/nrow(tmp_df),
    "mean" = mean(tmp_df$power),
    effect_size = i)
  
  all_crispr_prop = bind_rows(all_crispr_prop,crispr_genes)
  
}

all_crispr_prop