rm(list=ls())
library(data.table)
library(tidyverse)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

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
  mutate(pair_id = paste0(SS_coord,Gene)) %>% 
  mutate(gwas = as.character(gwas), data = "stingseq") %>%
  distinct(pair_id, .keep_all=T) %>% 
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`, pvalue = `Q-value (1 Mb)`,significant, ensembl_id = `Ensembl ID`, gene_name = Gene,
         tss_distance = `TSS Distance`, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+")))

ss_finemap_bcx = fread("00.intersect_data/sting_seq_credible_sets/BCX/merged_credset.txt") %>%
  rename(SS_coord = SNP_coord) %>%
  left_join(sting_seq,"SS_coord") %>%
  mutate(pair_id = paste0(SS_coord,Gene), data = "stingseq") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(significant = ifelse(`Q-value (1 Mb)`<0.1,1,0)) %>%
  dplyr::select(grna_target = SS_coord, logfc = `Log2 fold-change`,  pvalue = `Q-value (1 Mb)`, significant, ensembl_id = `Ensembl ID`, gene_name = Gene,
                tss_distance = `TSS Distance`, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data) %>%
  mutate(grna_pos = as.integer(str_extract(grna_target, "(?<=:)[0-9]+")))

missing_stingseq = sting_seq[!SS_coord %in% c(ss_finemap$grna_target,ss_finemap_bcx$grna_target),]
unique(missing_stingseq$Population)
length(unique(missing_stingseq$SS_coord))
missing_stingseq %>% group_by(Population) %>% summarise(n = n())
# EAS
missing_ea = missing_stingseq[Population == "EA",]

gas_finemap = fread("00.intersect_data/gasperini_credible_sets/merged_gasperini_credset_closest.txt") %>%
  left_join(resample_results, "target_site") %>%
  mutate(gwas = as.character(gwas), data= "gasperini") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS,
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

fwrite(gas_finemap, "Gasperini_gRNAs_intersecting_GWAS.txt", quote = F, row.names = F, sep = ",")

gas_finemap_bcx = fread("00.intersect_data/gasperini_credible_sets/BCX/merged_credset.txt") %>%
  left_join(resample_results, "target_site") %>%
  distinct(pair_id, .keep_all=T) %>%
  mutate(grna_target = paste0(chr.x, ":", lower_grna_target_site,"-",upper_grna_target_site),
         significant = as.integer(rejected),
         tss_distance = upper_grna_target_site-TSS, data = "gasperini",
         grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target, logfc = xi, pvalue = p_value, significant, ensembl_id = gene_id, gene_name = gene_short_name,
                tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)

# Write list of target sites for Jasper
#fwrite(list(c(unique(gas_finemap$target_site), gas_finemap_bcx$target_site)),"gasperini_target_sites_Jasper.txt")

encode_finemap = fread("00.intersect_data/encode_credible_sets/merged_encode_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(name, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS) %>%
  mutate(gwas = as.character(gwas), data = "encode", grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
  dplyr::select(grna_target = target_site, logfc = Effectsize, pvalue = pValue, significant = Significant, ensembl_id, 
                gene_name,tss_distance, finemap_snp_intersect_grna, sentinel_snp = finemap_snp, gwas, data, grna_pos)
  
encode_finemap_bcx = fread("00.intersect_data/encode_credible_sets/BCX/merged_credset.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  left_join(encode, "target_site") %>%
  distinct(name, .keep_all = T) %>%
  mutate(tss_distance = chromEnd_GRC38-startTSS, gwas = as.character(gwas), 
         data = "encode", grna_pos = round((lower_grna_target_site+upper_grna_target_site)/2)) %>%
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
  mutate(target_gene = paste0(grna_target,"_",ensembl_id)) %>%
  distinct(target_gene, .keep_all=T)

cat(length(unique(cres_w_grnas$grna_target)), "CREs after merging")

# Calculate dis rank
dis_rank = cres_w_grnas %>%
  distinct(target_gene, .keep_all=T) %>% 
  group_by(grna_target) %>% 
  mutate(dis_rank = rank(abs(tss_distance)))

####################################
# eQTL catalogue coloc results.    #
####################################
eqtl_w_ss_cat = fread("00.intersect_data/stingseq_eqtl_overlap/eqtl_catalogue_coloc_stingseq.txt") %>%
  mutate(gwas = as.character(gwas)) %>%
  dplyr::select(grna_target = SS_coord, eQTL_variant = eqtl_hit_hg19, PP.H4.abf,
                ensembl_id = molecular_id, finemap_snp_intersect_grna,
                sentinel_snp = GWAS_variant, gwas, eqtl = eqtl_name, gwas_hit_hg38, eqtl_hit_hg38,
                study_id, dataset_id, gwas_pval) %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id") %>% # Join with GRC38
  dplyr::select(grna_target, eQTL_variant, PP.H4.abf,
                ensembl_id, gene_name, finemap_snp_intersect_grna,
                sentinel_snp , gwas, eqtl , gwas_hit_hg38, eqtl_hit_hg38,
                study_id, dataset_id, gwas_pval)

eqtl_w_gas_cat = fread("00.intersect_data/gasperini_eqtl_overlap/eqtl_catalogue_coloc_gasperini.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant = eqtl_hit_hg19, PP.H4.abf,
                ensembl_id = molecular_id, finemap_snp_intersect_grna, sentinel_snp = GWAS_variant, 
                gwas = gwas_name, eqtl = eqtl_name, gwas_hit_hg38, eqtl_hit_hg38, upper_grna_target_site,
                study_id, dataset_id, gwas_pval) %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id") %>%
  dplyr::select(grna_target , eQTL_variant, PP.H4.abf,
                ensembl_id, gene_name, finemap_snp_intersect_grna, sentinel_snp,
                gwas , eqtl , gwas_hit_hg38, eqtl_hit_hg38,
                study_id, dataset_id, gwas_pval)

eqtl_w_encode_cat = fread("00.intersect_data/encode_eqtl_overlap/eqtl_catalogue_coloc_crispri.txt") %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant = eqtl_hit_hg19, PP.H4.abf,
                ensembl_id = molecular_id, finemap_snp_intersect_grna, sentinel_snp = GWAS_variant,
                gwas = gwas_name, eqtl = eqtl_name, gwas_hit_hg38, eqtl_hit_hg38,
                study_id, dataset_id, upper_grna_target_site, gwas_pval) %>% unique() %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id") %>%
  dplyr::select(grna_target , eQTL_variant ,PP.H4.abf,
                ensembl_id, gene_name, finemap_snp_intersect_grna, sentinel_snp,
                gwas , eqtl , gwas_hit_hg38, eqtl_hit_hg38,
                study_id, dataset_id, gwas_pval)

split1 = bind_rows(eqtl_w_ss_cat, eqtl_w_gas_cat, eqtl_w_encode_cat) %>%
  unique() %>%
  mutate(unique_id = paste(grna_target,ensembl_id,dataset_id, sep = "_")) %>%
  filter(is.na(gene_name) == F) # split df by missing gene names

split2 = bind_rows(eqtl_w_ss_cat, eqtl_w_gas_cat, eqtl_w_encode_cat) %>%
  unique() %>%
  mutate(unique_id = paste(grna_target,ensembl_id,dataset_id, sep = "_")) %>%
  filter(is.na(gene_name)) %>%# split df by missing gene names
  left_join(annot.GRC38.v33[,c("ensembl_id","gene_name")], "ensembl_id") %>%
  dplyr::select(-gene_name.x) %>% rename(gene_name = gene_name.y)

eqtl_cat_comb = bind_rows(split1,split2) %>% filter(eqtl != "fibroblast")

####################################
# GTEx coloc results.             ##
####################################
eqtl_w_ss_gtex =  fread("00.intersect_data/stingseq_eqtl_overlap/GTEx_coloc_stingseq.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>% 
  mutate(gwas = as.character(gwas)) %>%
  dplyr::select(grna_target = SS_coord, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                ensembl_id = phenotype_id, tss_distance, finemap_snp_intersect_grna, sentinel_snp = GWAS_variant, 
                gwas, gwas_pval = pval) %>%
  mutate(eqtl = "GTEx", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1]) %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id")

eqtl_w_gas_gtex = fread("00.intersect_data/gasperini_eqtl_overlap/GTEx_coloc_crispri.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                ensembl_id = phenotype_id, tss_distance, finemap_snp_intersect_grna, 
                sentinel_snp = GWAS_variant, gwas = gwas_name, gwas_pval = pval) %>%
  mutate(eqtl = "GTEx", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1]) %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id")

eqtl_w_encode_gtex = fread("00.intersect_data/encode_eqtl_overlap/GTEx_coloc_crispri.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>% 
  unique() %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>% 
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                ensembl_id = phenotype_id, tss_distance, finemap_snp_intersect_grna, sentinel_snp = GWAS_variant,
                gwas = gwas_name, gwas_pval = pval) %>%
  mutate(eqtl = "GTEx", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1]) %>%
  left_join(annot.GRC38[,c("ensembl_id","gene_name")], "ensembl_id")

eqtl_gtex_comb = bind_rows(eqtl_w_ss_gtex, eqtl_w_gas_gtex, eqtl_w_encode_gtex) %>%
  unique() %>% 
  mutate(coord = str_extract(sentinel_snp, "[0-9]\\:[0-9]+(?=\\:[A-Z])"),
         unique_id = paste0(grna_target,ensembl_id))

####################################
# MAGE coloc results.              #
####################################
eqtl_w_ss_mage =  fread("00.intersect_data/stingseq_eqtl_overlap/MAGE_coloc_stingseq.txt") %>%
  filter(`p-value` < 1e-05 & pval_nominal < 1e-03) %>% 
  mutate(gwas = as.character(gwas)) %>%
  dplyr::select(grna_target = SS_coord, eQTL_variant, beta = beta, pval_nominal, PP.H4.abf,
                gene_name = geneSymbol, ensembl_id = ensemblID, tss_distance, finemap_snp_intersect_grna, 
                sentinel_snp = GWAS_variant, gwas, gwas_pval = `p-value`) %>%
  mutate(eqtl = "MAGE", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1])

eqtl_w_gas_mage = fread("00.intersect_data/gasperini_eqtl_overlap/MAGE_coloc_gasperini.txt") %>%
  filter(`p-value` < 1e-05 & pval_nominal < 1e-03) %>% 
  mutate(gwas = as.character(gwas)) %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta, pval_nominal, PP.H4.abf,
                gene_name = geneSymbol, ensembl_id = ensemblID, tss_distance, finemap_snp_intersect_grna, 
                sentinel_snp = GWAS_variant, gwas, gwas_pval = `p-value`) %>%
  mutate(eqtl = "MAGE", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1])

eqtl_w_encode_mage = fread("00.intersect_data/encode_eqtl_overlap/MAGE_coloc_crispri.txt") %>%
  filter(`p-value` < 1e-05 & pval_nominal < 1e-03) %>%
  unique() %>% mutate(gwas = as.character(gwas)) %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta, pval_nominal, PP.H4.abf,
                gene_name = geneSymbol, ensembl_id = ensemblID, tss_distance, finemap_snp_intersect_grna, 
                sentinel_snp = GWAS_variant, gwas, gwas_pval = `p-value`) %>%
  mutate(eqtl = "MAGE", ensembl_id = str_split_fixed(ensembl_id, "\\.", 2)[,1])

eqtl_mage_comb = bind_rows(eqtl_w_ss_mage, eqtl_w_gas_mage, eqtl_w_encode_mage) %>%
  unique() %>% 
  mutate(coord = str_extract(eQTL_variant, "chr[0-9]+\\:[0-9]+(?=\\_[A-Z])"),
         unique_id = paste0(coord,ensembl_id))

####################################
# OneK1K coloc results.    #
####################################
eqtl_w_ss_1k = fread("00.intersect_data/stingseq_eqtl_overlap/All_coloc_stingseq_overlap.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>% 
  mutate(gwas = as.character(gwas)) %>%
  dplyr::select(grna_target = SS_coord, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                gene_name = phenotype_id, tss_distance = start_distance, finemap_snp_intersect_grna,
                sentinel_snp = GWAS_variant, gwas = gwas_name, eqtl, gwas_pval = pval)

eqtl_w_gas_1k = fread("00.intersect_data/gasperini_eqtl_overlap/All_coloc_crispri_overlap.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                gene_name = phenotype_id, tss_distance = start_distance, finemap_snp_intersect_grna, 
                sentinel_snp = GWAS_variant, gwas = gwas_name, eqtl, gwas_pval = pval)

eqtl_w_encode_1k = fread("00.intersect_data/encode_eqtl_overlap/All_coloc_crispri_overlap.txt") %>%
  filter(pval < 1e-05 & pval_nominal < 1e-03) %>%
  mutate(target_site = paste0(chr,":",lower_grna_target_site,"-",upper_grna_target_site)) %>%
  dplyr::select(grna_target = target_site, eQTL_variant, beta = beta.y, pval_nominal, PP.H4.abf,
                gene_name = phenotype_id, tss_distance = start_distance, finemap_snp_intersect_grna,
                sentinel_snp = GWAS_variant, gwas = gwas_name, eqtl, gwas_pval = pval)

eqtl_onek_comb = bind_rows(eqtl_w_ss_1k, eqtl_w_gas_1k, eqtl_w_encode_1k) %>% unique() %>%
  mutate(coord = str_extract(eQTL_variant, "[0-9]+\\:[0-9]+(?=\\:[A-Z])"),
         eqtl = ifelse(eqtl == "DC_mean", "DC_cells", eqtl),
         eqtl = ifelse(eqtl == "Other_cells", "other_cells", eqtl),
         eqtl = ifelse(eqtl == "Other_T_cells", "other_T_cells", eqtl),
         unique_id = paste0(grna_target,gene_name,eqtl)) %>%
  unique() %>% left_join(annot.GRC38, "gene_name")

cres_w_grnas_egene = bind_rows(eqtl_cat_comb,eqtl_gtex_comb,eqtl_mage_comb,eqtl_onek_comb) %>%
  dplyr::select("grna_target", "eQTL_variant", "PP.H4.abf", "ensembl_id", "gene_name",
         "tss_distance", "finemap_snp_intersect_grna", "sentinel_snp", "gwas", 
         "eqtl", "eqtl_hit_hg38", "study_id", 
         "dataset_id", "beta", "pval_nominal", "gwas_pval") %>% 
  left_join(merged_sites, "grna_target") %>%
  dplyr::select(-grna_target) %>% 
  rename(grna_target = grna_target_merge) %>%
  filter(is.na(grna_target) == F) %>%
  mutate(target_gene = paste0(grna_target,"_", ensembl_id)) %>%
  mutate(tss_distance = ifelse(tss_distance < -1000000, -1000000, tss_distance),
         tss_distance = ifelse(tss_distance > 1000000, 1000000, tss_distance))

# Final check before writing files
if(sum(cres_w_grnas_egene$grna_target %in% cres_w_grnas$grna_target) == nrow(cres_w_grnas_egene)){
  # write dataframes
  fwrite(cres_w_grnas, "cres_with_grnas.txt", quote = F, row.names = F)
  fwrite(cres_w_grnas_egene, "cres_with_grna_eqtls.txt", quote = F, row.names = F)
} else{
  warning("Something is wrong double check the CREs")
}

## Run scripts to calculate distance ranks
# Submit the SLURM job
job_id <- system("sbatch 06.Run_gene_distances.sh", intern = TRUE)

# Extract the job ID from the submission response
job_id <- sub("Submitted batch job ", "", job_id)

# Function to check if the job is still running or pending
is_job_running <- function(job_id) {
  # Use squeue to check the job status
  result <- system(paste("squeue -j", job_id), intern = TRUE)
  # If the result contains more than just the header, the job is still running/pending
  return(length(result) > 1)
}

# Wait until the job is done
while (is_job_running(job_id)) {
  cat("Job", job_id, "is still running...\n")
  Sys.sleep(10) # Wait for 10 seconds before checking again
}

cat("Job", job_id, "has completed.\n")

# Join calculated distance ranks and distance to main files
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

name = c("B", "Mono", "NK","DC","CD4_T","CD8_T", "other", "other_T","GTEx","MAGE")

df1 = cres_w_grnas_egene %>%
  mutate(unique_id = paste0(grna_target,ensembl_id,eqtl)) %>%
  filter(eqtl %in% c(paste0(name,"_cells"), "GTEx", "MAGE"))

df2 = cres_w_grnas_egene %>%
  mutate(unique_id = paste0(grna_target,ensembl_id,dataset_id)) %>%
  filter(!eqtl %in% c(paste0(name,"_cells"), "GTEx", "MAGE")) %>%
  dplyr::select(-tss_distance, -beta, -pval_nominal)

## Onek1k gene distances and ranks
onek_gene_distances = fread("eQTL_gene_distances/Onek1k_gene_distance.txt") %>%
  mutate(unique_id = paste0(grna_target,ensembl_id,cell_type)) %>%
  dplyr::select(unique_id, dis_rank)

# Load GTEx gene distances
gtex_gene_distances = fread("eQTL_gene_distances/GTEx_gene_distance.txt") %>%
  mutate(unique_id = paste0(grna_target,ensembl_id,"GTEx")) %>%
  dplyr::select(unique_id, dis_rank)

mage_gene_distances = fread("eQTL_gene_distances/MAGE_gene_distance.txt") %>%
  mutate(unique_id = paste0(grna_target,ensembl_id,"MAGE")) %>%
  dplyr::select(unique_id, dis_rank)

gene_distances = bind_rows(onek_gene_distances,gtex_gene_distances, mage_gene_distances)

gtex_onek_mage = df1 %>% left_join(gene_distances,"unique_id")

View(gtex_onek_mage[is.na(dis_rank),])

# Load cis-genes and calculate distance rank
cat_sumstats = fread("eQTL_gene_distances/eqtl_cat_gene_distances.txt") %>%
  mutate(unique_id = paste0(grna_target,ensembl_id, dataset_id),
         grna_target_dataset = paste0(grna_target,"_",dataset_id), 
         tss_distance = position-start) %>%
  arrange(tss_distance) %>% 
  distinct(unique_id, .keep_all = T) %>% # remove duplicates 
  group_by(grna_target_dataset) %>% 
  distinct(unique_id, .keep_all =T) %>% 
  mutate(dis_rank = rank(abs(tss_distance)))

eqtl_cat = df2 %>% 
  left_join(cat_sumstats[,c("unique_id","beta","pvalue","median_tpm","tss_distance","dis_rank")], "unique_id") %>%
  rename(pval_nominal = pvalue) %>%
  filter(pval_nominal < 1e-03 | is.na(pval_nominal) == T) %>%
  mutate(tss_distance = ifelse(tss_distance < -1000000, -1000000, tss_distance),
         tss_distance = ifelse(tss_distance > 1000000, 1000000, tss_distance))

cres_w_grnas_egene = bind_rows(gtex_onek_mage,eqtl_cat)

if(sum(cres_w_grnas_egene$grna_target %in% cres_w_grnas$grna_target) == nrow(cres_w_grnas_egene)){
  # write dataframes
  fwrite(cres_w_grnas_egene, "cres_with_grna_eqtls.txt", quote = F, row.names = F)
} else{
  warning("Something is wrong double check the CREs")
}

##############
## Analysis ##
##############
View(cres_w_grnas_egene[gene_name == ""])
View(cres_w_grnas_egene[ensembl_id == ""])

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

overlapping_cres_cgene = unique(cres_w_grnas$target_gene[cres_w_grnas$grna_target %in% overlapping_cres & cres_w_grnas$significant == 1])
cat("Number of target cgenes overlapping CRES =",length(overlapping_cres_cgene),"\n")

overlapping_cres_egene = unique(cres_w_grnas_egene$target_gene[cres_w_grnas_egene$grna_target %in% overlapping_cres])
cat("Number of target egenes overlapping CRES =",length(overlapping_cres_egene),"\n")

overlapping_target_genes = overlapping_cres_cgene[overlapping_cres_cgene %in% overlapping_cres_egene]
cat("Number of overlapping genes =",length(overlapping_target_genes),"\n")

overlapping_cres_egene = unique(cres_w_grnas_egene$target_gene[cres_w_grnas_egene$grna_target %in% overlapping_cres])
length(overlapping_cres_egene)

overlapping_target_genes = sum(overlapping_cres_cgene %in% overlapping_cres_egene)
