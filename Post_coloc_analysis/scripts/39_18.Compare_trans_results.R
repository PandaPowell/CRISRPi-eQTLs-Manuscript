rm(list=ls())
options(bitmapType="cairo")
library(data.table)
#library(arrow)
library(tidyverse)
library(UpSetR)
library(cowplot)
library(ggthemes)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

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

# Load gasperini cis-genes with trans effects
gas_df = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/genes_cis_w_trans.txt")
gas_cis_genes = unlist(gas_df$V1)

# Filter cis-genes for those with GWAS targets
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

gas_cis_genes = gas_cis_genes[gas_cis_genes %in% c(cres_w_grnas[significant == 1, gene_name], cres_w_grnas_egene$gene_name)]
morrs_cis_genes = c("GFI1B", "NFE2", "IKZF1", "HHEX", "RUNX1")

# Load MetaLCL trans-results
metalcl = fread("MetaLCL_extracted_transQTLs.txt") %>%
  left_join(annot.GRC38, by = c("molecular_trait_id" = "ensembl_id")) %>%
  mutate(tss_distance = position - start) %>% 
  filter((abs(tss_distance) > 5e+06 & chromosome == chr_number) | chromosome != chr_number, 
         cis_gene %in% gas_cis_genes)

p.11 = metalcl %>%
  filter(nlog10p > 1.3) %>%
  mutate(Cell.type = "LCL", pval = 10^(-nlog10p), af = ac/(2*sample_size)) %>%
  dplyr::select(Cell.type, variant_id = variant_kgpID, phenotype_id_trans = gene_name, pval, b = beta, b_se = se, af, phenotype_id_cis = cis_gene)

# Summarize and arrange in descending order
plot_data <- p.11 %>%
  group_by(phenotype_id_cis) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(cis_gene = factor(phenotype_id_cis, levels = phenotype_id_cis))  # Convert to factor for ordered plotting

png(filename = "plots/metalcl_barplot_p05.png", width = 10, height = 8, units = "in", res = 300)
# Create bar plot
ggplot(plot_data, aes(x = phenotype_id_cis, y = n)) +
  geom_col() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Number of Significant (p<0.05) \n trans-genes per cis-gene",
       x = "cis-gene",
       y = "n trans-genes")
dev.off()

# Load Onek1k sc-eQTL results
onek_dir = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/"

#### trans results #####
all_res_trans <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0(onek_dir,'analysis/02_QTL_calling//lvl1Res/', cells, '_mean_mx/'), pattern = '*trans_qtl_pairs.allpairs*')
  if (length(files)==0) {
    all_res_trans[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0(onek_dir,'analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    res_df$Cell.type <- cells
    res_df <- res_df[res_df$pval<0.05,]
    all_res_trans[[cells]] <- res_df
  }
}

all_res_trans <- do.call('rbind.data.frame', all_res_trans)
#table(all_res_trans$variant_id, all_res_trans$Cell.type)

## get cis gene ####
all_res <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0(onek_dir,'analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/'), pattern = paste0(cells,'\\.independent.*Filt*'))
  if (length(files)==0) {
    #all_res[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0(onek_dir,'analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    if (nrow(res_df)<1) {
      #all_res[[cells]] <- NA
    } else {
      res_df$Cell.type <- cells
      all_res[[cells]] <- res_df
    }
  }
}

all_res <- do.call('rbind.data.frame', all_res)
all_res = all_res[all_res$phenotype_id %in% gas_cis_genes,]

##analyze sharing #####
all_res_trans <- merge(all_res_trans,all_res[,c("phenotype_id","variant_id","pval_perm","slope","Cell.type")], by = c('Cell.type','variant_id'), suffixes = c("_trans","_cis"))
all_res_trans$cisQTL <- paste0(all_res_trans$phenotype_id_cis,'_',all_res_trans$variant_id)

all_res_trans = bind_rows(all_res_trans, p.11)

fwrite(all_res_trans, "Supplementary_Table9_trans_eGenes.txt", sep = ",", row.names = F, quote = F)

### upset plots #####
library(ComplexHeatmap)
## at the gene level

cisQTLs <- (sapply(c('LCL','B','CD4_T','CD8_T','Mono','NK'), function(cell)
  unique(all_res_trans[all_res_trans$Cell.type == cell,'phenotype_id_cis'])))

m1 = make_comb_mat(cisQTLs)

png("plots/LCL_onek1k_w_trans005.png",width = 6, height = 4, units = "in", res = 300)

UpSet(m1, set_order = c('LCL','CD4_T','CD8_T','NK','B','Mono'),
      comb_order = order(comb_size(m1), decreasing = T), 
      top_annotation = upset_top_annotation(m1, add_numbers = TRUE))

dev.off()

# Load Morris trans-gene results
trans_network_john <- readxl::read_xlsx("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/science.adh7699_table_s4.xlsx", skip = 2) %>%
  filter(`Q-value` < 0.1) %>%
  filter(`Network CRE` %in% c("GFI1B (CRE-1)", "HHEX","IKZF1","NFE2 (CRE-1)", "RUNX1")) # remove miRNAs and take most significant cre when multiple are present

# Barplot of number of trans-genes per cis-gene
n_trans_assoc <- trans_network_john %>%
  group_by(`Network CRE`) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(study = "Morris") %>%
  rename(cis_gene = `Network CRE`)

n_trans_assoc$cis_gene = gsub(" \\(CRE-1\\)", "", n_trans_assoc$cis_gene)

# Load Gasperini tran effects
trans_network_gas <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt") %>%
  filter(qvalue_trans < 0.10)

### Remove trans-results from genes with +/-5mbp of gRNAs
trans_network_gas$TSS = ifelse(trans_network_gas$strand == 1, trans_network_gas$start_position, trans_network_gas$end_position)
trans_network_gas$tss_distance = with(trans_network_gas, ifelse(gsub("chr","",gRNA_chr) == chromosome_name, as.integer(gRNA_end)-TSS, NA))
trans_network_gas = trans_network_gas[abs(tss_distance) > 5e+06 | is.na(tss_distance) == T,]

# Load gasperini data linking grnas to cis-target genes
gas_cis_data <- fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt")

# Make barplot of number of associated genes
# link cis-gene to gRNA target in trans data
gas_cis_data <- gas_cis_data[hgnc_symbol_closest %in% c(gas_cis_genes, morrs_cis_genes),]
trans_network_gas <- trans_network_gas[gRNA_target %in% gas_cis_data$gRNA_target,] %>%
  left_join(gas_cis_data[,c("gRNA_target","hgnc_symbol_closest")], by = "gRNA_target")

# Barplot of number of trans-genes per cis-gene
n_trans_assoc_gas <- trans_network_gas %>%
  group_by(hgnc_symbol_closest) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(study = "Gasperini") %>%
  rename(cis_gene = hgnc_symbol_closest)

svg(filename = "plots/supp_figs/trans_cgene_barplot.svg", width = 6, height = 6)

colors <- c(
  "Gasperini" = "#7c260b",
  "Morris" = "#ed111a"
)

# Create bar plot
ggplot(bind_rows(n_trans_assoc, n_trans_assoc_gas), aes(x = cis_gene, y = n, fill = study)) +
  geom_col(
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.8,
    alpha = 0.8
  ) +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  labs(title = "",
       x = "cis-gene",
       y = "n trans-genes") + coord_flip() +
  theme(legend.position = "top")


dev.off()

# Calculate correlation between crispri datasets
overlapping_genes= morrs_cis_genes[morrs_cis_genes %in% gas_cis_genes]

gene = "GFI1B"

print(gene)

# Filter crispri associations to cis-gene
trans_j <- trans_network_john[grep(gene, trans_network_john$`Network CRE`),]
trans_g = trans_network_gas[hgnc_symbol_closest == gene,]

comb_networks = trans_j %>%
  left_join(trans_g, by = c("Gene" = "hgnc_symbol")) %>%
  filter(is.na(log2fc) == F)

r = cor(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)
p = cor.test(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)$p.value
p = ifelse(p<2e-16, 2e-16, p)

svg(filename = "plots/supp_figs/cor_plots_GFI1B.svg", width = 4, height = 4)

ggplot(comb_networks, aes(x=`Log2 Fold Change`, y=log2fc)) + 
  geom_point(color = "black", alpha = 0.6) +
  annotate("text", x = -1, y = 1.5, label = paste0("r = ",round(r,2),", P < ", formatC(p, format = "e", digits = 1)), size = 3, hjust = 0.5, family = "sans") +
  ylab("Trans-gene LogFC - Gasperini dataset") +
  xlab("Trans-gene LogFC - Morris dataset") +
  theme_cowplot()

dev.off()

gene = "HHEX"

print(gene)

# Filter crispri associations to cis-gene
trans_j <- trans_network_john[grep(gene, trans_network_john$`Network CRE`),]
trans_g = trans_network_gas[hgnc_symbol_closest == gene,]

comb_networks = trans_j %>%
  left_join(trans_g, by = c("Gene" = "hgnc_symbol")) %>%
  filter(is.na(log2fc) == F)

r = cor(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)
p = cor.test(comb_networks$`Log2 Fold Change`, comb_networks$log2fc)$p.value
p = ifelse(p<2e-16, 2e-16, p)

svg(filename = "plots/supp_figs/cor_plot_HHEX.svg", width = 4, height = 4)

ggplot(comb_networks, aes(x=`Log2 Fold Change`, y=log2fc)) + 
  geom_point(color = "black", alpha = 0.6) +
  annotate("text", x = -0.4, y = 1, label = paste0("r = ",round(r,2),", P = ", formatC(p, format = "e", digits = 1)), size = 3, hjust = 0.5, family = "sans") +
  ylab("Trans-gene LogFC - Gasperini dataset") +
  xlab("Trans-gene LogFC - Morris dataset") +
  theme_cowplot()

dev.off()

# Make supplementary tables
fwrite(trans_network_john, "supp8_trans_genes_morris.txt", sep = ",", row.names = F, quote = F)

fwrite(trans_network_gas, "supp9_trans_genes_gas.txt", sep = ",", row.names = F, quote = F)

# load k562 genes
# Load gasperini data, filter to targeted enhancers
resample_results = fread("CRISPR_data/resampling_results.txt") %>% unique() %>% 
  filter(site_type != "NTC") #%>% distinct(gene_short_name, .keep_all = T)

#### compare full table trans qtl calling with John's data #####
# enrichment function
calculate_enrichment <- function(gene_list_a,gene_list_b,gene_sample) {
  n_iterations <- 10000
  set.seed(123)
  
  gene_overlap_counts <- numeric(n_iterations)
  
  for (i in seq_len(n_iterations)) {
    random_genes <- sample(gene_sample, size = length(gene_list_b), replace = FALSE)
    gene_overlap_counts[i] <- sum(random_genes %in% gene_list_a)
  }
  
  gene_random_counts <- median(gene_overlap_counts)
  observed_overlap <- sum(gene_list_b %in% gene_list_a)

  gene_table <- matrix(c(observed_overlap, length(gene_list_a) - observed_overlap,
                          gene_random_counts, length(gene_list_a) - gene_random_counts), nrow = 2)
  
  test=fisher.test(gene_table, alternative = 'two.sided', conf.int = T)
  
  or=test$estimate
  p=test$p.value
  conf.int=test$conf.int
  return(list(odds_ratio = or, p_value = p, conf.int = conf.int))
}

# Loop through all cis-genes and cell-types and calculate correlations
sting.seq.df = data.frame()

for (gene in c("GFI1B", "HHEX","IKZF1","NFE2", "RUNX1")){
  
  print(gene)
  
  # Filter crispri associations to cis-gene
  trans_g <- trans_network_john[grep(gene, trans_network_john$`Network CRE`),]
  # Filter eQTL associations to cis-gene
  cis.gene.trans = all_res_trans[all_res_trans$phenotype_id_cis == gene,]
  
  if(nrow(cis.gene.trans)<1){
    print("No significant trans-eQTLs")
    next
  }
  
  gene_df = data.frame()
  
  for (eqtl in unique(cis.gene.trans$Cell.type)){
    
    cell.type.trans = cis.gene.trans[cis.gene.trans$Cell.type == eqtl,]
    
    # Flip trans-associations to decreasing allele
    if(is.na(unique(cell.type.trans$slope)) == F && unique(cell.type.trans$slope) > 0){
      print("Flipping effect sizes")
      cell.type.trans$b = cis.gene.trans$b*-1
    }
    
    all <- merge(trans_g[,c("Gene","Log2 Fold Change","P-value")], cell.type.trans, by.x='Gene', by.y='phenotype_id_trans')
    # Calculate enrichment
    result <- calculate_enrichment(unique(cell.type.trans$phenotype_id_trans), 
                                   unique(trans_g$Gene), 
                                   unique(resample_results$gene_short_name))

    r <- round(cor(all$`Log2 Fold Change`, all$b), 4)
    p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
    conf.int = cor.test(all$`Log2 Fold Change`, all$b, conf.int = T)$conf.int
    cell_df = data.frame(cell.type = eqtl,
                         cor = r,
                         cor_pval = p,
                         cor_lower = conf.int[1],
                         cor_upper = conf.int[2],
                         n_cgenes = nrow(trans_g),
                         n_egenes = nrow(cell.type.trans),
                         n_overlap = nrow(all),
                         enrichment_OR = result$odds_ratio,
                         enrichment_p = result$p_value,
                         enrichment_lower_ci = result$conf.int[1],
                         enrichment_upper_ci = result$conf.int[2])
    gene_df = rbind(gene_df, cell_df)
  }
  
  gene_df$cis_gene = gene
  sting.seq.df = rbind(sting.seq.df, gene_df)
  
}

sting.seq.df$study = "morris"

#### compare full table trans qtl calling with Gasperini's data #####
# Load morris cis-genes with trans hits

gasperini.df = data.frame()

# Loop through cis-genes with trans-effects
for (gene in gas_cis_genes){
  
  print(gene)
  
  # link cis-gene to gRNA target in trans data
  cis_data_g <- gas_cis_data[gas_cis_data$hgnc_symbol_closest == gene,]
  cis_gene_effect = resample_results[grna_group %in% cis_data_g$gRNA_target] %>% 
    arrange(p_value) %>% 
    distinct(grna_group, .keep_all=T) %>%
    select(grna_group, cis_xi = xi)
  trans_g <- trans_network_gas[trans_network_gas$gRNA_target %in% cis_data_g$gRNA_target,] %>%
    left_join(cis_gene_effect, by = c("gRNA_target" = "grna_group"))
  
  if(unique(is.na(trans_g$cis_xi)) == F){
    
    print("Flipping CRISPRi associations to negative cis-effects")
    
    pos_idx <- trans_g$cis_xi > 0
    trans_g$cis_xi[pos_idx] <- -trans_g$cis_xi[pos_idx]
    trans_g$xi[pos_idx] <- -trans_g$xi[pos_idx]
    
  }
  
  # trans-eQTL data for cis-gene
  cis.gene.trans = all_res_trans[all_res_trans$phenotype_id_cis == gene,]

  if(nrow(cis.gene.trans)<1){
    print("No significant trans-eQTLs")
    next
  }
  
  print(nrow(cis.gene.trans))

  gene_df = data.frame()

    for (eqtl in unique(cis.gene.trans$Cell.type)){
      
      cis.gene.trans.2 = cis.gene.trans[cis.gene.trans$Cell.type == eqtl,]
      
      if(unique(is.na(cis.gene.trans.2$slope)) == F){
        
        print("Flipping eQTL associations to negative cis-effects")
        
        pos_idx <- cis.gene.trans.2$slope > 0
        cis.gene.trans.2$slope[pos_idx] <- -cis.gene.trans.2$slope[pos_idx]
        cis.gene.trans.2$b[pos_idx] <- -cis.gene.trans.2$b[pos_idx]
        
      }
      
      all <- merge(trans_g[,c("hgnc_symbol","log2fc","pvalue")], cis.gene.trans.2, by.x='hgnc_symbol', by.y='phenotype_id_trans')
        
        if(nrow(all) < 10){
          print("Too few overlaps")
          next
        }
        
        # Calculate enrichment
        result <- calculate_enrichment(unique(cell.type.trans$phenotype_id_trans), 
                                       unique(trans_g$hgnc_symbol), 
                                       unique(resample_results$gene_short_name))
        r <- round(cor(all$log2fc, all$b), 4)
        p <- cor.test(all$log2fc, all$b)$p.value
        conf.int = cor.test(all$log2fc, all$b, conf.int = T)$conf.int
        cell_df = data.frame(cell.type = eqtl,
                             cor = r,
                             cor_pval = p,
                             cor_lower = conf.int[1],
                             cor_upper = conf.int[2],
                             n_cgenes = nrow(trans_g),
                             n_egenes = nrow(cis.gene.trans.2),
                             n_overlap = nrow(all),
                             enrichment_OR = result$odds_ratio,
                             enrichment_p = result$p_value,
                             enrichment_lower_ci = result$conf.int[1],
                             enrichment_upper_ci = result$conf.int[2])
        gene_df = rbind(gene_df, cell_df)
      
      
    }
  if (nrow(gene_df) > 0) {
    gene_df$cis_gene = gene
    gasperini.df = rbind(gasperini.df, gene_df)
  } else {
    print(paste("Skipping", gene, "due to no valid overlaps"))
  }
  
}

gasperini.df$study = "gasperini"

all_res = bind_rows(sting.seq.df,gasperini.df)

# Forest plots:
library(ggplot2)
library(dplyr)

# Data for correlation estimates and confidence intervals
cor_data <- data.frame(
  Gene = c("GFI1B", "HHEX", "HHEX", "HHEX", "MRPL4", "MRPL4", "MRPL4", "GFI1B",
           "ZNHIT1", "ZNHIT1", "ZNHIT1", "ZNHIT1", "MAX", "SNRNP40", "SNRNP40"),
  Cell_Type = c("Mono", "CD8_T", "NK", "LCL", "CD4_T", "CD8_T", "LCL", "Mono",
                "CD4_T", "CD8_T", "NK", "LCL", "LCL", "B", "LCL"),
  Estimate = c(0.1609, 0.3921, 0.2402, -0.3394, 0.0410, 0.4054, 0.3170, 0.1984,
               -0.2802, -0.0967, -0.0910, -0.0364, -0.0837, 0.0093, -0.3474),
  P_value = c(0.0767427672, 0.2965876892, 0.3702429713, 0.2566283583,
              0.8043304015, 0.0682425004, 0.0006606825, 0.0837428382,
              0.1952450713, 0.5525978092, 0.5616351174, 0.8306203446,
              0.7494925636, 0.9759587006, 0.2684660852),
  Lower_CI = c(-0.01740814, -0.36777946, -0.29005631, -0.75009537,
               -0.27812577, -0.03181826, 0.13965411, -0.02680891,
               -0.62074379, -0.39631203, -0.38094125, -0.35622095,
               -0.54250629, -0.54447902, -0.76817515),
  Upper_CI = c(0.3291987, 0.8380097, 0.6575952, 0.2602877,
               0.3519536, 0.7124328, 0.4746331, 0.4043740,
               0.1491865, 0.2214364, 0.2152165, 0.2910483,
               0.4135984, 0.5574253, 0.2828545)
)

# Flip signs of LCL rows to harmonize direction
cor_data <- cor_data %>%
  mutate(
    flip = Cell_Type == "LCL" & !(Gene == "MRPL4" & abs(Estimate - 0.317) < 1e-6),
    Estimate = ifelse(flip, -Estimate, Estimate),
    Lower_CI = ifelse(flip, -Lower_CI, Lower_CI),
    Upper_CI = ifelse(flip, -Upper_CI, Upper_CI)
  ) %>%
  select(-flip)


# Sorting the data by Gene and Estimate for better visualization
cor_data <- cor_data %>% arrange(Gene, Estimate)
cor_data$Gene <- factor(cor_data$Gene, levels = rev(c("MRPL4", "GFI1B", "HHEX", "ZNHIT1", "MAX", "SNRNP40")))

cell_type_colors <- c(
  "CRISPRi" = "#d35e60",  # a complementary deep red
  "LCL"     = "#6794a7",  # light blue-gray
  "Mono"    = "#014d64",  # dark teal
  "CD8_T"   = "#01a2d9",  # bright blue
  "NK"      = "#76c0c1",  # soft teal
  "CD4_T"   = "#8e6c8a",  # muted purple
  "B"       = "#c85200"   # burnt orange
)

svg("plots/supp_figs/corr_foresplot.svg",width = 6, height = 4)

# Create the forest plot
ggplot(cor_data, aes(x = Estimate, y = Gene, color = Cell_Type)) +
  geom_point(size = 1, position = position_dodge(width = 0.6)) +  # Multiple points per gene
  #geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, position = position_dodge(width = 0.6), size=2, alpha=0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at 0
  labs(x = "Correlation Estimate", y = "", title = "",
       color = "Cell Type") +
  theme_cowplot() +
  scale_color_manual(values = cell_type_colors) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9)
  )

dev.off()

# Data for enrichment odds ratio and confidence intervals
enrichment_data <- data.frame(
  Gene = c("HHEX", "HHEX", "HHEX", "MRPL4", "MRPL4", "MRPL4", "GFI1B", "ZNHIT1",
           "ZNHIT1", "ZNHIT1", "ZNHIT1", "MAX", "SNRNP40", "SNRNP40"),
  Cell_Type = c("CD8_T", "NK", "LCL", "CD4_T", "CD8_T", "LCL", "Mono", "CD4_T",
                "CD8_T", "NK", "LCL", "LCL", "B", "LCL"),
  OR = c( 1.000, 1.464, 0.927, 1.088, 0.665, 1.227, 1.190, 0.784, 1.118, 1.286, 
         0.922, 1.000, 1.639, 0.629),
  Lower_CI = c( 0.3489962, 0.6331417, 0.3997391, 0.6645482, 0.3587289, 
               0.9177746, 0.8116576, 0.4272635, 0.6842802, 0.7872940, 0.5680369, 
               0.4784488, 0.6239715, 0.2779634),
  Upper_CI = c( 2.865361, 3.516484, 2.138152, 1.786404, 1.212700, 1.646479, 
               1.749397, 1.424647, 1.834401, 2.118019, 1.493634, 2.090088, 4.602003, 1.369979)
)
# Sorting the data by Gene and OR
enrichment_data <- enrichment_data %>% arrange(Gene, OR)
enrichment_data$Gene <- factor(enrichment_data$Gene, levels = rev(c("MRPL4", "GFI1B", "HHEX", "ZNHIT1", "MAX", "SNRNP40")))

svg("plots/supp_figs/enrichment_foresplot.svg",width = 6, height = 4)

# Create the forest plot for enrichment OR
ggplot(enrichment_data, aes(x = OR, y = Gene, color = Cell_Type)) +
  geom_point(size = 1, position = position_dodge(width = 0.6)) +  # Multiple points per gene
  #geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, position = position_dodge(width = 0.6), size=2, alpha=0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +  # Reference line at 0
  scale_x_log10() +  # Log scale for odds ratios
  labs(x = "Enrichment Odds Ratio", y = "Gene", title = "",
       color = "Cell Type") +
  theme_cowplot() +
  scale_color_manual(values = cell_type_colors) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9))

dev.off()

# adjust p-values for multiple testing
all_res$cor_fdr <- p.adjust(all_res$cor_pval, method = "BH")
all_res$enrich_fdr <- p.adjust(all_res$enrichment_p, method = "BH")

# Barplots of number of associated egenes and cgenes
all_res$cis_gene <- factor(all_res$cis_gene, levels = rev(c("MRPL4", "GFI1B", "HHEX", "ZNHIT1", "MAX", "SNRNP40")))

cell_type_colors <- c(
  "CRISPRi" = "#d35e60",  # a complementary deep red
  "LCL"     = "#6794a7",  # light blue-gray
  "Mono"    = "#014d64",  # dark teal
  "CD8_T"   = "#01a2d9",  # bright blue
  "NK"      = "#76c0c1",  # soft teal
  "CD4_T"   = "#8e6c8a",  # muted purple
  "B"       = "#c85200"   # burnt orange
)

# Set CRISPRi as the last level so it's drawn last (on top)
ordered_cell_types <- ordered_cell_types <- rev(c("CRISPRi","NK","LCL","CD8_T","CD4_T","Mono","B"))  # example order

df_plot = all_res %>%
  distinct(cell.type, cis_gene, .keep_all = T) %>%
  select(cis_gene, cell.type, n_egenes, n_overlap) %>%
  mutate(cell.type = factor(cell.type, levels = ordered_cell_types))  # order here matters!

#png("plots/figure_plots/trans_egenes_barplots.png",width = 3, height = 4, units = "in", res = 300)
svg("plots/figure_plots/trans_egenes_barplots.svg", width = 6, height = 6)

ggplot(df_plot, aes(x = cis_gene, y = n_overlap, fill = cell.type)) +
  geom_col(
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.5,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = cell_type_colors,
    name = NULL,
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "",
    x = "cis-gene",
    y = "n overlapping trans-genes"
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.05, 0.2))) +
  coord_flip() +
  theme_cowplot() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.line.y = element_blank()
  )

dev.off()

fwrite(all_res, "trans_comparison_results.txt", sep = ",", row.names = F, quote = F)

# All trans-eGenes plot
svg("plots/supp_figs/trans_egenes_barplots.svg", width = 6, height = 6)

ggplot(df_plot, aes(x = cis_gene, y = n_egenes, fill = cell.type)) +
  geom_col(
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.8,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = cell_type_colors,
    name = NULL,
    guide = guide_legend(nrow = 1)
  ) +  labs(
    title = "",
    x = "cis-gene",
    y = "n trans-genes"
  ) +
  coord_flip() +
  theme_cowplot() +
  theme(legend.position = "top")

dev.off()

fwrite(all_res, "trans_comparison_results.txt", sep = ",", row.names = F, quote = F)
