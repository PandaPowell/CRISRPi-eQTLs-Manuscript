suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(gridExtra)
  library(data.table)
  library(biomaRt)
})

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "chromosome_name", "transcription_start_site"),
                   filters = "biotype",
                   values = "protein_coding",
                   mart = ensembl)
protein_coding_genes <- gene_list[,c("ensembl_gene_id","hgnc_symbol")] %>% distinct(hgnc_symbol, .keep_all = T)

calculate_enrichment_gold <- function(cgenes, egenes, mendel_genes, protein_coding_genes, n_iterations = 10000, seed = 123) {
  
  set.seed(seed)
  
  # Simulate null overlaps
  cgene_null <- replicate(n_iterations, {
    sampled <- sample(protein_coding_genes, length(cgenes), replace = FALSE)
    sum(sampled %in% mendel_genes)
  })
  
  egene_null <- replicate(n_iterations, {
    sampled <- sample(protein_coding_genes, length(egenes), replace = FALSE)
    sum(sampled %in% mendel_genes)
  })
  
  # Observed overlaps
  cgene_obs <- sum(cgenes %in% mendel_genes)
  egene_obs <- sum(egenes %in% mendel_genes)
  
  # Median expected overlaps
  cgene_exp <- median(cgene_null)
  egene_exp <- median(egene_null)
  
  # Fisher test tables
  cgene_mat <- matrix(c(cgene_obs, length(cgenes) - cgene_obs, cgene_exp, length(cgenes) - cgene_exp), nrow = 2)
  egene_mat <- matrix(c(egene_obs, length(egenes) - egene_obs, egene_exp, length(egenes) - egene_exp), nrow = 2)
  compare_mat <- matrix(c(cgene_obs, length(cgenes) - cgene_obs, egene_obs, length(egenes) - egene_obs), nrow = 2)
  
  # Fisher tests
  cgene_fisher <- fisher.test(cgene_mat, alternative = "greater")
  egene_fisher <- fisher.test(egene_mat, alternative = "greater")
  compare_fisher <- fisher.test(compare_mat)
  
  # Stats table with p-values added
  stats <- data.frame(
    group = c("cgenes", "egenes"),
    observed = c(cgene_obs, egene_obs),
    expected = c(cgene_exp, egene_exp),
    n_total = c(length(cgenes), length(egenes)),
    p_value = c(cgene_fisher$p.value, egene_fisher$p.value)
  )
  
  return(list(
    stats = stats,
    compare_p = compare_fisher$p.value
  ))
}

# Load gasperini cis-genes with trans effects
gas_df = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/genes_cis_w_trans.txt")
gas_cis_genes = unlist(gas_df$V1)

# Filter cis-genes for those with GWAS targets
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

gas_cis_genes = gas_cis_genes[gas_cis_genes %in% c(cres_w_grnas[significant == 1, gene_name], cres_w_grnas_egene$gene_name)]
morrs_cis_genes = c("GFI1B", "NFE2", "IKZF1", "HHEX", "RUNX1")

# Load gasperini data linking grnas to cis-target genes
gas_cis_data <- fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt")

# Make barplot of number of associated genes
# link cis-gene to gRNA target in trans data
gas_cis_data <- gas_cis_data[hgnc_symbol_closest %in% c(gas_cis_genes, morrs_cis_genes),]

# Load trans morris data
trans_network_john <- readxl::read_xlsx("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/science.adh7699_table_s4.xlsx", skip = 2) %>%
  filter(`Q-value` < 0.1) %>%
  filter(`Network CRE` %in% c("GFI1B (CRE-1)", "HHEX","IKZF1","NFE2 (CRE-1)", "RUNX1")) # remove miRNAs and take most significant cre when multiple are present

# Load Gasperini tran effects
trans_network_gas <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt") %>%
  filter(qvalue_trans < 0.10)

### Remove trans-results from genes with +/-5mbp of gRNAs
trans_network_gas$TSS = ifelse(trans_network_gas$strand == 1, trans_network_gas$start_position, trans_network_gas$end_position)
trans_network_gas$tss_distance = with(trans_network_gas, ifelse(gsub("chr","",gRNA_chr) == chromosome_name, as.integer(gRNA_end)-TSS, NA))
trans_network_gas = trans_network_gas[abs(tss_distance) > 5e+06 | is.na(tss_distance) == T,]

# link cis-gene to gRNA target in trans data
gas_cis_data <- gas_cis_data[hgnc_symbol_closest %in% c(gas_cis_genes, morrs_cis_genes),]
trans_network_gas <- trans_network_gas[gRNA_target %in% gas_cis_data$gRNA_target,] %>%
  left_join(gas_cis_data[,c("gRNA_target","hgnc_symbol_closest")], by = "gRNA_target")

# Load trans-egenes
all_res_trans = fread("Supplementary_Table9_trans_eGenes.txt")

# Load gold-genes
gw_gold_genes = fread("gw_gold_genes.txt")

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
  
  a = sum(unique(trans_g$Gene) %in% gw_gold_genes$V1)
  b = length(unique(trans_g$Gene))-a
  
  gene_df = data.frame()
  
  for (eqtl in unique(cis.gene.trans$Cell.type)){
    
    cell.type.trans = cis.gene.trans[cis.gene.trans$Cell.type == eqtl,]
    all <- merge(trans_g[,c("Gene","Log2 Fold Change","P-value")], cell.type.trans, by.x='Gene', by.y='phenotype_id_trans')
    
    cell_df = calculate_enrichment_gold(mendel_genes = gw_gold_genes$V1,
                              protein_coding_genes = protein_coding_genes$hgnc_symbol,
                              cgenes = unique(trans_g$Gene), 
                              egenes = unique(cell.type.trans$phenotype_id_trans))
    
    cell_df = cell_df$stats
    cell_df$cell = eqtl
    
    gene_df = bind_rows(gene_df, cell_df)
    
  }
  
  gene_df$cis_gene = gene
  sting.seq.df = rbind(sting.seq.df, gene_df)
  
}

sting.seq.df$study = "morris"

gasperini.df = data.frame()

# Loop through cis-genes with trans-effects
for (gene in gas_cis_genes){
  
  print(gene)
  
  # link cis-gene to gRNA target in trans data
  cis_data_g <- gas_cis_data[gas_cis_data$hgnc_symbol_closest == gene,]
  trans_g <- trans_network_gas[trans_network_gas$gRNA_target %in% cis_data_g$gRNA_target,]
  # trans-eQTL data for cis-gene
  cis.gene.trans = all_res_trans[all_res_trans$phenotype_id_cis == gene,]
  
  if(nrow(cis.gene.trans)<1){
    print("No significant trans-eQTLs")
    next
  }
  
  print(nrow(cis.gene.trans))
  
  gene_df = data.frame()
  
  for (eqtl in unique(cis.gene.trans$Cell.type)){
    
    cell.type.trans = cis.gene.trans[cis.gene.trans$Cell.type == eqtl,]
    all <- merge(trans_g[,c("hgnc_symbol","log2fc","pvalue")], cell.type.trans, by.x='hgnc_symbol', by.y='phenotype_id_trans')

    
    cell_df = calculate_enrichment_gold(mendel_genes = gw_gold_genes$V1,
                                        protein_coding_genes = protein_coding_genes$hgnc_symbol,
                                        cgenes = unique(trans_g$hgnc_symbol), 
                                        egenes = unique(cell.type.trans$phenotype_id_trans))
    
    cell_df = cell_df$stats
    cell_df$cell = eqtl
    
    gene_df = bind_rows(gene_df, cell_df)

  }
  
  gene_df$cis_gene = gene
  gasperini.df = rbind(gasperini.df, gene_df)
  
}

gasperini.df$study = "gasperini"

all_res = bind_rows(sting.seq.df,gasperini.df)

all_res$fdr <- p.adjust(all_res$p_value, method = "BH")

fwrite(all_res, "supplementary_table_13.txt", quote = F, row.names = F, sep = ",")
