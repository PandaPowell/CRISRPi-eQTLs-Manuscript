### Identify gens with trans networks to test ####
library(data.table)
library(ggplot2)
library(dplyr)

setwd("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K")

### read Data ####
trans_network <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt") 
head(trans_network)
length(table(trans_network$gRNA_target))
## 895 gRNAs

### Read cis genes ####
cis <- fread('analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt')
head(cis)

## Remove trans-results from genes with +/-5mbp of gRNAs
trans_network$TSS = ifelse(trans_network$strand == 1, trans_network$start_position, trans_network$end_position)
trans_network$tss_distance = with(trans_network, ifelse(gsub("chr","",gRNA_chr) == chromosome_name, as.integer(gRNA_end)-TSS, NA))
trans_network = trans_network[abs(tss_distance) > 5e+06 | is.na(tss_distance) == T,]

### Get gRNAs with 5 or more trans #### 
trans_gRNA <- as.data.frame(table(trans_network$gRNA_target[trans_network$sig_trans==1]))
trans_gRNA_5 <- trans_gRNA[trans_gRNA$Freq>5,]
head(trans_gRNA_5)

cis_genes <- cis$hgnc_symbol_closest[cis$gRNA_target %in% trans_gRNA_5$Var1]
head(cis_genes)

cis_genes_chr <- trans_network[trans_network$hgnc_symbol %in% cis_genes, c('hgnc_symbol', 'chromosome_name')] %>% distinct()
head(cis_genes_chr)

write.table(cis_genes_chr, 'analysis/02_QTL_calling/genes_cis_w_trans.txt', quote = F, sep = '\t', col.names = F, row.names = F)
