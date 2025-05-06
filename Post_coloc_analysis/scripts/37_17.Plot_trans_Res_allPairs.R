############################
# Author: Winona Oliveros
# Date: 2024/01/11 
# Analize trans results
############################
library(data.table)

setwd("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K")

#### trans results #####
all_res_trans <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0('analysis/02_QTL_calling//lvl1Res/', cells, '_mean_mx/'), pattern = '*trans_qtl_pairs.allpairs*')
  if (length(files)==0) {
    all_res_trans[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    res_df$Cell.type <- cells
    res_df <- res_df[res_df$pval<0.05,]
    all_res_trans[[cells]] <- res_df
  }
}

all_res_trans <- do.call('rbind.data.frame', all_res_trans)
table(all_res_trans$variant_id, all_res_trans$Cell.type)

## get cis gene ####
all_res <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/'), pattern = paste0(cells,'\\.independent.*Filt*'))
  if (length(files)==0) {
    #all_res[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
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

##Analize sharing #####
all_res_trans <- merge(all_res_trans,all_res[,c("phenotype_id","variant_id","pval_perm","slope","Cell.type")], by = c('Cell.type','variant_id'), suffixes = c("_trans","_cis"))
head(all_res_trans)
all_res_trans$cisQTL <- paste0(all_res_trans$phenotype_id_cis,'_',all_res_trans$variant_id)

### upset plots #####
library(ComplexHeatmap)

cisQTLs <- (sapply(c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), function(cell) 
  unique(all_res_trans[all_res_trans$Cell.type == cell,'cisQTL'])
)
)

m1 = make_comb_mat(cisQTLs)
m1

"/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/plots/Azimuth_lvl1_upset_cisqtls_w_trans005.pdf"

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl1_upset_cisqtls_w_trans005.pdf",
    width = 6, height = 4)
UpSet(m1, set_order = c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), 
      comb_order = order(comb_size(m1), decreasing = T), top_annotation = upset_top_annotation(m1, add_numbers = TRUE))
dev.off()

## at the gene level 
cisQTLs <- (sapply(c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), function(cell) 
  unique(all_res_trans[all_res_trans$Cell.type == cell,'phenotype_id_cis'])
)
)

m1 = make_comb_mat(cisQTLs)
m1

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl1_upset_cisgenes_w_trans005.pdf",
    width = 6, height = 4)
UpSet(m1, set_order = c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), 
      comb_order = order(comb_size(m1), decreasing = T), top_annotation = upset_top_annotation(m1, add_numbers = TRUE))
dev.off()

### cis upset #####
## at the gene level 
cisQTLs <- (sapply(c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), function(cell) 
  unique(all_res[all_res$Cell.type == cell,'phenotype_id'])
)
)

m1 = make_comb_mat(cisQTLs)
m1

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl1_upset_cisgenes_w_network.pdf",
    width = 6, height = 4)
UpSet(m1, set_order = c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), 
      comb_order = order(comb_size(m1), decreasing = T), top_annotation = upset_top_annotation(m1, add_numbers = TRUE))
dev.off()

## Level 2
cell_types <- read.delim('analysis/02_QTL_calling/cell_types_lvl2.txt', header = F)
all_res_trans_lvl2 <- list()

for (cells in cell_types$V1) {
  files <- list.files(path = paste0('analysis/02_QTL_calling//lvl2Res/', cells, '_mean_mx/'), pattern = '*trans_qtl_pairs.allpairs*')
  if (length(files)==0) {
    #all_res_trans_lvl2[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl2Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    res_df$Cell.type <- cells
    res_df <- res_df[res_df$pval<0.05,]
    all_res_trans_lvl2[[cells]] <- res_df
  }
}

all_res_trans_lvl2 <- do.call('rbind.data.frame', all_res_trans_lvl2)
table(all_res_trans_lvl2$variant_id, all_res_trans_lvl2$Cell.type)

## get cis gene level 2####
all_res_level2 <- list()

for (cells in cell_types$V1) {
  files <- list.files(path = paste0('analysis/02_QTL_calling/lvl2Res/', cells, '_mean_mx/'), pattern = paste0(cells,'\\.independent.*Filt*'))
  if (length(files)==0) {
    #all_res_level2[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl2Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    if (nrow(res_df)<1) {
      #all_res_level2[[cells]] <- NA
    } else {
      res_df$Cell.type <- cells
      all_res_level2[[cells]] <- res_df
    }
  }
}

all_res_level2 <- do.call('rbind.data.frame', all_res_level2)

##Analize sharing ####
all_res_trans_lvl2 <- merge(all_res_trans_lvl2,all_res_level2[,c("phenotype_id","variant_id","pval_perm","slope","Cell.type")], by = c('Cell.type','variant_id'), suffixes = c("_trans","_cis"))
head(all_res_trans_lvl2)  
all_res_trans_lvl2$cisQTL <- paste0(all_res_trans_lvl2$phenotype_id_cis,'_',all_res_trans_lvl2$variant_id)

### upset plots #####
library(ComplexHeatmap)

cisQTLs <- (sapply(cell_types$V1, function(cell) 
  unique(all_res_trans_lvl2[all_res_trans_lvl2$Cell.type == cell,'cisQTL'])
)
)

m1 = make_comb_mat(cisQTLs)
m1

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl2_upset_cisqtls_w_trans005.pdf",
    width = 6, height = 6)
UpSet(m1[comb_degree(m1) > 0], #set_order = cell_types$V1, 
      comb_order = order(comb_size(m1[comb_degree(m1) > 0]), decreasing = T), top_annotation = upset_top_annotation(m1[comb_degree(m1) > 0], add_numbers = TRUE))
dev.off()

## at the gene level #####
cisQTLs <- (sapply(cell_types$V1, function(cell) 
  unique(all_res_trans_lvl2[all_res_trans_lvl2$Cell.type == cell,'phenotype_id_cis'])
)
)

m1 = make_comb_mat(cisQTLs)
m1

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl2_upset_cisgenes_w_trans005.pdf",
    width = 6, height = 6)
UpSet(m1[comb_degree(m1) > 0], #set_order = c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), 
      comb_order = order(comb_size(m1[comb_degree(m1) > 0]), decreasing = T), top_annotation = upset_top_annotation(m1[comb_degree(m1) > 0], add_numbers = TRUE))
dev.off()


### plot number of trans genes at p<0.05  #####
trans_genes_summary <- as.data.frame(table(all_res_trans$Cell.type, all_res_trans$cisQTL))
trans_genes_summary <- trans_genes_summary[trans_genes_summary$Freq>0,]
head(trans_genes_summary)

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl1_density_n_transgenes.pdf",
    width = 4, height = 4)
ggplot(trans_genes_summary, aes(x=Freq)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Nº trans-Genes / cisQTL") +
  theme_classic()
dev.off()

pdf("analysis/02_QTL_calling/plots/Azimuth_lvl1_density_n_transgenes_ct.pdf",
    width = 4, height = 4)
ggplot(trans_genes_summary, aes(x=Freq,group=Var1, fill=Var1)) +
  geom_density(adjust=1.5, alpha=.4) +
  ggtitle("Nº trans-Genes / cisQTL") +
  theme_classic()
dev.off()


#### compare full table trans qtl calling with gasperini and/or John #####
trans_network_john <- readxl::read_xlsx('data/science.adh7699_table_s4.xlsx', skip = 2)

# NK HHEX_10:94429467:C:T 719 #####
trans_g <- trans_network_john[grep('HHEX', trans_network_john$`Network CRE`),]
trans_g <- as.data.frame(trans_g)
table(trans_g$`SNP ID`)
snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% c('10:94429467:C:T'),]
#trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
# for (variant in snps_trans_crispr) {
#   for (cell_type in unique(trans_genes_lvl2$cell_type)) {
    #trans_crispr <- trans_g[trans_g$`SNP ID` == gsub(':','-',variant),]
trans_all <- read.delim("analysis/02_QTL_calling/lvl1Res/NK_mean_mx/trans_res/OneK1K_NK.HHEX.['10:94429467:C:T'].trans_qtl_pairs.allpairs.csv")
trans_all <- read.delim("analysis/02_QTL_calling/lvl2Res/NK_mean_mx/trans_res/OneK1K_NK.HHEX.['10:94429467:C:T'].trans_qtl_pairs.allpairs.csv")
    all <- merge(trans_g[,c("Gene","Log2 Fold Change","P-value")], trans_all, by.x='Gene', by.y='phenotype_id')
    
    r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
    p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
    
    pdf(paste0('analysis/02_QTL_calling/plots/HHEX_','10:94429467:C:T','_','NKlvl2','_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
    print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    all <- all[all$Gene!='HHEX',]
    r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], (all$b[all$`P-value`<0.05])), 2)
    p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/HHEX_','10:94429467:C:T','_','NKlvl2','_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
    print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], (all$b[all$`P-value`<0.05 & all$pval<0.1])), 2)
    p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], all$b[all$`P-value`<0.05 & all$pval<0.1])$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/HHEX_','10:94429467:C:T','_','NKlvl2','_cor_CRISPR_tensor_trans01.p005.pdf'), width = 5, height = 4)
    print(  ggplot(all[all$`P-value`<0.05 & all$pval<0.1,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
#   }
#   
# }
    
# Mono GFI1B_9:135879542:C:T 510 #####
    
    trans_g <- trans_network_john[grep('GFI1B', trans_network_john$`Network CRE`),]
    trans_g <- as.data.frame(trans_g)
    table(trans_g$`SNP ID`)
    snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))
    
    #trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% c('10:94429467:C:T'),]
    #trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]
    
    # trans_genes_lvl1$phenotype_id %in% trans_g$Gene
    # trans_genes_lvl2$phenotype_id %in% trans_g$Gene
    
    ### plot 
    # for (variant in snps_trans_crispr) {
    #   for (cell_type in unique(trans_genes_lvl2$cell_type)) {
    trans_crispr <- trans_g[trans_g$`SNP ID` == '9-135879542-C-T',]
    trans_all <- read.delim("analysis/02_QTL_calling/lvl2Res/CD16_Mono_mean_mx///trans_res/OneK1K_CD16_Mono.GFI1B.['9:135879542:C:T'].trans_qtl_pairs.allpairs.csv")
    all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_all, by.x='Gene', by.y='phenotype_id')

    r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
    p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value

    pdf(paste0('analysis/02_QTL_calling/plots/GFI1B_','9:135879542:C:T','_','CD16_Mono','_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
    print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) +
              geom_point() +
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()

    all <- all[all$Gene!='GFI1B',]
    r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], (all$b[all$`P-value`<0.05])), 2)
    p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value

    pdf(paste0('analysis/02_QTL_calling//plots/GFI1B_','9:135879542:C:T','_','CD16_Mono','_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
    print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) +
              geom_point() +
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()



    r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], (all$b[all$`P-value`<0.05 & all$pval<0.1])), 2)
    p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], all$b[all$`P-value`<0.05 & all$pval<0.1])$p.value

    pdf(paste0('analysis/02_QTL_calling//plots/GFI1B_','9:135879542:C:T','_','CD16_Mono','_cor_CRISPR_tensor_trans01.p005.pdf'), width = 5, height = 4)
    print(  ggplot(all[all$`P-value`<0.05 & all$pval<0.1,], aes(y=`Log2 Fold Change`, x=b)) +
              geom_point() +
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    ## read CHIP 
    chip_john <- readxl::read_xlsx('data/John.science.adh7699_tables_s1_to_s4/science.adh7699_table_s4.xlsx', skip = 2,range = 'Table S4B!A3:S20330')
    
    filt_chip <- all[all$Gene %in% chip_john$Gene[chip_john$GFI1B_ChIP == 1],]
    r <- round(cor(filt_chip$`Log2 Fold Change`, (filt_chip$b)), 2)
    p <- cor.test(filt_chip$`Log2 Fold Change`, filt_chip$b)$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/GFI1B_','9:135879542:C:T','_','CD16_Mono','_cor_CRISPR_tensor_transCHIP.pdf'), width = 5, height = 4)
    print(  ggplot(filt_chip, aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    r <- round(cor(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], (filt_chip$b[filt_chip$`P-value`<0.05])), 2)
    p <- cor.test(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], filt_chip$b[filt_chip$`P-value`<0.05])$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/GFI1B_','9:135879542:C:T','_','Mono','_cor_CRISPR_tensor_transCHIP.p005.pdf'), width = 5, height = 4)
    print(  ggplot(filt_chip[filt_chip$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    r <- round(cor(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1], (filt_chip$b[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1])), 2)
    p <- cor.test(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1], filt_chip$b[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1])$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/GFI1B_','9:135879542:C:T','_','CD16_Mono','_cor_CRISPR_tensor_transCHIP01.p005.pdf'), width = 5, height = 4)
    print(  ggplot(filt_chip[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()

    ## IKZF1 #####
    
    trans_g <- trans_network_john[grep('IKZF1', trans_network_john$`Network CRE`),]
    trans_g <- as.data.frame(trans_g)
    table(trans_g$`SNP ID`)
    snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))
    
    #trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% c('10:94429467:C:T'),]
    #trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]
    
    # trans_genes_lvl1$phenotype_id %in% trans_g$Gene
    # trans_genes_lvl2$phenotype_id %in% trans_g$Gene
    
    ### plot 
    snps <- paste0(trans_genes_summary$Var1[grep('IKZF1',trans_genes_summary$Var2)],'-',trans_genes_summary$Var2[grep('IKZF1',trans_genes_summary$Var2)])
    # for (variant in snps_trans_crispr) {
       for (var in snps) {
    trans_crispr <- trans_g[trans_g$`SNP ID` == '7-50427982-A-G',]
    variant <- gsub('.*_','',var)
    cell_type <- gsub('-.*','',var)
    trans_all <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/",cell_type,"_mean_mx///trans_res/OneK1K_",cell_type,".IKZF1.['",variant,"'].trans_qtl_pairs.allpairs.csv"))
    all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_all, by.x='Gene', by.y='phenotype_id')
    
    # r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
    # p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
    # print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(all$`Log2 Fold Change`)), hjust=0) +
    #           theme_classic() )
    # dev.off()
    
    all <- all[all$Gene!='IKZF1',]    
    # r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], (all$b[all$`P-value`<0.05])), 2)
    # p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
    # print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(all$`Log2 Fold Change`[all$`P-value`<0.05])), hjust=0) +
    #           theme_classic() )
    # dev.off()
    # 
    # 
    # 
    # r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], (all$b[all$`P-value`<0.05 & all$pval<0.1])), 2)
    # p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1], all$b[all$`P-value`<0.05 & all$pval<0.1])$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans01.p005.pdf'), width = 5, height = 4)
    # print(  ggplot(all[all$`P-value`<0.05 & all$pval<0.1,], aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(all$`Log2 Fold Change`[all$`P-value`<0.05 & all$pval<0.1])), hjust=0) +
    #           theme_classic() )
    # dev.off()
    
    ## read CHIP 
    #chip_john <- readxl::read_xlsx('data/John.science.adh7699_tables_s1_to_s4/science.adh7699_table_s4.xlsx', skip = 2,range = 'Table S4B!A3:S20330')
    filt_chip <- all[all$Gene %in% chip_john$Gene[chip_john$IKZF1_ChIP == 1],]
    # r <- round(cor(filt_chip$`Log2 Fold Change`, (filt_chip$b)), 2)
    # p <- cor.test(filt_chip$`Log2 Fold Change`, filt_chip$b)$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_transCHIP.pdf'), width = 5, height = 4)
    # print(  ggplot(filt_chip, aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(filt_chip$`Log2 Fold Change`)), hjust=0) +
    #           theme_classic() )
    # dev.off()
    # 
    # r <- round(cor(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], (filt_chip$b[filt_chip$`P-value`<0.05])), 2)
    # p <- cor.test(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], filt_chip$b[filt_chip$`P-value`<0.05])$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_transCHIP.p005.pdf'), width = 5, height = 4)
    # print(  ggplot(filt_chip[filt_chip$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05])), hjust=0) +
    #           theme_classic() )
    # dev.off()
    
    r <- round(cor(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1], (filt_chip$b[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1])), 2)
    p <- cor.test(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1], filt_chip$b[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1])$p.value
    
    pdf(paste0('analysis/02_QTL_calling//plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_transCHIP01.p005.pdf'), width = 5, height = 4)
    print(  ggplot(filt_chip[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              annotate("text", x=max(filt_chip$b)-0.10, y=min(filt_chip$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05 & filt_chip$pval<0.1])), hjust=0) +
              theme_classic() )
    dev.off()
    
       }

    ## RPL23A #####
    
trans_network <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt")
head(trans_network)

cis_data <- fread('analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt')
head(cis_data)

cis_data_g <- cis_data[cis_data$hgnc_symbol_closest == 'RPL23A',]
cis_data_g

trans_g <- trans_network[trans_network$gRNA_target %in% cis_data_g$gRNA_target,]
#trans_g <- trans_g[trans_g$hgnc_symbol != 'RPL23A',]

#trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% c('10:94429467:C:T'),]
#trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
snps <- paste0(trans_genes_summary$Var1[grep('RPL23A',trans_genes_summary$Var2)],'-',trans_genes_summary$Var2[grep('RPL23A',trans_genes_summary$Var2)])
# for (variant in snps_trans_crispr) {
for (var in snps) {
  trans_crispr <- trans_g
  variant <- gsub('.*_','',var)
  cell_type <- gsub('-.*','',var)
  trans_all <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/",cell_type,"_mean_mx///trans_res/OneK1K_",cell_type,".RPL23A.['",variant,"'].trans_qtl_pairs.allpairs.csv"))
  all <- merge(trans_crispr[,c("hgnc_symbol","log2fc","qvalue_trans")], trans_all, by.x='hgnc_symbol', by.y='phenotype_id')
  
  r <- round(with(subset(all, (is.finite(log2fc)) & (is.finite(b))), cor(log2fc, b, method = "pearson", use = "complete.obs")), 2)
  p <- with(subset(all, (is.finite(log2fc)) & (is.finite(b))), cor.test(log2fc, b, method = "pearson", use = "complete.obs"))$p.value
  
  pdf(paste0('analysis/02_QTL_calling/plots/RPL23A_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
  print(  ggplot(all, aes(y=`log2fc`, x=b)) + 
            geom_point() + 
            geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
            annotate("text", x=max(all$b)-0.15, y=max(all$`log2fc`), label=paste0("r = ", r), hjust=0) +
            annotate("text", x=max(all$b)-0.15, y=max(all$`log2fc`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
            annotate("text", x=max(all$b)-0.15, y=min(all$`log2fc`)+0.25, label=paste0("n = ", length(all$`log2fc`)), hjust=0) +
            theme_classic() )
  dev.off()
  
  all <- all[all$hgnc_symbol!='RPL23A',]    
  r <- round(cor(all$`log2fc`[all$`qvalue_trans`<0.05], (all$b[all$`qvalue_trans`<0.05])), 2)
  p <- cor.test(all$`log2fc`[all$`qvalue_trans`<0.05], all$b[all$`qvalue_trans`<0.05])$p.value
  
  pdf(paste0('analysis/02_QTL_calling/plots/RPL23A_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
  print(  ggplot(all[all$`qvalue_trans`<0.05,], aes(y=`log2fc`, x=b)) + 
            geom_point() + 
            geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.05]), label=paste0("r = ", r), hjust=0) +
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.05])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05])-0.15, y=min(all$`log2fc`[all$`qvalue_trans`<0.05])+0.25, label=paste0("n = ", length(all$`log2fc`[all$`qvalue_trans`<0.05])), hjust=0) +
            theme_classic() )
  dev.off()
  
  
  
  r <- round(cor(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1], (all$b[all$`qvalue_trans`<0.05 & all$pval<0.1])), 2)
  p <- cor.test(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1], all$b[all$`qvalue_trans`<0.05 & all$pval<0.1])$p.value
  
  pdf(paste0('analysis/02_QTL_calling/plots/RPL23A_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans01.p005.pdf'), width = 5, height = 4)
  print(  ggplot(all[all$`qvalue_trans`<0.05 & all$pval<0.1,], aes(y=`log2fc`, x=b)) + 
            geom_point() + 
            geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05 & all$pval<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1]), label=paste0("r = ", r), hjust=0) +
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05 & all$pval<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
            annotate("text", x=max(all$b[all$`qvalue_trans`<0.05 & all$pval<0.1])-0.15, y=min(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1])+0.25, label=paste0("n = ", length(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1])), hjust=0) +
            theme_classic() )
  dev.off()
  
}    

### try with gasperini res #######

# trans_network <- fread('/gpfs/commons/groups/lappalainen_lab/jmorris/230525_OneK1K/data/231011_GasperiniTransNetwork.txt')
# head(trans_network)
# 
# cis_data <- fread('analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt')
# head(cis_data)
# 
# cis_data_g <- cis_data[cis_data$hgnc_symbol_closest == 'HHEX',]
# cis_data_g
# 
# trans_g <- trans_network[trans_network$gRNA_target %in% cis_data_g$gRNA_target,]
# #trans_g <- trans_g[trans_g$hgnc_symbol != 'RPL23A',]
# table(trans_g$gRNA_target)
# snps_trans_crispr <- gsub('-',':',unique(trans_g$gRNA_target))
# 
# #trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% c('10:94429467:C:T'),]
# #trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]
# 
# # trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# # trans_genes_lvl2$phenotype_id %in% trans_g$Gene
# 
# ### plot 
# snps <- paste0(trans_genes_summary$Var1[grep('HHEX',trans_genes_summary$Var2)],'-',trans_genes_summary$Var2[grep('HHEX',trans_genes_summary$Var2)])
# for (variant2 in snps_trans_crispr) {
#   for (var in snps) {
#     trans_crispr <- trans_g[trans_g$gRNA_target == variant2,]
#     variant <- gsub('.*_','',var)
#     cell_type <- gsub('-.*','',var)
#     trans_all <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/",cell_type,"_mean_mx///trans_res/OneK1K_",cell_type,".HHEX.['",variant,"'].trans_qtl_pairs.allpairs.csv"))
#     all <- merge(trans_crispr[,c("hgnc_symbol","log2fc","qvalue_trans")], trans_all, by.x='hgnc_symbol', by.y='phenotype_id')
#     
#     r <- round(with(subset(all, (is.finite(log2fc)) & (is.finite(b))), cor(log2fc, b, method = "pearson", use = "complete.obs")), 2)
#     p <- with(subset(all, (is.finite(log2fc)) & (is.finite(b))), cor.test(log2fc, b, method = "pearson", use = "complete.obs"))$p.value
#     
#     pdf(paste0('analysis/02_QTL_calling/plots/HHEX_',variant,'_',cell_type,'_',variant2,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
#     print(  ggplot(all, aes(y=`log2fc`, x=b)) + 
#               geom_point() + 
#               geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
#               annotate("text", x=max(all$b)-0.15, y=max(all$`log2fc`), label=paste0("r = ", r), hjust=0) +
#               annotate("text", x=max(all$b)-0.15, y=max(all$`log2fc`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
#               annotate("text", x=max(all$b)-0.15, y=min(all$`log2fc`)+0.25, label=paste0("n = ", length(all$`log2fc`)), hjust=0) +
#               theme_classic() )
#     dev.off()
#     
#     all <- all[all$hgnc_symbol!='HHEX',]    
#     r <- round(cor(all$`log2fc`[all$`qvalue_trans`<0.1], (all$b[all$`qvalue_trans`<0.1])), 2)
#     p <- cor.test(all$`log2fc`[all$`qvalue_trans`<0.1], all$b[all$`qvalue_trans`<0.1])$p.value
#     
#     pdf(paste0('analysis/02_QTL_calling/plots/HHEX_',variant,'_',cell_type,variant2,'_cor_CRISPR_tensor_trans.p01.pdf'), width = 5, height = 4)
#     print(  ggplot(all[all$`qvalue_trans`<0.1,], aes(y=`log2fc`, x=b)) + 
#               geom_point() + 
#               geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.1]), label=paste0("r = ", r), hjust=0) +
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.1])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1])-0.15, y=min(all$`log2fc`[all$`qvalue_trans`<0.1])+0.25, label=paste0("n = ", length(all$`log2fc`[all$`qvalue_trans`<0.05])), hjust=0) +
#               theme_classic() )
#     dev.off()
#     
#     
#     
#     r <- round(cor(all$`log2fc`[all$`qvalue_trans`<0.1 & all$pval<0.1], (all$b[all$`qvalue_trans`<0.1 & all$pval<0.1])), 2)
#     p <- cor.test(all$`log2fc`[all$`qvalue_trans`<0.1 & all$pval<0.1], all$b[all$`qvalue_trans`<0.1 & all$pval<0.1])$p.value
#     
#     pdf(paste0('analysis/02_QTL_calling/plots/HHEX_',variant,'_',cell_type,variant2,'_cor_CRISPR_tensor_trans01.p01.pdf'), width = 5, height = 4)
#     print(  ggplot(all[all$`qvalue_trans`<0.1 & all$pval<0.1,], aes(y=`log2fc`, x=b)) + 
#               geom_point() + 
#               geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1 & all$pval<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.1 & all$pval<0.1]), label=paste0("r = ", r), hjust=0) +
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1 & all$pval<0.1])-0.15, y=max(all$`log2fc`[all$`qvalue_trans`<0.1 & all$pval<0.1])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
#               annotate("text", x=max(all$b[all$`qvalue_trans`<0.1 & all$pval<0.1])-0.15, y=min(all$`log2fc`[all$`qvalue_trans`<0.1 & all$pval<0.1])+0.25, label=paste0("n = ", length(all$`log2fc`[all$`qvalue_trans`<0.05 & all$pval<0.1])), hjust=0) +
#               theme_classic() )
#     dev.off()
    
    ## read CHIP 
    #chip_john <- readxl::read_xlsx('data/John.science.adh7699_tables_s1_to_s4/science.adh7699_table_s4.xlsx', skip = 2,range = 'Table S4B!A3:S20330')
    # filt_chip <- all[all$Gene %in% chip_john$Gene[chip_john$IKZF1_ChIP == 1],]
    # r <- round(cor(filt_chip$`Log2 Fold Change`, (filt_chip$b)), 2)
    # p <- cor.test(filt_chip$`Log2 Fold Change`, filt_chip$b)$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_transCHIP.pdf'), width = 5, height = 4)
    # print(  ggplot(filt_chip, aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(filt_chip$`Log2 Fold Change`)), hjust=0) +
    #           theme_classic() )
    # dev.off()
    # 
    # r <- round(cor(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], (filt_chip$b[filt_chip$`P-value`<0.05])), 2)
    # p <- cor.test(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05], filt_chip$b[filt_chip$`P-value`<0.05])$p.value
    # 
    # pdf(paste0('analysis/02_QTL_calling/plots/IKZF1_',variant,'_',cell_type,'_cor_CRISPR_tensor_transCHIP.p005.pdf'), width = 5, height = 4)
    # print(  ggplot(filt_chip[filt_chip$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
    #           geom_point() + 
    #           geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
    #           annotate("text", x=max(filt_chip$b)-0.15, y=max(filt_chip$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
    #           annotate("text", x=max(all$b)-0.15, y=min(all$`Log2 Fold Change`)+0.25, label=paste0("n = ", length(filt_chip$`Log2 Fold Change`[filt_chip$`P-value`<0.05])), hjust=0) +
    #           theme_classic() )
    # dev.off()
    
 #  }    
 # }

#### Get distribution pvals and b from genes with chip info ####
### No differences #####

library(ggpubr)

genes <- c('GFI1B','IKZF1')
for (gene in genes) {
  snps <- paste0(trans_genes_summary$Var1[grep(gene,trans_genes_summary$Var2)],'-',trans_genes_summary$Var2[grep(gene,trans_genes_summary$Var2)])
  for (var in snps) {
    variant <- gsub('.*_','',var)
    cell_type <- gsub('-.*','',var)
    trans_all <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/",cell_type,"_mean_mx///trans_res/OneK1K_",cell_type,".",gene,".['",variant,"'].trans_qtl_pairs.allpairs.csv"))
    
    chip_john <- readxl::read_xlsx('data/John.science.adh7699_tables_s1_to_s4/science.adh7699_table_s4.xlsx', skip = 2,range = 'Table S4B!A3:S20330')
    coln <- paste0(gene,'_ChIP')
    trans_all$ChIP <- 'No'
    trans_all$ChIP[trans_all$phenotype_id %in% chip_john$Gene[chip_john[,coln] == 1]] <- 'Yes'
    
    pdf(paste0('analysis/02_QTL_calling/plots/',gene,'_',variant,'_',cell_type,'_violin_b_ChIP.pdf'), width = 5, height = 4)
    print(  ggplot(trans_all, aes(y=abs(b), x=ChIP, fill=ChIP)) + 
              #geom_density(adjust=1.5, alpha=.4) +
              geom_violin(outlier.shape = NA)+
              geom_boxplot(width=.1)+
              annotate("text", y=max(abs(trans_all$b))-0.01, x=1, label=paste0("n = ", sum(trans_all$ChIP == 'No')), hjust=0, col='red') +
              annotate("text", y=max(abs(trans_all$b))-0.15, x=1, label=paste0("n = ", sum(trans_all$ChIP == 'Yes')), hjust=0,col='blue') +
              stat_compare_means()+
              theme_classic() )
    dev.off()
    
  }
}

### Correlate trans networks same variant different cell-types ######


# load libraries ggplot2 and ggally 
library(ggplot2) 
library(GGally) 

#### load data #####
gene <- 'RPL23A'
snps <- paste0(trans_genes_summary$Var1[grep(gene,trans_genes_summary$Var2)],'-',trans_genes_summary$Var2[grep(gene,trans_genes_summary$Var2)])

trans_all_B <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/B_mean_mx//trans_res/OneK1K_B.RPL23A.['17:26936608:C:T'].trans_qtl_pairs.allpairs.csv"))
trans_all_cd4 <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/CD4_T_mean_mx///trans_res/OneK1K_CD4_T.RPL23A.['17:26936608:C:T'].trans_qtl_pairs.allpairs.csv"))
trans_all_nk <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/NK_mean_mx///trans_res/OneK1K_NK.RPL23A.['17:26936608:C:T'].trans_qtl_pairs.allpairs.csv"))

trans_all_cd4 <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/Mono_mean_mx////trans_res/OneK1K_Mono.RPL23A.['17:27065211:T:C'].trans_qtl_pairs.allpairs.csv"))
trans_all_nk <- read.delim(paste0("analysis/02_QTL_calling/lvl1Res/other_T_mean_mx////trans_res/OneK1K_other_T.RPL23A.['17:27065211:T:C'].trans_qtl_pairs.allpairs.csv"))

genes_common <- trans_all_B$phenotype_id[trans_all_B$phenotype_id %in% trans_all_cd4$phenotype_id]
genes_common <- trans_all_nk$phenotype_id[trans_all_nk$phenotype_id %in% trans_all_cd4$phenotype_id]
genes_common <-  genes_common[genes_common %in% trans_all_nk$phenotype_id & genes_common != 'RPL23A']

trans_all_B_s <- trans_all_B[trans_all_B$pval<0.1,]
trans_all_cd4_s <- trans_all_cd4[trans_all_cd4$pval<0.1,]
trans_all_nk_s <- trans_all_nk[trans_all_nk$pval<0.1,]
genes_common_sig <- trans_all_B_s$phenotype_id[trans_all_B_s$phenotype_id %in% trans_all_cd4_s$phenotype_id & trans_all_B_s$phenotype_id %in% trans_all_nk_s$phenotype_id]
genes_common_sig <- trans_all_cd4_s$phenotype_id[trans_all_cd4_s$phenotype_id %in% trans_all_nk_s$phenotype_id]
genes_union_sig <-unique(c(trans_all_B_s$phenotype_id,trans_all_cd4_s$phenotype_id,trans_all_nk_s$phenotype_id))
genes_union_sig <-unique(c(trans_all_cd4_s$phenotype_id,trans_all_nk_s$phenotype_id))
genes_union_sig <- genes_union_sig[genes_union_sig %in% genes_common]
genes_common_sig <- genes_common_sig[genes_common_sig %in% genes_common]

rownames(trans_all_B) <- trans_all_B$phenotype_id
rownames(trans_all_cd4) <- trans_all_cd4$phenotype_id
rownames(trans_all_nk) <- trans_all_nk$phenotype_id

sample_data <- data.frame(trans_all_B[genes_common,'b'], trans_all_cd4[genes_common,'b'], trans_all_nk[genes_common,'b']) 
sample_data <- data.frame(trans_all_B[genes_common_sig,'b'], trans_all_cd4[genes_common_sig,'b'], trans_all_nk[genes_common_sig,'b']) 
sample_data <- data.frame(trans_all_B[genes_union_sig,'b'], trans_all_cd4[genes_union_sig,'b'], trans_all_nk[genes_union_sig,'b']) 
colnames(sample_data) <- c('B','CD4','NK')

sample_data <- data.frame(trans_all_cd4[genes_common,'b'], trans_all_nk[genes_common,'b']) 
sample_data <- data.frame(trans_all_cd4[genes_common_sig,'b'], trans_all_nk[genes_common_sig,'b']) 
sample_data <- data.frame(trans_all_cd4[genes_union_sig,'b'], trans_all_nk[genes_union_sig,'b']) 
colnames(sample_data) <- c('Mono','Other T')

# create pairs plot 
pdf(paste0('analysis/02_QTL_calling/plots/',gene,'_','17:27065211:T:C','_','comp_celltypesTensor_union01.pdf'), width = 5, height = 4)
ggpairs( sample_data )
dev.off()

