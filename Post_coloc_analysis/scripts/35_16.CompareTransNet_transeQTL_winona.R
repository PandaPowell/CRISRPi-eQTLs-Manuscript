############################
# Author: Winona Oliveros
# Date: 2024/01/11 
# Compare trans network with trans-eQTLs 
############################
library(data.table)

setwd("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K")

## Read Gasperini ##
cis_data <- fread('analysis/01_ExploringData/Cis_genes_trans_network_Gasperini.txt')
head(cis_data)

##trans networks ##
trans_network <- fread("/gpfs/commons/groups/lappalainen_lab/jmorris/230525_GasperiniTransNetworks/data/231011_GasperiniTransNetwork.txt") 
head(trans_network)
trans_network <- trans_network[trans_network$sig_trans==1,]

## trans results ###
all_res_trans <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/'), pattern = '*trans*')
  if (length(files)==0) {
    all_res_trans[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    #res_df$Cell.type <- cells
    all_res_trans[[cells]] <- res_df
  }
}

all_res_trans <- do.call('rbind.data.frame', all_res_trans)

all_res_lvl2_trans <- list()

for (cells in cell_types$V1) {
  files <- list.files(path = paste0('analysis/02_cisQTL/lvl2Res/', cells, '_mean_mx/'), pattern='*trans*')
  if (length(files)==0) {
    all_res_lvl2[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_cisQTL/lvl2Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    #res_df$Cell.type <- cells
    all_res_lvl2_trans[[cells]] <- res_df
  }
}

all_res_lvl2_trans <- do.call('rbind.data.frame', all_res_lvl2_trans)


## cis results ## 
#### now get results for specific calling
all_res <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/'), pattern = paste0(cells,'\\.independent.*Filt*'))
  if (length(files)==0) {
    all_res[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    if (nrow(res_df)<1) {
      all_res[[cells]] <- NA
    } else {
    res_df$Cell.type <- cells
    all_res[[cells]] <- res_df
    }
  }
}

all_res <- do.call('rbind.data.frame', all_res)

# number_of_cisQTLs <- t(sapply(c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T'), function(cell) 
#   if (is.null(nrow(all_res[[cell]]))) {
#     0
#   }else {nrow(all_res[[cell]])}
# )
# )

## lvl2 
cell_types <- read.delim('analysis/02_QTL_calling//cell_types_lvl2.txt', header = F)

all_res_lvl2 <- list()

for (cells in cell_types$V1) {
  files <- list.files(path = paste0('analysis/02_QTL_calling/lvl2Res/', cells, '_mean_mx/'), paste0(cells,'\\.independent.*Filt*'))
  if (length(files)==0) {
    all_res_lvl2[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_QTL_calling/lvl2Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    if (nrow(res_df)<1) {
      all_res_lvl2[[cells]] <- NA
    } else {
      res_df$Cell.type <- cells
      all_res_lvl2[[cells]] <- res_df
    }
  }
}

all_res_lvl2 <- do.call('rbind.data.frame', all_res_lvl2)

## gene 1 HHEX
cis_data_g <- cis_data[cis_data$hgnc_symbol_closest == 'HHEX',]
cis_data_g

trans_g <- trans_network[trans_network$gRNA_target %in% cis_data_g$gRNA_target,]
trans_g <- trans_g[trans_g$hgnc_symbol != 'HHEX',]

cis_variants_lvl1 <- all_res[all_res$phenotype_id=='HHEX',]
cis_variants_lvl2 <- all_res_lvl2[all_res_lvl2$phenotype_id=='HHEX',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

trans_genes_lvl1$phenotype_id %in% trans_g$hgnc_symbol
trans_genes_lvl2$phenotype_id %in% trans_g$hgnc_symbol

### plot genes
eqtls_interest <- c('10:94429467:C:T:NK:MGAT4A','10:94429467:C:T:NK:ATP6V1F',
                    '10:94429467:C:T:NK:LTBP1','10:94429467:C:T:NK:ADD2',
                    '10:94429467:C:T:NK:H2AFV', '10:94429467:C:T:NK:GMPPA',
                    '10:94429467:C:T:NK:GYPE', '10:94429467:C:T:NK:GSTO1',
                    '10:94429467:C:T:NK:UBTF', '10:94429467:C:T:NK:RPL21')

eqtls_interest <- c('10:94429467:C:T:NK:C1orf162','10:94429467:C:T:NK:CD2',
                    '10:94429467:C:T:NK:CD160','10:94429467:C:T:NK:TMEM123',
                    '10:94429467:C:T:NK:ESD', '10:94429467:C:T:NK:CMC1',
                    '10:94429467:C:T:NK:SORBS2', '10:94429467:C:T:NK:DOCK5',
                    '10:94429467:C:T:NK:GLIPR2')
for (eqtl in eqtls_interest) {
  
  info <- unlist(str_split(eqtl, ":"))
  
  #expression_bed <- paste0("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_cisQTL/AverageExprBed/Paper/OneK1K_980_samples_",info[[5]],".bed.gz")
  expression_bed <- paste0("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_cisQTL/Residuals/Granular/",info[[5]],".exprs_residuals.tsv")
  genotypes <- paste0('data/FinalGenotypePLINK/chr',info[[1]],'_001.bed')
  file_bim <- paste0('data/FinalGenotypePLINK/chr',info[[1]],'_001.bim')
  file_fam <- paste0('data/FinalGenotypePLINK/chr',info[[1]],'_001.fam')
  
  # read annotation tables
  bim <- read_bim(file_bim)
  fam <- read_fam(file_fam)
  
  gen <- read_bed(genotypes, bim$id, gsub('-','_',fam$id))  
  head(gen)  
  
  # read expr
  expr <- read.delim(expression_bed)
  head(colnames(expr))
  
  ## get snp 
  snp_id <- paste0(info[[1]],':',info[[2]],':',info[[3]],':',info[[4]])
  snp <- as.data.frame((gen[snp_id,]))
  colnames(snp) <- c('variant')
  snp$ind <- rownames(snp)
  
  ## get expr 
  #expr_g <- as.data.frame(t(expr[expr$gene_name == info[[6]], c(5:ncol(expr))]))
  expr_g <- as.data.frame(t(expr[info[[6]],]))
  colnames(expr_g) <- c('expression')
  expr_g$ind <- gsub('X','',rownames(expr_g))
  
  ## merge
  all <- merge(snp, expr_g)
  head(all)
  all$variant <- as.factor(all$variant)
  
  Means <- all %>% group_by(variant) %>% 
    dplyr::summarize(Avg = median(expression))
  
  library(dplyr)
  # compute lower and upper whiskers for each group
  ylims <- all %>%
    group_by(as.factor(variant)) %>%
    summarise(q1 = quantile(expression, 1/4), q3 = quantile(expression, 3/4)) %>%
    ungroup() %>%
    #get lowest q1 and highest q3
    summarise(lowq1 = min(q1), highq3 = max(q3))
  
  
  pdf(paste0('analysis/02_cisQTL/plots/',info[[6]],'_',info[[2]],'_eQTL_res.transegenes.HHEX.pdf'), width = 6, height = 4)
  print(ggplot(all, aes(x=as.factor(variant), y=expression, fill=as.factor(variant))) + 
          geom_violin(outlier.shape = NA)+
          geom_boxplot(width=.1)+
          # geom_boxplot(width=.7,outlier.shape = NA,notch = T,
          #            notchwidth = 0.10) + 
          scale_fill_brewer(palette="RdBu") + ggtitle(paste0(snp_id,' ',info[[6]]))+
          geom_line(data = Means, mapping = aes(x = variant, y = Avg, group = 1)) +
          theme_classic()) + coord_cartesian(ylim = as.numeric(ylims)*1.5)
  dev.off()
  
}

### get all trans results #####
all_res_trans <- list()

for (cells in c('B','CD4_T','CD8_T','DC','Mono','NK','other', 'other_T')) {
  files <- list.files(path = paste0('analysis/02_cisQTL/lvl1Res/', cells, '_mean_mx/'), pattern = '*trans.*allpairs*')
  if (length(files)==0) {
    all_res[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_cisQTL/lvl1Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    #res_df$Cell.type <- cells
    all_res_trans[[cells]] <- res_df
  }
}

all_res_trans <- do.call('rbind.data.frame', all_res_trans)
all_res_trans$cell_type <- gsub('\\..*','',rownames(all_res_trans))

all_res_lvl2_trans <- list()

for (cells in cell_types$V1) {
  files <- list.files(path = paste0('analysis/02_cisQTL/lvl2Res/', cells, '_mean_mx/'), pattern='*trans.*allpairs*')
  if (length(files)==0) {
    all_res_lvl2[[cells]] <- NA
  } else {
    res <- lapply(files, function(file) read.csv(paste0('analysis/02_cisQTL/lvl2Res/', cells, '_mean_mx/', file), sep = '\t'))
    res_df <- do.call('rbind.data.frame', res)
    #res_df$Cell.type <- cells
    all_res_lvl2_trans[[cells]] <- res_df
  }
}

all_res_lvl2_trans <- do.call('rbind.data.frame', all_res_lvl2_trans)
all_res_lvl2_trans$cell_type <- gsub('\\..*','',rownames(all_res_lvl2_trans))

## Read Johns trans network results ###
trans_network_john <- readxl::read_xlsx('data/science.adh7699_tables_s1_to_s4/science.adh7699_table_s4.xlsx', skip = 2)

## GFI1B ##

trans_g <- trans_network_john[grep('GFI1B', trans_network_john$`Network CRE`),]
trans_g <- as.data.frame(trans_g)
table(trans_g$`SNP ID`)
snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))

cis_variants_lvl1 <- all_res[!is.na(all_res$X) & all_res$phenotype_id=='GFI1B',]
cis_variants_lvl2 <- all_res_lvl2[!is.na(all_res_lvl2$X) & all_res_lvl2$phenotype_id=='GFI1B',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
for (variant in snps_trans_crispr) {
  for (cell_type in unique(trans_genes_lvl2$cell_type)) {
    trans_crispr <- trans_g[trans_g$`SNP ID` == gsub(':','-',variant),]
    all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_genes_lvl2[trans_genes_lvl2$cell_type==cell_type,], by.x='Gene', by.y='phenotype_id')
    
    r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05]), 2)
    p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
    
    pdf(paste0('analysis/02_cisQTL/plots/GFI1B_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
    print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
    r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
    p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
    
    pdf(paste0('analysis/02_cisQTL/plots/GFI1B_',variant,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
    print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
              geom_point() + 
              geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
              annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
              theme_classic() )
    dev.off()
    
  }
  
}

## IKZF1 ##

trans_g <- trans_network_john[grep('IKZF1', trans_network_john$`Network CRE`),]
trans_g <- as.data.frame(trans_g)
table(trans_g$`SNP ID`)
snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))

cis_variants_lvl1 <- all_res[!is.na(all_res$X) & all_res$phenotype_id=='IKZF1',]
cis_variants_lvl2 <- all_res_lvl2[!is.na(all_res_lvl2$X) & all_res_lvl2$phenotype_id=='IKZF1',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
for (variant in snps_trans_crispr) {
  for (cell_type in unique(trans_genes_lvl1$cell_type)) {
    for (variant2 in unique(trans_genes_lvl1$variant_id[trans_genes_lvl1$cell_type == cell_type])) {
      trans_crispr <- trans_g[trans_g$`SNP ID` == gsub(':','-',variant),]
      all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_genes_lvl1[trans_genes_lvl1$cell_type==cell_type & trans_genes_lvl1$variant_id==variant2,], by.x='Gene', by.y='phenotype_id')
      
      r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05]), 2)
      p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/IKZF1_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
      print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05]), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
      r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
      p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/IKZF1_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
      print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
    
    }
  }
  
}

## HHEX ##

trans_g <- trans_network_john[grep('HHEX', trans_network_john$`Network CRE`),]
trans_g <- as.data.frame(trans_g)
table(trans_g$`SNP ID`)
snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))

cis_variants_lvl1 <- all_res[!is.na(all_res$X) & all_res$phenotype_id=='HHEX',]
cis_variants_lvl2 <- all_res_lvl2[!is.na(all_res_lvl2$X) & all_res_lvl2$phenotype_id=='HHEX',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
for (variant in snps_trans_crispr) {
  for (cell_type in unique(trans_genes_lvl1$cell_type)) {
    for (variant2 in unique(trans_genes_lvl1$variant_id[trans_genes_lvl1$cell_type == cell_type])) {
      trans_crispr <- trans_g[trans_g$`SNP ID` == gsub(':','-',variant),]
      all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_genes_lvl1[trans_genes_lvl1$cell_type==cell_type & trans_genes_lvl1$variant_id==variant2,], by.x='Gene', by.y='phenotype_id')
      
      r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05]), 2)
      p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/HHEX_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
      print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05]), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
      r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
      p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/HHEX_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
      print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
    }
  }
  
}


## RUNX1 ##

trans_g <- trans_network_john[grep('RUNX1', trans_network_john$`Network CRE`),]
trans_g <- as.data.frame(trans_g)
table(trans_g$`SNP ID`)
snps_trans_crispr <- gsub('-',':',unique(trans_g$`SNP ID`))

cis_variants_lvl1 <- all_res[!is.na(all_res$X) & all_res$phenotype_id=='RUNX1',]
cis_variants_lvl2 <- all_res_lvl2[!is.na(all_res_lvl2$X) & all_res_lvl2$phenotype_id=='RUNX1',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
for (variant in snps_trans_crispr) {
  for (cell_type in unique(trans_genes_lvl1$cell_type)) {
    for (variant2 in unique(trans_genes_lvl1$variant_id[trans_genes_lvl1$cell_type == cell_type])) {
      trans_crispr <- trans_g[trans_g$`SNP ID` == gsub(':','-',variant),]
      all <- merge(trans_crispr[,c("Gene","Log2 Fold Change","P-value")], trans_genes_lvl1[trans_genes_lvl1$cell_type==cell_type & trans_genes_lvl1$variant_id==variant2,], by.x='Gene', by.y='phenotype_id')
      
      r <- round(cor(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05]), 2)
      p <- cor.test(all$`Log2 Fold Change`[all$`P-value`<0.05], all$b[all$`P-value`<0.05])$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/HHEX_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
      print(  ggplot(all[all$`P-value`<0.05,], aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05]), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b[all$`P-value`<0.05])-0.05, y=max(all$`Log2 Fold Change`[all$`P-value`<0.05])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
      r <- round(cor(all$`Log2 Fold Change`, all$b), 2)
      p <- cor.test(all$`Log2 Fold Change`, all$b)$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/HHEX_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
      print(  ggplot(all, aes(y=`Log2 Fold Change`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b)-0.15, y=max(all$`Log2 Fold Change`)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
    }
  }
  
}

## RPL23A 17:27065211:T:C ##

cis_data_g <- cis_data[cis_data$hgnc_symbol_closest == 'RPL23A',]
cis_data_g

trans_network <- fread('/gpfs/commons/groups/lappalainen_lab/jmorris/230525_OneK1K/data/231011_GasperiniTransNetwork.txt') 
head(trans_network)
trans_g <- trans_network[trans_network$gRNA_target %in% cis_data_g$gRNA_target,]
trans_g <- trans_g[trans_g$hgnc_symbol != 'RPL23A',]

table(trans_g$gRNA_target)
snps_trans_crispr <- gsub('-',':',unique(trans_g$gRNA_target))

cis_variants_lvl1 <- all_res[!is.na(all_res$X) & all_res$phenotype_id=='RPL23A',]
cis_variants_lvl2 <- all_res_lvl2[!is.na(all_res_lvl2$X) & all_res_lvl2$phenotype_id=='RPL23A',]

trans_genes_lvl1 <- all_res_trans[all_res_trans$variant_id %in% cis_variants_lvl1$variant_id,]
trans_genes_lvl2 <- all_res_lvl2_trans[all_res_lvl2_trans$variant_id %in% cis_variants_lvl2$variant_id,]

# trans_genes_lvl1$phenotype_id %in% trans_g$Gene
# trans_genes_lvl2$phenotype_id %in% trans_g$Gene

### plot 
for (variant in snps_trans_crispr) {
  for (cell_type in unique(trans_genes_lvl2$cell_type)) {
    for (variant2 in unique(trans_genes_lvl2$variant_id[trans_genes_lvl2$cell_type == cell_type])) {
      trans_crispr <- trans_g[trans_g$gRNA_target == variant,]
      all <- merge(trans_crispr[,c("hgnc_symbol","log2fc","qvalue_trans")], trans_genes_lvl2[trans_genes_lvl2$cell_type==cell_type & trans_genes_lvl2$variant_id==variant2,], by.x='hgnc_symbol', by.y='phenotype_id')
      
      r <- round(cor(all$`log2fc`[all$qvalue_trans<0.05], all$b[all$qvalue_trans<0.05], method = 'spearman'), 2)
      p <- cor.test(all$`log2fc`[all$qvalue_trans<0.05], all$b[all$qvalue_trans<0.05], method = 'spearman')$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/RPL23A_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.p005.pdf'), width = 5, height = 4)
      print(  ggplot(all[all$qvalue_trans<0.05,], aes(y=`log2fc`, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b[all$qvalue_trans<0.05])-0.05, y=max(all$`log2fc`[all$qvalue_trans<0.05]), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b[all$qvalue_trans<0.05])-0.05, y=max(all$`log2fc`[all$qvalue_trans<0.05])-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
      r <- round(cor(as.numeric(all$log2fc), all$b, method = 'spearman'), 2)
      p <- cor.test(as.numeric(all$log2fc), all$b, use = "complete.obs" , method = 'spearman')$p.value
      
      pdf(paste0('analysis/02_cisQTL/plots/RPL23A_',variant,'_',variant2,'_',cell_type,'_cor_CRISPR_tensor_trans.pdf'), width = 5, height = 4)
      print(  ggplot(all, aes(y=log2fc, x=b)) + 
                geom_point() + 
                geom_smooth(method="lm", color="red", fill="#69b3a2", se=TRUE) + 
                annotate("text", x=max(all$b)-0.15, y=max(all$log2fc), label=paste0("r = ", r), hjust=0) +
                annotate("text", x=max(all$b)-0.15, y=max(all$log2fc)-0.25, label=paste0("p = ", round(p, 3)), hjust=0) +
                theme_classic() )
      dev.off()
      
    }
  }
  
}



