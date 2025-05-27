rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(scales)
library(tidyverse)
library(cowplot)
library(gridExtra)

setDTthreads(8)
setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load cres
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

# Load gold standard causal gene list
gold_genes = fread("cis_gold_genes.txt", header = T)
colnames(gold_genes)[2] = "rare_genes"

# filter eQTLs and crispri overlapping gold standard gene list
cgenes = cres_w_grnas %>% 
  filter(ensembl_id %in% gold_genes$rare_genes, significant == T) %>% 
  distinct(grna_target, .keep_all = T)

egenes = cres_w_grnas_egene %>% 
  filter(ensembl_id %in% gold_genes$rare_genes) %>% 
  distinct(grna_target, .keep_all = T)

# Load eqtl distances
eqtl.cat.distances = fread("eQTL_gene_distances/eqtl_cat_gene_distances.txt")
sc.eqtl.distances = fread("eQTL_gene_distances/Onek1k_gene_distance.txt")
gtex.distance = fread("eQTL_gene_distances/GTEx_gene_distance.txt")
mage.distance = fread("eQTL_gene_distances/MAGE_gene_distance.txt")

# x = egenes$grna_target[5]
# y = egenes$gene_name[5]

plot.heatmap = function(x, y){
  
  cat("Plotting heatmap for", x , "\n")
  
  # filter distance dataframe to snp
  eqtl.cat.region = eqtl.cat.distances %>% 
    filter(grna_target == x) %>% 
    mutate(tss_distance = position-start) %>%
    dplyr::select(grna_target, study = eqtl, gene_name, pvalue, tss_distance) %>%
    filter(study != "fibroblast")
  
  eqtl.cat.region$study[eqtl.cat.region$study == "blood"] = "Twins UK whl blood"
  
  sc.eqtl.region = sc.eqtl.distances %>% 
    filter(grna_target == x, !cell_type %in% c("Other_T_cells", "Other_cells")) %>%
    dplyr::select(grna_target, study = cell_type, gene_name, pvalue = pval_nominal, tss_distance = distance)
  
  gtex.region = gtex.distance %>% 
    filter(grna_target == x) %>%
    mutate(study = "GTEx whl blood") %>%
    dplyr::select(grna_target, study,gene_name, pvalue = pval_nominal, tss_distance)
  
  mage.region = mage.distance %>%
    filter(grna_target == x) %>%
    mutate(study = "MAGE LCL") %>%
    dplyr::select(grna_target, study, gene_name=geneSymbol, pvalue = pval_nominal, tss_distance)
  
  eqtl.region = rbind(eqtl.cat.region,sc.eqtl.region,gtex.region,mage.region) %>%
    mutate(pvalue = (pvalue/1e-03)*0.10) %>% # scale p-value to a line with crispr
    mutate(pvalue = ifelse(pvalue > 1, 1, pvalue),
           id = paste0(study,"_",gene_name)) %>% 
    arrange(pvalue) %>%
    distinct(id, .keep_all = T) %>% dplyr::select(-id)
  
  crispr.region = cres_w_grnas %>% 
    filter(grna_target == x) %>%
    mutate(study = "CRISPRi") %>%
    dplyr::select(grna_target, study, gene_name, pvalue, tss_distance) %>%
    distinct(gene_name, .keep_all = T)
  
  if(x %in% c("7:50427982", "chr11:5268102-5269244")){
    # Load gene positions
    # Load Gencode annotation data
    annot_file <- "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/Gencode/gencode.v33lift37.GRCh38.genes.gtf"
    # Read the annotation file
    annot <- read.table(annot_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    # Keep only genes from chr1-22
    annot <- annot[annot$V1 %in% paste0("chr", 1:22), ]
    annot <- annot[annot$V3 == "gene", ]
    # Extract Ensembl gene ID
    annot$ensembl_id <- sapply(annot$V9, function(x) {
      unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
    })
    # Extract gene name
    annot$gene_name <- sapply(annot$V9, function(x) {
      unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[4]), "[.]"))[1]
    })
    # Add start (TSS -1) and end (TSS) based on strand
    annot$start <- ifelse(annot$V7 == "+", annot$V4 - 1, annot$V5 - 1)
    annot$end <- ifelse(annot$V7 == "+", annot$V4, annot$V5)
    # Extract chromosome number
    annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
    # Sort and select relevant columns
    annot.GRC37 <- as.data.table(annot[order(annot$chr_number, annot$start), c("chr_number", "start", "ensembl_id", "gene_name")])
    K562_bulk = fread("/gpfs/commons/groups/lappalainen_lab/jmorris/210205_STINGseq-v2/data/K562_scRNAseq_bulkRNAseq.txt")
    
    if(x %in% c("chr11:5268102-5269244")){
      
      gene_region = annot.GRC37[chr_number == 11 & inrange(start, 5268102-1e+06, 5268102+1e+06),]
      gene_region = gene_region[gene_name %in% K562_bulk$GENE,]
      gene_region$tss_distance = gene_region$start - 5268102
      gene_region$study = "CRISPRi"
      gene_region$grna_target = "chr11:5268102-5269244"
      gene_region$pvalue = 1
      gene_region = gene_region[,c("grna_target","study","gene_name", "pvalue", "tss_distance")]
      colnames(gene_region) = c("grna_target","study","gene_name","pvalue", "tss_distance")
      crispr.region = rbind(crispr.region, gene_region[!gene_name %in% crispr.region$gene_name,])
    }
    
  }
  
  dis_rank = rbind(eqtl.region,crispr.region) %>%
    filter(is.na(tss_distance) == F) %>%
    arrange(pvalue) %>%
    distinct(gene_name,.keep_all=T) %>%
    mutate(dis_rank = rank(tss_distance*-1),
           zero_rank = rank(abs(tss_distance)))
  
  closest.gene = dis_rank$gene_name[dis_rank$zero_rank == 1]
  
  heatmap_data <- rbind(eqtl.region,crispr.region) %>%
    mutate(log_pval = -log(pvalue)) %>%
    dplyr::select(study, gene_name, log_pval) %>%
    filter(gene_name != "") %>%
    arrange(desc(log_pval)) %>%
    #distinct(study, .keep_all=T) %>%
    pivot_wider(names_from = gene_name, values_from = log_pval) %>%
    column_to_rownames("study")
  
  # remove rows with no sifnificant pvalues
  heatmap_data = heatmap_data[apply(heatmap_data, 1, max, na.rm=T) > -log(0.20),]
  
  # reduce size of heatmap by smallest p-value or <20 genes
  if(nrow(dis_rank) > 20){
    
    min.p.up = min(dis_rank$pvalue[dis_rank$zero_rank>20 & dis_rank$tss_distance > 0], na.rm = T)
    min.p.down = min(dis_rank$pvalue[dis_rank$zero_rank>20 & dis_rank$tss_distance < 0], na.rm = T)
    
    if(min.p.up <= 0.10 | min.p.down <= 0.10){
      
      if(sum(dis_rank$pvalue < 0.20) > 19){ # are there at least 20 values with p<0.20 
        
        heatmap_data = heatmap_data[, colnames(heatmap_data) %in% c(dis_rank$gene[dis_rank$pvalue < 0.20], y, closest.gene) ]
        
      } else {
        
        keep.genes = dis_rank$gene_name[dis_rank$pvalue < 0.20]
        heatmap_data = heatmap_data[,colnames(heatmap_data) %in% c(dis_rank$gene[dis_rank$zero_rank<15], keep.genes, y, closest.gene)]
        
      }
      
    } else{
      
      heatmap_data = heatmap_data[,colnames(heatmap_data) %in% c(dis_rank$gene[dis_rank$zero_rank<=20], y, closest.gene) ]
    
      }

  }
  
  #If CRISPRi not present make empty row
  if(!"CRISPRi" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "CRISPRi"
    # Load gene expression
    K562_bulk = fread("/gpfs/commons/groups/lappalainen_lab/jmorris/210205_STINGseq-v2/data/K562_scRNAseq_bulkRNAseq.txt")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% K562_bulk$GENE)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  #If GTEx not present make empty row
  if(!"GTEx whl blood" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "GTEx whl blood"
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  #If GTEx not present make empty row
  if(!"NK_cells" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "NK_cells"
    # Load gene expression
    gene_ex = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/AverageExprBed/Original/OneK1K_980_samples_NK.bed.gz")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% gene_ex$gene_name)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  #If GTEx not present make empty row
  if(!"CD8_T_cells" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "CD8_T_cells"
    # Load gene expression
    gene_ex = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/AverageExprBed/Original/OneK1K_980_samples_CD8_T.bed.gz")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% gene_ex$gene_name)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  if(!"CD4_T_cells" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "CD4_T_cells"
    # Load gene expression
    gene_ex = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/AverageExprBed/Original/OneK1K_980_samples_CD4_T.bed.gz")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% gene_ex$gene_name)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  if(!"Mono_cells" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "Mono_cells"
    # Load gene expression
    gene_ex = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/AverageExprBed/Original/OneK1K_980_samples_Mono.bed.gz")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% gene_ex$gene_name)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  if(!"B_cells" %in% rownames(heatmap_data)){
    crispri_row <- as.data.frame(t(rep(0, ncol(heatmap_data))))
    colnames(crispri_row) <- colnames(heatmap_data)
    rownames(crispri_row) = "B_cells"
    # Load gene expression
    gene_ex = fread("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/AverageExprBed/Original/OneK1K_980_samples_B.bed.gz")
    # non-expressed genes NA
    crispri_row[, !(colnames(crispri_row) %in% gene_ex$gene_name)] <- NA
    heatmap_data = rbind(heatmap_data, crispri_row)
  }
  
  rank.list= dis_rank[gene_name %in% colnames(heatmap_data)] %>% arrange(dis_rank)
  
  heatmap_long <- heatmap_data %>%
    rownames_to_column("study") %>%
    pivot_longer(cols = -study, names_to = "gene_name", values_to = "pvalue") %>%
    filter(gene_name %in% dis_rank$gene_name) %>%
    left_join(dis_rank[,c("gene_name","tss_distance")], "gene_name")
  
  heatmap_long$pval_filled <- as.numeric(ifelse(is.na(heatmap_long$pvalue), -Inf, heatmap_long$pvalue))
  heatmap_long$gene_name = factor(heatmap_long$gene_name, levels = rank.list$gene_name) # set order of x-axis
  
  base.studies = c("CRISPRi","GTEx whl blood", "CD4_T_cells", "CD8_T_cells", "B_cells","Mono_cells", "NK_cells")
  
  # Set CRISPRi & GTEX to the bottom of the heatmap
  heatmap_long$study <- factor(heatmap_long$study, 
                               levels = c("CRISPRi","GTEx whl blood", "CD4_T_cells", "CD8_T_cells", "B_cells","Mono_cells","NK_cells",
                                          unique(heatmap_long$study)[!unique(heatmap_long$study) %in% base.studies]))
  
  # set center of plot for title
  center_position <- which(levels(heatmap_long$gene_name) == closest.gene)/length(rank.list$gene_name)
  half.square = (1/length(rank.list$gene_name))/2 # size of half a square
  
  if (!dir.exists("plots/heatmaps")) {
    
    print("No directory")
    
    dir.create("plots/heatmaps", recursive = TRUE)
  }
  
  if(length(unique(heatmap_long$study)) > 9){
    
    # Plotting the heatmap
    svg(paste0("plots/heatmaps/", gsub(":","_",x), ".svg"), width = 14, height = length(unique(heatmap_long$study))*0.5, bg = "transparent")
    
    plot = ggplot(heatmap_long, aes(x = gene_name, y = study, fill = pval_filled)) +
      geom_tile(aes(color = is.na(pvalue)), show.legend = TRUE) + # Show legend for z-scores
      scale_fill_gradientn(colors = c("white", "white", "red"),
                           values = rescale(c(0, 1.4, 10)), na.value = "gainsboro", oob = squish, limits = c(0, 10)) +  
      scale_color_manual(values = c("TRUE" = "gainsboro", "FALSE" = NA), guide = "none") + # Exclude color legend for NA tiles
      theme_cowplot() +
      labs(title = gsub("chr","",x), x = "Gene", y = "", fill = "-log10 p-value") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, color = ifelse(levels(heatmap_long$gene_name) == y, "gold", 
                                                                         ifelse(levels(heatmap_long$gene_name) == dis_rank$gene_name[dis_rank$zero_rank == 1],"blue",
                                                                                "black"))),
        plot.title = element_text(hjust = center_position-half.square, size = 18, margin = margin(b = 25)),# Center the title
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16)
      )
    
  } else{
    
    # Plotting the heatmap
    svg(paste0("plots/heatmaps/", gsub(":","_",x), ".svg"), width = 14, height = ifelse(length(unique(heatmap_long$study)) < 4, 4, length(unique(heatmap_long$study))), bg = "transparent")
    
    plot = ggplot(heatmap_long, aes(x = gene_name, y = study, fill = pval_filled)) +
      geom_tile(aes(color = is.na(pvalue)), show.legend = TRUE) + # Show legend for z-scores
      scale_fill_gradientn(colors = c("white", "white", "red"),
                           values = rescale(c(0, 1.4, 10)), na.value = "gainsboro", oob = squish, limits = c(0, 10)) +  
      scale_color_manual(values = c("TRUE" = "gainsboro", "FALSE" = NA), guide = "none") + # Exclude color legend for NA tiles
      theme_cowplot() +
      labs(title = gsub("chr","",x), x = "Gene", y = "", fill = "-log10 p-value") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = ifelse(levels(heatmap_long$gene_name) == y, "gold", 
                                                                         ifelse(levels(heatmap_long$gene_name) == dis_rank$gene_name[dis_rank$zero_rank == 1],"blue",
                                                                                "black"))),
        plot.title = element_text(hjust = center_position-half.square, size = 18, margin = margin(b = 25)),# Center the title
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16)
      )
    
  }
  
  print(plot)
  dev.off()
  
  gold = unique(heatmap_long$gene_name[heatmap_long$gene_name == y])
  closest = unique(dis_rank$gene_name[dis_rank$zero_rank == 1])
  
  if (gold == closest) {
    print("gold is closest")
    return(as.character(closest))  # return the gene name
  } else {
    return(NA)  # or NULL if you prefer
  }
}

#plot.heatmap(egenes$grna_target[3], egenes$gene_name[3])

library(purrr)
results <- pmap_chr(
  list(egenes$grna_target, egenes$gene_name),
  plot.heatmap
)

results <- na.omit(gold_genes)
length(unique(gold_genes))

results2 <- pmap_chr(
  list(cgenes$grna_target, cgenes$gene_name),
  plot.heatmap
)

results2 <- na.omit(gold_genes)
length(unique(gold_genes))/23


plots <- list()

for (i in seq_along(cgenes$grna_target)) {
  plots[[i]] <- plot.heatmap(cgenes$grna_target[i], cgenes$gene_name[i])
}

row_counts <- sapply(plots, function(plot) {
  data <- ggplot_build(plot)$data[[1]]
  length(unique(data$y))
})

ordered_indices <- order(row_counts, decreasing = TRUE)  # Use decreasing = FALSE for ascending order
plots <- plots[ordered_indices]

for (i in seq(1, length(plots), by = 5)) {
  # Select 5 plots (or fewer if at the end)
  group <- plots[i:min(i + 4, length(plots))]
  
  # Determine row counts
  row_counts <- sapply(group, function(plot) {
    data <- ggplot_build(plot)$data[[1]]
    length(unique(data$y))
  })
  
  # Set dynamic height
  row_sum = sum(row_counts, na.rm = T)
  dynamic_height <- ifelse(row_sum < 31, row_sum + 10, row_sum * 0.5)
  
  # Combine plots
  combined_plot <- plot_grid(plotlist = group, ncol = 1)
  
  # Save combined plot
  ggsave(
    filename = paste0("plots/heatmaps/cgenes", ceiling(i / 5), ".pdf"),
    plot = combined_plot,
    width = 8.27,
    height = dynamic_height,
    units = "in",
    dpi = 300, limitsize = F
  )
}

#plot.heatmap(egenes$grna_target[35], egenes$gene_name[35])
