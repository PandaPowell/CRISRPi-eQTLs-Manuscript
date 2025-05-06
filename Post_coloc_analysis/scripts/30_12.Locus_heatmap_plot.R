rm(list=ls())

library(data.table)
library(scales)
library(arrow)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(tidyplots)

setDTthreads(10)
setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

# Load cres
cres_w_grnas = fread("cres_with_grnas.txt")
cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt")

# Load gold standard causal gene list
gold_genes = fread("rare_genes_list.txt")

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

# load power calculations
gtex = fread("power_results/GTEx.power.results.txt")

x = egenes$grna_target[9]
y = egenes$gene_name[9]

plot.heatmap = function(x, y){
  
  cat("Plotting heatmap for", x , "\n")
  
  # filter distance dataframe to snp
  eqtl.cat.region = eqtl.cat.distances %>% 
    filter(grna_target == x) %>% 
    mutate(tss_distance = position-start) %>%
    dplyr::select(grna_target, study = eqtl, gene_name, pvalue, tss_distance) %>%
    filter(study != "fibroblast")
  
  sc.eqtl.region = sc.eqtl.distances %>% 
    filter(grna_target == x, !cell_type %in% c("Other_T_cells", "Other_cells")) %>%
    dplyr::select(grna_target, study = cell_type, gene_name, pvalue = pval_nominal, tss_distance = distance)
  
  gtex.region = gtex.distance %>% 
    filter(grna_target == x) %>%
    mutate(study = "GTEx") %>%
    dplyr::select(grna_target, study,gene_name, pvalue = pval_nominal, tss_distance)
  
  mage.region = mage.distance %>%
    filter(grna_target == x) %>%
    mutate(study = "MAGE") %>%
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
    dplyr::select(grna_target, study, gene_name, pvalue, tss_distance)
  
  dis_rank = rbind(eqtl.region,crispr.region) %>%
    filter(is.na(tss_distance) == F) %>%
    arrange(pvalue) %>%
    distinct(gene_name,.keep_all=T) %>%
    mutate(dis_rank = rank(tss_distance*-1),
           zero_rank = rank(abs(tss_distance)))
  
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
    
    
    if(min.p.up <= 0.10 & min.p.down <= 0.10){
      
      heatmap_data = heatmap_data[,colnames(heatmap_data) %in% 
                                    c(dis_rank$gene_name[ dis_rank$tss_distance >= dis_rank$tss_distance[dis_rank$pvalue == min.p.down] & 
                                                     dis_rank$tss_distance <= dis_rank$tss_distance[dis_rank$pvalue == min.p.up]],y) ]
    } else if(min.p.up <= 0.10) {
      
      # list of genes upstream up until min p.val
      down.stream.genes = dis_rank$gene_name[dis_rank$tss_distance < 0 & dis_rank$pvalue < 0.10]
      up.stream.genes = dis_rank$gene_name[dis_rank$tss_distance > 0 & dis_rank$tss_distance <= dis_rank$tss_distance[dis_rank$pvalue == min.p.up]]
      heatmap_data = heatmap_data[,colnames(heatmap_data) %in% c(up.stream.genes,down.stream.genes,y)]
      
    } else if(min.p.down <= 0.10){
      
      # list of genes upstream up until min p.val
      down.stream.genes = dis_rank$gene_name[dis_rank$tss_distance >= dis_rank$tss_distance[dis_rank$pvalue == min.p.down]]
      up.stream.genes = dis_rank$gene_name[dis_rank$tss_distance > 0 & dis_rank$pvalue < 0.10]
      heatmap_data = heatmap_data[,colnames(heatmap_data) %in% c(up.stream.genes,down.stream.genes,y)]
      
    } else{
      
      heatmap_data = heatmap_data[,colnames(heatmap_data) %in% c(dis_rank$gene[dis_rank$zero_rank<=20],y)]
    
      }

  }
  
  rank.list= dis_rank[gene_name %in% colnames(heatmap_data)] %>% arrange(dis_rank)
  
  heatmap_long <- heatmap_data %>%
    rownames_to_column("study") %>%
    pivot_longer(cols = -study, names_to = "gene_name", values_to = "pvalue") %>%
    filter(gene_name %in% dis_rank$gene_name) %>%
    left_join(dis_rank[,c("gene_name","tss_distance")], "gene_name")
  
  heatmap_long$pval_filled <- ifelse(is.na(heatmap_long$pvalue), -Inf, heatmap_long$pvalue)
  heatmap_long$gene_name = factor(heatmap_long$gene_name, levels = rank.list$gene_name) # set order of x-axis
  
  if("CRISPRi" %in% heatmap_long$study){
    heatmap_long$study <- factor(heatmap_long$study, levels = c("CRISPRi",unique(heatmap_long$study)[unique(heatmap_long$study) != "CRISPRi"]))
  }
  
  # set center of plot for title
  center_position <- which(levels(heatmap_long$gene_name) == rank.list$gene_name[rank.list$zero_rank == 1])/length(rank.list$gene_name)
  half.square = (1/length(rank.list$gene_name))/2 # size of half a square
  
  if (!dir.exists("plots/heatmaps")) {
    
    print("No directory")
    
    dir.create("plots/heatmaps", recursive = TRUE)
  }
  
  # Plotting the heatmap
  png(paste0("plots/heatmaps/", x, ".png"), width = 14, height = 4, units = "in", res = 300)
  
  plot = ggplot(heatmap_long, aes(x = gene_name, y = study, fill = pval_filled)) +
    geom_tile(aes(color = is.na(pvalue)), show.legend = TRUE) + # Show legend for z-scores
    scale_fill_gradientn(colors = c("white", "white", "red"),
                         values = rescale(c(0, 2, 10)), na.value = "gainsboro", oob = squish, limits = c(0, 10)) +  
    scale_color_manual(values = c("TRUE" = "gainsboro", "FALSE" = NA), guide = "none") + # Exclude color legend for NA tiles
    theme_cowplot() +
    labs(title = x, x = "Gene", y = "", fill = "-log10 p-value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = ifelse(levels(heatmap_long$gene_name) == y, "gold", 
                                                                       ifelse(levels(heatmap_long$gene_name) == dis_rank$gene_name[dis_rank$zero_rank == 1],"blue",
                                                                              "black"))),
      plot.title = element_text(hjust = center_position-half.square, size = 12, margin = margin(b = 25)) # Center the title
    )
  
  print(plot)
  
  dev.off()
  
  print("Done")
  
}

plot.heatmap(egenes$grna_target[10], egenes$gene_name[10])

lapply(1:10, function(x) plot.heatmap(egenes$grna_target[x], egenes$gene_name[x]))
