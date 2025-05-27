rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(ggtranscript)
library(purrr)
library(ggrepel)

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

x = cgenes$grna_target[1]
y = cgenes$gene_name[1]

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
  
  rank.list= dis_rank[gene_name %in% colnames(heatmap_data)] %>% arrange(dis_rank)
  
  heatmap_long <- heatmap_data %>%
    rownames_to_column("study") %>%
    pivot_longer(cols = -study, names_to = "gene_name", values_to = "pvalue") %>%
    filter(gene_name %in% dis_rank$gene_name) %>%
    left_join(dis_rank[,c("grna_target","gene_name","tss_distance")], "gene_name") %>%
    rowwise() %>%
    mutate(variant_pos = as.numeric(strsplit(gsub(".*?:", "", grna_target), "-")[[1]][1]),
           gene_start = variant_pos - tss_distance) %>%
    ungroup()
  
  heatmap_long$pval_filled <- as.numeric(ifelse(is.na(heatmap_long$pvalue), 0, heatmap_long$pvalue))
  
  gold = unique(heatmap_long$gene_name[heatmap_long$gene_name == y])
  closest = unique(dis_rank$gene_name[dis_rank$zero_rank == 1])

  variant_pos <- unique(heatmap_long$variant_pos)
  
  genes <- heatmap_long %>%
    distinct(gene_name, .keep_all = T) %>%
    select(start = gene_start, gene_name) %>%
    mutate(
      direction = case_when(
        start < variant_pos ~ "left",
        start > variant_pos ~ "right",
        TRUE ~ "center"
      )
    )
  
  genes <- genes %>%
    group_by(direction) %>%
    mutate(distance = abs(variant_pos-start)) %>%
    arrange(distance, .by_group = TRUE) %>%
    mutate(
      layout_pos = case_when(
        direction == "left" ~ -row_number(),
        direction == "right" ~ row_number(),
        TRUE ~ 0
      )
    ) %>%
    ungroup()
  
  svg(paste0("plots/figure_plots/",variant_pos,"_",gold,".svg"), width = 4, height = 0.5)
  
  # Wrap ggplot in print()!
  print(
    ggplot(genes, aes(x = layout_pos)) +
      # Gene blocks
      geom_range(
        aes(xstart = layout_pos - 0.2, xend = layout_pos + 0.2, y = 0),
        height = 0.005,
        fill = "black"
      ) +
      # Variant at center
      geom_point(aes(x = 0, y = 0), shape = 21, fill = "red", size = 1) +
      # Gene labels below
      geom_text(
        aes(x = layout_pos - 0.4, y = -0.05, label = gene_name),
        angle = 45,
        vjust = 0,
        hjust = 0.5,
        size = 2
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      theme_void() +
      ggtitle("")
  )
  
  dev.off()
  
  if (gold == closest) {
    print("gold is closest")
    return(as.character(closest))  # return the gene name
  } else {
    return(NA)  # or NULL if you prefer
  }
}

results <- pmap_chr(
  list(cgenes$grna_target, cgenes$gene_name),
  plot.heatmap
)

genes = data.frame(seqnames = rep(1,5),
                start = seq(1,1000000,100000),
                end = seq(1,1000000,100000) + 10000,
                strand = "+",
                type = "exons",
                gene_name = "Gene",
                transcript_name = "transcript_name",
                transcript_biotype = "transcript_biotype")

eqtl_links = data.frame(seqnames = rep(1,10),
                      start = rep(500000,10), # variant position
                      end = genes$start, # tss position
                      strand = "+",
                      n_studies = c(2,5,10, 30, 1,4, 5,1 ,1 ,1)) # thickness of lines/ p-value

crispri_links = data.frame(seqnames = 1,
                             start = 500000, # variant position
                             end = 700000, # tss position
                             strand = "+",
                             n_studies = c(1)) # thickness of lines/ p-value

genes %>% ggplot(aes(xstart = start,xend = end, y = "")) + 
  geom_range(fill = "black", height = 0.025) +
  geom_junction(
    data = eqtl_links,
    aes(size = n_studies),
    junction.orientation = "alternating",
    junction.y.max = 0.2,
    angle = 90,
    ncp = 100, 
    colour = "#00BFC4") +
  geom_junction(
    data = crispri_links,
    aes(size = n_studies),
    junction.orientation = "alternating",
    junction.y.max = 0.2,
    angle = 90,
    ncp = 100, 
    colour = "#F8766D") +
  ylab("") + xlab("Position on chr 1") +
  scale_size_continuous(range = c(0.2, 2), guide = "none") +
  theme_minimal()


data.frame(seqnames = rep(1,5),
           start = seq(1,1000000,100000),
           end = seq(1,1000000,100000) + 10000,
           strand = "+",
           type = "exons",
           gene_name = "Gene",
           transcript_name = "transcript_name",
           transcript_biotype = "transcript_biotype")


library(ggrepel)

genes = heatmap_long %>%
  select(start = gene_start, gene_name) %>%
  mutate(end = start + 10000, 
         strand = "+", 
         type = "exons", 
         transcript_name = "transcript_name",
         transcript_biotype = "transcript_biotype")

eqtl_links = heatmap_long %>%
  filter(pvalue > 3, is.na(pvalue) == F) %>%
  mutate(strand = "+") %>%
  rename(start = variant_pos, end = gene_start) %>%
  mutate(line_thickness = pvalue/10)
eqtl_links$seqnames = rep(1,nrow(eqtl_links))

genes %>% distinct(gene_name, .keep_all = T) %>%
  ggplot(aes(xstart = start, xend =end, y = "")) + 
  geom_range(fill = "black", height = 0.025) +
  geom_text_repel(
    aes(x = start, label = gene_name),
    y = 0,
    size = 3,
    direction = "y",
    nudge_y = 0.05,
    segment.color = "black"
  )

  geom_junction(
    data = eqtl_links,
    aes(xstart = start, xend = end, y = "", size = line_thickness, colour = study),
    junction.orientation = "alternating",
    junction.y.max = 0.1,
    angle = 90,
    ncp = 50) +
  ylab("") + xlab("Position on chr 1") +
  scale_size_continuous(range = c(0.2, 2), guide = "none") +
  theme_minimal()
  