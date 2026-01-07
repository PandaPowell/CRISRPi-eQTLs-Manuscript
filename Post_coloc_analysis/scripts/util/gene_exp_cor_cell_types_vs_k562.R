rm(list = ls())
library(data.table)
library(tidyverse)
library(Matrix)
library(gridExtra)

# Load gencode annotations
gencode = fread("data/gencode.v33lift37.annotation.gtf", skip = 5)
gencode = gencode[V3 == "gene"]
#Extract `gene_id` and `gene_name` into new columns
gencode[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
gencode[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
gencode[, tss := ifelse(V7 == "+",V4,V5)]
gencode = gencode[,.(V1,tss,gene_id,gene_name)]
gencode[, ensembl_id := str_split_fixed(gencode$gene_id,"\\.",2)[,1]]
colnames(gencode)[1] = "chr"

# Load K562 gene expression count matrix
K562_bulk = fread("/gpfs/commons/groups/lappalainen_lab/jmorris/210205_STINGseq-v2/data/K562_scRNAseq_bulkRNAseq.txt") %>%
  mutate(ensembl_id = str_split_fixed(ENSG, "\\.",2)[,1]) %>% left_join(gencode, "ensembl_id")

names=c("B","Mono","NK","DC","CD4_T","CD8_T")

plot_list <- list()  # Create an empty list to store plots

for (i in seq_along(names)){
  
  cell.name <- names[i]  # Access the current cell type name
  
  print(cell.name)
  
  #Load example count matrices
  dir = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/cell_specific_rds/lvl1/"
  count.matrix = readRDS(paste0(dir,cell.name,"_raw_counts.RDS"))
  #Load the associated annotation data frame with the necessary meta data columns 
  #for cell type and donor information
  annot.df = readRDS(paste0(dir,cell.name,"_meta_data.rds"))
  # How many cells contain the gene
  cells.with.molecules = Matrix::rowSums(count.matrix > 0)
  # total molecules
  total.molecules<-Matrix::rowSums(count.matrix)
  # Make df
  count.df = data.frame(gene_name = names(cells.with.molecules),
                        cells.with.molecules = cells.with.molecules,
                        total.molecules = total.molecules,
                        molecules.per.cell = total.molecules/cells.with.molecules)
  
  count.df$molecules.per.cell[count.df$molecules.per.cell == "NaN"] = 0

  # Add ensemble ids
  count.df = count.df %>% left_join(gencode, "gene_name")
  
  merge_df = K562_bulk %>% left_join(count.df,"ensembl_id")
  merge_df$molecules.per.cell[is.na(merge_df$molecules.per.cell)] = 0
  
  spearmans = print(cor(merge_df$MoleculesPerCell, merge_df$molecules.per.cell, method = "spearman"))
  
  # Create the plot
  plot <- ggplot(merge_df, aes(x = MoleculesPerCell, y = molecules.per.cell)) +
    geom_point(color = "grey") +
    theme_cowplot() +
    labs(title = cell.name, x = "K562 mean count per cell", y = paste0(cell.name," mean count per cell")) +
    annotate("text", x = 50, y = 60, label = paste0("cor: ", round(spearmans, 2)), size = 5, hjust = 0.5)
  # Store the plot in the plot list
  plot_list[[i]] <- plot  # Use the index i to store plots
  
}

png("plots/gene_exp_cor_scatter.png", width = 15, height = 10,units = "in", res = 300)
# Use grid.arrange to print all plots together
do.call(grid.arrange, c(plot_list, nrow = 2))
dev.off()
