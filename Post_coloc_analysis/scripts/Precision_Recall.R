rm(list=ls())
options(bitmapType="cairo")
# Load libraries
library(ggplot2)
library(pROC)
library(dplyr)
library(data.table)
library(tidyr)
library(cowplot)
library(ggplot2)
library(biomaRt)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis")

# Get Ensembl v99 protein coding genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "chromosome_name", "transcription_start_site"),
                   filters = "biotype",
                   values = "protein_coding",
                   mart = ensembl)
protein_coding_genes <- gene_list$ensembl_gene_id

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
annot.GRC37 <- annot[order(annot$chr_number, annot$start), c("chr_number", "start", "ensembl_id", "gene_name")]

# Filter for protein-coding genes (ensure protein_coding_genes is defined)
# Example: protein_coding_genes <- c("ENSG00000139618", "ENSG00000157764") # Replace with actual list
annot.GRC37 <- annot.GRC37[annot.GRC37$ensembl_id %in% protein_coding_genes, ]

# Load Gold Standard Genes
gold_genes <- fread("cis_gold_genes.txt")
# Split grna_target into separate rows for better filtering
gold_genes <- gold_genes[, .(grna = unlist(strsplit(grna_target, ","))), by = gold_gene] %>%
  mutate(cre_gene = paste(grna, gold_gene, sep = "_"))

# Load cgenes
cgenes = fread("processed_data/cres_with_grnas.txt") %>%
  filter(significant == 1, ensembl_id %in% protein_coding_genes)

# Load egenes
egenes = fread("processed_data/cres_with_grna_eqtls_interval.txt") %>%
  filter(ensembl_id %in% protein_coding_genes) %>%
  distinct(target_gene, .keep_all = T)

egenes_pph4 = egenes %>% 
  filter(PP.H4.abf > 0.9)

# Load ABC gene targets
all_genes_grc38 <- fread("data/ensembl_gencode_gene_list.txt")
cres_w_grnas_abc <- fread("ABC/cres_w_grnas_ABC_intersections.txt")[, -5] %>%
  left_join(all_genes_grc38[, c("ensembl_id", "gene_name")], by = c("TargetGene" = "gene_name")) %>%
  filter(!is.na(ensembl_id), ensembl_id %in% protein_coding_genes) %>%
  arrange(desc(ABC.Score)) %>%
  distinct(target_site, .keep_all = TRUE) %>%
  mutate(target_gene = paste(target_site, ensembl_id, sep = "_"))

# load Hi-C gene targets
# --- Evaluate ABC predictions ---
hic = fread("Hi_C/K562/K562.hg19.AllInteractions.SP4.FDR0.1.txt")
hic_cres_w_grnas = fread("Hi_C/K562/cres_w_grnas_HiC_interactions.bed") %>%
  left_join(hic, "InteractorID") %>% 
  left_join(gene_list, by = c("RefSeqName" = "external_gene_name")) %>%
  filter(ensembl_gene_id %in% protein_coding_genes) %>%
  mutate(target_gene = paste0(grna_target, "_", ensembl_gene_id)) %>%
  distinct(target_gene, .keep_all = T)

# Establish all CRE-gene links within 1MBP
cres <- fread("processed_data/cres_with_grnas.txt") %>%
  dplyr::select(chr,grna_target,grna_pos) %>%
  distinct() %>% 
  filter(grna_target %in% gold_genes$grna)

library(readxl)

twas = read_excel("data/Rowland_TWAS_supplemental_tables.xlsx", sheet = 2) %>%
  left_join(gene_list[,c("ensembl_gene_id", "external_gene_name")], by = c("gene_name" = "external_gene_name")) %>%
  filter(ensembl_gene_id %in% protein_coding_genes) %>% distinct()

# Find cis TWAS genes and take most significant
find.cis.genes = function(chromosome, pos, grna_target){
  
  cis.genes = twas %>%
    filter(chr == chromosome) %>%
    filter(start_pos >= (pos-1e6) & start_pos <= (pos+1e6) | end_pos >= (pos-1e6) & end_pos <= (pos+1e6)) %>%
    mutate(grna_target = grna_target)
  
  if(nrow(cis.genes) < 1){
    return(tibble())  # Return empty tibble
  } else{
    cis.max = cis.genes[cis.genes$log10_regenie_p == max(cis.genes$log10_regenie_p),] # max twas in region
    return(cis.max)
  }
  
}

cis.twas.tmp = lapply(1:nrow(cres), function(x) {
  find.cis.genes(chromosome = cres$chr[x], 
                 pos = cres$grna_pos[x],
                 grna_target = cres$grna_target[x])
})

cis.twas.genes = bind_rows(cis.twas.tmp) %>%
  filter(marginal_significant == 1) %>%
  mutate(target_gene = paste(grna_target,ensembl_gene_id, sep = "_"))

# Find closest gene to all GWAS CREs
find.cis.genes = function(chr, pos, grna_target){
  
  cis.genes = annot.GRC37 %>%
    filter(chr_number == chr & start >= (pos - 1e6) & start <= (pos + 1e6)) %>%
    mutate(distance = abs(start-pos)) %>%
    mutate(grna_target = grna_target)
  
  return(cis.genes)
}

cis.genes_tmp = lapply(1:nrow(cres), function(x) {
  find.cis.genes(chr = cres$chr[x], 
                    pos = cres$grna_pos[x],
                    grna_target = cres$grna_target[x])
})

cis.genes = bind_rows(cis.genes_tmp) %>%
  mutate(cre_gene = paste(grna_target, ensembl_id, sep = "_")) %>%
  group_by(grna_target) %>%
  mutate("Closest gene" = ifelse(distance == min(distance),1,0)) %>% # assign closest gene
  ungroup() %>%
  mutate(gold_gene = ifelse(cre_gene %in% gold_genes$cre_gene,1,0),
         cGenes = ifelse(cre_gene %in% cgenes$target_gene,1,0),
         "eGenes (H4>0.5)" = ifelse(cre_gene %in% egenes$target_gene,1,0),
         "cGenes or eGenes" = ifelse(cre_gene %in% cgenes$target_gene | cre_gene %in% egenes_pph4$target_gene,1,0),
         "eGenes (H4>0.9)" = ifelse(cre_gene %in% egenes_pph4$target_gene,1,0),
         "ABC-Max" = ifelse(cre_gene %in% cres_w_grnas_abc$target_gene,1,0),
         "Hi-C" = ifelse(cre_gene %in% hic_cres_w_grnas$target_gene,1,0),
         "TWAS" = ifelse(cre_gene %in% cis.twas.genes$target_gene,1,0))

# Function to compute performance metrics
evaluate_method <- function(method_name) {
  
  cat(sum(cis.genes[method_name]), method_name)
  
  conf_matrix <- table(
    gold = cis.genes$gold_gene, 
    pred = cis.genes[[method_name]]
  )
  
  print(conf_matrix)
  
  # Handle missing cases
  TP <- ifelse(!is.na(conf_matrix["1", "1"]), conf_matrix["1", "1"], 0)
  FP <- ifelse(!is.na(conf_matrix["0", "1"]), conf_matrix["0", "1"], 0)
  FN <- ifelse(!is.na(conf_matrix["1", "0"]), conf_matrix["1", "0"], 0)
  TN <- ifelse(!is.na(conf_matrix["0", "0"]), conf_matrix["0", "0"], 0)
  
  precision <- TP / (TP + FP)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  cat(sprintf("[%s]\nPrecision: %.3f\nRecall: %.3f\nSpecificity: %.3f\n\n",
              method_name, precision, sensitivity, specificity))
  
  return(data.frame(Method = method_name,
                    Precision = precision,
                    Recall = sensitivity,
                    Specificity = specificity))
}

# Initialize results data frame
method_metrics <- data.frame()
method_metrics <- rbind(method_metrics, evaluate_method("cGenes"))
method_metrics <- rbind(method_metrics, evaluate_method("eGenes (H4>0.5)"))
method_metrics <- rbind(method_metrics, evaluate_method("eGenes (H4>0.9)"))
method_metrics <- rbind(method_metrics, evaluate_method("cGenes or eGenes"))
method_metrics <- rbind(method_metrics, evaluate_method("ABC-Max"))
method_metrics <- rbind(method_metrics, evaluate_method("Hi-C"))
method_metrics <- rbind(method_metrics, evaluate_method("Closest gene"))
method_metrics <- rbind(method_metrics, evaluate_method("TWAS"))
method_metrics$f1 = f1 = 2 * (method_metrics$Precision * method_metrics$Recall) / (method_metrics$Precision + method_metrics$Recall)

library(ggrepel)
library(colorspace)

colors <- c(
  "cGenes" = "#d35e60",
  "eGenes (H4>0.5)" = "#76c0c1",
  "eGenes (H4>0.9)" = "#76c0c1",
  "cGenes or eGenes" = "#ed111a",
  "Hi-C" = "#ebebeb",
  "Closest gene" = '#c9c9c9',
  "ABC-Max" = '#adadad',
  "TWAS" = "#a18376"
)

svg("plots/interval/gold_standard_genes_hic_barplots.svg", width = 6, height = 6)

ggplot(method_metrics, aes(x = Recall, y = Precision)) +
  geom_point(
    aes(color = Method, fill = Method),
    size = 7, alpha = 0.8, 
    shape = 21 # This is a dot with both border (color) and fill.
  ) +
  # Add auto-positioned text
  geom_text_repel(
    aes(label = Method),
    color = "black",
    size = 16/.pt, # font size 9 pt
    point.padding = 0.1, 
    box.padding = 0.6,
    min.segment.length = 0,
    max.overlaps = 1000,
    seed = 7654 # For reproducibility reasons
  ) +
  scale_color_manual(
    name = NULL, # it's one way to omit the legend title
    values = darken(colors, 0.3) # dot borders are a darker than the fill
  ) +
  scale_fill_manual(
    name = NULL,
    values = colors
  ) +
  # Add labels and customize axes
  scale_x_continuous(
    name = "Recall",
    limits = c(0, 0.10),
    breaks = c(0.02, 0.04, 0.06, 0.08,0.10),
    expand = c(0, 0) # This removes the default padding on the ends of the axis
  ) + 
  scale_y_continuous(
    name = "Precision",
    limits = c(0, 0.25),
    expand = c(0, 0)
  ) +
  theme_cowplot() +  # Base font size increased
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

dev.off()

svg("plots/interval/precision_recall_plot.svg", width = 6, height = 6)

ggplot(method_metrics, aes(x = Precision, y = Recall)) +
  geom_point(
    aes(color = Method, fill = Method),
    size = 7, alpha = 0.8,
    shape = 21
  ) +
  geom_text_repel(
    aes(label = Method),
    color = "black",
    size = 16/.pt,
    point.padding = 0.1,
    box.padding = 0.6,
    min.segment.length = 0,
    max.overlaps = 1000,
    seed = 7654
  ) +
  scale_color_manual(name = NULL, values = darken(colors, 0.3)) +
  scale_fill_manual(name = NULL,  values = colors) +
  # Precision on x
  scale_x_continuous(
    name   = "Precision",
    limits = c(0, 0.25),
    expand = c(0, 0)
  ) +
  # Recall on y
  scale_y_continuous(
    name   = "Recall",
    limits = c(0, 0.10),
    breaks = c(0.02, 0.04, 0.06, 0.08, 0.10),
    expand = c(0, 0)
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16)
  )
dev.off()
