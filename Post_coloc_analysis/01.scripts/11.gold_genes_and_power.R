rm(list=ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(readxl)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")
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

# Download Mendelian genes file
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867420309995-mmc7.xlsx", "mendelian_genes.xlsx", mode = "wb")

# Load and process Mendelian genes data
curated_genes <- read_xlsx("mendelian_genes.xlsx") %>%
  rename(gene_name = Gene_Symbol_HGNC) %>%
  left_join(annot.GRC37, by = "gene_name") %>%
  drop_na() %>%
  distinct(ensembl_id, .keep_all = TRUE) %>%
  dplyr::select(gene_name, chr_number, start, ensembl_id)

# Load UK Biobank WES blood cell trait genes
ukbb_genes <- fread("data/genebass/blood_trait_burden_test_genes.csv", sep = ",")

# Filter for pLoF annotations, remove duplicates, and join with gene annotations
ukbb_genes <- ukbb_genes %>%
  filter(grepl("pLoF", annotation)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  rename(ensembl_id = gene_id) %>%
  left_join(annot.GRC37, by = "ensembl_id") %>%
  drop_na() %>%
  dplyr::select(gene_name, chr_number, start, ensembl_id)

#-----------------------------#
# Combine UKBB and Mendelian curated genes
#-----------------------------#
gw.rare.genes <- bind_rows(ukbb_genes, curated_genes)

length(gw.rare.genes$ensembl_id)

fwrite(gw.rare.genes, "gw_gold_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
#-----------------------------#
# Load CRISPR eQTL and GRNA target data
#-----------------------------#
cres_w_grnas_egene <- fread("cres_with_grna_eqtls.txt")
cres_w_grnas <- fread("cres_with_grnas.txt")

# Significant gRNA targets and genes
cgene.targets <- cres_w_grnas %>%
  filter(significant == 1) %>%
  dplyr::select(grna_target, cgene = ensembl_id)

egene.targets = cres_w_grnas_egene %>%
  distinct(target_gene, .keep_all = T) %>%
  dplyr::select(grna_target, egene = ensembl_id)
  
all.grna.targets = cres_w_grnas %>%
  distinct(grna_target, .keep_all = TRUE)

# Find closest gene to all GWAS CREs
find.closest.gene = function(chr, pos, grna_target){
  
  cis.genes = annot.GRC37 %>%
    filter(chr_number == chr & start >= (pos - 1e6) & start <= (pos + 1e6)) %>%
    mutate(distance = abs(start-pos)) %>%
    mutate(grna_target = grna_target)
  
  closest.gene = cis.genes[cis.genes$distance == min(cis.genes$distance),]
  
  return(closest.gene)
}

closest.genes = lapply(1:nrow(all.grna.targets), function(x) {
  find.closest.gene(chr = all.grna.targets$chr[x], 
                    pos = all.grna.targets$grna_pos[x],
                    grna_target = all.grna.targets$grna_target[x])
})

closest.genes = bind_rows(closest.genes)

fwrite(closest.genes, "closest_gene_to_GWAS_CRE.txt", row.names = F, quote = F)

# Find all annotated genes with 1 Mb of gRNAs
find.cis <- function(chr, pos, grna_target) {
  cre.rare.genes <- gw.rare.genes %>%
    filter(chr_number == chr & start >= (pos - 1e6) & start <= (pos + 1e6)) %>%
    mutate(grna_target = grna_target)
  return(cre.rare.genes)
}

# Apply the function across all significant targets
results <- lapply(1:nrow(all.grna.targets), function(x) {
  find.cis(chr = all.grna.targets$chr[x], 
           pos = all.grna.targets$grna_pos[x],
           grna_target = all.grna.targets$grna_target[x])
})

# Combine results into one dataframe
# Consolidate duplicate genes and combine grna targets
gold_cis_genes <- bind_rows(results) %>%
  left_join(cgene.targets, by = "grna_target") %>%
  left_join(egene.targets, by = "grna_target") %>%
  rename(gold_gene = ensembl_id) %>%
  group_by(gene_name, chr_number, start, gold_gene) %>%
  summarise(
    grna_target = paste(unique(grna_target), collapse = ","),
    cgene = paste(unique(cgene), collapse = ","),
    egene = paste(unique(egene), collapse = ","),
    .groups = "drop"
  ) %>%
  # Calculate cgene_multiple excluding NA, empty strings, and genes equal to gold_gene
  mutate(
    cgene_multiple = mapply(function(cgenes, gold) {
      sum(str_trim(str_split(cgenes, ",")[[1]]) != "NA" &
            str_trim(str_split(cgenes, ",")[[1]]) != "" &
            str_trim(str_split(cgenes, ",")[[1]]) != gold)
    }, cgene, gold_gene),
    
    # Calculate egene_multiple with similar conditions
    egene_multiple = mapply(function(egenes, gold) {
      sum(str_trim(str_split(egenes, ",")[[1]]) != "NA" &
            str_trim(str_split(egenes, ",")[[1]]) != "" &
            str_trim(str_split(egenes, ",")[[1]]) != gold)
    }, egene, gold_gene)
  )

fwrite(gold_cis_genes[,c("grna_target", "gold_gene")], "cis_gold_genes.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#------------------------------------------#
# Filter cGenes and eGenes within 1 Mb of 
#------------------------------------------#
# Load Gold Standard Genes
gold_genes <- fread("cis_gold_genes.txt")

# Split grna_target into separate rows for better filtering
gold_genes <- gold_genes[, .(grna = unlist(strsplit(grna_target, ","))), by = gold_gene]

cgenes = cres_w_grnas %>%
  filter(grna_target %in% gold_genes$grna, significant == 1, ensembl_id %in% protein_coding_genes) %>%
  distinct(ensembl_id)

cgenes = unlist(cgenes, use.names = F)

egenes = cres_w_grnas_egene %>%
  filter(grna_target %in% gold_genes$grna, ensembl_id %in% protein_coding_genes) %>%
  distinct(ensembl_id)

egenes = unlist(egenes, use.names = F)

## Load Hi-c genes
hic = fread("Hi_C/K562/K562.hg19.AllInteractions.SP4.FDR0.1.txt")

intersect = fread("Hi_C/K562/cres_w_grnas_HiC_interactions.bed") %>%
  left_join(hic, "InteractorID") %>% 
  left_join(gene_list, by = c("RefSeqName" = "external_gene_name")) %>%
  filter(grna_target %in% gold_genes$grna) %>%
  distinct(ensembl_gene_id)

hic_genes = unlist(intersect, use.names = F)
#-----------------------------#
# Function: Calculate enrichment
#-----------------------------#
calculate_enrichment = function(mendel_genes){
  # Set the number of iterations
  n_iterations <- 10000
  
  # Initialize a vector to store the overlaps
  cgene_overlap_counts <- numeric(n_iterations)
  egene_overlap_counts <- numeric(n_iterations)
  
  # Calculate the overlap in each iteration
  set.seed(123) # for reproducibility
  for (i in 1:n_iterations) {
    # Randomly sample a set of genes with the same size as your GWAS list
    crandom_genes <- sample(protein_coding_genes, size = length(cgenes), replace = FALSE)
    erandom_genes = sample(protein_coding_genes, size = length(egenes), replace = FALSE)
    # Calculate the number of overlapping genes
    cgene_overlap_counts[i] <- sum(crandom_genes %in% mendel_genes)
    egene_overlap_counts[i] <- sum(erandom_genes %in% mendel_genes)
  }
  
  # Take the median of the random gene overlaps
  cgene_random_counts = median(cgene_overlap_counts)
  egene_random_counts = median(egene_overlap_counts)
  
  # Observed overlap
  cgene_observed_overlap <- sum(cgenes %in% mendel_genes)
  egene_observed_overlap <- sum(egenes %in% mendel_genes)
  
  a = cgene_observed_overlap # cgene in mendel set
  b = length(cgenes) - cgene_observed_overlap # cgenes not in mendel set
  c = cgene_random_counts
  d = length(crandom_genes) - cgene_random_counts
  
  cgene_table = matrix(c(a,b,c,d), nrow = 2)
  print(cgene_table)
  
  print(fisher.test(cgene_table, alternative = "greater"))
  
  a1 = egene_observed_overlap
  b1 = length(egenes) - egene_observed_overlap
  c1 = egene_random_counts
  d1 = length(erandom_genes) - egene_random_counts
  
  egene_table = matrix(c(a1,b1,c1,d1), nrow = 2)
  print(egene_table)
  
  print(fisher.test(egene_table, alternative = 'greater'))
  
  compare_table = matrix(c(a,b,a1,b1), nrow = 2)
  print(compare_table)
  
  print(fisher.test(compare_table))
}

calculate_enrichment(gold_cis_genes$gold_gene)

#-----------------------------#
# Function: Calculate enrichment (with hic_genes)
#-----------------------------#
calculate_enrichment = function(mendel_genes){
  # Set the number of iterations
  n_iterations <- 10000
  
  # Initialize vectors to store the overlaps
  cgene_overlap_counts <- numeric(n_iterations)
  egene_overlap_counts <- numeric(n_iterations)
  hicgene_overlap_counts <- numeric(n_iterations)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Calculate the overlap in each iteration
  for (i in 1:n_iterations) {
    # Randomly sample a set of genes with the same size as your gene lists
    crandom_genes <- sample(protein_coding_genes, size = length(cgenes), replace = FALSE)
    erandom_genes <- sample(protein_coding_genes, size = length(egenes), replace = FALSE)
    hicrandom_genes <- sample(protein_coding_genes, size = length(hic_genes), replace = FALSE)
    
    # Calculate the number of overlapping genes with the Mendelian set
    cgene_overlap_counts[i] <- sum(crandom_genes %in% mendel_genes)
    egene_overlap_counts[i] <- sum(erandom_genes %in% mendel_genes)
    hicgene_overlap_counts[i] <- sum(hicrandom_genes %in% mendel_genes)
  }
  
  # Take the median of the random gene overlaps
  cgene_random_counts <- median(cgene_overlap_counts)
  egene_random_counts <- median(egene_overlap_counts)
  hicgene_random_counts <- median(hicgene_overlap_counts)
  
  # Observed overlaps
  cgene_observed_overlap <- sum(cgenes %in% mendel_genes)
  egene_observed_overlap <- sum(egenes %in% mendel_genes)
  hicgene_observed_overlap <- sum(hic_genes %in% mendel_genes)
  
  # Construct contingency tables
  # cgenes
  a <- cgene_observed_overlap
  b <- length(cgenes) - a
  c <- cgene_random_counts
  d <- length(crandom_genes) - c
  cgene_table <- matrix(c(a, b, c, d), nrow = 2)
  print("cgene_table:")
  print(cgene_table)
  print(fisher.test(cgene_table, alternative = "greater"))
  
  # egenes
  a1 <- egene_observed_overlap
  b1 <- length(egenes) - a1
  c1 <- egene_random_counts
  d1 <- length(erandom_genes) - c1
  egene_table <- matrix(c(a1, b1, c1, d1), nrow = 2)
  print("egene_table:")
  print(egene_table)
  print(fisher.test(egene_table, alternative = "greater"))
  
  # hic_genes
  a2 <- hicgene_observed_overlap
  b2 <- length(hic_genes) - a2
  c2 <- hicgene_random_counts
  d2 <- length(hicrandom_genes) - c2
  hicgene_table <- matrix(c(a2, b2, c2, d2), nrow = 2)
  print("hicgene_table:")
  print(hicgene_table)
  print(fisher.test(hicgene_table, alternative = "greater"))
  
  # Optional: Compare observed overlaps between the groups
  compare_table <- matrix(c(a, b, a1, b1, a2, b2), nrow = 2)
  print("Observed overlap comparison table:")
  print(compare_table)
  print("Fisher test comparing cgenes vs egenes:")
  print(fisher.test(matrix(c(a, b, a1, b1), nrow = 2)))
  print("Fisher test comparing cgenes vs hic_genes:")
  print(fisher.test(matrix(c(a, b, a2, b2), nrow = 2)))
  print("Fisher test comparing egenes vs hic_genes:")
  print(fisher.test(matrix(c(a1, b1, a2, b2), nrow = 2)))
}

calculate_enrichment(gold_cis_genes$gold_gene)

#-----------------------------#
# Create and save bar plot
#-----------------------------#
data <- data.frame(
  Gene_type = factor(c("cGenes", "eGenes", "Hi-C"), levels = c("cGenes", "eGenes", "Hi-C")),
  Proportion = c(23 / length(cgenes), 53 / length(egenes), 34 / length(hic_genes))
)

colors <- c(
  "cGenes" = "#d35e60",
  "eGenes" = "#76c0c1",
  "Hi-C" = "#014d64"
)

svg("plots/figure_plots/gold_standard_genes_hic_barplots.svg", width = 4, height = 6)

ggplot(data, aes(x = Gene_type, y = Proportion, fill = Gene_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.5, position = position_dodge(width = 0.7)) +
  labs(title = "", x = "", y = "Proportion of genes") +
  scale_fill_manual(values = colors) +
  annotate("text", x = 1.5, y = 0.25, label = "P = 6.9e-03", size = 5, hjust = 0.5, family = "sans") +
  annotate("text", x = 2.5, y = 0.15, label = "P = 0.57", size = 5, hjust = 0.5, family = "sans") +
  theme_cowplot() +
  theme(legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(l = 20)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16)
  )

dev.off()

#-----------------------------#
# Load and process CRISPR power data
#-----------------------------#
stingseq.p1 <- fread("power_results/crispri_power/morris_large_power_per_pair.csv") %>%
  filter(effect_size == 0.85)
stingseq.p2 <- fread("power_results/crispri_power/morris_small_power_per_pair.csv") %>%
  filter(effect_size == 0.85)
stingseq.p <- bind_rows(stingseq.p1, stingseq.p2)

gasperini.p <- fread("power_results/crispri_power/gasperini_power_per_pair.csv") %>%
  filter(effect_size == 0.85) %>%
  mutate(cre_gene = paste0(cre_pert, "_", gene))

#-----------------------------#
# Load GRC38 annotations
#-----------------------------#
annot_file <- "data/gencode.v47.annotation.gtf.gz"
annot <- read.table(annot_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

annot <- annot %>%
  filter(V1 %in% paste0("chr", 1:22), V3 == "gene") %>%
  mutate(
    ensembl_id = sapply(V9, function(x) unlist(strsplit(sub("gene_id ", "", strsplit(x, ";")[[1]][1]), "[.]"))[1]),
    gene_name = sapply(V9, function(x) unlist(strsplit(sub(".*gene_name ", "", strsplit(x, ";")[[1]][3]), "[.]"))[1]),
    start = ifelse(V7 == "+", V4 - 1, V5 - 1),
    end = ifelse(V7 == "+", V4, V5),
    chr_number = as.numeric(sub("chr", "", V1))
  )

annot.GRC38 <- annot %>%
  arrange(chr_number, start) %>%
  dplyr::select(chr_number, start, ensembl_id, gene_name)

#-----------------------------#
# Process CRISPR encode data
#-----------------------------#
encode <- fread("CRISPR_data/NoGasperini_crispri_data.tsv") %>%
  mutate(target_site = paste0(chrom_GRC37, ":", chromStart_GRC37, "-", chromEnd_GRC37)) %>%
  rename(gene_name = measuredGeneSymbol) %>%
  left_join(annot.GRC38[, c("ensembl_id", "gene_name")], by = "gene_name") %>%
  dplyr::select(gene = ensembl_id, fraction_sig = PowerAtEffectSize15)

crispri.power <- bind_rows(encode, gasperini.p, stingseq.p) %>%
  arrange(desc(fraction_sig)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(gene %in% gold_cis_genes$gold_gene, !gene %in% cgenes, fraction_sig < 0.8)

## Number of gold-genes underpowered and not expressed in k562
sting_seq = fread("CRISPR_data/All_STING_seq_CREs.csv") %>%
  rename(SS_coord = `SNP Coordinates (hg19)`)

# Load gasperini data, filter to targeted enhancers
resample_results = fread("CRISPR_data/resampling_results.txt") %>% unique() %>% 
  filter(site_type == "DHS" & quality_rank_grna == "top_two", is.na(target_site.start) == F)

k562_genes = unique(c(sting_seq$`Ensembl ID`,resample_results$gene_id))

sum(gold_cis_genes$gold_gene %in% k562_genes) # expressed genes
not_expressed = gold_cis_genes$gold_gene[!gold_cis_genes$gold_gene %in% k562_genes] # not expressed genes
sum(gold_cis_genes$gold_gene %in% crispri.power$gene)# underpowered

not_expressed %in% crispri.power$gene
sum(!gold_cis_genes$gold_gene %in% k562_genes & !gold_cis_genes$gold_gene %in% crispri.power$gene)# underpowered and not expressed

#--------------------------------------#
# Process & combine eQTL power data    #
#--------------------------------------#
power.files = system("ls power_results/power.results.beta0.5*", intern = T)
eqtl.power = data.table()

for (i in 1:length(power.files)){
  
  name = str_split_fixed(pattern = "\\.",string = power.files[i], 5)[,4]
  name = gsub(pattern = "5", replacement = "",x = name)
  tmp = fread(power.files[i])
  tmp$study = name
  tmp$effect_size = 0.5
  eqtl.power = bind_rows(eqtl.power,tmp)
  
}

eqtl.power = eqtl.power %>% 
  arrange(desc(eQTL.power)) %>% 
  mutate(cre_gene = paste0(snp, "_",ensembl_id)) %>%
  distinct(cre_gene, .keep_all=T)

fwrite(eqtl.power, "power_results/eqtl.power.combined.txt", quote = F, row.names = F, col.names = T)

# sort by power and remove duplicated
eqtl.power = eqtl.power %>% 
  arrange(desc(eQTL.power)) %>% 
  distinct(ensembl_id, .keep_all=T) %>%
  filter(eQTL.power < 0.8)

## Determine the number of associated genes per CRE
cgenes_per_cre = cres_w_grnas %>% 
  filter(significant ==1) %>% 
  distinct(target_gene, .keep_all=T) %>% 
  group_by(grna_target) %>% 
  summarise(n = n())

egenes_per_cre = cres_w_grnas_egene %>% 
  distinct(target_gene, .keep_all=T) %>% 
  group_by(grna_target) %>% 
  summarise(n = n())

#-----------------------------#
# Create grid data and assign categories
#-----------------------------#
# Create a gene grid
gene_grid <- gold_cis_genes %>%
  mutate(
    egene = gold_gene %in% egenes,
    cgene = gold_gene %in% cgenes,
    crispri_underpowered = gold_gene %in% crispri.power$gene,
    eqtl_underpowered = gold_gene %in% eqtl.power$ensembl_id
  )
# Define the group names with priority order
group_names <- c("eQTL underpowered", "egene", "CRISPRi underpowered", "cgene")

# Assign hierarchical categories
data <- gene_grid %>%
  rowwise() %>%
  mutate(
    eqtl_underpowered = eqtl_underpowered & !egene,          # Override eQTL underpowered if egene is present
    crispri_underpowered = crispri_underpowered & !cgene,    # Override CRISPRi underpowered if cgene is present
    active_groups = list(group_names[c(eqtl_underpowered, egene, crispri_underpowered, cgene)]), # Active groups
    category = ifelse(length(active_groups) > 0, paste(active_groups, collapse = " & "), "None"), # Concatenate groups
  ) %>%
  mutate(gene_exp_crispri = ifelse(gold_gene %in% cres_w_grnas$ensembl_id, 1, 0),
         gene_exp_eqtl = ifelse(gold_gene %in% cres_w_grnas_egene$ensembl_id, 1, 0))

table = data %>%
  filter(cgene == T | egene == T) %>%
  group_by(category) %>%
  summarise(
    cgene_multiple_gt_0 = sum(cgene_multiple > 0 & egene_multiple == 0),
    egene_multiple_gt_0 = sum(egene_multiple > 0 & cgene_multiple == 0),
    both_gt_0 = sum(cgene_multiple > 0 & egene_multiple > 0),
    both_zero = sum(cgene_multiple == 0 & egene_multiple == 0),
    .groups = "drop"
  )

table

# CRISPRi Venn diagram
cgene_venn = data %>% 
  filter(cgene == F)

# another gene in locus & powered
another_cgene = sum(cgene_venn$cgene_multiple > 0 & cgene_venn$crispri_underpowered == F)
# underpowered/not expressed
egene_unpow = sum((cgene_venn$crispri_underpowered | cgene_venn$gene_exp_crispri == 0) & cgene_venn$cgene_multiple == 0)
# another gene in locus & underpowered
intersection = sum(cgene_venn$cgene_multiple > 0 & cgene_venn$crispri_underpowered == T)

another_cgene
egene_unpow
intersection

# CRISPRi Venn diagram
egene_venn = data %>% 
  filter(egene == F)

# another gene in locus & powered
another_egene = sum(egene_venn$egene_multiple > 0 & egene_venn$eqtl_underpowered == F)
# underpowered/not expressed
egene_unpow = sum((egene_venn$crispri_underpowered | egene_venn$gene_exp_crispri == 0) & egene_venn$cgene_multiple == 0)
# another gene in locus & underpowered
intersection = sum(cgene_venn$cgene_multiple > 0 & cgene_venn$crispri_underpowered == T)

another_cgene
egene_unpow
intersection

# genes not discovered & not expressed
sum(cgene_venn$gene_exp_crispri == 0)

sum(285,136,113)
sum(data$crispri_underpowered & data$gene_exp_crispri == 0)

## These numbers are used to create the Venn diagram

cat("Number of gold-genes underpowered in eQTL analysis", sum(data$eqtl_underpowered))

cat("Number of gold-genes underpowered in CRISPRi analysis", sum(data$crispri_underpowered))


cat("Overlap between underpowered & another egene", sum(gene_categories[grepl("eQTL underpowered",gene_categories$category),c(3:4)]))
cat("Overlap between underpowered & another egene", sum(gene_categories[grepl("CRISPRi underpowered",gene_categories$category),c(2,4)]))

# Number of cgenes
sum(gene_categories$both_gt_0) + sum(gene_categories$cgene_multiple_gt_0)

# Number of egenes
sum(gene_categories$both_gt_0) + sum(gene_categories$egene_multiple_gt_0)

cres_w_grnas = fread("cres_with_grnas.txt") %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id))

cres_w_grnas_egene = fread("cres_with_grna_eqtls.txt") %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id))

total_cres = unique(cres_w_grnas$grna_target)
grna_cres_w_cgenes = unique(cres_w_grnas$grna_target[cres_w_grnas$significant == 1])
target_cgenes = unique(cres_w_grnas$target_gene[cres_w_grnas$significant == 1])
cgenes = unique(cres_w_grnas$ensembl_id[cres_w_grnas$significant == 1])
grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
target_egenes = unique(cres_w_grnas_egene$target_gene)
egenes = unique(cres_w_grnas_egene$ensembl_id)
cres_no_target = unique(cres_w_grnas$grna_target[!cres_w_grnas$grna_target %in% c(grna_cres_w_cgenes,grna_cres_w_egenes)])
overlapping_cres = grna_cres_w_cgenes[grna_cres_w_cgenes %in% grna_cres_w_egenes]
overlapping_cres_cgene = unique(cres_w_grnas$target_gene[cres_w_grnas$grna_target %in% overlapping_cres & cres_w_grnas$significant == 1])
overlapping_cres_egene = unique(cres_w_grnas_egene$target_gene[cres_w_grnas_egene$grna_target %in% overlapping_cres])
overlapping_target_genes = overlapping_cres_cgene[overlapping_cres_cgene %in% overlapping_cres_egene]

# Load power
crispr.power = fread("cres_with_grnas_power.txt") %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id)) %>%
  filter(effect_size == 0.85)

eqtl.power = fread("power_results/eqtl.power.combined.txt") %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id)) 

sum(eqtl.power$eQTL.power > 0.79)/nrow(eqtl.power)

# Filter to overlapping significant CREs
crispr.power.overlap = crispr.power %>%
  distinct(target_gene, .keep_all=T) %>%
  filter(target_gene %in% c(overlapping_cres_egene,overlapping_cres_cgene)) %>%
  mutate(intersect = as.factor(ifelse(target_gene %in% overlapping_target_genes, "Intersection", "No intersection")))

#Fisher exact test
a = sum(crispr.power.overlap$power > 0.79 & crispr.power.overlap$intersect == "Intersection")
b = sum(crispr.power.overlap$power > 0.79 & crispr.power.overlap$intersect == "No intersection")
c = sum(crispr.power.overlap$power < 0.8 & crispr.power.overlap$intersect == "Intersection")
d = sum(crispr.power.overlap$power < 0.8 & crispr.power.overlap$intersect == "No intersection")
print(matrix(c(a,b,c,d),nrow=2,ncol=2))
test = fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))

eqtl.power.overlap = eqtl.power %>%
  distinct(target_gene, .keep_all=T) %>%
  filter(target_gene %in% c(overlapping_cres_egene,overlapping_cres_cgene)) %>%
  mutate(intersect = as.factor(ifelse(target_gene %in% overlapping_target_genes, "Intersection", "No intersection")))

sum(eqtl.power.overlap$eQTL.power > 0.79 & eqtl.power.overlap$intersect == "Intersection")

#Fisher exact test
a = sum(eqtl.power.overlap$eQTL.power > 0.79 & eqtl.power.overlap$intersect == "Intersection")
b = sum(eqtl.power.overlap$eQTL.power > 0.79 & eqtl.power.overlap$intersect == "No intersection")
c = sum(eqtl.power.overlap$eQTL.power < 0.8 & eqtl.power.overlap$intersect == "Intersection")
d = sum(eqtl.power.overlap$eQTL.power < 0.8 & eqtl.power.overlap$intersect == "No intersection")
print(matrix(c(a,b,c,d),nrow=2,ncol=2))
test = fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))

eqtl.power = fread("power_results/eqtl.power.combined.txt") %>%
  mutate(target_gene = paste0(grna_target,"_",ensembl_id)) 

# Load sumstats to obtain MAFs
# Load GWAS summary stats
sumstats = fread("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/UKBB_sumstats/30000_formatted.tsv")
sumstats$snp_pos = paste0(sumstats$Chr,":",sumstats$Pos)
sumstats_filtered = sumstats %>% filter(snp_pos %in% eqtl.power$grna_target)

plot_df <- eqtl.power %>%
  left_join(sumstats_filtered[, c("snp_pos", "minor_AF")], by = c("grna_target" = "snp_pos"))

test = wilcox.test(plot_df$minor_AF[plot_df$unique_CREs == "cgenes"], plot_df$minor_AF[plot_df$unique_CREs == "egenes"])
print(test)