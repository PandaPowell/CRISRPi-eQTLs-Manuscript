install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

.libPaths("/gpfs/commons/home/sghatan/R/x86_64-pc-linux-gnu-library/4.4")
print(.libPaths())

# List of required packages
packages <- c("data.table", "tidyverse", "scPower", "Matrix")

# Install and load each package
lapply(packages, install_if_missing)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("supply eQTL name \n")
} else {
  cell_name = args[1]
}

#cell_name = "DC"

# Obtain chr, pos of genes in GRC37
# Load gencode data
annot_file = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/Gencode/gencode.v33lift37.GRCh38.genes.gtf"
# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$ensembl_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[4]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot.GRC37 <- annot[order(annot$chr_number, annot$start),c("chr_number","start","ensembl_id", "gene_name")]

# Load count matrices
dir = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/cell_specific_rds/lvl1/"
count.matrix <- readRDS(paste0(dir, cell_name, "_raw_counts.RDS"))

# Load metadata (must contain 'individual' column)
metadata <- readRDS(paste0(dir, cell_name, "_meta_data.rds"))

# Ensure metadata matches filtered cells
metadata <- metadata[rownames(metadata) %in% colnames(count.matrix),]

# Filter to genes with annot
count.matrix = count.matrix[rownames(count.matrix) %in% annot.GRC37$gene_name,]

# Convert counts to presence/absence (1 if expressed, 0 otherwise)
binary_count_matrix <- count.matrix > 0

# Count the number of unique genes expressed per cell (column-wise sum)
unique_gene_counts <- Matrix::colSums(binary_count_matrix)

# Filter cells with unique gene counts between 100 and 2500
filtered_cells <- unique_gene_counts > 100 & unique_gene_counts < 2500

# Subset the count matrix based on this filtering
count.matrix <- count.matrix[, filtered_cells]

# Remove genes that are not expressed in at least 7 cells
count.matrix <- count.matrix[Matrix::rowSums(count.matrix) > 6, ]

# Convert counts to presence/absence (1 if expressed, 0 otherwise)
binary_count_matrix <- count.matrix > 0

# Create a mapping from cell to individual
cell_to_individual <- setNames(metadata$individual, rownames(metadata))

# Convert binary_count_matrix column names from cells to individual IDs
individual_labels <- cell_to_individual[colnames(binary_count_matrix)]

# Get unique individual IDs and their corresponding column indices
unique_individuals <- unique(individual_labels)
num_individuals <- length(unique_individuals)

# Create an empty sparse matrix (genes x individuals)
gene_individual_matrix <- Matrix(0, nrow = nrow(binary_count_matrix), ncol = num_individuals, sparse = TRUE)
rownames(gene_individual_matrix) <- rownames(binary_count_matrix)
colnames(gene_individual_matrix) <- unique_individuals

# Aggregate expression per individual (1 if any cell from an individual expresses the gene)
for (ind in unique_individuals) {
  # Get all cells belonging to the individual
  cell_indices <- which(individual_labels == ind)
  
  # Perform element-wise OR across all cells of this individual
  gene_individual_matrix[, ind] <- Matrix::rowSums(binary_count_matrix[, cell_indices, drop = FALSE]) > 0
}

# Count unique individuals expressing each gene
gene_individual_counts <- Matrix::rowSums(gene_individual_matrix)

# Compute 20% threshold
threshold <- 0.2 * num_individuals

# Keep only genes expressed in at least 20% of individuals
filtered_genes <- rownames(gene_individual_matrix)[gene_individual_counts >= threshold]

# Subset the count matrix to retain only these genes
count_matrix <- count.matrix[rownames(count.matrix) %in% filtered_genes, ]

# Print final number of genes
cat("Final number of genes after filtering:", nrow(count_matrix), "\n")

# Take the mean of these counts
#count.mean.sum.tmp = Matrix::rowSums(count.matrix)/length(unique(annot.df$individual))
# Calculate gene mean per individual
count.sum = Matrix::rowSums(count.matrix)
count.mean.sum = Matrix::rowSums(count.matrix)/length(unique(metadata$individual))
# calculate mean per cell
gene.means.per.cell<-Matrix::rowSums(count.matrix)/ncol(count.matrix)
# Make df
count.df = data.frame(gene_name = names(count.mean.sum),
                      count.mean = gene.means.per.cell,
                      count.mean.sum = count.mean.sum)

# Add ensemble ids
count.df = count.df %>%
  left_join(annot.GRC37, "gene_name") %>%
  filter(is.na(ensembl_id) == F)

# To speed up computation only test the power of cis genes
# Load targeted GWAS variants
cres_w_grnas = fread("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/cres_with_grnas.txt") %>% distinct(grna_target, .keep_all = T)
# flip alleles to match gwas snp ids
temp = gsub("_", ":", cres_w_grnas$finemap_snp_intersect_grna)
temp = str_split_fixed(temp,":",4)[,2]
cres_w_grnas$snp_pos = paste0(cres_w_grnas$chr, ":", temp)
cres_w_grnas$pos = temp

# Load GWAS summary stats to obtain MAFs
sumstats = fread("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/UKBB_sumstats/30000_formatted.tsv")
sumstats$snp_pos = paste0(sumstats$Chr,":",sumstats$Pos)
sumstats_filtered = sumstats %>% filter(snp_pos %in% cres_w_grnas$snp_pos)

power.per.cre = function(chr, pos, snp, grna_target, beta){
  
  print(snp)
  
  # Calculate heritability/variance explained by MAF and effect size
  maf = sumstats_filtered$minor_AF[sumstats_filtered$snp_pos == snp]
  
  if(length(maf) == 0){
    maf = mean(sumstats_filtered$minor_AF)
  }
  
  q2 = 2*maf*(1-maf)*(beta^2)
  
  # Obtain genes within 1mbp of SNP
  cis.count.df <- count.df[count.df$chr_number == chr & 
                             count.df$start >= (pos - 1000000) & 
                             count.df$start <= (pos + 1000000), ]
  

  if(nrow(cis.count.df) < 1){
    stop("No cis genes in the region")
  }
  
  #Skip power calculation for not expressed genes (<0.01)
  power_results <- lapply(1:nrow(cis.count.df), function(x) {
    if (cis.count.df$count.mean.sum[x] < 0.01) {
      return(0) # Set power to 0 if count.mean < 0.01
    } else {
      return(scPower:::power.eqtl(
        count.mean = cis.count.df$count.mean.sum[x], # Expression mean in the pseudobulk
        heritability = q2,                      # Heritability
        sig.level = 1e-03,                      # Significance threshold
        nSamples = length(unique(metadata$individual)) # Sample size
      ))
    }
  })
  
  power_results_df = data.frame(snp = snp,
                                genes = cis.count.df$gene_name,
                                ensembl_id = cis.count.df$ensembl_id,
                                count.mean = cis.count.df$count.mean,
                                count.mean.sum = cis.count.df$count.mean.sum,
                                eQTL.power = unlist(power_results),
                                grna_target = grna_target)
  
  return(power_results_df)
  
}

power.per.cre = purrr::possibly(power.per.cre, otherwise = NA, quiet = F)

results = lapply(1:nrow(cres_w_grnas), function(x) power.per.cre(chr = cres_w_grnas$chr[x],
                                                                 pos = as.numeric(cres_w_grnas$pos[x]),
                                                                 snp = cres_w_grnas$snp_pos[x],
                                                                 grna_target = cres_w_grnas$grna_target[x],
                                                                 beta = 0.5))

power_results_df = do.call(rbind, results) %>%
  mutate(id = paste0(snp,genes))

results2 = lapply(1:nrow(cres_w_grnas), function(x) power.per.cre(chr = cres_w_grnas$chr[x],
                                                                  pos = as.numeric(cres_w_grnas$pos[x]),
                                                                  snp = cres_w_grnas$snp_pos[x],
                                                                  grna_target = cres_w_grnas$grna_target[x],
                                                                  beta = 0.25))

power_results_df2 = do.call(rbind, results2) %>%
  mutate(id = paste0(snp,genes))

#Take cell type frequencies estimated from BioRad (see above)
cell.types.biorad<-data.frame(ct=c("B", "Mono", "NK","DC","CD4_T","CD8_T"),
                              alt.ct.name=c("B cells", "CD14+ Monocytes", "NK cells","Dendritic cells","CD4 T cells","CD8 T cells"),
                              freq=c(13,13,7,1.5,30,13)/100)

# Calculate per gene expression probabilites
#Get the fraction of cell type cells
ctCells<-round(ncol(count.matrix)/length(unique(metadata$individual)))
ct = cell.types.biorad$alt.ct.name[cell.types.biorad$ct == cell_name]

#Fit dispersion parameter dependent on mean parameter
disp.fun<-disp.fun.param[disp.fun.param$ct==ct,]

calculate_gene_prob = function(df, beta){
  
  # calculate dispersion
  df$disp<-sample.disp.values(df$count.mean,disp.fun)
  df$disp.sum<-df$disp/ctCells
  
  #Sort simulated genes by mean expression
  df<-df[order(df$count.mean, decreasing = TRUE),]
  
  #Calculate for each gene the expression probability
  df$exp.probs<-estimate.exp.prob.values(df$count.mean, # mean per cell
                                                       1/df$disp, # disp per cell
                                                       ctCells, # num cells per cell type
                                                       nSamples=length(unique(annot.df$individual)),
                                                       min.counts=6,
                                                       perc.indiv.expr=0.2)
  
  summary(df$exp.probs)
  
  df$overall_power = df$eQTL.power*df$exp.probs
  
  print(summary(df$overall_power))
  
  system("mkdir -p power_results/")
  
  fwrite(df, paste0("power_results/",cell_name,".power.results.",beta,".txt"), row.names = F, quote = F)
  
  
}

calculate_gene_prob(power_results_df, 0.5)
calculate_gene_prob(power_results_df2, 0.25)
