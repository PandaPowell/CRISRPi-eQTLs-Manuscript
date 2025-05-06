libraries <- c("data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "ggplot2",
               "biomaRt",
                "cowplot")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/")

cell_types= c("NK_cells", "B_cells" ,"CD4_T_cells" ,"CD8_T_cells", "DC_mean", "Other_T_cells", "Mono_cells", "Other_cells")

for (cell in cell_types){
  print(cell)
  # Count number of colocalised indepedent variants
  coloc_res = fread(paste0("coloc_results_combined/",cell,"_coloc_results.txt"))
  # Filter our snps with low pvalues
  coloc_res = coloc_res[pval_nominal < 1e-05 & pval < 1e-05,]
  print(nrow(coloc_res))
  # Assign single cv coloc credible sets to 1
  coloc_res[is.na(gwas_credible_set) == T, "gwas_credible_set"] = 1
  coloc_res[is.na(eqtl_credible_set) == T, "eqtl_credible_set"] = 1
  # create a credible set id for which to count colocolized independent variants
  coloc_res[,eqtl_cred_set := paste0(chr.x, ":",region, ":", eqtl_credible_set)]
  coloc_res = coloc_res %>% distinct(eqtl_cred_set, .keep_all = T)
  print(length(unique(coloc_res$eqtl_cred_set)))
}

# Make lollipop graph of results
plot_table = data.frame("Cell_type" = cell_types,
                        "Number of indepedent variants" = c(239,187,368,237,45,45,124,22),
                        "Percentage of colocalised variants" = c(10,12,6,9,25,25,16,45))

png("plots/no_coloc_snps_cell_by_cell_type.png", units="in", width=10, height=6, res=300)

plot_table %>%
  arrange(Number.of.indepedent.variants) %>%
  ggplot( aes(x=Cell_type, y=Number.of.indepedent.variants) ) +
  geom_segment( aes(xend=Cell_type, yend=0), ) +
  geom_segment( aes(xend=Cell_type, yend=0), ) +
  geom_point( size=3) +
  coord_flip()+
  theme_cowplot() +
  xlab("Cell types") +
  ylab("Number of colocalized causal variants")
  
dev.off()

png("plots/percent_coloc_snps_cell_by_cell_type.png", units="in", width=10, height=6, res=300)

plot_table %>%
  arrange(Percentage.of.colocalised.variants) %>%
  ggplot( aes(x=Cell_type, y=Percentage.of.colocalised.variants) ) +
  geom_segment( aes(xend=Cell_type, yend=0), ) +
  geom_segment( aes(xend=Cell_type, yend=0), ) +
  geom_point( size=3) +
  coord_flip()+
  theme_cowplot() +
  xlab("Cell types") +
  ylab("Percentage of colocalized causal variants")

dev.off()

### # For each GWAS count the number of colocalized variants
combined_coloc = fread("coloc_results_combined/All_cell_coloc_results.txt")[gwas != "",]

# Make dot plot of eQTL and GWAS
data = as.data.frame(table(combined_coloc$gwas, combined_coloc$eqtl))

png("plots/gwas_by_eqtl_dotplot.png", units="in", width=10, height=8, res=300)

# Create a scatter plot
ggplot(data, aes(x = Var2, y = Var1, size = Freq, color = Freq)) +
  geom_point(alpha = 0.8) +  # Adjust point transparency
  scale_color_gradient2(low = "white", mid = "red", high = "black", midpoint = 85) +
  scale_size(range = c(0, 8), breaks = c(10, 50, 100, 150), labels = c("10", "50", "100", "150")) +
  theme_cowplot() +
  labs(title = "",
       x = "eQTL cell type",
       y = "Blood cell GWAS",
       size = "Number of colocalizations",
       color = "Number of colocalizations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

dev.off()


gwas_names = unique(combined_coloc$gwas)

# Initialize vectors to store results
gwas_list = vector()
gene_trait_pairs = vector()
colocalized_cred_sets = vector()

for (i in gwas_names) {
  print(i)
  gwas = combined_coloc[combined_coloc$gwas == i, ]
  gwas = gwas[pval_nominal < 1e-05 & pval < 1e-05, ]
  cat("Number of gene-trait pairs", nrow(gwas), "\n")
  gene_trait_pairs = c(gene_trait_pairs, nrow(gwas))
  
  # Assign single cv coloc credible sets to 1
  gwas[is.na(gwas_credible_set) == T, "gwas_credible_set"] = 1
  gwas[is.na(eqtl_credible_set) == T, "eqtl_credible_set"] = 1
  
  # Create a credible set id for which to count colocalized independent variants
  gwas[, gwas_cred_set := paste0(chr.x, ":", region, ":", gwas_credible_set)]
  gwas = gwas %>% distinct(gwas_cred_set, .keep_all = TRUE)
  cat("Number of colocalized credible sets", length(unique(gwas$gwas_cred_set)), "\n")
  colocalized_cred_sets = c(colocalized_cred_sets, length(unique(gwas$gwas_cred_set)))
  
  gwas_list = c(gwas_list, i)
}

# Create a DataFrame with the results
results_df = data.frame(
  GWAS = gwas_list,
  Gene_Trait_Pairs = gene_trait_pairs,
  Colocalized_Credible_Sets = colocalized_cred_sets,
  Percentage = c(
    19.6181, 13.9738, 13.8943, 14.5969, 17.0963, 14.8784, 16.3265, 16.6934,
    16.0134, 15.3946, 14.3322, 15.9226, 20.1709, 14.8607, 19.2872, 15.5722,
    19.5122, 22.0884, 18.4211, 22.2222, 17.0437, 19.8630, 14.0449, 13.6276,
    18.2266, 17.7852, 17.4263, 13.6525, 14.7425),
  Group = c(
    "White blood cells", "Red blood cells", "Red blood cells", "Red blood cells",
    "Red blood cells", "Red blood cells", "Red blood cells", "Red blood cells",
    "Platelets", "Platelets", "Platelets", "Platelets", "White blood cells",
    "White blood cells", "White blood cells", "White blood cells", "White blood cells",
    "White blood cells", "White blood cells", "White blood cells", "White blood cells",
    "White blood cells", "Red blood cells", "Red blood cells", "Red blood cells",
    "Red blood cells", "Red blood cells", "Red blood cells", "Red blood cells")
)

# Make plot

png("plots/no_credsets_by_gwas.png", units="in", width=10, height=6, res=300)

results_df %>%
  arrange(Percentage) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Trait=factor(GWAS, levels=GWAS)) %>%   # First sort by val. This sort the dataframe but NOT the factor levels
  ggplot( aes(x=Trait, y=Percentage, color = Group) ) +
  geom_segment( aes(xend=Trait, yend=0), ) +
  geom_segment( aes(xend=Trait, yend=0), ) +
  geom_point( size=3) +
  coord_flip()+
  theme_cowplot() +
  scale_color_manual(values = c("thistle", "tomato1", "lightblue2")) +
  xlab("Traits") +
  ylab("Percentage (%) of colocalized causal variants")

dev.off()

# Intersect STING-seq variants with eQTLs
sting_seq = fread("Post_coloc_analysis/All_STING_seq_CREs.csv")
# Load finemapped sting-seq credible sets
sting_seq_cred_sets = fread("data/Stingseq_data/Finemapped_variants_STINGseq_UKBB.csv")
sting_seq_cred_sets[,SNP_coords := paste0(Chromosome,":",`Position (hg19)`)]
sting_seq_cred_sets = sting_seq_cred_sets[SNP_coords %in% sting_seq$`SNP Coordinates (hg19)`,] %>%
  arrange(desc(PIP)) %>%
  distinct(`SNP ID`, .keep_all=T)

determine_overlap = function(x){
  # Loop through each variant and extract the all the snps in that credible set
  gwas_id = sting_seq_cred_sets$`GWAS ID`[x]
  chr = sting_seq_cred_sets$Chromosome[x]
  ss_variant = sting_seq_cred_sets$`SNP ID`[x]
  lower = sting_seq_cred_sets$`Position (hg19)`[x] - 1000
  upper = sting_seq_cred_sets$`Position (hg19)`[x] + 1000
  # Load appropiate gwas
  finemap = fread(paste0("results/UKBB_SuSiE_finemap/",gwas_id,"/",gwas_id,"_chr",chr,"_finemap_results.txt"))
  # Check if sting-seq varaint with 10kb of credible set
  cred_set = finemap[inrange(pos, lower, upper), cond_indep]
  cred_set_snps = finemap[cred_set == cond_indep, variant]
  # Check if any of the snps in this eqtl over lap with a colocalised variant
  name = "Mono_cells"
  eqtl = fread(paste0("results/Coloc_results_V2/",name,"/",gwas_id,"/","chr",chr,"_coloc_results.txt"))
  # Create snp coords column
  eqtl[,SNP_coords:=paste0(Chr,":",Pos)]
  eqtl_overlap = eqtl[SNP_coords %in% cred_set_snps,]
  print(eqtl_overlap)
  # Then I want the output just to be the sting seq variant and the colocalised eqtl varaint
}

determine_overlap2 = purrr::possibly(determine_overlap, otherwise = NA, quiet = F)
lapply(1:nrow(sting_seq_cred_sets), function(x) determine_overlap2(x))


cell_names = c("NK_cells", "B_cells","CD4_T_cells",  "CD8_T_cells" , "DC_mean" , "Mono_cells" , "Other_cells",  "Other_T_cells")

single_cell_results = data.table()

for (name in cell_names){
  
  cat("Loading",name,"\n")
  # Load coloc results for each cell type
  cell_type = fread(paste0("results/Coloc_results_V2/",name,"/", name,"_combined_coloc_results.txt"), fill=T)[,c(1:7,12:15,29,31,47,55:57)]
  cell_type[,unique_id := paste0(Chr.x,":",cond_indep,"_", trait.x)]
  # Filter for sig coloc unique ids
  sig_coloc = cell_type[pval < 1e-05 & pval_nominal < 1e-03 & PP.H4.abf > 0.8,]
  pph4 = unique(sig_coloc$unique_id)
  genes = unique(sig_coloc$trait.y)
  
  cat(length(pph4),"coloc independent variant for", name, "\n")
  cat(length(genes),"genes", name, "\n")
  
  # Rejoin all snps in credible set back to sig colocs
  #snps_pips = cell_type[,c("SNP","Chr.x", "region","cs" ,"trait.x", "pip","pval","Pos","pval_nominal", "trait.y", "cond_indep","gwas","unique_id")]
  snps_pips = cell_type
  sig_cell_type = snps_pips %>% left_join(sig_coloc[,c("PP.H4.abf","pval_nominal","trait.y","unique_id")],"unique_id") %>%
    filter(is.na(PP.H4.abf.y) == F)
  sig_cell_type = sig_cell_type[,c(1:7,10:13,16:21)]
  sig_cell_type[,cell_type := name]
  single_cell_results = rbind(single_cell_results, sig_cell_type)
}

# Load GTEx results
gtex = fread("results/Coloc_results/GTEx_all_coloc_results.txt")
colnames(gtex) = c(colnames(gtex)[2:53],"gwas")
gtex[,unique_id := paste0(Chr.x,":",cond_indep,"_", gwas)]
gtex_sig = gtex[PP.H4.abf > 0.8 & pval < 1e-05 & pvalue < 1e-03, ]
pph4 = unique(gtex_sig$unique_id)
genes = unique(gtex_sig$trait.y)

cat(length(pph4),"coloc independent variant for", name, "\n")

cat(length(genes),"genes", name, "\n")

gtex_results = gtex[unique_id %in% pph4,]

# Convert ensemble ids
convert_ensembl_to_gene_name <- function(ensembl_ids) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = ensembl_ids,
                      mart = ensembl)
  return(gene_names)
}

# Example usage
ensembl_ids <- genes  # Add your IDs here
gtex_genes <- convert_ensembl_to_gene_name(ensembl_ids)
colnames(gtex_genes)[1] = c("trait.y")

gtex_results = gtex_results %>% left_join(gtex_genes, "trait.y")

sig_sc_results = single_cell_results[pval < 1e-05 & pval_nominal.y < 1e-03 & PP.H4.abf.y > 0.8,]
sc_genes = unique(sig_sc_results$trait.y)
all_genes = unique(c(gtex_genes$external_gene_name, sc_genes))

grep("GFI1B|NFE2|IKZF1|HHEX|RUNX1", unique(c(gtex_genes$external_gene_name, sc_genes)))

# Create SNP coords column in eQTL data
gtex_results[,SNP_coords := paste0(Chr.y,":",Pos)]
# Determine which cond indep coloc
gtex_overlap = gtex_results[SNP_coords %in% sting_seq$SNP_coords,c("SNP_coords", "cond_indep")] %>% distinct(SNP_coords,.keep_all = T)
# Make new column
gtex_overlap[ ,stingseq_variants := SNP_coords]
# Get all variants in those credible sets
gtex_overlap_cs = gtex_results[cond_indep %in% gtex_overlap$cond_indep,] %>% filter(is.na(external_gene_name) == F, external_gene_name!="") %>%
  mutate(temp = paste0(SNP_coords,external_gene_name)) %>% distinct(temp, .keep_all = T) %>% 
  left_join(gtex_overlap, "cond_indep")


single_cell_results[,SNP_coords := paste0(Chr.x,":",Pos)]
snp_overlap = single_cell_results[SNP_coords %in% sting_seq$SNP_coords,c("SNP_coords","trait.y.y")]
colnames(snp_overlap) = c("SNP_coords","gene")




one_table = sting_seq %>% left_join(snp_overlap[is.na(gene)==F,],"SNP_coords") %>%
  left_join(gtex_overlap[is.na(gene)==F,],"SNP_coords") %>% filter(is.na(gene.y) == F | is.na(gene) ==F)

cat("Number of overlaping variants between eqtl and ss", length(unique(one_table$SNP_coords)))

ov = one_table[one_table$gene.x == one_table$gene.y | one_table$gene.x == one_table$gene]

# Overlap of multi ethnic sting-seq hits
multi_coloc = fread("New_multi_coloc_pipeline/multi_ancestry_coloc_res.txt") %>%
  mutate(ss_SNP_coords = gsub("_.*", "",sting_seq_variant))
BCX_ss = sting_seq %>% filter(SNP_coords %in% multi_coloc$ss_SNP_coords) %>% mutate(ss_SNP_coords = SNP_coords) %>%
  left_join(multi_coloc[,c("ss_SNP_coords","gwas","geneSymbol","tss_distance","slope","pval_nominal")], "ss_SNP_coords")

sum(BCX_ss$gene == BCX_ss$geneSymbol)