# Function to check and install missing packages
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

.libPaths("/gpfs/commons/home/sghatan/R/x86_64-pc-linux-gnu-library/4.4")
print(.libPaths())  # Check if the path is correctly set

# List of required packages
packages <- c("data.table", "tidyverse", "scPower")

# Install and load each package
lapply(packages, install_if_missing)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap")

data.sets = fread("eQTL_catalogue/eQTL_catalogue_datasets.txt")
files = system("ls eQTL_catalogue/count_matrices/*.gz", intern = T)
meta.files = system("ls eQTL_catalogue/count_matrices/*.tsv", intern = T)

calculate_means = function(x){
  
  print(x)
  
  count = fread(files[x])
  meta.data = fread(meta.files[x])
  name = str_split_fixed(str_split_fixed(meta.files[x], "/",4)[,3], "\\.",2)[,1]
  
  if(name == "GTExV8"){
    conditions = "blood"
  } else{
    conditions = unique(meta.data$qtl_group)
  }
  
  # Filter samples by condition
  for (i in 1:length(conditions)){
    
    samples = meta.data$pseudo_id[meta.data$qtl_group == conditions[i]]
    count.filtered =  count[, .SD, .SDcols = colnames(count) %in% c("phenotype_id", samples)]
    count.filtered$count_mean = rowSums(count.filtered[,2:ncol(count.filtered)])/(ncol(count.filtered) - 1)
    count.filtered$n = length(samples)
    
    system("mkdir -p eQTL_catalogue/count_means")
    
    fwrite(count.filtered[, c("phenotype_id","count_mean","n")], 
           paste0("eQTL_catalogue/count_means/",name,"_",conditions[i], ".csv"),
           quote = F, 
           row.names = F)
  }
  
}

# Calcylate count means and write to files
lapply(1:length(files), function(x) calculate_means(x))

# Remove non-blood related datasets
system("rm eQTL_catalogue/count_means/HipSci_iPSC.csv")
system("rm eQTL_catalogue/count_means/Schwartzentruber_2018_sensory_neuron.csv")
system("rm eQTL_catalogue/count_means/van_de_Bunt_2015_pancreatic_islet.csv")
system("rm eQTL_catalogue/count_means/BrainSeq_brain.csv")
system("rm eQTL_catalogue/count_means/FUSION_adipose_naive.csv")
system("rm eQTL_catalogue/count_means/FUSION_muscle_naive.csv")
system("rm eQTL_catalogue/count_means/GENCORD_fibroblast.csv")
# Remove duplicate studies with low counts
system("rm eQTL_catalogue/count_means/BLUEPRINT_PE_monocyte.csv")
system("rm eQTL_catalogue/count_means/BLUEPRINT_PE_neutrophil.csv")

# obtain list of count files
count.files = system("ls eQTL_catalogue/count_means/*.csv", intern = T)
# Load gencode with GRC37 positions
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

# To speed up computation only test the power of cis genes
# Load targeted GWAS variants
cres_w_grnas = fread("Post_coloc_analysis/cres_with_grnas.txt") %>% distinct(grna_target, .keep_all = T)
# flip alleles to match gwas snp ids
temp = gsub("_", ":", cres_w_grnas$finemap_snp_intersect_grna)
temp = str_split_fixed(temp,":",4)[,2]
cres_w_grnas$snp_pos = paste0(cres_w_grnas$chr, ":", temp)
cres_w_grnas$pos = temp

# Load GWAS summary stats to obtain MAFs
sumstats = fread("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/UKBB_sumstats/30000_formatted.tsv")
sumstats$snp_pos = paste0(sumstats$Chr,":",sumstats$Pos)
sumstats_filtered = sumstats %>% filter(snp_pos %in% cres_w_grnas$snp_pos)

power.per.cre = function(count.df, chr, pos, snp, grna_target, beta){
  
  print(snp)
  
  # Calculate heritability/variance explained by MAF and effect size
  maf = sumstats_filtered$minor_AF[sumstats_filtered$snp_pos == snp]
  
  if(length(maf) == 0){
    maf = mean(sumstats_filtered$minor_AF)
  }
  
  q2 = 2*maf*(1-maf)*(beta^2)
  
  # Obtain genes within 1mb of snp add extra as snp is in GRC37
  cis.count.df = count.df[chr == chr_number & inrange(start, pos-1000000, pos+1000000),]
  
  if(nrow(cis.count.df) < 1){
    stop("No cis genes in the region")
  }
  
  #Skip power calculation for not expressed genes (<0.01)
  power_results <- lapply(1:nrow(cis.count.df), function(x) {
    if (cis.count.df$count.mean[x] < 0.01) {
      return(0) # Set power to 0 if count.mean < 0.01
    } else {
      return(scPower:::power.eqtl(
        count.mean = cis.count.df$count.mean[x], # Expression mean in the pseudobulk
        heritability = q2,                      # Heritability
        sig.level = 1e-03,                      # Significance threshold
        nSamples = cis.count.df$n[1]) # Sample size
      )
    }
  })
  
  power_results_df = data.frame(snp = snp,
                                genes = cis.count.df$gene_name,
                                ensembl_id = cis.count.df$ensembl_id,
                                count.mean = cis.count.df$count.mean,
                                eQTL.power = unlist(power_results),
                                grna_target = grna_target)
  
  return(power_results_df)
  
}

power.per.cre = purrr::possibly(power.per.cre, otherwise = NA, quiet = F)

# Load and merge into one file
for (x in 1:length(count.files)){
  
  print(count.files[x])
  
  #Load count matrices
  count.df = fread(count.files[x])
  name = str_split_fixed(str_split_fixed(count.files[x], "/",4)[,3], "\\.",2)[,1]
  count.df$ensembl_id = str_split_fixed(count.df$phenotype_id,"\\.",2)[,1]
  colnames(count.df)[2] = c("count.mean")
  # Add ensemble ids
  count.df = count.df %>% left_join(annot.GRC37, "ensembl_id")
  
  eqtl.name = str_split_fixed(count.files[x], "/", 3)[,3]
  
  results = lapply(1:nrow(cres_w_grnas), function(x) power.per.cre(count.df = count.df,
                                                                   chr = cres_w_grnas$chr[x],
                                                                   pos = as.numeric(cres_w_grnas$pos[x]),
                                                                   snp = cres_w_grnas$snp_pos[x],
                                                                   grna_target = cres_w_grnas$grna_target[x],
                                                                   beta=0.5))
  power_results_df = do.call(rbind, results)
  
  system("mkdir -p power_results/")
  
  fwrite(power_results_df, paste0("Post_coloc_analysis/power_results/power.results.beta0.5", eqtl.name), row.names = F, quote = F)
  
  results = lapply(1:nrow(cres_w_grnas), function(x) power.per.cre(count.df = count.df,
                                                                   chr = cres_w_grnas$chr[x],
                                                                   pos = as.numeric(cres_w_grnas$pos[x]),
                                                                   snp = cres_w_grnas$snp_pos[x],
                                                                   grna_target = cres_w_grnas$grna_target[x],
                                                                   beta=0.25))
  power_results_df = do.call(rbind, results)
  
  system("mkdir -p power_results/")
  
  fwrite(power_results_df, paste0("Post_coloc_analysis/power_results/power.results.beta0.25", eqtl.name), row.names = F, quote = F)
  
}
