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
cres_w_grnas = fread("Post_coloc_analysis/processed_data/cres_with_grnas.txt") %>% 
  distinct(grna_target, .keep_all = T)
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
        nSamples = 4732) # Sample size
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

count.df = fread("eQTL_catalogue/count_means/INTERVAL_salmon_tximport_gene_gene_mean_counts.csv.gz")
colnames(count.df) = c("ensembl_id","count.mean")
count.df$n = 4732
name = "Interval"

# Add ensemble ids
count.df = count.df %>% left_join(annot.GRC37, "ensembl_id")
eqtl.name = "Interval"

results = lapply(1:nrow(cres_w_grnas), function(x) power.per.cre(count.df = count.df,
                                                                 chr = cres_w_grnas$chr[x],
                                                                 pos = as.numeric(cres_w_grnas$pos[x]),
                                                                 snp = cres_w_grnas$snp_pos[x],
                                                                 grna_target = cres_w_grnas$grna_target[x],
                                                                 beta=0.27))
power_results_df = do.call(rbind, results)

system("mkdir -p power_results/")

fwrite(power_results_df, paste0("Post_coloc_analysis/power_results/power.results.beta0.27", eqtl.name), row.names = F, quote = F)
