libraries <- c("coloc",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "ggplot2",
               "biomaRt")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Supply necessary files")
} else {
  gwas_name = args[1]
}

print(paste("GWAS name:", gwas_name))

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/")

# Download finemapping results from BCX consortium
# dest = "data/BCX_sumstats/BCX_finemapping_results.zip"
# source = "http://www.mhi-humangenetics.org/dataset/BCX2_CR95_VEP.zip"
# download.file(source, dest)
# system("mkdir -p data/BCX_sumstats/BCX_finemapping_results")
# system("unzip data/BCX_sumstats/BCX_finemapping_results.zip -d data/BCX_sumstats/BCX_finemapping_results/")

# Load trans finemapping results
finemap = fread(paste0("data/BCX_sumstats/BCX_finemapping_results/", gwas_name,"_Trans_Credible_Sets_VEP.txt")) %>%
  separate(MarkerName, into = c("chr", "pos_tmp"), remove = F, sep = ":") %>%
  separate(pos_tmp, into=c("pos","alleles"), sep = "_") %>% mutate(pos = as.numeric(pos)) %>%
  dplyr::select(MarkerName, chr, pos, Sentinel_MarkerName) %>% 
  mutate(lower = pos-250000, upper = pos+250000, lower_bed = pos-1, chr_bed = paste0("chr",chr))
finemap = as.data.table(finemap)

# Make bed file finemap
fwrite(finemap[,c("chr_bed","lower_bed","pos","MarkerName","Sentinel_MarkerName")], paste0("results/MAGE_finemap/",gwas_name,"_finemap_results.bed"), sep = "\t", row.names = F, col.names = F, quote = F)

# Define regions for coloc based on finemapping results
n_region = 0
iteration = 0
keep_going = T

while(keep_going){
  
  finemap_regions = c()
  options(warn=1)
  
  for (i in 1:22){
  
    # Filter by chromosome & pvalue
    gwas_chr = finemap[chr == i,]
    gwas_chr = gwas_chr[order(lower,upper),]
    # Calculate differences between regions
    gwas_chr[, difference := lower - shift(upper, 1, type = "lag"),]
    gwas_chr$difference[is.na(gwas_chr$difference) == T] = 1
    # Group variants with overlaping regions
    gwas_chr[, group := cumsum(difference > 0) + 1]
    # Calculate lower and upper bounds of grouped region
    gwas_chr[, `:=` (regions = paste0(i, ":",min(lower),"-", max(upper))), by = group]
    # Get list of regions
    regions = gwas_chr %>% distinct(regions, .keep_all=T)
    finemap_regions = c(finemap_regions, unique(gwas_chr$regions))
  }
  
  combined_df = data.table(chr = as.numeric(sub(":.*$", "",finemap_regions)),
                           lower = as.numeric(str_split_fixed(gsub(".*:","",finemap_regions), "-",2)[,1]),
                           upper = as.numeric(str_split_fixed(gsub("*:","",finemap_regions), "-",2)[,2]),
                           finemap_regions = finemap_regions)
  
  if (n_region != nrow(combined_df)){
    keep_going = T
    n_region = nrow(combined_df)
    iteration = iteration + 1
    cat("Regions are not yet converged. Interation", iteration,"\n")
  } else {
    print("Regions have converged!")
    keep_going = F
  }
  
}

# Function to harmonize alleles without merging
harmonize_alleles <- function(dt1, dt2) {
  # Set SNP as key for quick reference lookup
  setkey(dt1, SNP_coords)
  setkey(dt2, SNP_coords)
  
  # Merge datasets
  dt_merged <- merge(dt1, dt2, by = "SNP_coords", suffixes = c(".1", ".2"))
  # Check for inconsistencies in reference alleles
  dt_merged[, match_ref := (variantRef != reference_allele & variantRef == other_allele)]
  # Create the new columns in dt2_common before performing the conditional operation
  dt_merged[, `:=` (ref_h = variantRef, alt_h = variantAlt, raf = maf)]
  # Flip betas where reference alleles do not match
  dt_merged[match_ref == TRUE, `:=` (
    ref_h = variantAlt,
    alt_h = variantRef,
    slope = -slope,
    raf = 1-maf
  )]
  
  dt_merged[,variant := paste0(SNP_coords,"_",ref_h, "_",alt_h)]
  
  # Split data.tables again
  dt1_h = dt_merged[,1:20] %>% distinct(rs_number, .keep_all=T) %>%
    mutate(variant = paste0(SNP_coords,"_",reference_allele, "_",other_allele))
  dt2_h = dt_merged[,c(1,21:47)]
  
  # Drop temporary columns
  dt2_h[, c("match_ref", "maf", "variantRef","variantAlt") := NULL]
  
  # Return both datasets as a list if further joint analysis is needed
  return(list(dt1 = dt1_h, dt2 = dt2_h))
}

# Load appropriate sumstats
sumstats = fread(paste0("data/BCX_sumstats/",gwas_name, ".tsv.gz"))
sumstats[,chr := str_split_fixed(sumstats$rs_number,":",2)[,1]]
sumstats[,pos := str_split_fixed(sumstats$rs_number,":",2)[,2]]
sumstats[,pos := as.numeric(gsub("_.*", "", pos))]

# chromosome = combined_df$chr[1]
# lower = combined_df$lower[1]
# upper = combined_df$upper[1]

run_coloc = function(chromosome,lower,upper){
  
  cat(paste0("Preparing to run COLOC for region ",chromosome,lower,"-",upper,"\n"))

  gwas_region = sumstats[chr == chromosome & inrange(pos, lower,upper),]
  gwas_region[,SNP_coords:=paste0("chr",chr,":",pos)]
  
  # Load MAGE sumstats
  mage = fread(paste0("data/MAGE/eQTL_FastQTL_results_chr",chromosome,".hg19.txt"))
  eqtl_region = mage[ inrange(upper_hg19, as.numeric(lower),as.numeric(upper))]
  eqtl_region[,variant := paste0(chr,":",upper_hg19,"_",variantAlt,"_",variantRef)]
  eqtl_region[,SNP_coords := paste0(chr,":",upper_hg19)]
  
  if (min(eqtl_region$pval_nominal, na.rm=T) > 1e-05){
    stop("No GWS significant variants in this region for this eQTL Lower p-value= ", min(eqtl_region$pval_nominal, na.rm=T))
  }
  
  # Loop through genes in eqtl region
  genes = unique(eqtl_region$geneSymbol)
  
  coloc_tr_snps = data.frame()
  
  for (name in genes){
    
    gene_region = eqtl_region[geneSymbol == name]
    
    if (min(gene_region$pval_nominal, na.rm=T) > 1e-05){
      warning("No GWS significant variants in this region for the genes ",name, " Lower p-value= ", min(gene_region$pval_nominal, na.rm=T))
      next
    } else{
      cat("significant variants found for gene",name,"\n")
    }
    
    harmonized_datasets <- harmonize_alleles(gwas_region, gene_region)
    gwas_region_h = harmonized_datasets$dt1
    gene_region_h = harmonized_datasets$dt2
    
    int = intersect(gwas_region_h$variant, gene_region_h$variant)
    
    gwas_region_h = gwas_region_h[variant %in% int,]
    gwas_region_h = unique(gwas_region_h, by = c("variant"))
    gwas_region_h[, varbeta := se^2]
    gwas_region_h[, z := beta / se]
    
    gene_region_h = gene_region_h[variant %in% int,]
    gene_region_h = unique(gene_region_h, by = c("variant"))
    gene_region_h[, varbeta := slope_se^2]
    gene_region_h[, z := slope / slope_se]
    
    # Format eQTL data
    b2 = c(gene_region_h$slope)
    names(b2) = gene_region_h$variant
    vb2 = c(gene_region_h$varbeta)
    names(vb2) = gene_region_h$variant
    maf2 = c(gene_region_h$raf)
    names(maf2) = c(gene_region_h$variant)
    
    D1 = list(gene_region_h$variant, gene_region_h$upper_hg19, b2, vb2, maf2, 731, "quant")
    names(D1) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type")
    
    # Format GWAS data
    b2 = c(gwas_region_h$beta)
    names(b2) = gwas_region_h$variant
    vb2 = c(gwas_region_h$varbeta)
    names(vb2) = gwas_region_h$variant
    maf2 = c(gwas_region_h$eaf)
    names(maf2) = c(gwas_region_h$variant)
    
    D2 = list(gwas_region_h$variant, gwas_region_h$pos, b2, vb2,maf2, 575453, "quant")
    names(D2) = c("snp", "position", "beta", "varbeta","MAF", "N", "type")
    
    res = coloc.abf(D1,D2)
    
    if(res$summary['PP.H4.abf'] < 0.5){
      
      cat("No significant colocalization detected \n")
      
    } else{
      
      coloc_snps = res$results[,c(1,12)]
      colnames(coloc_snps)[1] = "SNP"
      coloc_snps$hit2 = coloc_snps$SNP
      coloc_snps = coloc_snps[which(coloc_snps$SNP.PP.H4 == max(coloc_snps$SNP.PP.H4)),]
      coloc_res = as.data.frame(t(res$summary))
      coloc_res = coloc_res[rep(seq_len(nrow(coloc_res)),nrow(coloc_snps)), ]
      coloc_res = cbind(coloc_snps, coloc_res)
      coloc_res$gwas_credible_set = NA
      coloc_res$eqtl_credible_set = NA
      
      # order columns of dataframe
      coloc_res = coloc_res[,c("SNP","hit2","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","nsnps","gwas_credible_set","eqtl_credible_set","SNP.PP.H4")]
      
      df = filter(gene_region_h, variant %in% coloc_snps$SNP) %>%
        mutate(trait = name) %>% rename(SNP = variant)
      
      df2 = filter(gwas_region_h, variant %in% coloc_snps$SNP) %>% rename(SNP = variant) 
      
      EX_SNP = coloc_res %>% left_join(df2, "SNP") %>%
        left_join(df, "SNP") %>% mutate(region = paste0(lower,"-",upper)) %>%
        mutate( method = "coloc",gwas = gwas_name)
      
      coloc_tr_snps = bind_rows(coloc_tr_snps, EX_SNP)
    }
  }
  return(coloc_tr_snps)
}

options(warn=1)
### In case of an error return NA ### 
run_coloc2 = purrr::possibly(run_coloc, otherwise = NA, quiet = F)
testie = lapply(1:nrow(combined_df), function(x) run_coloc2(combined_df$chr[x],
                                                            combined_df$lower[x],
                                                            combined_df$upper[x]))

t2 = do.call(rbind, testie)

fwrite(t2, paste0("results/Coloc_results_V2/MAGE/MAGE_",gwas_name,"_coloc_res.txt"), sep = ",", row.names = F, quote=F)



