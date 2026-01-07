libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "gassocplot2")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("supply in order: \n
       1. GWAS summary statistics file path \n
       2. LD regions file path \n
       3. Path to folder containing LD matrices by chr (i.e data/UKBB_LDmatrices/9/)")
} else {
  gwas_path = args[1]
  ld_regions = args[2]
  ld_path = args[3]
}

# Uncomment the following variables to run interactively
# ld_regions = "New_multi_coloc_pipeline/lead_loci/chr1.regions"
# gwas_path = "data/UKBB_sumstats/30000_formatted.tsv"
# ld_path = "data/UKBB_LDmatrices/1/"

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap" # Working directory
setwd(wd)

cat("Loading GWAS summary statistcs\n")
# Loads GWAS summary statistcs input ising bash command
gwas <- fread(gwas_path, sep=" ")
gwas = gwas[,Pos:=as.integer(Pos)]

# Obtain name of GWAS summary statistics loaded
gwas_name <- gsub(".*\\/|_formatted.tsv$", "", gwas_path) # Simplify string operations

# Load regions file
cat("Loading primary regions\n")
sig_regions <- fread(file = ld_regions, header = F, stringsAsFactors = FALSE, fill=T)
setnames(sig_regions, c("CHR", "lower", "upper"))
sig_regions[,CHR:=as.numeric(CHR)]
sig_regions[,lower:=as.numeric(lower)]
sig_regions[,upper:=as.numeric(upper)]
sig_regions[,hg19_region:=paste("chr",CHR,":",lower,"-",upper,sep="")]

check_distance = function(x){
  
  CHR = sig_regions$CHR[x]
  LOWER = sig_regions$lower[x]
  UPPER = sig_regions$upper[x]
  REGION = paste(sig_regions$lower[x],sig_regions$upper[x],sep=".")
  
  cat("Checking the distance between region bounds and significant snps for",sig_regions$hg19_region[x], "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = gwas[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  if (min(snp_region$pval, na.rm=T) > 5e-08){
    stop("No GWS significant variants in this region for this GWAS. Lower p-value= ", min(snp_region$pval, na.rm=T))
  }
  
  snp_region[,lower_bound:=abs(Pos-LOWER)]
  snp_region[,upper_bound:=abs(Pos-UPPER)]
  
  sig_snps = snp_region[pval<5e-08,]
  
  if (max(sig_snps$lower_bound) < 100000 | max(sig_snps$upper_bound) < 100000){
    cat("Significant SNPs close to region boundary for", sig_regions$hg19_region[x], "\n", "Plotting region \n")
    
    cat("Loading GWAS LD and BIM files\n")
    
    # Define filenames
    LDfilename <- paste(ld_path,REGION,"/",REGION,sep = "")
    BIMfilename <- paste(ld_path,REGION,"/",REGION,sep = "")
    # Load LD and bim file data
    LD <- fread(paste(LDfilename, "ld", sep = "." ))
    BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
    # assign headers
    setnames(BIM, c("chr", "rsid", "dk", "pos", "alt", "ref"))
    BIM[, SNP := paste(chr, pos, alt, ref, sep = ":")]
    # assign SNP labels to LD matrix
    setnames(LD, BIM$SNP)
    LD[, SNP := BIM$SNP]
    # Remove columns from the LD matrix that are all NAs
    na_columns = colSums(is.na(LD)) == nrow(LD)
    nacol_names = names(na_columns[na_columns == TRUE])
    LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
    # Filter rows with complete cases
    LD = LD[complete.cases(LD)]
    # Filter snps in LD matrix for those in GWAS sum stats
    LD <- LD[LD$SNP %in% snp_region$variant,]
    # Make sure LD only contains columns in LD$SNP and "SNP"
    required_cols <- unique(c(LD$SNP,"SNP"))
    LD <- LD[, ..required_cols]
    # Remove duplicated SNPs
    LD_gwas <- unique(LD, by = "SNP")
    # Remove duplicates from BIM file
    BIM_gwas <- BIM[!duplicated(BIM, by = "SNP"), ]
    # Which SNPs intersect
    int1 = intersect(snp_region$variant, LD_gwas$SNP)
    
    # Remove NAs from GWAS LD and format into matrix
    LD_gwas2 = LD_gwas[LD_gwas$SNP %in% int1,]
    LD_gwas2 = LD_gwas2[, colnames(LD_gwas2) %in% int1, with=F]
    LD_gwas_coloc = as.matrix(LD_gwas2)
    dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]
    
    # Filter GWAS data to those snps present in LD matrix
    gwas_region <- snp_region[variant %in% colnames(LD_gwas2)]
    # Calculate varbeta and z
    gwas_region[, varbeta := se^2]
    gwas_region[, z := beta / se]
    setkey(gwas_region, variant)
    
    markers <- gwas_region %>%
      select(variant, Chr, Pos, z)
    colnames(markers) = c('marker','chr','pos',"z")
    
    snp = gwas_region$variant[which(gwas_region$pval == min(gwas_region$pval))]
    
    stack_plot <- assoc_plot(markers, LD_gwas_coloc ,title = gwas_name, top.marker = snp, legend=F, labels = "")
    
    dir_path <- paste0("results/Locus_plots/LDblock_check/", gwas_name)
    
    # Check if the directory exists before attempting to create it
    if (!file.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    ggplot2::ggsave(stack_plot, filename = paste0("results/Locus_plots/LDblock_check/",gwas_name,"/",REGION,".png"), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)
    dev.off()
    
  }
  
}
### In case of an error return NA ###
options(warn=1)
check_distance = purrr::possibly(check_distance, otherwise = NA, quiet = F)
results = lapply(1:nrow(sig_regions), function(x) check_distance(x))
