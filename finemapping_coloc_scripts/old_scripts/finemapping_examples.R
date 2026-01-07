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

# Argument 1 - Pathway to directory
# Argument 2 - SuSie coloc results file
# if (length(args) == 0) {
#   stop("Supply PATH file as well as SNPs file")
# } else {
#   DIR <- args[1]
#   RES <- args[2]
# }

DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/"

setwd(DIR)

finemap = fread("results/UKBB_SuSiE_finemap/30010_total_finemap_results.txt", colClasses = list(character = c('region')))
finemap[, region := as.character(region)]
sum = fread("data/UKBB_sumstats/30010_formatted.tsv")

find_non_sma = function(reg){
  
  finemap_region = finemap[region == reg, ]
  CHR = finemap_region$Chr.x[1]
  lower = as.integer(str_split(finemap_region$cond_indep[1],"\\.")[[1]][1])
  upper = as.integer(str_split(finemap_region$cond_indep[1],"\\.")[[1]][2])
  REGION = paste0(lower, ".", upper)
  
  snp_region = sum[Chr == CHR & inrange(Pos, lower, upper)]
  
  # Load LD and bim file data
  LDfilename <- paste("data/UKBB_LDmatrices/", CHR, "/",REGION,"/",REGION,sep = "") #
  BIMfilename <- paste("data/UKBB_LDmatrices/", CHR, "/",REGION,"/",REGION,sep = "") #
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
  
  int1 = intersect(snp_region$variant, LD_gwas$SNP)
  
  LD_gwas2 = LD_gwas[, colnames(LD_gwas) %in% int1, with=F]
  
  snp_region <- snp_region[variant %in% int1]
  
  sma = snp_region[which(snp_region$pval == min(snp_region$pval, na.rm = T))]
  cs = finemap_region[which(finemap_region$pip == max(finemap_region$pip, na.rm = T))]
  
  if ( (abs(sma$Pos - cs$Pos.x) > 2000) & (abs(sma$Pos - cs$Pos.x) < 50000) & cs$pip > 0.1 & (sma$pval/cs$pval) < 0.1 ){
    print(paste0(lower, ".",upper))
    cat("sma =", sma$variant, " the highest pip =", cs$variant, " ", abs(sma$Pos - cs$Pos.x),"\n")
  } else {
  }
}

find_non_sma = purrr::possibly(find_non_sma, otherwise = NA, quiet = F)
reg= "70855.121359"
res = lapply(unique(finemap$region), find_non_sma)
sig_snps = c("11:30777014:G:GA", "11:111196937:C:T", "12:121128699:A:G", "14:64726503:A:G", "14:65487694:T:C", 
             "15:76461656:G:A", "16:51188432:G:A", "16:72105965:T:C", "1:23851680:AC:A", "1:147290790:G:A")

res2 = lapply(unique(finemap$region[finemap$variant %in% sig_snps]), find_non_sma)
