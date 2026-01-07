options(bitmapType="cairo")

libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "Rfast",
               "geni.plots",
               "bigsnpr")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

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
  chromosome= args[4]
}

cat("The following column names are required within the GWAS file (Cap sensitive): \n
    variant \n
    Chr \n
    Pos \n
    minor_allele \n
    minor_AF \n
    n_complete_samples \n
    beta \n
    se \n
    pval \n")

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap" # Working directory
setwd(wd)

#Uncomment the following variables to run interactively
chromosome=20
ld_regions ="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/Interval_LDmatrices/approx_LD_blocks_hg38.bed"
eqtl_path = paste0("/gpfs/commons/datasets/controlled/INTERVAL/public_sumstats/cis/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_nominal_chr",chromosome,".tsv")
ld_path = "data/Interval_LDmatrices/"

check_headers <- function(df, required_headers) {
  # Check which required headers are missing
  missing_headers <- setdiff(required_headers, names(df))
  
  # If there are missing headers, stop the execution and show an error message
  if (length(missing_headers) > 0) {
    stop("Missing required headers: ", paste(missing_headers, collapse=", "))
  } else {
    message("All required headers are present.")
  }
}

eqtl_headers = c("phenotype_id", "variant_id", "af", "ma_samples","ma_count","pval_nominal","beta","slope_se")

cat("Loading eQTL summary statistcs\n")

eqtl = fread(eqtl_path)[,-c("V14","V16")]
if("maf" %in% colnames(eqtl) | "slope" %in% colnames(eqtl) ){
  colnames(eqtl)[colnames(eqtl) == 'maf'] = "af"
  colnames(eqtl)[colnames(eqtl) == 'slope'] = "beta"
}
check_headers(eqtl,eqtl_headers)
eqtl = eqtl[af > 0.005,]

# We need to get the ref/eff allele from the snp id, in tensor this is coded as chr:pos:ref:alt
eqtl[, `:=` (
  chr = as.integer(chr),
  pos = as.integer(pos_b37),
  a0 = effect_allele,
  a1 = other_allele)]

# Load primary lead variant loci and regions
cat("Loading primary regions\n")
sig_regions <- fread(file = ld_regions, header = F, stringsAsFactors = FALSE, fill=T)
setnames(sig_regions, c("CHR", "lower", "upper"))
sig_regions[,CHR:=gsub("chr","",sig_regions$CHR)]
sig_regions[,CHR:=as.numeric(CHR)]
sig_regions[,lower:=as.numeric(lower)]
sig_regions[,upper:=as.numeric(upper)]
sig_regions[,region:=paste("chr",CHR,":",lower,"-",upper,sep="")]
sig_regions = sig_regions[CHR == chromosome,]

load_LD = function(CHR,LOWER,UPPER,ld_path,sumstats){
  
  cat("Loading LD and BIM files\n")
  
  REGION = paste(LOWER,UPPER,sep=".")
  
  if(file.exists(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".bk"))){
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".rds"))
  } else {
    cat("Generating bk & rds files\n")
    # Read from bed/bim/fam, it generates .bk and .rds files.
    snp_readBed(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".bed"))
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".rds"))
  }
  
  # Match variant effects to allele LD was calcualted for
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a0", "a1"))
  sumstats_matched <- snp_match(sumstats, map, join_by_pos = T, remove_dups = T, strand_flip = F)
  # Make new variant id based on matched snps
  sumstats_matched$variant = with(sumstats_matched,paste(chr,pos,a0,a1,sep=":"))
  
  # Define filenames
  LDfilename <- paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION)
  BIMfilename <- paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION)
  # Load LD and bim file data
  LD <- fread(paste(LDfilename, "ld", sep = "." ))
  BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
  # assign headers
  setnames(BIM, c("chr", "rsid", "dk", "pos", "a0", "a1"))
  BIM[, SNP := paste(chr, pos, a0, a1, sep = ":")]
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
  LD <- LD[LD$SNP %in% sumstats_matched$variant,]
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$SNP,"SNP"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "SNP")
  # Remove duplicates from BIM file
  BIM_clean <- BIM[!duplicated(BIM, by = "SNP"), ]
  cat("Finished \n")
  return(list(LD_clean, sumstats_matched))
}

cat("Fine-mapping Interval", "\n")



run_coloc_susie = function(x){
  
  CHR = sig_regions$CHR[x]
  LOWER = sig_regions$lower[x]
  UPPER = sig_regions$upper[x]
  REGION = paste(sig_regions$lower[x],sig_regions$upper[x],sep=".")
  
  cat("Preparing to run SuSiE FINEMAP for region",sig_regions$hg19_region[x], "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = gwas[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  if (min(snp_region$pval, na.rm=T) > 6.6e-09){
    stop("No GWS significant variants in this region for this GWAS. Lower p-value= ", min(snp_region$pval, na.rm=T))
  }
  
  cat("Loading GWAS LD and BIM files\n")
  
  temp = load_LD(CHR = CHR,LOWER = LOWER,UPPER = UPPER,
                 ld_path = ld_path, sumstats = snp_region)
  
  LD_gwas2 = temp[[1]]
  snp_region = as.data.table(temp[[2]])
  
  # Filter GWAS data to those snps present in LD matrix
  LD_gwas2 = LD_gwas2[SNP %in% snp_region$variant,]
  LD_gwas2 = LD_gwas2[,colnames(LD_gwas2) %in% snp_region$variant, with=F]
  LD_gwas_coloc = as.matrix(LD_gwas2)
  dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]
  
  # Calculate varbeta and z
  gwas_region = snp_region
  gwas_region[, varbeta := se^2]
  gwas_region[, z := beta / se]
  setkey(gwas_region, variant)
  
  if (nrow(gwas_region) != nrow(LD_gwas2)){
    stop("Number of SNPs in GWAS regions and LD don't match")
  }
  
  # Order DF
  gwas_region = gwas_region[match(rownames(LD_gwas_coloc), gwas_region$variant),]
  
  if( !identical(gwas_region$variant, rownames(LD_gwas_coloc)) ){
    stop("SNP ordering between sum stats and LD don't match")
  }
  
  # Here we loop through all the molecular_trait_ids within the region and perform coloc
  cat("Running susie for region CHR",CHR,REGION,"\n")
  
  # Format GWAS data
  b2 = c(gwas_region$beta)
  names(b2) = gwas_region$variant
  vb2 = c(gwas_region$varbeta)
  names(vb2) = gwas_region$variant
  maf2 = c(gwas_region$minor_AF)
  names(maf2) = c(gwas_region$variant)
  
  D2 = list(gwas_region$variant, gwas_region$Pos, b2, vb2, maf2, gwas_region$n_complete_samples[1], "quant", LD_gwas_coloc)
  names(D2) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
  check_alignment(D2)
  #plot_dataset(D2, main = "GWAS")
  
  # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
  cat("running susie for GWAS summary statistic data \n")
  S2 = try(susie_rss(z = gwas_region$z,R = LD_gwas_coloc,gwas_region$n_complete_samples[1],
    estimate_residual_variance = FALSE,
    prior_variance = 50,
    check_prior = TRUE))
  sets = S2$sets
  
  cat("Extracting finemapped SNPs from GWAS summary statistics \n")
  
  fine_res = data.frame()
  
  if (is.null(sets$cs)) {
    
    warning("No credible sets found")
      
    pr = data.frame(variant = character(),lbf_1 = numeric(),lbf_2 = numeric(),lbf_3 = numeric(),
      lbf_4 = numeric(),lbf_5 = numeric(),lbf_6 = numeric(),lbf_7 = numeric(),lbf_8 = numeric(),
      lbf_9 = numeric(),lbf_10 = numeric(),Chr = double(),region = character(),stringsAsFactors = FALSE)
    fine_res = bind_rows(fine_res, pr)
    
    } else {
      
    pr = as.data.frame(t(S2$lbf_variable)) %>% rownames_to_column("SNP") %>%
      mutate(Chr = CHR, region = sig_regions$hg19_region[x])
    colnames(pr) = c("variant", paste("lbf_",seq(1,10), sep=""),"Chr", "region")
    fine_res = bind_rows(fine_res, pr)
    
    }
  
  fine_res_out = fine_res %>% left_join(gwas_region, "variant")
  
  return(fine_res_out[,c(1:17,27)])
}

### In case of an error return NA ###
options(warn=1)
run_coloc_susie2 = purrr::possibly(run_coloc_susie, otherwise = NA, quiet = F)
results = lapply(1:nrow(sig_regions), function(x) run_coloc_susie2(x))
t2 = do.call(rbind, results) %>% filter(is.na(variant) == F)

fwrite(t2, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap_lbf/",gwas_name,"/",gwas_name,"_chr",chromosome,"_finemap_results.txt", sep=""), sep = ",", row.names = F, quote=F)

# Produce bed file also
bed = data.frame(Chr = paste0("chr", t2$chr),
                lower = t2$pos-1,
                upper = t2$pos,
                snp = t2$variant)

fwrite(bed, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap_lbf/",gwas_name,"/",gwas_name,"_chr",chromosome,"_finemap_results.bed", sep=""), sep = "\t", row.names = F, quote=F, col.names = F)
