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
ld_regions ="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/gwas_regions/merged_blood_trait_regions.txt"
gwas_path = "data/UKBB_sumstats/30220_formatted.tsv"
ld_path = "data/UKBB_LDmatrices/"
chromosome=9

cat("Loading GWAS summary statistcs\n")
# Loads GWAS summary statistcs input ising bash command
gwas <- fread(gwas_path, sep=" ")
# Obtain name of GWAS summary statistics loaded
gwas_name <- gsub(".*\\/|_formatted.tsv$", "", gwas_path)
dict = c(
  '30000' = 'White blood cell (leukocyte) count',
  '30010' = 'Red blood cell (erythrocyte) count',
  '30020' = 'Haemoglobin concentration',
  '30030' = 'Haematocrit percentage',
  '30040' = 'Mean corpuscular volume',
  '30050' = 'Mean corpuscular haemoglobin',
  '30060' = 'Mean corpuscular haemoglobin concentration',
  '30070' = 'Red blood cell (erythrocyte) distribution width',
  '30080' = 'Platelet count',
  '30090' = 'Platelet crit',
  '30100' = 'Mean platelet (thrombocyte) volume',
  '30110' = 'Platelet distribution width',
  '30120' = 'Lymphocyte count',
  '30130' = 'Monocyte count',
  '30140' = 'Neutrophil count',
  '30150' = 'Eosinophil count',
  '30160' = 'Basophil count',
  '30180' = 'Lymphocyte percentage',
  '30190' = 'Monocyte percentage',
  '30200' = 'Neutrophil percentage',
  '30210' = 'Eosinophil percentage',
  '30220' = 'Basophil percentage',
  '30240' = 'Reticulocyte percentage',
  '30250' = 'Reticulocyte count',
  '30260' = 'Mean reticulocyte volume',
  '30270' = 'Mean sphered cell volume',
  '30280' = 'Immature reticulocyte fraction',
  '30290' = 'High light scatter reticulocyte percentage',
  '30300' = 'High light scatter reticulocyte count'
)
translated_gwas_name = dict[gwas_name]
gwas = gwas[minor_AF > 0.01,]
# Format GWAS Variant identifier in the form "chr:pos:ref:alt", where "ref" is aligned to the forward strand.
gwas[, `:=` (
  chr = as.integer(Chr),
  pos = as.integer(Pos),
  alleles = sub(".*[0-9]:", "\\1",variant))]
gwas[, `:=` (
  a0 = tstrsplit(alleles, ":", fixed = TRUE)[[1]],
  a1 = tstrsplit(alleles, ":", fixed = TRUE)[[2]])][, alleles := NULL]

# Load primary lead variant loci and regions
cat("Loading primary regions\n")
sig_regions <- fread(file = ld_regions, header = F, stringsAsFactors = FALSE, fill=T)
setnames(sig_regions, c("CHR", "lower", "upper"))
sig_regions[,CHR:=as.numeric(CHR)]
sig_regions[,lower:=as.numeric(lower)]
sig_regions[,upper:=as.numeric(upper)]
sig_regions[,hg19_region:=paste("chr",CHR,":",lower,"-",upper,sep="")]
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

cat("Analysing GWAS", gwas_name, "\n")

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
  
  S2 <- tryCatch(
    susie_rss(
      z = gwas_region$z,
      R = LD_gwas_coloc,
      n = gwas_region$n_complete_samples[1],
      L = 10,
      estimate_residual_variance = F,
      prior_variance = 50,
      check_prior = TRUE,
      max_iter = 100
    ),
    warning = function(w) {
      if (grepl("IBSS algorithm did not converge in 100 iterations!", w$message)) {
        message("Warning detected: IBSS algorithm did not converge!")
        return(NULL)  # You can change this to handle the warning appropriately
      }
    },
    error = function(e) {
      message("Error detected:", e$message)
      return(NULL)  # Handle error case
    }
  )
  
  sets = S2$sets
  
  cat("Extracting finemapped SNPs from GWAS summary statistics \n")
  
  fine_res = data.frame()
  
  if (is.null(S2)) {
    
    warning("No credible sets found")
    cat("Decreasing maximum number of causal variants in region from 10 to 1/n")
    
    S2 = try(susie_rss(z = gwas_region$z,
                       R = LD_gwas_coloc,
                       n = gwas_region$n_complete_samples[1],
                       estimate_residual_variance = F,
                       prior_variance = 50,
                       check_prior = TRUE,
                       L=1))
    sets = S2$sets
    
    
    
    if (is.null(sets$cs)) {
      
      warning("Still no credible sets found after lowering L=1")
      
      pr = as.data.frame(cbind(S2$pip[S2$sets$cs$L1], t(S2$lbf_variable))) %>%
        mutate(SNP = names(S2$pip), Chr = CHR, region = sig_regions$hg19_region[x], cs = NA)
      empt = data.frame(matrix(ncol=9,nrow=1))
      pr = bind_cols(pr, empt)
      colnames(pr) = c("pip", "lbf_cs1", "SNP", "Chr", "region", "cs", paste("lbf_cs",seq(2,10), sep=""))
      fine_res = bind_rows(fine_res, pr)
      
    } else {
      
      cat("Converged when L lowered to 1 /n")
      
      pr = as.data.frame(cbind(S2$pip, t(S2$lbf_variable))) %>%
          mutate(SNP = names(S2$pip), Chr = CHR, region = sig_regions$hg19_region[x], cs = 1)
      empt = data.frame(matrix(ncol=9,nrow=1))
      pr = bind_cols(pr, empt)
      colnames(pr) = c("pip", "lbf_cs1", "SNP", "Chr", "region", "cs", paste("lbf_cs",seq(2,10), sep=""))
      fine_res = bind_rows(fine_res, pr)
    }
    
  } else {
  
    # produce credible set plot
    dir_path <- paste0("results/credible_set_plots/", gwas_name)
    
    # Check if the directory exists before attempting to create it
    if (!file.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    png(paste0("results/credible_set_plots/",gwas_name,"/",REGION,".png"),width = 9, height = 3+3*1, res = 900, units = "in")
    susie_plot(S2,"PIP")
    dev.off()
    
    for (i in 1:length(sets$cs_index)) {
      pr = as.data.frame(cbind(S2$pip[sets$cs[[i]]], t(S2$lbf_variable[,sets$cs[[i]]]))) %>%
        mutate(SNP = names(S2$pip[sets$cs[[i]]]), Chr = CHR, region = sig_regions$hg19_region[x], cs = i)
      colnames(pr) = c("pip", paste("lbf_cs",seq(1,10), sep=""), "SNP", "Chr", "region", "cs")
      fine_res = bind_rows(fine_res, pr)
    
    }
  }
  
  cat("Ploting region \n")
  
  markers <- gwas_region %>%
    dplyr::select(variant, Chr, Pos, z) %>% 
    mutate(Chr = as.numeric(Chr))
  colnames(markers) = c('marker','chr','pos',"z")
  
  markers <- markers[match(colnames(LD_gwas2), markers$marker), ]
  
  snp = fine_res$SNP[which(fine_res$lbf_cs1 == max(fine_res$lbf_cs1))]
  
  # limit window size
  if(max(markers$pos)-min(markers$pos) > 5000000){
    lead_pos = markers$pos[markers$marker == snp]
    markers = markers[inrange(pos, lead_pos-2499999, lead_pos+2499999)]
    plot_LD = LD_gwas_coloc[rownames(LD_gwas_coloc) %in% markers$marker,]
    plot_LD = plot_LD[,colnames(plot_LD) %in% markers$marker]
    stack_plot <- fig_region(data = markers, corr = plot_LD ,title = translated_gwas_name, top_marker = snp, build =37)
  } else{
    stack_plot <- fig_region(data = markers, corr = LD_gwas_coloc ,title = translated_gwas_name, top_marker = snp, build =37)
  }
  
  dir_path <- paste0("results/Locus_plots/UKBB/", gwas_name)

  # Check if the directory exists before attempting to create it
  if (!file.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  ggplot2::ggsave(stack_plot, filename = paste0(dir_path,"/",REGION,".png"), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)

  fine_res_out = fine_res %>% rename(variant = SNP)
  fine_res_out = fine_res_out %>% left_join(gwas_region, "variant") %>%  
    mutate(cond_indep = paste(region,cs,sep="."))
  
  return(fine_res_out)
}

### In case of an error return NA ###
options(warn=1)
run_coloc_susie2 = purrr::possibly(run_coloc_susie, otherwise = NA, quiet = F)
results = lapply(1:nrow(sig_regions), function(x) run_coloc_susie2(x))
t2 = do.call(rbind, results) %>% filter(is.na(variant) == F)

fwrite(t2, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/",gwas_name,"/",gwas_name,"_chr",chromosome,"_finemap_results.txt", sep=""), sep = ",", row.names = F, quote=F)

# Produce bed file also
bed = data.frame(Chr = paste0("chr", t2$chr),
                lower = t2$pos-1,
                upper = t2$pos,
                snp = t2$variant,
                cond_indep = t2$cond_indep)

fwrite(bed, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/",gwas_name,"/",gwas_name,"_chr",chromosome,"_finemap_results.bed", sep=""), sep = "\t", row.names = F, quote=F, col.names = F)
