libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "Rfast",
               "geni.plots",
               "bigsnpr",
               "arrow")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Supply necessary files")
} else {
  eqtl_name = args[1]
  chromosome= args[2]
  gwas_path <- args[3]
  eqtl_path = args[4]
}

#Uncomment the following variables to run interactively
gwas_path = "data/UKBB_sumstats/30000_formatted.tsv"
chromosome = 22
eqtl_path = paste0("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/eQTL_sumstats/GTEx_whlbld.sumstats.chr",chromosome,".hg19.txt")
eqtl_name = "GTEx"

# Defined variables
wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap" # Working directory
hg19_regions_path = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/OneK1K_LDmatrices/merged_blood_trait_regions.txt"
output_dir = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/Coloc_results_V2/"
gwas_ld_path = "data/UKBB_LDmatrices/"
eqtl_ld_path = "data/GTEx_EU_LDmatrices/"
gwas_n=350473
gwas_type = "quant" # cc or quant
n_cases = NA
eqtl_n=573
NCORES=10

if(gwas_type == "cc" & is.na(n_cases) == T){
  stop("Provide number of cases")
}

# Loads GWAS summary statistcs input ising bash command
cat("This script is designed to work with output of tensorqtl\n")
# Define the required headers
gwas_headers <- c("variant", "Chr", "Pos", "minor_allele", "beta", "se", "pval")
eqtl_headers = c("phenotype_id", "variant_id", "af", "ma_samples","ma_count","pval_nominal","beta","slope_se")
# Function to check headers in a dataframe
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

setwd(wd)

cat("Loading GWAS summary statistcs\n")

# Loads GWAS summary statistcs input ising bash command
gwas <- fread(gwas_path)
check_headers(gwas,gwas_headers)
# Filter for common alleles
gwas = gwas[minor_AF > 0.01,]
# Format GWAS Variant identifier in the form "chr:pos:ref:alt", where "ref" is aligned to the forward strand.
gwas[, `:=` (
     chr = as.integer(Chr),
     pos = as.integer(Pos),
     alleles = sub(".*[0-9]:", "\\1",variant))]
gwas[, `:=` (
  a0 = tstrsplit(alleles, ":", fixed = TRUE)[[1]],
  a1 = tstrsplit(alleles, ":", fixed = TRUE)[[2]])][, alleles := NULL]

# Obtain name of GWAS summary statistics loaded
gwas_name <- gsub(".*\\/|_formatted.tsv$", "", gwas_path) # Simplify string operations
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
# Load primary lead variant loci and regions
cat("Loading primary regions\n")

hg19_regions <- fread(file = hg19_regions_path, header = T, stringsAsFactors = FALSE, fill=T) %>%
  rename(CHR = chr) %>% mutate_at(c("CHR","lower", "upper"), as.numeric) %>%
  mutate(hg19_region = paste0(CHR,":",lower,"-",upper)) %>% filter(CHR == chromosome)

cat("Analysing GWAS id", gwas_name, "\n")

hg19_regions = hg19_regions[order(hg19_regions$CHR,hg19_regions$lower), ]

################### Load eQTL summary statistics per chrosome ##########################################
cat("Loading eQTL summary statistics\n")
eqtl = fread(eqtl_path)
if("maf" %in% colnames(eqtl) | "slope" %in% colnames(eqtl) ){
  colnames(eqtl)[colnames(eqtl) == 'maf'] = "af"
  colnames(eqtl)[colnames(eqtl) == 'slope'] = "beta"
}
check_headers(eqtl,eqtl_headers)
eqtl = eqtl[af > 0.01,]

# We need to get the ref/eff allele from the snp id, in tensor this is coded as chr:pos:ref:alt
eqtl[, `:=` (
  chr = as.integer(gsub("chr","",str_split_fixed(variant_id,"_",5)[,1])),
  pos = as.integer(upper_hg19),
  a0 = str_split_fixed(variant_id,"_",5)[,3],
  a1 = str_split_fixed(variant_id,"_",5)[,4])]
# x = 4,5,6

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
  sumstats_matched <- snp_match(sumstats, map, join_by_pos = T, remove_dups = F)
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

CHR=chromosome
hg19_lower=hg19_regions$lower[3]
hg19_upper=hg19_regions$upper[3]

run_coloc_susie = function(CHR, hg19_lower, hg19_upper){
  
  hg19_region = paste(hg19_lower,hg19_upper,sep=".")
  
  cat("Preparing to run SuSiE COLOC for region",CHR,hg19_region,"\n")
  
  # Filter GWAS sumstats for SNPs in region
  gwas_region = gwas[chr == CHR & inrange(pos, hg19_lower, hg19_upper)]
  gwas_region$variant = gsub("_",":",gwas_region$variant)
  
  if (min(gwas_region$pval, na.rm=T) > 5e-08){
    stop("No GWS significant variants in this region for this sumstats Lower p-value= ", min(gwas_region$pval, na.rm=T))
  }
  
  temp = load_LD(CHR = CHR,LOWER = hg19_lower,UPPER = hg19_upper,
                    ld_path = gwas_ld_path, sumstats = gwas_region)
  
  gwas_LDmat = temp[[1]]
  gwas_sumstats = as.data.table(temp[[2]])
  
  # Which SNPs intersect
  int1 = intersect(gwas_sumstats$variant, gwas_LDmat$SNP)
  
  cat("Loading eQTL LD and BIM files\n")

  eqtl_sumstats = eqtl[chr == CHR & inrange(pos, hg19_lower, hg19_upper) ]
  
  if (min(eqtl_sumstats$pval_nominal, na.rm=T) > 1e-03){
    stop("No GWS significant variants in this region for this eQTL Lower p-value= ", min(eqtl_sumstats$pval_nominal, na.rm=T))
  }

  temp2 = load_LD(CHR = CHR,LOWER = hg19_lower,UPPER = hg19_upper,
                    ld_path = eqtl_ld_path, sumstats = eqtl_sumstats)
  
  eqtl_LDmat = temp2[[1]]
  eqtl_sumstats = as.data.table(temp2[[2]])
  
  genes = unique(eqtl_sumstats$phenotype_id)
  
  coloc_tr_snps = data.frame()
  
  # Loop through eqtl file names/trait names that we want to test for coloc
  for (name in genes){
    
    ##### Here we would loop through the genes in the unique gene list #####
    snp_region_gene = eqtl_sumstats[phenotype_id == name,]
    
    if (min(snp_region_gene$pval_nominal, na.rm=T) > 1e-03){
      warning("No GWS significant variants in this region for the genes ",name, " Lower p-value= ", min(snp_region_gene$pval_nominal, na.rm=T))
      next
    } else{
      cat("significant variants found for gene",name,"\n")
    }
    
    # Intersect with GWAS SNPs
    int2 = intersect(int1, eqtl_LDmat$SNP)
    
    # Remove NAs from LD and format into matrix
    eqtl_LD = eqtl_LDmat[SNP %in% int2 & SNP %in% snp_region_gene$variant,]
    eqtl_LD = eqtl_LD[, colnames(eqtl_LD) %in% int2 & colnames(eqtl_LD) %in% snp_region_gene$variant, with=F]
    LD_eqtl_coloc = as.matrix(eqtl_LD)
    dimnames(LD_eqtl_coloc)[[1]] <- dimnames(LD_eqtl_coloc)[[2]]
    
    # Filter eqtl data to those snps present in LD matrix
    eqtl_region <- snp_region_gene[variant %in% colnames(eqtl_LD)]
    # Remove duplicate SNPs while keeping all columns
    eqtl_region <- unique(eqtl_region, by = "variant")
    # Calculate varbeta and z
    eqtl_region[, varbeta := slope_se^2]
    eqtl_region[, z := beta / slope_se]
    # Arrange by SNP_GRCh37 - for data.table, setkey() can be used for arranging and faster subsequent joins
    setkey(eqtl_region, variant)
    
    if (nrow(eqtl_region) != nrow(LD_eqtl_coloc)){
      stop("Number of SNPs in eQTL regions and LD don't match")
    }
    
    # Remove NAs from GWAS LD and format into matrix
    LD_gwas2 = gwas_LDmat[SNP %in% int2 & SNP %in% eqtl_region$variant,]
    LD_gwas2 = LD_gwas2[, colnames(LD_gwas2) %in% int2 & colnames(LD_gwas2) %in% eqtl_region$variant, with=F]
    LD_gwas_coloc = as.matrix(LD_gwas2)
    dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]
    
    # Filter GWAS data to those snps present in LD matrix
    gwas_region <- gwas_sumstats[variant %in% colnames(LD_gwas2)]
    # Calculate varbeta and z
    gwas_region[, varbeta := se^2]
    gwas_region[, z := beta / se]
    setkey(gwas_region, variant)
    
    if (nrow(gwas_region) != nrow(LD_gwas_coloc)){
      stop("Number of SNPs in GWAS regions and LD don't match")
    }
    
    if (nrow(gwas_region) != nrow(eqtl_region)){
      stop("Number of SNPs in GWAS and eQTL region")
    }
    
    # Allign to LD matrix
    gwas_region = gwas_region[match(rownames(LD_gwas_coloc), gwas_region$variant),]
    eqtl_region = eqtl_region[match(rownames(LD_eqtl_coloc), eqtl_region$variant),]
    
    cat("Colocalising region CHR",CHR,hg19_lower, hg19_upper,"and gene",name,"\n")
    
    # Format eQTL data
    b2 = c(eqtl_region$beta)
    names(b2) = eqtl_region$variant
    vb2 = c(eqtl_region$varbeta)
    names(vb2) = eqtl_region$variant
    maf2 = c(eqtl_region$af)
    names(maf2) = c(eqtl_region$variant)
    
    D1 = list(eqtl_region$variant, eqtl_region$pos, b2, vb2, maf2, eqtl_n, "quant", LD_eqtl_coloc)
    names(D1) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
    #check_dataset(D1, req= "LD", warn.minp = 5e-03)
    #check_alignment(D1)
    #plot_dataset(D1, main = "eQTL")
    cat("running susie finemap for eQTL summary statistic data \n")
    S1 = suppressWarnings(try(runsusie(D1,
                                       estimate_residual_variance = TRUE)))
    
    # Format GWAS data
    b2 = c(gwas_region$beta)
    names(b2) = gwas_region$variant
    vb2 = c(gwas_region$varbeta)
    names(vb2) = gwas_region$variant
    maf2 = c(gwas_region$minor_AF)
    names(maf2) = c(gwas_region$variant)
    
    if (gwas_type == "cc"){
      prop_cases = n_cases/gwas_n
      D2 = list(gwas_region$variant, gwas_region$Pos, b2, vb2, gwas_n, "cc", prop_cases, LD_gwas_coloc)
      names(D2) = c("snp", "position", "beta", "varbeta", "N", "type", "s","LD")
    } else {
      D2 = list(gwas_region$variant, gwas_region$Pos, b2, vb2,maf2, gwas_n, "quant", LD_gwas_coloc)
      names(D2) = c("snp", "position", "beta", "varbeta","MAF", "N", "type","LD")
    }
    #check_dataset(D2, req= "LD",warn.minp = 5e-05)
    #check_alignment(D2)
    #plot_dataset(D2, main = "GWAS")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie for GWAS summary statistic data \n")
    S2 = suppressWarnings(try(runsusie(D2, 
                                       estimate_residual_variance = F)))
    #print(summary(S2))
    sets = S2$sets
    eqtl_sets = S1$sets
    
    # if no credible sets lower coverage
    if(is.null(S2$sets$cs)){
      S2 = suppressWarnings(try(runsusie(D2,
                                         coverage=0.1,
                                         estimate_residual_variance = F)))
      sets = S2$sets
    }
    
    if(is.null(S1$sets$cs)){
      S1 = suppressWarnings(try(runsusie(D1,
                                         coverage=0.1, 
                                         estimate_residual_variance = TRUE)))
      eqtl_sets = S1$sets
    }
    
    # If for what ever reason finemapping fails run coloc
    if(is.null(sets$cs) | is.null(eqtl_sets$cs) | class(S1) == "try-error" | class(S2) == "try-error") {
      
      cat("no credible sets identified via SuSiE COLOC, running COLOC \n")
      res = suppressWarnings(coloc.abf(D1,D2))
      
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
        
        df = filter(eqtl_region, variant %in% coloc_snps$SNP) %>%
          mutate(trait = name) %>% rename(SNP = variant)
        
        df2 = filter(gwas_region, variant %in% coloc_snps$SNP) %>% rename(SNP = variant) 
        
        EX_SNP = coloc_res %>% left_join(df2, "SNP") %>%
          left_join(df, "SNP") %>% mutate(region = paste0(hg19_lower,"-",hg19_upper)) %>%
          mutate( method = "coloc")
        
        coloc_tr_snps = bind_rows(coloc_tr_snps, EX_SNP)
        
        cat("Sucessful Colocalization Ploting region",hg19_region,name,"\n")
        
        identical(gwas_region$variant, eqtl_region$variant)
        
        markers <- data.frame(marker = gwas_region$variant,
                              chr = CHR,
                              pos = gwas_region$Pos,
                              z_1 = gwas_region$z,
                              z_2 = eqtl_region$z)
        
        markers <- markers[match(colnames(LD_gwas2), markers$marker), ]
        
        snp = coloc_tr_snps$SNP[which(coloc_tr_snps$PP.H4.abf == max(coloc_tr_snps$PP.H4.abf, na.rm = T))]
        
        stack_plot <- fig_region_stack(data = markers, trait = c(translated_gwas_name, name), corr = LD_gwas_coloc,
                                       top_marker = snp, build =37, title_center = TRUE)
        
        dir_path <- paste0(output_dir,"Locus_plots/", eqtl_name, "/",gwas_name)
        
        # Check if the directory exists before attempting to create it
        if (!file.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
        }
        
        ggplot2::ggsave(stack_plot, filename = paste0(dir_path,"/",CHR,"_",hg19_region,"_",name,".png"), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)
        
      }
      
    } else {
      
      cat("running SuSiE COLOC \n")
      
      res = coloc.susie(S2,S1)
      
      if(max(res$summary[,"PP.H4.abf"]) > 0.5){
        
        cat("Sucessfull coloc",hg19_region,name,gwas_name, "\n")
        
        coloc_snps = res$results
        coloc_res = res$summary
        colnames(coloc_res)[c(2,9,10)] = c("SNP","gwas_credible_set","eqtl_credible_set")
        coloc_snps = coloc_snps[snp %in% coloc_res$SNP,]
        colnames(coloc_snps)[1] = "SNP"
        
        coloc_res = coloc_res[PP.H4.abf > 0.5,] %>% left_join(coloc_snps[,1:2], "SNP")
        colnames(coloc_res)[11] = "SNP.PP.H4"
        
        # Order columns of dataframe
        coloc_res = coloc_res[,c("SNP","hit2","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","nsnps","gwas_credible_set","eqtl_credible_set","SNP.PP.H4")]

        df2 = filter(gwas_region, variant %in% coloc_res$SNP) %>% rename(SNP = variant)
        
        # Add if statement here for palindromic snps
        # Need to only filter out those pallindromic with intermediate allele frequencies
        
        #find_prox(df2$hit2)
        # sapply(df2$hit2, find_prox)
        # ld_prox = sapply(fil_LD, max, na.rm=T)
        # data.frame(palidromic)
        # sapply(colnames(fil_LD), function(x) rownames(fil_LD)[max_ld[names(max_ld) == x]])
        
        df = filter(eqtl_region, variant %in% coloc_res$hit2) %>%
          mutate(trait = name) %>% rename(hit2 = variant)
        
        EX_SNP = coloc_res %>%  
          left_join(df2, "SNP") %>% left_join(df, "hit2") %>% 
          mutate(region = paste0(hg19_lower,"-",hg19_upper), method = "susie_coloc")
        
        coloc_tr_snps = bind_rows(coloc_tr_snps, EX_SNP)
        
        cat("Sucessful Colocalization Ploting region",hg19_region,name,"\n")
        
        identical(gwas_region$variant, eqtl_region$variant)
        
        markers <- data.frame(marker = gwas_region$variant,
                              chr = CHR,
                              pos = gwas_region$Pos,
                              z_1 = gwas_region$z,
                              z_2 = eqtl_region$z)
        
        markers <- markers[match(colnames(LD_gwas2), markers$marker), ]
        
        snp = coloc_res$SNP
        
        stack_plot <- fig_region_stack(data = markers, trait = c(gwas_name, name), corr = LD_gwas_coloc,
                                       top_marker = snp, build =37, title_center = TRUE)
        
        dir_path <- paste0(output_dir,"Locus_plots/", eqtl_name, "/",gwas_name)
        
        # Check if the directory exists before attempting to create it
        if (!file.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
        }
        
        ggplot2::ggsave(stack_plot, filename = paste0(dir_path,"/",CHR,"_",hg19_region,"_",name,".png"), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)
        
      }

    }
    
  }
  
   return(coloc_tr_snps)
  
}

options(warn=1)
### In case of an error return NA ### 
run_coloc_susie2 = purrr::possibly(run_coloc_susie, otherwise = NA, quiet = F)
testie = lapply(1:nrow(hg19_regions), function(x) run_coloc_susie2(chromosome, 
                                                                   hg19_regions$lower[x], 
                                                                   hg19_regions$upper[x]))
t2 = do.call(rbind, testie)
t3 = filter(t2, is.na(SNP) == F) %>%
  mutate(gwas = translated_gwas_name)

# Write results
res_dir = paste0(output_dir,eqtl_name,"/", gwas_name)
system(paste0("bash -c 'mkdir -p ", res_dir, "'"))

fwrite(t3, paste0(res_dir,"/chr",chromosome,"_coloc_results.txt"), sep = ",", row.names = F, quote=F)

# # Make df on conditionally independent variants
# con_indp = t3 %>% filter(PP.H4.abf >= 0.8)
# 
# cat(nrow(con_indp), "causal variants, withtin", length(unique(con_indp$region)),"region colocalised for the GWAS",gwas_name,"\n")
# 
# cat("Outputting regions that did not converge")
# 
# res = fread(paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"_coloc_results.txt", sep=""))
# res$region = gsub("\\.", "-",res$region)
# res$region = paste(res$Chr, res$region, sep=":")
# 
# no_converge = sig_regions$region[!sig_regions$region %in% res$region]
# nc = which(!sig_regions$region %in% res$Region)
# 
# fwrite(as.list(no_converge), paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"non_converged_loci.txt",sep=""), quote = F, row.names = F, sep=",")
