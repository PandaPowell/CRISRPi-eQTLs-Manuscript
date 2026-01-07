libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr",
               "Rfast",
               "geni.plots")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

# if (length(args) == 0) {
#   stop("Supply GWAS summary statistics file")
# } else {
#   base_gwas <- args[1]
#   ld_regions = args[2]
# }

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap" # Working directory
setwd(wd)

# Uncomment the following variables to run interactively
ld_regions = "New_multi_coloc_pipeline/lead_loci/chr11.regions"
base_gwas = "data/UKBB_sumstats/30040_formatted.tsv"

cat("Loading GWAS summary statistcs\n")
# Loads GWAS summary statistcs input ising bash command
gwas <- fread(base_gwas, sep=" ")
gwas = gwas[,Pos:=as.integer(Pos)]

# Obtain name of GWAS summary statistics loaded
gwas_name <- gsub(".*\\/|_formatted.tsv$", "", base_gwas) # Simplify string operations

# Load primary lead variant loci and regions
cat("Loading primary regions\n")

hg19_regions <- fread(file = ld_regions, header = F, stringsAsFactors = FALSE, fill=T) %>%
  rename(CHR = V1, lower = V2, upper = V3) %>% mutate_at(c("CHR","lower", "upper"), as.numeric) %>%
  mutate(hg19_region = paste("chr",CHR,":",lower,"-",upper,sep=""))

cat("Analysing GWAS id", gwas_name, "\n")

# Load hg38 converted sig regions, we extracted the hg38 regions and then lifted over for compute efficency,
# So the names of the files will be the hg38 regions which is why we need the hg38 region file
hg38_regions = fread("data/GTEx_EU_LDmatrices/merged_hg38_blood_regions.txt")
colnames(hg38_regions) = c("CHR","lower","upper","hg19_region")
hg38_regions[,hg19_region:=paste0("chr",hg19_region)]
hg38_regions = hg38_regions %>% filter(CHR == paste("chr",hg19_regions$CHR[1], sep=""))
# filter only overlapping regions
hg19_regions = hg19_regions %>% filter(hg19_region %in% hg38_regions$hg19_region)
# Order the same
hg38_regions <- hg38_regions[order(hg38_regions$hg19_region), ]
hg19_regions = hg19_regions[order(hg19_regions$hg19_region), ]

if (identical(hg19_regions$hg19_region, hg38_regions$hg19_region) == F){
  stop("hg19 and hg38 region files don't match")
}

################### Load eQTL summary statistics per chrosome ##########################################
cat("Loading eQTL summary statistics\n")
eqtl = fread(paste0("data/eQTL_sumstats/GTEx_whlbld.sumstats.chr",hg19_regions$CHR[1],".hg19.txt"))

### Start of SUSIE COLOC function ###
# x = row number of regions file containing: range, lower, upper, chr.
# x = 1209 there should be a significant coloc here with 30,050
# x = 43 there should be a significant multiple hits here with 30,000

for (x in 1:10){
  
  CHR = hg19_regions$CHR[x]
  LOWER = hg19_regions$lower[x]
  UPPER = hg19_regions$upper[x]
  REGION = paste(hg19_regions$lower[x],hg19_regions$upper[x],sep=".")
  
  cat("Preparing to run SuSiE COLOC for region",CHR,REGION, "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = gwas[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  if (min(snp_region$pval, na.rm=T) > 6.6e-09){
    cat("No GWS significant variants in this region for this GWAS. Lower p-value= ", min(snp_region$pval, na.rm=T))
  }
}

# x = 4,5,6

run_coloc_susie = function(x){
  
  CHR = hg19_regions$CHR[x]
  LOWER = hg19_regions$lower[x]
  UPPER = hg19_regions$upper[x]
  REGION = paste(hg19_regions$lower[x],hg19_regions$upper[x],sep=".")
  
  cat("Preparing to run SuSiE COLOC for region",CHR,REGION, "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = gwas[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  if (min(snp_region$pval, na.rm=T) > 6.6e-09){
    stop("No GWS significant variants in this region for this GWAS. Lower p-value= ", min(snp_region$pval, na.rm=T))
  }
  
  cat("Loading GWAS LD and BIM files\n")
  
  # Define filenames
  LDfilename <- paste0("data/UKBB_LDmatrices/", "chr",CHR, "/",REGION,"/",REGION)
  BIMfilename <- paste0("data/UKBB_LDmatrices/", "chr",CHR, "/",REGION,"/",REGION)
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
  
  cat("Loading eQTL LD and BIM files\n")
  
  # Filter eQTL sumstats for SNPs in region
  snp_region0 = eqtl[chr == CHR & inrange(lower_hg19, LOWER, UPPER) ]
  
  if (min(snp_region0$pvalue, na.rm=T) > 5e-03){
    stop("No GWS significant variants in this region for this eQTL Lower p-value= ", min(snp_region0$pvalue, na.rm=T))
  }
  
  genes = unique(snp_region0$gene_id)
  
  coloc_tr_snps = data.frame()
  #  ENSG00000123405 = NFE2
  # Loop through eqtl file names/trait names that we want to test for coloclaisation
  for (name in genes){
    
    cat(name, "\n")
    
    ##### Here we would loop through the genes in the unique gene list #####
    snp_region_gene = snp_region0[gene_id == name]
    
    if (min(snp_region_gene$pvalue, na.rm=T) > 5e-03){
      warning("No GWS significant variants in this region for this eQTL Lower p-value= ", min(snp_region_gene$pvalue, na.rm=T))
      next
    }
    
    # Define hg38 converted positions
    
    LOWER2 = sig_regions2$lower[x]
    UPPER2 = sig_regions2$upper[x]
    REGION2 = paste(sig_regions2$lower[x],sig_regions2$upper[x],sep=".")
  
    # Remove SNPs not in sum stats, duplicates and NAs
    # Load LD and bim file data
    LDfilename <- paste("data/GTEX_ADMIX_LD_matrices/", CHR, "/",REGION2,"/",REGION2,sep = "")
    BIMfilename <- paste("data/GTEX_ADMIX_LD_matrices/", CHR, "/",REGION2,"/",REGION2,sep = "")
    LD <- fread(paste(LDfilename, "ld", sep = "." ))
    BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
    
    # assign headers
    setnames(BIM, c("chr", "rsid", "dk", "pos", "alt", "ref"))
    BIM[, SNP := paste(chr, pos, ref, alt, sep = ":")]
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
    LD <- LD[LD$SNP %in% snp_region_gene$SNP_GRCh37,]
    # Make sure LD only contains columns in LD$SNP and "SNP"
    required_cols <- unique(c(LD$SNP,"SNP"))
    LD <- LD[, ..required_cols]
    # remove duplicates
    LD_eqtl <- unique(LD, by = "SNP")
    # Remove duplicates from BIM file
    BIM_eqtl  = BIM[!duplicated(BIM, by = "SNP"), ]
    # Intersect with GWAS SNPs
    int2 = intersect(int1, LD_eqtl$SNP)
    
    # Remove NAs from LD and format into matrix
    LD_eqtl = LD_eqtl[LD_eqtl$SNP %in% int2,]
    LD_eqtl = LD_eqtl[, colnames(LD_eqtl) %in% int2, with=F]
    LD_eqtl_coloc = as.matrix(LD_eqtl)
    dimnames(LD_eqtl_coloc)[[1]] <- dimnames(LD_eqtl_coloc)[[2]]
    
    # Filter eqtl data to those snps present in LD matrix
    eqtl_region <- snp_region_gene[SNP_GRCh37 %in% colnames(LD_eqtl)] 
    # Remove duplicate SNPs while keeping all columns
    eqtl_region <- unique(eqtl_region, by = "SNP_GRCh37")
    # Calculate varbeta and z
    eqtl_region[, varbeta := se^2]
    eqtl_region[, z := beta / se]
    # Arrange by SNP_GRCh37 - for data.table, setkey() can be used for arranging and faster subsequent joins
    setkey(eqtl_region, SNP_GRCh37)
    
    if (nrow(eqtl_region) != nrow(LD_eqtl_coloc)){
      stop("Number of SNPs in eQTL regions and LD don't match")
    }
    
    # Remove NAs from GWAS LD and format into matrix
    LD_gwas2 = LD_gwas[LD_gwas$SNP %in% int2,]
    LD_gwas2 = LD_gwas2[, colnames(LD_gwas2) %in% int2, with=F]
    LD_gwas_coloc = as.matrix(LD_gwas2)
    dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]
    
    # Filter GWAS data to those snps present in LD matrix
    gwas_region <- snp_region[variant %in% colnames(LD_gwas2)]
    # Calculate varbeta and z
    gwas_region[, varbeta := se^2]
    gwas_region[, z := beta / se]
    setkey(gwas_region, variant)
    
    if (nrow(gwas_region) != nrow(LD_gwas2)){
      stop("Number of SNPs in GWAS regions and LD don't match")
    }
    
    if (nrow(gwas_region) != nrow(eqtl_region)){
      stop("Number of SNPs in GWAS and eQTL region")
    }
    
    cat("Colocalising region CHR",CHR,REGION,"\n")
    
    # Format eQTL data
    b2 = c(eqtl_region$beta)
    names(b2) = eqtl_region$SNP_GRCh37
    vb2 = c(eqtl_region$varbeta)
    names(vb2) = eqtl_region$SNP_GRCh37
    maf2 = c(eqtl_region$maf)
    names(maf2) = c(eqtl_region$SNP_GRCh37)
    
    D1 = list(eqtl_region$SNP_GRCh37, eqtl_region$Pos_GRCh37, b2, vb2, maf2, 838, "cc", LD_eqtl_coloc)
    names(D1) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
    #check_dataset(D1, req= "LD", warn.minp = 5e-03)
    #check_alignment(D1)
    #plot_dataset(D1, main = "eQTL")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie finemap for eQTL summary statistic data \n")
    S1 = try(runsusie(D1))
    #print(summary(S1))
    
    # Format GWAS data
    b2 = c(gwas_region$beta)
    names(b2) = gwas_region$variant
    vb2 = c(gwas_region$varbeta)
    names(vb2) = gwas_region$variant
    maf2 = c(gwas_region$minor_AF)
    names(maf2) = c(gwas_region$variant)
    
    D2 = list(gwas_region$variant, gwas_region$Pos, b2, vb2, maf2, gwas_region$n_complete_samples[1], "cc", LD_gwas_coloc)
    names(D2) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
    #check_dataset(D2, req= "LD",warn.minp = 5e-05)
    #check_alignment(D2)
    #plot_dataset(D2, main = "GWAS")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie for GWAS summary statistic data \n")
    S2 = try(runsusie(D2))
    #print(summary(S2))
    
    sets = S2$sets
    
    fine_res = data.frame()
    
    # Check if credible sets were found
    if (is.null(sets$cs)) {
      
      warning("No credible sets found")
      
      pr = data.frame("pip"=0, "SNP"=0, "Chr"=0, "region"=0, "cs"=0, "trait"= as.character(0))
      fine_res = bind_rows(fine_res, pr[,c(2:6,1)])
      
    } else {
      
      cat("Extracting base GWAS finemapped SNPs from summary statistics \n")
      
      for (i in 1:length(sets$cs_index)) {
        
        pr = as.data.frame(S2$pip[sets$cs[[i]]]) %>%
          mutate(SNP = names(S2$pip[sets$cs[[i]]]), Chr = CHR, region = as.character(REGION), cs = i, trait = gwas_name)
        
        colnames(pr) = c("pip", "SNP", "Chr", "region", "cs", "trait")
        fine_res = bind_rows(fine_res, pr[,c(2:6,1)])
      }
    }
    
    
    if( class(S1) == "try-error" | class(S2) == "try-error" ){
      test_res=NA
    } else{
      test_res = coloc.susie(S1,S2)
    }
    
    # Check if any credible sets colocalised
    if(is.na(test_res[1]) == F){
      tr2 = test_res$summary[PP.H4.abf >= 0.5]
    } else {
      tr2 = data.table()
    }
    
    # If susie cant identify credible sets or doesnt detect coloclaisation then run coloc
    if(is.null(sets$cs)) {
      
      cat("no colocalisation or credible sets identified via SuSiE COLOC, running COLOC \n")
      res = coloc.abf(D1,D2)
      
      coloc_snps = res$results
      colnames(coloc_snps)[1] = "SNP"
      coloc_snps$hit1 = coloc_snps$SNP
      coloc_snps = coloc_snps[which.max(res$results$SNP.PP.H4),c(1,13)]
      coloc_res = as.data.frame(t(res$summary))
      coloc_res = coloc_res[rep(seq_len(nrow(coloc_res)),nrow(coloc_snps)), ]
      coloc_res = cbind(coloc_snps, coloc_res)
      
      df = filter(eqtl_region, SNP_GRCh37 %in% coloc_snps$SNP) %>%
        mutate(trait = name) %>% rename(SNP = SNP_GRCh37)
      
      df2 = filter(gwas_region, variant %in% coloc_snps$SNP) %>% rename(SNP = variant) 
      
      fine_res$SNP = coloc_res$SNP
      
      EX_SNP = fine_res %>% left_join(coloc_res, "SNP") %>% left_join(df2, "SNP") %>%
        left_join(df, "SNP") %>% mutate(region = as.character(REGION)) %>% mutate(idx1 = 1, idx2 = 1)
      
      coloc_tr_snps = bind_rows(EX_SNP, coloc_tr_snps)
      
    } else if(is.na(test_res[1]) |  nrow(tr2) == 0){
      
      cat("no colocalisation identified via SuSiE COLOC, running COLOC \n")
      res = coloc.abf(D1,D2)
      
      coloc_snps = res$results[res$results[,1] %in% fine_res$SNP,c(1,12)]
      colnames(coloc_snps)[1] = "SNP"
      coloc_snps$hit1 = coloc_snps$SNP
      coloc_snps = coloc_snps[which.max(coloc_snps$SNP.PP.H4),c(1,3)]
      coloc_res = as.data.frame(t(res$summary))
      coloc_res = coloc_res[rep(seq_len(nrow(coloc_res)),nrow(coloc_snps)), ]
      coloc_res = cbind(coloc_snps, coloc_res)

      df = filter(eqtl_region, SNP_GRCh37 %in% coloc_snps$SNP) %>%
        mutate(trait = name) %>% rename(SNP = SNP_GRCh37)
      
      df2 = filter(gwas_region, variant %in% coloc_snps$SNP) %>% rename(SNP = variant) 
      
      EX_SNP = fine_res %>% left_join(coloc_res, "SNP") %>% left_join(df2, "SNP") %>%
        left_join(df, "SNP") %>% mutate(region = as.character(REGION)) %>% mutate(idx1 = 1, idx2 = 1)
      
      coloc_tr_snps = bind_rows(EX_SNP, coloc_tr_snps)

  
    } else {
      
      cat("running SuSiE COLOC \n")
      
      res = coloc.susie(S1,S2)
      
      # Check to make sure we have results to work with
      # colnames(res$summary)
      # print(res$summary)
      # print(res$summary[PP.H4.abf > 0.6])
      
      coloc_res = res$summary
      colnames(coloc_res)[3] = "SNP"
      
      cat("Extracting base GWAS colocalised SNPs from summary statistics of both traits \n")
      
      df2 = filter(gwas_region, variant %in% fine_res$SNP) %>% rename(SNP = variant)
      
      # Add if statement here for palindromic snps
      # Need to only filter out those pallindromic with intermediate allele frequencies
      
      # Obtain snps in LD
      find_prox = function(pal){
        
        # Transpose the subset of LD_coloc and filter
        fil_LD <- transpose(LD_coloc[rownames(LD_coloc) == pal, .SD, .SDcols = colnames(LD_coloc) != pal])
        fil_LD <- as.data.table(fil_LD)
        setnames(fil_LD, colnames(LD_coloc)[colnames(LD_coloc) != pal])
        fil_LD <- fil_LD[, SNP := .I][rowSums(.SD > 0.7) > 0, .(SNP)]
        
        # Filter d2 for SNPs in fil_LD and perform transformations
        d2_filtered <- d2[SNP %in% fil_LD$SNP]
        d2_filtered[, c("effect_allele", "other_allele", "MAF", "beta", "N") := .(EA, NEA, EAF, Beta, n)]
        d2_filtered[, trait := "T2D"]
        d2_filtered[, palindromic := fifelse(effect_allele %in% c("A", "T") & other_allele %in% c("A", "T") | 
                                               effect_allele %in% c("C", "G") & other_allele %in% c("C", "G"), 1, 0)]
        d2_filtered <- d2_filtered[palindromic == 0]
        
        # Join fil_LD with d2_filtered to get SNPs
        result <- d2_filtered[fil_LD, on = "SNP"]
        result[, hit2 := SNP]
        
        # Extract the row with the maximum value in the 15th column (assumed here based on your original code)
        # Ensure you adjust the column index if necessary
        extr_prx <- result[which.max(unlist(result[[15]])), .SD, .SDcols = -"palindromic"]
        extr_prx[, original_snp := pal]
        
        return(extr_prx)
        
      }
      
      #find_prox(df2$hit2)
      # sapply(df2$hit2, find_prox)
      # ld_prox = sapply(fil_LD, max, na.rm=T)
      # data.frame(palidromic)
      # sapply(colnames(fil_LD), function(x) rownames(fil_LD)[max_ld[names(max_ld) == x]])
      
      df = filter(eqtl_region, SNP_GRCh37 %in% coloc_res$hit1) %>%
        mutate(trait = name) %>% rename(hit1 = SNP_GRCh37)
      
      EX_SNP = fine_res %>% left_join(coloc_res, "SNP") %>%  
        left_join(df2, "SNP") %>% left_join(df, "hit1") %>% 
        mutate(region = as.character(REGION))
      
      coloc_tr_snps = bind_rows(EX_SNP, coloc_tr_snps)
      
      }
  
    }
  
  return(coloc_tr_snps)
}

options(warn=1)
### In case of an error return NA ### 
run_coloc_susie2 = purrr::possibly(run_coloc_susie, otherwise = NA, quiet = F)

print(system.time({
  run_coloc_susie2(43)
}))

# testie = lapply(1:nrow(sig_regions), function(x) run_coloc_susie2(x))
# t2 = do.call(rbind, testie)
# t2$Pos = str_split_fixed(t2$SNP,":",4)[,2]
# t3 = t2 %>%
#   mutate(cond_indep = paste(region,cs,sep="."))
# 
# fwrite(t3, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"_coloc_results.txt", sep=""), sep = ",", row.names = F, quote=F)
# 
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
# 
