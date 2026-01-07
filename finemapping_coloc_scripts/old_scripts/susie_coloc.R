libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "foreign",
               "purrr")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Supply GWAS summary statistics file")
} else {
  base_gwas <- args[1]
  ld_regions = args[2]
}

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap" # Working directory
setwd(wd)

# Uncomment the following variables to run interactively
ld_regions = "New_multi_coloc_pipeline/lead_loci/chr17.regions"
base_gwas = "UKBB_sumstats/30290_formatted.tsv"

cat("Loading GWAS summary statistcs\n")
# Loads GWAS summary statistcs input ising bash command
maha <- fread(base_gwas, sep=" ")
maha = maha[,Pos:=as.integer(Pos)]

# Obtain name of GWAS summary statistics loaded
gwas_name <- gsub(".*\\/|_formatted.tsv$", "", base_gwas) # Simplify string operations

# Load primary lead variant loci and regions
cat("Loading primary regions\n")
sig_regions <- fread(file = ld_regions, header = F, stringsAsFactors = FALSE, fill=T) %>%
  rename(CHR = V1, lower = V2, upper = V3) %>% mutate_at(c("CHR","lower", "upper"), as.numeric) %>%
  mutate(hg19_region = paste("chr",CHR,":",lower,"-",upper,sep=""))

cat("Analysing GWAS id", gwas_name, "\n")

# Load hg38 converted sig regions, we extracted the hg38 regions and then lifted over for compute efficency
sig_regions2 = fread("GTEX_LD_matrices/approx_LDblocks_hg38.bed")
colnames(sig_regions2) = c("CHR","lower","upper","hg19_region")
sig_regions2 = sig_regions2 %>% filter(CHR == paste("chr",sig_regions$CHR[1], sep=""))

# filter only overlapping regions
sig_regions = sig_regions %>% filter(hg19_region %in% sig_regions2$hg19_region)

if (identical(sig_regions$hg19_region, sig_regions2$hg19_region) == F){
  stop("GWAS and eQTL regions don't match")
}

################### Load Primary GWAS summary statistics ##########################################
cat("Loading eQTL summary statistics\n")
eqtl = fread(paste("eqtl_sumstats/split_by_chrm/GTEx_blood_GRCh37_chr",sig_regions$CHR[1],".txt",sep=""))[,-c(2,4,11,20)]
colnames(eqtl) = c("SNP_GRCh37", "Pos_GRCh37", "SNP_GRCh38","molecular_trait_id", "Chr","Pos_GRCh38", 
                   "ref", "alt",  "ma_samples", "maf", "pvalue", "beta", "se","type", "ac", "an",
                   "molecular_trait_object_id", "gene_id", "median_tpm",  "rsid")

### Start of SUSIE COLOC function ###
# x = row number of file containing: range, lower, upper, chr.
# x = 1209 there should be a significant coloc here with 30,050
# x = 43 there should be a significant multiple hits here with 30,000
run_coloc_susie = function(x){
  
  CHR = sig_regions$CHR[x]
  LOWER = sig_regions$lower[x]
  UPPER = sig_regions$upper[x]
  REGION = paste(sig_regions$lower[x],sig_regions$upper[x],sep=".")
  
  cat("Preparing to run SuSiE COLOC for region",CHR,REGION, "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = maha[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  if (min(snp_region$pval, na.rm=T) > 5e-05){
    stop("No GWS significant variants in this region for this GWAS. Lower p-value= ", min(snp_region$pval, na.rm=T))
  }
  
  cat("Loading GWAS LD and BIM files\n")
  
  # Load LD and bim file data
  LDfilename <- paste("./UKBB_LDmatrices/", CHR, "/",REGION,"/",REGION,sep = "") #
  BIMfilename <- paste("./UKBB_LDmatrices/", CHR, "/",REGION,"/",REGION,sep = "") #
  LD <- fread(paste(LDfilename, "ld", sep = "." ))
  BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
  colnames(BIM) = c("chr", "rsid", "dk", "pos", "alt", "ref")
  BIM$SNP <- paste(BIM$chr, BIM$pos, BIM$alt, BIM$ref, sep = ":")
  
  # Remove SNPs not in sum stats, duplicates and NAs
  colnames(LD) <- BIM$SNP
  LD$SNP <- BIM$SNP
  LD  = LD[apply(LD, 1, function(x) sum(is.na(x))) < nrow(LD)]
  LD = LD[, colnames(LD) %in% c(LD$SNP,"SNP"), with=F]
  LD = LD[complete.cases(LD)]
  LD = LD[, colnames(LD) %in% c(LD$SNP,"SNP"), with=F]
  LD <- LD[LD$SNP %in% snp_region$variant,]
  LD_gwas <- LD[!duplicated(LD$SNP),]
  
  # Remove duplicates from BIM file
  BIM_gwas  = distinct(BIM, SNP, .keep_all = T)
  
  int1 = intersect(snp_region$variant, LD_gwas$SNP)
  
  cat("Loading eQTL LD and BIM files\n")
  
  # Filter sumstats for SNPs in region
  snp_region0 = eqtl[Chr == CHR & inrange(Pos_GRCh37, LOWER, UPPER)]
  
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
    LDfilename <- paste("./GTEX_ADMIX_LD_matrices/", CHR, "/",REGION2,"/",REGION2,sep = "") #x
    BIMfilename <- paste("./GTEX_ADMIX_LD_matrices/", CHR, "/",REGION2,"/",REGION2,sep = "") #
    LD <- fread(paste(LDfilename, "ld", sep = "." ))
    BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
    colnames(BIM) = c("chr", "SNP", "dk", "pos", "alt", "ref")
    BIM$SNP <- paste(BIM$chr, BIM$pos, BIM$ref, BIM$alt, sep = ":")
    colnames(LD) <- BIM$SNP
    LD$SNP <- BIM$SNP
    LD  = LD[apply(LD, 1, function(x) sum(is.na(x))) < nrow(LD)]
    LD = LD[, colnames(LD) %in% c(LD$SNP,"SNP"), with=F]
    LD = LD[complete.cases(LD)]
    LD = LD[, colnames(LD) %in% c(LD$SNP,"SNP"), with=F]
    LD <- LD[LD$SNP %in% snp_region_gene$SNP_GRCh37,]
    LD_eqtl <- LD[!duplicated(LD$SNP),]
    
    # Remove duplicates from BIM file
    BIM_eqtl  = distinct(BIM, SNP, .keep_all = T)
    
    # Intersect with GWAS SNPs
    int2 = intersect(int1, LD_eqtl$SNP)
    
    # Remove NAs from LD and format into matrix
    LD_eqtl = LD_eqtl[LD_eqtl$SNP %in% int2,]
    LD_eqtl = LD_eqtl[, colnames(LD_eqtl) %in% int2, with=F]
    LD_eqtl_coloc = as.matrix(LD_eqtl)
    dimnames(LD_eqtl_coloc)[[1]] <- dimnames(LD_eqtl_coloc)[[2]]
    
    # Filter eqtl data to those snps present in LD matrix
    eqtl_region <- snp_region_gene[SNP_GRCh37 %in% colnames(LD_eqtl)] %>% distinct(SNP_GRCh37, .keep_all=T) %>%
      mutate(varbeta = se^2, z = beta/se) %>%
      arrange(SNP_GRCh37)
    
    if (nrow(eqtl_region) != nrow(LD_eqtl_coloc)){
      stop("Number of SNPs in eQTL regions and LD don't match")
    }
    
    # Remove NAs from GWAS LD and format into matrix
    LD_gwas2 = LD_gwas[LD_gwas$SNP %in% int2,]
    LD_gwas2 = LD_gwas2[, colnames(LD_gwas2) %in% int2, with=F]
    LD_gwas_coloc = as.matrix(LD_gwas2)
    dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]
    
    # Filter GWAS data to those snps present in LD matrix
    gwas_region <- snp_region[variant %in% colnames(LD_gwas2)] %>%
      mutate(varbeta = se^2, z = tstat) %>%
      arrange(variant)
    
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
    check_dataset(D1, req= "LD", warn.minp = 5e-03)
    check_alignment(D1)
    #plot_dataset(D1, main = "eQTL")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie for eQTL summary statistic data \n")
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
    check_dataset(D2, req= "LD",warn.minp = 5e-05)
    check_alignment(D2)
    #plot_dataset(D2, main = "GWAS")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie for GWAS summary statistic data \n")
    S2 = try(runsusie(D2))
    #print(summary(S2))
    
    sets = S2$sets
    fine_res = data.frame()
    
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
    
    if(is.na(test_res[1]) == F){
      tr2 = test_res$summary[PP.H4.abf >= 0.5]
    } else {tr2 = data.table()
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
        
        fil_LD = as.data.frame(t(as.data.frame(LD_coloc)[rownames(LD_coloc) == pal,colnames(LD_coloc) != pal])) %>%
          filter(if_any(everything(), ~. > 0.7)) %>% rownames_to_column("SNP")
        
        df2.2 = filter.(d2, SNP %in% fil_LD$SNP) %>%
          rename.(effect_allele = EA, other_allele = NEA , MAF = EAF, beta = Beta, N = n) %>% 
          mutate.(trait = "T2D") %>%
          mutate(palindromic = if_else(effect_allele %in% c("A","T") & other_allele %in% c("A","T") |
                                         effect_allele %in% c("C","G") & other_allele %in% c("C","G"),1,0 )) %>%
          filter(palindromic == 0) %>%
          left_join(fil_LD, "SNP") %>% 
          rename(hit2 = SNP)
        
        extr_prx = df2.2[which.max(unlist(df2.2[,15])),] %>% select(-palindromic) %>% mutate(original_snp = pal)
        
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
testie = lapply(1:nrow(sig_regions), function(x) run_coloc_susie2(x))
t2 = do.call(rbind, testie)
t2$Pos = str_split_fixed(t2$SNP,":",4)[,2]
t3 = t2 %>%
  mutate(cond_indep = paste(region,cs,sep="."))

fwrite(t3, paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"_coloc_results.txt", sep=""), sep = ",", row.names = F, quote=F)

# Make df on conditionally independent variants
con_indp = t3 %>% filter(PP.H4.abf >= 0.8)

cat(nrow(con_indp), "causal variants, withtin", length(unique(con_indp$region)),"region colocalised for the GWAS",gwas_name,"\n")

cat("Outputting regions that did not converge")

res = fread(paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"_coloc_results.txt", sep=""))
res$region = gsub("\\.", "-",res$region)
res$region = paste(res$Chr, res$region, sep=":")

no_converge = sig_regions$region[!sig_regions$region %in% res$region]
nc = which(!sig_regions$region %in% res$Region)

fwrite(as.list(no_converge), paste("/gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results_V2/",gwas_name,"/",gwas_name,"_chr",sig_regions$CHR[1],"non_converged_loci.txt",sep=""), quote = F, row.names = F, sep=",")

