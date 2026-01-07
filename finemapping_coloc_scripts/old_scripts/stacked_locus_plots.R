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
if (length(args) == 0) {
  stop("Supply PATH file as well as SNPs file")
} else {
  DIR <- args[1]
  RES <- args[2]
}

DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/"

setwd(DIR)

# Load results from susie coloc
RES="results/Coloc_results/ALL_GWASs_coloc_results.txt"
coloc_res = fread(file = RES, header = T, stringsAsFactors = FALSE) %>% 
  mutate_at(c("PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf","pval","pvalue"), as.numeric) %>%
  filter(PP.H4.abf > 0.799) #, pval < 6.6e-09, pvalue < 1e-05)

#View(coloc_res[,c("SNP","pip","region","cs","idx1","idx2","gene_id")])

length(unique(coloc_res$molecular_trait_id))
length(unique(paste(coloc_res$region, coloc_res$idx1, sep = ":")))
length(unique(coloc_res$hit1))

# We need the:
gwas_traits = c("30260","30020","30030")
region = "100196651.101199253"
mol_trait = c("ENSG00000146830")
snp = "7:100240296:A:G" 
snp2 = "7:100309544:A:G"

plotRegion = function(gwas_traits, mol_trait, region, snp, snp2){
  
  CHR = as.integer(str_split(snp,":")[[1]][1])
  LOWER = as.integer(str_split(region, "\\.")[[1]][1])
  UPPER = as.integer(str_split(region, "\\.")[[1]][2])
  REGION = region
  
  # Filter GWAS sumstats for SNPs in region
  df_names = paste("data/UKBB_sumstats/",gwas_traits,"_formatted.tsv", sep = "")
  testie = lapply(df_names, fread)
  
  for (i in seq_along(testie)) {
    testie[[i]][["trait"]] <- gwas_traits[i]
  }

  gwas = do.call(rbind, testie) 
  
  snp_region = gwas[Chr == CHR & inrange(Pos, LOWER, UPPER)]
  
  cat("Loading GWAS LD and BIM files\n")
  
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
  
  cat("Loading eQTL summary statistics\n")
  
  eqtl = fread(paste("data/eqtl_sumstats/split_by_chrm/GTEx_blood_GRCh37_chr",CHR,".txt",sep=""))[,-c(2,4,11,20)]
  colnames(eqtl) = c("SNP_GRCh37", "Pos_GRCh37", "SNP_GRCh38","molecular_trait_id", "Chr","Pos_GRCh38", 
                     "ref", "alt",  "ma_samples", "maf", "pvalue", "beta", "se","type", "ac", "an",
                     "molecular_trait_object_id", "gene_id", "median_tpm",  "rsid")
  
  # Filter sumstats for SNPs in region
  snp_region0 = eqtl[Chr == CHR & inrange(Pos_GRCh37, LOWER, UPPER)]
  
  snp_region_gene = snp_region0[gene_id %in% mol_trait]
  
  LD_gwas = LD_gwas[LD_gwas$SNP %in% snp_region_gene$SNP_GRCh37,]
  
  int2 = intersect(int1, LD_gwas$SNP)
  
  LD_gwas2 = LD_gwas[, colnames(LD_gwas) %in% int2, with=F]
  
  eqtl_region <- snp_region_gene[SNP_GRCh37 %in% int2] %>%
    mutate(varbeta = se^2, z = beta/se, uniq_id = paste(SNP_GRCh37, molecular_trait_id,sep=":")) %>%
    distinct(uniq_id, .keep_all = T) %>%
    arrange(SNP_GRCh37) %>% select(SNP_GRCh37, Chr, Pos_GRCh37, z, molecular_trait_id) %>% 
    pivot_wider(names_from = "molecular_trait_id", values_from = "z")
  
  gwas_region <- snp_region[variant %in% int2] %>%
    mutate(varbeta = se^2, z = tstat) %>%
    arrange(variant) %>% select(SNP_GRCh37 = variant, Chr, Pos_GRCh37 = Pos, z, trait) %>% 
    pivot_wider(names_from = "trait", values_from = "z")
  
  markers = eqtl_region[,1:3]
  colnames(markers) = c('marker','chr','pos')
  
  z_scores = merge(eqtl_region[,c(-2,-3)],gwas_region[,c(-2,-3)],  "SNP_GRCh37") %>% 
    dplyr::arrange(SNP_GRCh37) %>% distinct(SNP_GRCh37,.keep_all = T) %>%
    tibble::column_to_rownames(var = "SNP_GRCh37")

  # filter out bad SNPs from LD matrix
  # LD = LD %>% filter.(SNP %in% markers$SNP)
  # LD <- LD[, c(markers$SNP,"SNP"), with=F ]
  # # remove SNP column from LD dataframe
  # LD_coloc = as.matrix( LD[, 1:(ncol(LD)-1)] )
  
  LD_coloc = as.matrix(LD_gwas2)
  
  n_traits <- length(colnames(z_scores))
  # adjust the width of the plot area in the R studio
  stack_plot <- stack_assoc_plot(markers, z_scores, LD_coloc ,traits = colnames(z_scores), top.marker = snp, legend=F, labels = snp2)
  # stack_assoc_plot_save(x = stack_plot, file = paste("./susie_coloc_results/",coloc_range[x],".png",sep=""), n_traits = length(colnames(z_scores)))
  ggplot2::ggsave(stack_plot, filename = paste("results/Locus_plots/",paste(mol_trait, collapse = "_"),".png",sep=""), width = 9, height = 3+3*n_traits, dpi = 300, units = "in", limitsize = F)
  dev.off()
}

plotRegion(gwas_traits = "30290",mol_trait = "ENSG00000132591",region = "21290357.27334244",snp = "17:27180784:C:T", snp2="")

plotRegion(gwas_traits = "30040",mol_trait = "ENSG00000107290",region = "135298842.137041122",snp = "9:136128000:G:C", snp2="")

plotRegion(gwas_traits = "30040",mol_trait = "ENSG00000107290",region = "135298842.137041122",snp = "9:135864436:C:G", snp2="9:135864167:A:G")

stack_plot <- stack_assoc_plot(markers[1800:2799,], z_scores[1800:2799,], LD_coloc[1800:2799,1800:2799] ,traits = colnames(z_scores), top.marker = snp, legend=F, labels = snp2)

# TFR2, MOSPD3, GIGYF1
plotRegion(gwas_traits = c("30260","30020","30030"), mol_trait = "ENSG00000146830",region = "100196651.101199253",snp = "7:100240296:A:G")

plotRegion(gwas_traits = c("30260","30020","30030"), mol_trait = "ENSG00000146830",region = "100196651.101199253",snp = "7:100240296:A:G")

plotRegion(gwas_traits = c("30280","30010","30110"), mol_trait = c("ENSG00000106330","ENSG00000121716"), 
           region = "100196651.101199253",snp = "7:100214762:G:A") # 7:100309171:C:T

# CD52
plotRegion(gwas_traits = c("30190"), mol_trait = c("ENSG00000169442"), 
           region = "25516845.27401867",snp = "1:26647949:G:C", snp2 = "1:26648034:C:T")

# IKZF1
plotRegion(gwas_traits = c("30190"), mol_trait = c("ENSG00000185811"), 
           region = "49212278.51675322",snp = "7:50361683:G:A", snp2 = "")

plotRegion(gwas_traits = c("30100"), mol_trait = c("ENSG00000123405"), 
           region = "53039004.54778823",snp = "12:54685880:C:T", snp2 = "")



