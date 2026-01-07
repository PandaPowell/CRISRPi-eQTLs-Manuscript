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
  
  if ( (abs(sma$Pos - cs$Pos.x) > 2000) & (abs(sma$Pos - cs$Pos.x) < 50000) & cs$pip > 0.8 ){
    print(paste0(lower, ".",upper))
    cat("sma =", sma$variant, " the highest pip =", cs$variant, "\n" )
  } else {
    print("sma is in credible set")
  }
}

find_non_sma = purrr::possibly(find_non_sma, otherwise = NA, quiet = F)
reg= "70855.121359"
res = lapply(unique(finemap$region), find_non_sma)

plot = ggplot(snp_region, aes(x = Pos, y = -log10(pval))) +
  geom_point() + # Add points
  theme_minimal() + # Use a minimal theme
  labs(x = "Position", y = "-log10(p-value)", title = "P-values vs. Positions Plot") # Labeling
ggsave("results/Locus_plots/preLD_plot.pdf", plot, width = 11, height = 8.5) # Save the plot to a PDF file

plotRegion = function(gwas_traits, region, snp, snp2){
  
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
  
  LD_gwas2 = LD_gwas[, colnames(LD_gwas) %in% int1, with=F]
  
  gwas_region <- snp_region[variant %in% int1] %>%
    mutate(varbeta = se^2, z = tstat) %>%
    arrange(variant) %>% select(SNP_GRCh37 = variant, Chr, Pos_GRCh37 = Pos, z, trait) %>% 
    pivot_wider(names_from = "trait", values_from = "z")
  
  markers = gwas_region
  colnames(markers) = c('marker','chr','pos',"z")
  
  LD_coloc = as.matrix(LD_gwas2)
  
  # adjust the width of the plot area in the R studio
  stack_plot <- assoc_plot(markers, LD_coloc ,title = gwas_traits, top.marker = snp, legend=F, labels = snp2)
  # stack_assoc_plot_save(x = stack_plot, file = paste("./susie_coloc_results/",coloc_range[x],".png",sep=""), n_traits = length(colnames(z_scores)))
  ggplot2::ggsave(stack_plot, filename = paste("results/Locus_plots/",paste(gwas_traits,region, sep="_"),".png",sep=""), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)
  dev.off()
}

plotRegion(gwas_traits = "30010",region = "119754110.122007651",snp = "12:121128699:A:G", snp2="")

gwas_traits = "30010"
region = "119754110.122007651"
snp = "12:121128699:A:G"

snp2="12:121128699:A:G"

stack_plot <- assoc_plot(markers[3500:5500,], LD_coloc[3500:5500,3500:5500] ,title = gwas_traits, top.marker = snp, legend=F, labels = "")
