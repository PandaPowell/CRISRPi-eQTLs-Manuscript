sumstats = snp_region
ld_path = gwas_ld_path
# Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".bed"))
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".rds"))
# set sumstats sample size
if(gwas_type=="cc"){
  sumstats$n_eff <- 4 / (1 / n_cases + 1 / (gwas_n-n_cases))
  sumstats$n_case <- sumstats$n_control <- NULL
} else {
  sumstats$n_eff <- gwas_n
}
# Match variants between genotype and sumstats
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map)

# Calculate correlations
LD <- snp_cor(obj.bigSNP$genotypes,ncores = 1)
corr0 <- runonce::save_run(
  snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr], size = 3 / 1000, ncores = NCORES),
  file = paste0("data/corr_hm3_altpop/chr", chr, ".rds")
)