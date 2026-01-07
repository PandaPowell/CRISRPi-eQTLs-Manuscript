library(bigsnpr)
library(data.table)
library(susieR)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap")

# Inputs
plink_prefix <- "data/Interval_LDmatrices/chr1/100725569.102121621/"   # genotypes.{bed,bim,fam}
pheno_file   <- "data/pheno.tsv"   # tab/space-delimited with IID, y, covariates
chr  <- 7
start <- 55e6
end   <- 56e6
L <- 10

# One-time: convert PLINK bed to bigsnpr .rds/.bk (skip if already done)
if (!file.exists(paste0(plink_prefix, ".rds"))) {
  snp_readBed(paste0(plink_prefix, ".bed"))
}
obj <- snp_attach(paste0(plink_prefix, ".rds"))
G   <- obj$genotypes
map <- obj$map
# map columns: chr, rsid, genetic.dist, physical.pos, allele1, allele2
names(map) <- c("chr", "rsid", "pos_cm", "pos_bp", "a1", "a0")

# Keep variants in the interval
ind_region <- which(map$chr == chr & map$pos_bp >= start & map$pos_bp <= end)

# MAF filter
maf <- snp_MAF(G, ind.col = ind_region)
ind_keep <- ind_region[maf >= 0.01 & maf <= 0.99]

# Simple mean imputation of missing genotypes (in place)
snp_fastImputeSimple(G, ind.col = ind_keep)

# Extract region to an ordinary R matrix (ok for locus-sized p)
X <- G[, ind_keep]         # copies to RAM
colnames(X) <- map$rsid[ind_keep]

# Phenotype and covariates; align to PLINK samples
pheno <- fread(pheno_file)
# bigsnpr sample IDs (often from .fam) are in obj$fam$sample.ID
ids <- obj$fam$sample.ID
pheno <- pheno[match(ids, pheno$IID), ]   # reorder to match genotypes

y <- pheno$y
# Example covariates: sex, age, PC1..PC10 if present
covar_cols <- intersect(c("sex", "age", paste0("PC", 1:10)), names(pheno))
C <- if (length(covar_cols)) as.matrix(pheno[, ..covar_cols]) else NULL

# Residualize and run SuSiE
residualize <- function(M, C) {
  if (is.null(C)) return(M)
  qrC <- qr(cbind(1, C))
  qr.resid(qrC, M)
}
y_res <- residualize(y, C)
X_res <- residualize(X, C)

fit <- susie(X = X_res, y = y_res,
             L = L, coverage = 0.95,
             estimate_residual_variance = TRUE,
             estimate_prior_variance = TRUE)

pip <- fit$pip
cs  <- susie_get_cs(fit, X = X_res, coverage = 0.95, min_abs_corr = 0.3)

result <- data.frame(
  chr = chr,
  pos = map$pos_bp[ind_keep],
  rsid = colnames(X_res),
  pip = pip,
  cs_id = NA_integer_
)
if (length(cs$cs) > 0) {
  for (j in seq_along(cs$cs)) result$cs_id[cs$cs[[j]]] <- j
}
result <- result[order(-result$pip), ]
fwrite(result, sprintf("susie_chr%s_%d_%d.csv", chr, start, end))
