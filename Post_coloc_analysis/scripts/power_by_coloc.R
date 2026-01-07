options(bitmapType="cairo")
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

.libPaths("/gpfs/commons/home/sghatan/R/x86_64-pc-linux-gnu-library/4.4")
print(.libPaths())  # Check if the path is correctly set

# List of required packages
packages <- c("data.table", "tidyverse", "scPower", "cowplot", "ggplot2", "readxl")

# Install and load each package
lapply(packages, install_if_missing)

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/")

gene_list <- read_excel("Post_coloc_analysis/data/interval_genes.xlsx",sheet = 1) %>%
  filter(gene_biotype == "protein_coding")

power.files = system("ls Post_coloc_analysis/power_results/power.results.beta0.25*", intern = T)
eqtl.power = data.table()

for (i in 1:length(power.files)){
  
  name = str_split_fixed(pattern = "\\.",string = power.files[i], 5)[,4]
  name = gsub(pattern = "25", replacement = "",x = name)
  tmp = fread(power.files[i])
  tmp$study = name
  tmp$effect_size = 0.25
  eqtl.power = bind_rows(eqtl.power,tmp)
  
}

just_top = eqtl.power %>%
  arrange(desc(eQTL.power)) %>%
  distinct(grna_target,ensembl_id, .keep_all = T)

fwrite(just_top, "Post_coloc_analysis/power_results/eqtl.power.combined_just_top.txt", quote = F, row.names = F)

eqtl.power = eqtl.power[ensembl_id %in% gene_list$feature_id]

# Loop through studies and calculate proportion of CRE-genes that have power >0.79
proportions = data.frame()

for (study_name in unique(eqtl.power$study)){
  tmp = eqtl.power[study == study_name,]
  N_80 = tmp[eQTL.power > 0.8,]
  if(study_name == "Interval"){
    n_tmp = 4732
  } else{
    n_tmp = unlist(fread(paste0("eQTL_catalogue/count_means/",study_name,".csv"))[1,3])
  }
  tmp_df = data.frame(study = study_name,
                      n_above_80 = nrow(N_80)/nrow(tmp),
                      sample_size = n_tmp)
  proportions = rbind(proportions, tmp_df)
}

# assuming your dataframe is called proportions
ggplot(proportions, aes(x = sample_size, y = n_above_80)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred") +
  labs(
    x = "Sample size",
    y = "Proportion above 80% power",
    title = "Trend of Power vs. Sample Size across Studies"
  ) +
  theme_minimal()

### Means ###
means = data.frame()

for (study_name in unique(eqtl.power$study)){
  tmp = eqtl.power[study == study_name,]
  mean = mean(tmp$eQTL.power, na.rm = T)
  if(study_name == "Interval"){
    n_tmp = 4732
  } else{
    n_tmp = unlist(fread(paste0("eQTL_catalogue/count_means/",study_name,".csv"))[1,3])
  }
  tmp_df = data.frame(study = study_name,
                      mean_power = mean,
                      sample_size = n_tmp)
  means = rbind(means, tmp_df)
}

png("Post_coloc_analysis/plots/interval/power_vs_samplesize.png", width = 6, height = 4, res = 300, units = "in")
# assuming your dataframe is called proportions
ggplot(means, aes(x = sample_size, y = mean_power)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred") +
  labs(
    x = "Sample size",
    y = "Mean eQTL power",
    title = "Trend of Power vs. Sample Size across Studies"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 22, margin = margin(t = 8)),
    axis.title.y = element_text(size = 22, margin = margin(r = 8)),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    plot.title   = element_text(size = 20)
  )
dev.off()

gtex = data.table()

for (id in seq(30000,30300, 10)){
  
  for (i in c(1:22)){
    
    if(file.exists(paste0("results/Coloc_results_V2/GTEx/",id,"/chr",i,"_coloc_results.txt"))){
      
      gtex_chr = fread(paste0("results/Coloc_results_V2/GTEx/",id,"/chr",i,"_coloc_results.txt"))
      gtex = rbind(gtex, gtex_chr)
      
    } else{
      next()
    }
    
  }
  
}

gtex = gtex %>%
  filter(pval_nominal < 1e-03, pval < 1e-05)

mean(abs(gtex$beta.y))

interval = data.table()

for (id in seq(30000,30300, 10)){
  
  for (i in c(1:22)){
    
    if(file.exists(paste0("results/Coloc_results_V2/Interval/",id,"/chr",i,"_coloc_results.txt"))){
      
      interval_chr = fread(paste0("results/Coloc_results_V2/Interval/",id,"/chr",i,"_coloc_results.txt"))
      interval = rbind(interval, interval_chr, fill=T)
      
    } else{
      next()
    }
    
  }
  
}

interval = interval %>%
  filter(pval_nominal < 1e-03, pval < 1e-05)

cat("the ",mean(abs(interval$beta.y)))

cat("the ",mean( c(abs(interval$beta.y), abs(gtex$beta.y)) ))

cat("the difference in the number of colocalizations is ",nrow(interval)/nrow(gtex),"x")

cat("the difference in the sample size is ",4732/573,"x")

# eQTL catalogue
eqtl_cat = data.table()

for (id in seq(30000,30300, 10)){
    
    if(file.exists(paste0("eQTL_catalogue/coloc_results/",id,"_coloc_results.txt"))){
      
      eqtl_tmp = fread(paste0("eQTL_catalogue/coloc_results/",id,"_coloc_results.txt"))
      eqtl_cat = rbind(eqtl_cat, eqtl_tmp, fill=T)
      
    } else{
      next()
    }
    
}

# load eqtl cat meta data
meta = fread("eQTL_catalogue/eQTL_catalogue_datasets.txt")

# Add on sumstats to get eqtl p-values
eqtl_cat_table = eqtl_cat %>% 
  filter(gwas_pval < 1e-08) %>%
  group_by(dataset_id) %>%
  summarise(n_coloc = n()) %>%
  left_join(meta, "dataset_id")

other_table = data.frame(n_coloc = c(nrow(interval),nrow(gtex)),
                         sample_size = c(4732,573))

plot_table = bind_rows(eqtl_cat_table[,c("n_coloc","sample_size")], other_table)

png("Post_coloc_analysis/plots/interval/sample_size_vs_ncoloc.png", width = 6, height = 6, res = 300, units = "in")

ggplot(plot_table, aes(x = sample_size, y = n_coloc)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  labs(
    x = "Sample size",
    y = "Number of colocalized variants",
    title = "Colocalizations vs. sample size across studies"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 22, margin = margin(t = 8)),
    axis.title.y = element_text(size = 22, margin = margin(r = 8)),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    plot.title   = element_text(size = 20)
  )

dev.off()

library(dplyr)
library(mgcv)

df <- plot_table %>% arrange(sample_size)
idx_max <- which.max(df$sample_size)
test_pt <- df[idx_max, , drop = FALSE]
train   <- df[-idx_max, , drop = FALSE]

# --- OPTION A: linear model on training points ---
lm_fit <- lm(n_coloc ~ sample_size, data = train)
pi_lm  <- predict(lm_fit, newdata = test_pt, interval = "prediction", level = 0.95)
outside_lm <- with(test_pt, n_coloc < pi_lm[1,"lwr"] || n_coloc > pi_lm[1,"upr"])

# Influence of the max point in the full linear fit
lm_all <- lm(n_coloc ~ sample_size, data = df)
t_out   <- rstudent(lm_all)[idx_max]
p_out   <- 2*pt(abs(t_out), df = df.residual(lm_all), lower.tail = FALSE)
cook_out<- cooks.distance(lm_all)[idx_max]
lev_out <- hatvalues(lm_all)[idx_max]

cat(sprintf("[LM] Last point outside 95%% PI? %s  (studentized resid = %.2f, p = %.3g, Cook's D = %.3g, leverage = %.3g)\n",
            outside_lm, t_out, p_out, cook_out, lev_out))

# Plot (LM): LOOCV line + 95% PI, highlight last point
# rebuild newx if needed
# assume df, train, test_pt, lm_fit already defined
x_start <- min(train$sample_size)
x_end   <- test_pt$sample_size  # extend to the max point's x

newx <- data.frame(sample_size = seq(x_start, x_end, length.out = 300))
pred  <- predict(lm_fit, newdata = newx, se.fit = TRUE)
sigma <- summary(lm_fit)$sigma
tcrit <- qt(0.975, df = df.residual(lm_fit))

newx$fit  <- pred$fit
newx$pi_l <- pred$fit - tcrit * sqrt(pred$se.fit^2 + sigma^2)
newx$pi_u <- pred$fit + tcrit * sqrt(pred$se.fit^2 + sigma^2)

# predicted value at the last x
y_pred_last <- as.numeric(predict(lm_fit, newdata = test_pt))

png("Post_coloc_analysis/plots/interval/leave_one_out_regression.png", width = 10, height = 10, res = 300, units = "in")

p_lm <- ggplot(train, aes(x = sample_size, y = n_coloc)) +
  geom_point(size = 3,alpha = 0.5) +
  geom_ribbon(data = newx,
              aes(x = sample_size, ymin = pi_l, ymax = pi_u),
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = newx,
            aes(x = sample_size, y = fit),
            inherit.aes = FALSE, color = "darkred", linewidth = 1) +
  # predicted value at max-N (open circle)
  geom_point(data = data.frame(sample_size = x_end, n_coloc = y_pred_last),
             aes(x = sample_size, y = n_coloc),
             inherit.aes = FALSE, shape = 21, size = 4, stroke = 1.1,
             fill = "white", color = "darkred") +
  # observed max-N point (solid red)
  geom_point(data = test_pt,
             aes(x = sample_size, y = n_coloc),
             inherit.aes = FALSE, color = "red", size = 4) +
  # dashed residual showing deviation at max-N
  geom_segment(aes(x = x_end, xend = x_end, y = y_pred_last, yend = test_pt$n_coloc),
               inherit.aes = FALSE, linetype = "dashed", color = "red") +
  labs(title = "",
       x = "Sample size", y = "Number of colocalized variants") +
  theme_cowplot(font_size = 24)

p_lm
dev.off()

# compare power of CREs without a cGene
cres_w_grnas = fread("Post_coloc_analysis/processed_data/cres_with_grnas.txt")
cres_w_grnas$target_gene = paste0(cres_w_grnas$grna_target,"_",cres_w_grnas$ensembl_id)
cres_w_grnas_egene = fread("Post_coloc_analysis/processed_data/cres_with_grna_eqtls_interval.txt")
cres_w_grnas_egene$target_gene = paste0(cres_w_grnas_egene$grna_target,"_",cres_w_grnas_egene$ensembl_id)

total_cres = unique(cres_w_grnas$grna_target)
cat("Total cres tested =",length(total_cres),"\n")

grna_cres_w_egenes = unique(cres_w_grnas_egene$grna_target)
cat("Total number CREs with egenes =",length(grna_cres_w_egenes),"\n")

target_egenes = unique(cres_w_grnas_egene$target_gene)
cat("Total number of target-egene pairs =",length(target_egenes),"\n")

egenes = unique(cres_w_grnas_egene$ensembl_id)
cat("Total number of egenes =",length(egenes),"\n")

cres_no_target = unique(cres_w_grnas$grna_target[!cres_w_grnas$grna_target %in% c(grna_cres_w_egenes)])
cat("CREs without target genes =",length(cres_no_target),"\n")

# Filter eQTL power to CREs without a target
no_target_snps = unique(cres_w_grnas$grna_target[!cres_w_grnas$grna_target %in% c(grna_cres_w_egenes)])

no_egene_power = eqtl.power %>%
  filter(grna_target %in% no_target_snps) %>%
  arrange(desc(eQTL.power)) %>%
  distinct(grna_target,ensembl_id, .keep_all = T)

egene_power = eqtl.power %>%
  filter(!grna_target %in% no_target_snps) %>%
  arrange(desc(eQTL.power)) %>%
  distinct(grna_target,ensembl_id, .keep_all = T)

mean(no_egene_power$eQTL.power)
mean(egene_power$eQTL.power)

test = wilcox.test(no_egene_power$eQTL.power, 
                   egene_power$eQTL.power)
print(test)

# Calculate the sample size needed to detect the mean interval coloc beta
cat("the mean beta",mean(abs(interval$beta.y)))

cres_no_egene = cres_w_grnas[!grna_target %in% c(grna_cres_w_egenes)] %>%
  distinct(grna_target, .keep_all = T)

# flip alleles to match gwas snp ids
temp = gsub("_", ":", cres_no_egene$finemap_snp_intersect_grna)
temp = str_split_fixed(temp,":",4)[,2]
cres_no_egene$snp_pos = paste0(cres_no_egene$chr, ":", temp)
cres_no_egene$pos = temp

# Load GWAS summary stats to obtain MAFs
sumstats = fread("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/UKBB_sumstats/30000_formatted.tsv")
sumstats$snp_pos = paste0(sumstats$Chr,":",sumstats$Pos)
sumstats_filtered = sumstats %>% 
  filter(snp_pos %in% cres_no_egene$snp_pos)

count.df = fread("eQTL_catalogue/count_means/INTERVAL_salmon_tximport_gene_gene_mean_counts.csv.gz")
colnames(count.df) = c("ensembl_id","count.mean")

gene_list <- read_excel("Post_coloc_analysis/data/interval_genes.xlsx",sheet = 1) %>%
  filter(gene_biotype == "protein_coding")

count.df = count.df[ensembl_id %in% gene_list$feature_id]

# Add ensemble ids
# Load gencode with GRC37 positions
annot_file = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/Gencode/gencode.v33lift37.GRCh38.genes.gtf"
# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$ensembl_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[4]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot.GRC37 <- annot[order(annot$chr_number, annot$start),c("chr_number","start","ensembl_id", "gene_name")]

count.df = count.df %>% 
  left_join(annot.GRC37, "ensembl_id")

power.per.cre = function(count.df, chr, pos, snp, grna_target, beta, n){
  
  print(snp)
  
  # Calculate heritability/variance explained by MAF and effect size
  maf = sumstats_filtered$minor_AF[sumstats_filtered$snp_pos == snp]
  
  if(length(maf) == 0){
    maf = mean(sumstats_filtered$minor_AF)
  }
  
  q2 = 2*maf*(1-maf)*(beta^2)
  
  # Obtain genes within 1mb of snp add extra as snp is in GRC37
  cis.count.df = count.df[chr == chr_number & inrange(start, pos-1000000, pos+1000000),]
  
  
  if(nrow(cis.count.df) < 1){
    stop("No cis genes in the region")
  }
  
  #Skip power calculation for not expressed genes (<0.01)
  power_results <- lapply(1:nrow(cis.count.df), function(x) {
    if (cis.count.df$count.mean[x] < 0.01) {
      return(0) # Set power to 0 if count.mean < 0.01
    } else {
      return(scPower:::power.eqtl(
        count.mean = cis.count.df$count.mean[x], # Expression mean in the pseudobulk
        heritability = q2,                      # Heritability
        sig.level = 1e-05,                      # Significance threshold
        nSamples = n) # Sample size
      )
    }
  })
  
  power_results_df = data.frame(snp = snp,
                                genes = cis.count.df$gene_name,
                                ensembl_id = cis.count.df$ensembl_id,
                                count.mean = cis.count.df$count.mean,
                                eQTL.power = unlist(power_results),
                                grna_target = grna_target,
                                MAF = maf)
  
  return(power_results_df)
  
}

power.per.cre = purrr::possibly(power.per.cre, otherwise = NA, quiet = F)

results = lapply(1:nrow(cres_no_egene), function(x) power.per.cre(count.df = count.df,
                                                                 chr = cres_no_egene$chr[x],
                                                                 pos = as.numeric(cres_no_egene$pos[x]),
                                                                 snp = cres_no_egene$snp_pos[x],
                                                                 grna_target = cres_no_egene$grna_target[x],
                                                                 beta = 0.05,
                                                                 n=100000))

power_results_beta0.05 = do.call(rbind, results)

mean(power_results_beta0.05$eQTL.power)
sum(power_results_beta0.05$eQTL.power > 0.8)/nrow(power_results_beta0.05)

results = lapply(1:nrow(cres_no_egene), function(x) power.per.cre(count.df = count.df,
                                                                  chr = cres_no_egene$chr[x],
                                                                  pos = as.numeric(cres_no_egene$pos[x]),
                                                                  snp = cres_no_egene$snp_pos[x],
                                                                  grna_target = cres_no_egene$grna_target[x],
                                                                  beta = 0.1,
                                                                  n=25000))

power_results_beta0.1 = do.call(rbind, results)

mean(power_results_beta0.1$eQTL.power)
sum(power_results_beta0.1$eQTL.power > 0.8)/nrow(power_results_beta0.1)

results = lapply(1:nrow(cres_no_egene), function(x) power.per.cre(count.df = count.df,
                                                                  chr = cres_no_egene$chr[x],
                                                                  pos = as.numeric(cres_no_egene$pos[x]),
                                                                  snp = cres_no_egene$snp_pos[x],
                                                                  grna_target = cres_no_egene$grna_target[x],
                                                                  beta = 0.25,
                                                                  n=4732))

power_results_beta0.25 = do.call(rbind, results)

Ns     <- c(5000, 10000, 25000, 50000, 100000)
Betas  <- c(0.05, 0.10, 0.25)

# ----------------------------
# Helper: run power for a given (n, beta)
# ----------------------------
run_power_for <- function(n, beta) {
  results <- lapply(seq_len(nrow(cres_no_egene)), function(x)
    power.per.cre(
      count.df    = count.df,
      chr         = cres_no_egene$chr[x],
      pos         = as.numeric(cres_no_egene$pos[x]),
      snp         = cres_no_egene$snp_pos[x],
      grna_target = cres_no_egene$grna_target[x],
      beta        = beta,
      n           = n
    )
  )
  df <- do.call(rbind, results)
  tibble(
    sample_size = n,
    beta        = beta,
    mean_power  = mean(df$eQTL.power, na.rm = TRUE),
    prop_ge_80  = mean(df$eQTL.power >= 0.8, na.rm = TRUE)
  )
}

# ----------------------------
# Run all (beta, N) combos
# ----------------------------
summary_all <- map_dfr(Betas, function(b)
  map_dfr(Ns, function(n) run_power_for(n, b))
)

# Optional: inspect
# print(summary_all)

# Nice labels for legend
summary_all <- summary_all %>%
  mutate(beta_f = factor(beta, levels = Betas,
                         labels = paste0("\u03B2 = ", format(Betas, trim = TRUE))))

fwrite(summary_all, "Post_coloc_analysis/processed_data/power_vs_N_all_betas.txt", quote = F, row.names = F)

summary_all = fread("Post_coloc_analysis/processed_data/power_vs_N_all_betas.txt")

Ns     <- c(5000, 10000, 25000, 50000, 100000)
Betas  <- c(0.05, 0.10, 0.25)

# ----------------------------
# Plot: Proportion of loci with power ≥ 0.8 vs sample size
# ----------------------------
png("Post_coloc_analysis/plots/interval/power_vs_N_all_betas.png",
    width = 10, height = 10, units = "in", res = 300)

ggplot(summary_all, aes(x = sample_size, y = prop_ge_80, color = beta_f, group = beta_f)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(trans = "log10", breaks = Ns) +
  labs(
    x = "Sample size (log scale)",
    y = "Proportion of loci with power \u2265 0.8",
    color = "Effect size",
    title = ""
  ) +
  theme_cowplot(font_size = 24) +
  theme(legend.position = "top")

dev.off()

# ----------------------------
# (Optional) Also save the summary table
# ----------------------------
readr::write_csv(summary_all, "Post_coloc_analysis/plots/interval/power_summary_all_betas.csv")
