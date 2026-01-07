# collapse_transcripts_to_genes.R
# Run: Rscript collapse_transcripts_to_genes.R

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
})

# ---- Inputs (edit only if paths change) ----
rds_path  <- "/gpfs/commons/datasets/controlled/INTERVAL_sanger/rna_processed/INTERVAL_FilteredSamples_salmon_tximport.rds"
gtf       <- "/gpfs/commons/groups/lappalainen_lab/sghatan/Homo_sapiens.GRCh38.99.gtf.gz"
out_prefix <- "INTERVAL_salmon_tximport_gene"

stopifnot(file.exists(rds_path), file.exists(gtf))

# ---- Load transcript-level tximport object ----
txi_tx <- readRDS(rds_path)

# ---- Build tx2gene from Ensembl v99 GTF ----
# Lightweight parse via data.table + shell gunzip
gtf_cmd <- sprintf("gunzip -c %s | grep -v '^#'", shQuote(gtf))
gtf_dt <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE,
                col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))

tx_attr <- gtf_dt[feature == "transcript", .(attribute)]

get_attr <- function(x, key) {
  # pull e.g. gene_id / transcript_id from the attributes column
  # handles Ensembl-style attributes; keeps everything before quotes
  sub(sprintf('.*%s "([^"]+)".*', key), "\\1", x)
}

# strip version suffix (".1", ".2", …) to match Salmon’s tx IDs if needed
strip_ver <- function(ids) sub("\\..*$", "", ids)

tx2gene <- unique(data.frame(
  transcript = strip_ver(get_attr(tx_attr$attribute, "transcript_id")),
  gene       = strip_ver(get_attr(tx_attr$attribute, "gene_id")),
  stringsAsFactors = FALSE
))

# ---- Ensure rownames in txi match tx2gene$transcript (strip version if needed) ----
fix_rownames <- function(M) { if (!is.null(M)) { rownames(M) <- strip_ver(rownames(M)); M } else M }
txi_tx$counts    <- fix_rownames(txi_tx$counts)
txi_tx$abundance <- fix_rownames(txi_tx$abundance)
txi_tx$length    <- fix_rownames(txi_tx$length)

# ---- Summarize transcripts -> genes ----
txi_gene <- summarizeToGene(txi_tx, tx2gene = tx2gene)

mapped_frac <- mean(rownames(txi_tx$counts) %in% tx2gene$transcript)
cat(sprintf("Transcripts with gene mapping (Ensembl v99): %.2f%%\n", 100 * mapped_frac))
if (!is.null(txi_gene$countsFromAbundance)) {
  cat("countsFromAbundance: ", txi_gene$countsFromAbundance, "\n", sep = "")
}

# ---- Fast writers (gzipped CSVs) ----
setDTthreads(8)

fast_write_matrix <- function(mat, path_base) {
  df <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  setDT(df)
  out <- paste0(path_base, ".csv.gz")
  fwrite(df, file = out, sep = ",", quote = FALSE,
         nThread = getDTthreads(), compress = "gzip", showProgress = TRUE)
  out
}

counts_path <- fast_write_matrix(txi_gene$counts,  paste0("eQTL_catalogue/count_matrices/",out_prefix, "_counts"))

# ---- Per-gene and per-sample means ----
gene_mean_counts <- rowMeans(txi_gene$counts, na.rm = TRUE)
fwrite(data.table(gene_id = rownames(txi_gene$counts),
                  mean_count = as.numeric(gene_mean_counts)),
       file = paste0("eQTL_catalogue/count_means/",out_prefix, "_gene_mean_counts.csv.gz"),
       sep = ",", quote = FALSE, nThread = getDTthreads(), compress = "gzip", showProgress = TRUE)


# ---- Save an RDS for quick re-use ----
saveRDS(txi_gene, paste0(out_prefix, "_txi_gene.rds"))

cat("Done.\n",
    "Counts:    ", counts_path, "\n",
    "Means per gene:     ", paste0(out_prefix, "_gene_mean_counts.csv.gz"), "\n",
    "Means per sample:   ", paste0(out_prefix, "_sample_mean_counts.csv.gz"), "\n", sep = "")


