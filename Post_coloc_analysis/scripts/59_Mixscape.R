rm(list=ls())
# Load packages.
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)

# Download dataset using SeuratData.
InstallData(ds = "thp1.eccite")

# Setup custom theme for plotting.
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

# Load object.
eccite <- LoadData(ds = "thp1.eccite")

View(eccite@meta.data)

# Normalize protein.
eccite <- NormalizeData(
  object = eccite, 
  assay = "ADT", 
  normalization.method = "CLR", 
  margin = 2)

# Prepare RNA assay for dimensionality reduction: 
# Normalize data, find variable features and scale data.
DefaultAssay(object = eccite) <- 'RNA'
eccite <- NormalizeData(object = eccite) %>% FindVariableFeatures() %>% ScaleData()

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite)

# Run Uniform Manifold Approximation and Projection (UMAP) to visualize clustering in 2-D.
eccite <- RunUMAP(object = eccite, dims = 1:40)

# Calculate perturbation signature (PRTB).
eccite<- CalcPerturbSig(
  object = eccite, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  split.by = "replicate", 
  new.assay.name = "PRTB")

# Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.
DefaultAssay(object = eccite) <- 'PRTB'

# Use variable features from RNA assay.
VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[["RNA"]])
eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
eccite <- RunUMAP(
  object = eccite, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

library(mixtools)

# Run mixscape.
eccite <- RunMixscape(
  object = eccite, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "gene", # metadata column with target gene labels
  nt.class.name = "NT", 
  min.de.genes = 5, 
  iter.num = 10, 
  de.assay = "RNA", 
  verbose = F,
  prtb.type = "KO")

# Calculate percentage of KO cells for all target gene classes.
df <- prop.table(table(eccite$mixscape_class.global, eccite$NT),2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("NT", "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "g")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), 
                           function(x) strsplit(x, split = "g")[[1]][2])
df3 <- df2[-c(which(df2$gene == "NT")),]


### Mixscale ####
# Download and format the data
# Create directories and link raw data (same logic as before)
system(paste0("mkdir -p ./joust"))
# Download barcodes
system("wget -O ./joust/barcodes.tsv.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/GSE132080%5F10X%5Fbarcodes%2Etsv%2Egz")
#Download features
system("wget -O ./joust/features.tsv.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/GSE132080%5F10X%5Fgenes%2Etsv%2Egz")
#Download matrix
system("wget -O ./joust/matrix.mtx.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/GSE132080%5F10X%5Fmatrix%2Emtx%2Egz")
# Download cell identidies
system("wget -O ./joust/cell_identities.csv ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/GSE132080%5Fcell%5Fidentities%2Ecsv%2Egz")

ct_mat = ReadMtx(mtx = "./joust/matrix.mtx.gz", 
                 cells = "./joust/barcodes.tsv.gz", 
                 features = "./joust/features.tsv.gz")
# load the meta_data
meta_data = read.csv("joust//cell_identities.csv")
rownames(meta_data) = meta_data$cell_barcode

# create a seurat object 
seurat_obj = CreateSeuratObject(counts = ct_mat, meta.data = meta_data)
rm(ct_mat, meta_data)

# retrieve the guide information for each cell
txt = seurat_obj$guide_identity
txt2 = str_extract(txt, "^[^_]+")
txt3 = gsub(pattern = "^[^_]+_", replacement = "", txt)
seurat_obj[['gene']] = txt2
seurat_obj[['gRNA_name']] = txt3
seurat_obj[['cell_type']] = "K562"
rm(txt, txt2, txt3)

head(seurat_obj@meta.data)

system(paste0("mkdir -p ./epapalexi"))

# -------------------------------
# Add GDO (Guide-DNA Oligonucleotide) data as an assay
# -------------------------------
# Create directories and link raw data (same logic as before)
system(paste0("mkdir -p ./epapalexi_gdo"))
# Download barcodes
system("wget -O ./epapalexi_gdo/barcodes.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4633nnn/GSM4633618/suppl/GSM4633618%5FECCITE%5FGDO%5FBarcodes%2Ecsv%2Egz")
# Download features
system("wget -O ./epapalexi_gdo/features.tsv.gz ")
# Download matrix
system("wget -O ./epapalexi_gdo/matrix.mtx.gz ")

gdo_data <- Read10X(data.dir = "./morris_v1_gdo/")