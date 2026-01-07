# Post-Colocalization Analysis Scripts

This directory contains scripts for analyzing the overlap between CRISPR-based functional genomics data (STING-seq, Gasperini, ENCODE) and eQTL studies, specifically focusing on colocalization with GWAS signals.

## Overview

The analysis pipeline identifies cis-regulatory elements (CREs) that intersect with fine-mapped GWAS variants and compares CRISPR-validated target genes (cGenes) with eQTL-identified genes (eGenes).

## Pipeline Workflow

### 1. CRISPR Data Processing

#### `01.determine_encode_targetseq.sh`
**Purpose**: Process ENCODE CRISPRi data and map gRNA spacer sequences to genome coordinates
- Downloads ENCODE harmonized CRISPRi data (GRCh38)
- Creates BLAST database from reference genome (GRCh38)
- Extracts gRNA spacer sequences and aligns them using BLAST
- Lifts over coordinates from GRCh38 to GRCh37
- Filters out Gasperini data to avoid duplication
- **SBATCH Resources**: 20G memory, 8 CPUs

#### `01.determine_gasperini_targetseq.sh`
**Purpose**: Process Gasperini et al. CRISPRi-seq data and map gRNA target sites
- Downloads and formats Gasperini supplementary data
- Creates fasta files from spacer sequences
- Aligns gRNA spacers to GRCh37.p13 using BLAST
- **SBATCH Resources**: 20G memory, 8 CPUs

### 2. Intersection Analysis

#### `02.run_intersections.sh`
**Purpose**: Master script that orchestrates all intersection analyses
- Executes multiple intersection scripts sequentially
- Intersects CRISPR data (ENCODE, Gasperini, STING-seq) with fine-mapped GWAS variants
- Intersects eQTL data (OneK1K, GTEx, eQTL Catalogue, MAGE) with CRISPR target sites
- Stops pipeline if any script fails (`set -e`)

**Sub-scripts executed:**
1. `01.Filter_finemap_pip.sh` - Filters fine-mapped variants by PIP threshold
2. `01.Intersect_ENCODE_w_finemap.sh` - ENCODE × GWAS fine-mapping
3. `01.Intersect_BCX.sh` - BCX dataset intersections
4. `01.Intersect_sting_seq_w_finemap.sh` - STING-seq × GWAS fine-mapping
5. `02.Intersect_gasperini_w_finemap_closest.sh` - Gasperini × GWAS fine-mapping
6. `03.Intersect_onek1k_gtex_with_stingseq.sh` - OneK1K/GTEx × STING-seq
7. `04.Intersect_onek1k_gtex_with_gasperini.sh` - OneK1K/GTEx × Gasperini
8. `05.Intersect_onek1k_gtex_with_encode.sh` - OneK1K/GTEx × ENCODE
9. `06.Intersect_eqtl_cat_with_all_crispr.sh` - eQTL Catalogue × all CRISPR
10. `07.Intersect_MAGE_with_all_crispr.sh` - MAGE × all CRISPR

### 3. Gene Distance Calculations

#### `03.Run_gene_distances.sh`
**Purpose**: Calculate gene distance ranks for eQTL datasets
- Extracts gene distance information for eQTL Catalogue
- Processes GTEx gene distances
- Processes OneK1K gene distances
- Processes MAGE gene distances
- **SBATCH Resources**: 100G memory, 8 CPUs

### 4. Data Formatting and Integration

#### `04.format_merge_data_interval.R`
**Purpose**: Main data integration and formatting script (549 lines)
- Loads and formats CRISPR data from three sources:
  - STING-seq CREs
  - Gasperini targeted enhancers
  - ENCODE CRISPRi (excluding Gasperini)
- Loads Gencode annotations (GRCh37 and GRCh38 versions)
- Merges CRISPR target sites with fine-mapped GWAS variants
- Integrates eQTL colocalization results from:
  - Interval (scRNA-seq eQTLs)
  - eQTL Catalogue
  - GTEx
  - MAGE
  - OneK1K
- Merges overlapping CREs (within 4kb)
- Calculates distance ranks for CRE-gene pairs
- Submits gene distance calculation job and waits for completion
- **Outputs**:
  - `cres_with_grnas.txt` - All CRE-gene pairs from CRISPR
  - `cres_with_grna_eqtls_interval.txt` - CRE-eGene pairs from eQTL colocalization
  - `Gasperini_gRNAs_intersecting_GWAS.txt` - Gasperini targets intersecting GWAS

### 5. Statistical Power Analysis

#### `05.Run_sceQTL_power.high.mem.sh`
**Purpose**: Calculate statistical power for single-cell eQTL detection (CD4/CD8 T cells)
- Array job for CD4_T and CD8_T cell types
- Calls R script `08.calculate_sceqtl_power.R`
- **SBATCH Resources**: 150G memory, 8 CPUs per task

#### `05.Run_sceQTL_power.low.mem.sh`
**Purpose**: Calculate statistical power for single-cell eQTL detection (other cell types)
- Array job for B, Mono, NK, DC cell types
- Calls R script `08.calculate_sceqtl_power.R`
- **SBATCH Resources**: 50G memory, 8 CPUs per task

#### `05.Run_bulk_eQTL_power.sh`
**Purpose**: Calculate statistical power for bulk eQTL detection
- Runs power calculations for bulk eQTL studies
- Processes Interval cohort power analysis
- **SBATCH Resources**: 80G memory, 1 CPU, 7-day time limit
- **Scripts executed**:
  - `08.calculate_bulk_eqtl_power.R`
  - `08.calculate_interval_power.R`

### 6. Power Results Processing

#### `06.process_format_power.R`
**Purpose**: Format and integrate statistical power results with association data (309 lines)
- Loads CRISPR power calculations for STING-seq, Gasperini, and ENCODE
- Merges power estimates with CRE-gene association results
- Handles overlapping CREs (within 4kb window)
- Calculates proportion of CRE-gene pairs with power ≥80% at different effect sizes
- Generates visualization of power by effect size
- **Outputs**:
  - `cres_with_grnas_power.txt` - CRE-gene pairs with power estimates
  - `plots/interval/proportion_by_effect_size.svg` - Power visualization

### 7. Comparison and Analysis

#### `07.compare_targets.R`
**Purpose**: Compare CRISPR target genes (cGenes) with eQTL genes (eGenes)
- Identifies overlapping CREs between CRISPR and eQTL datasets
- Calculates overlap statistics for genes and CRE-gene pairs
- Generates distance-based comparisons
- Creates visualizations comparing significant vs. non-significant targets

### 8. Trans-eQTL Network Analysis

#### `08.CompareTransNet_transeQTL.R`
**Purpose**: Compare trans networks from CRISPR screens with trans-eQTLs (405 lines)
- Loads Gasperini trans network data and OneK1K cis-eQTL results
- Identifies cis-variants for genes with trans-effects (HHEX, GFI1B, IKZF1, RUNX1, RPL23A)
- Loads trans-eQTL results from OneK1K across 8 cell types
- Merges trans-network data with trans-eQTL data by gene
- Creates correlation plots comparing CRISPR trans-gene effects vs. eQTL trans-gene effects
- Generates violin plots showing genotype-expression relationships for specific trans-eQTL associations
- Uses Morris et al. supplementary data (science.adh7699_table_s4.xlsx)

#### `09.IdentifyGeneswTransNetw.R`
**Purpose**: Identify cis-genes with significant trans networks (35 lines)
- Loads Gasperini trans network data
- Filters trans-effects removing genes within ±5Mbp of gRNAs (true trans, not proximal cis)
- Identifies gRNAs with >5 significant trans-effects (q<0.1)
- Links gRNA targets to their cis-target genes using Gasperini data
- **Output**: `genes_cis_w_trans.txt` - List of cis-genes with trans networks

#### `10.Plot_trans_Res_allPairs.R`
**Purpose**: Analyze and visualize trans-eQTL results across cell types (741 lines)
- Loads trans-eQTL results (all pairs) from OneK1K for 8 cell types
- Loads cis-eQTL results and merges with trans results by variant
- Creates UpSet plots showing sharing of cis-genes with trans-effects across cell types
- Generates density plots of number of trans-genes per cis-QTL
- Compares trans-eQTL results with Morris and Gasperini CRISPR trans networks
- Creates correlation plots (trans-eGene beta vs. trans-cGene LogFC) for specific genes:
  - HHEX, GFI1B, IKZF1, RUNX1, RPL23A
- Filters to ChIP-seq validated targets (for GFI1B and IKZF1)
- Generates pairs plots comparing trans-effects across cell types for the same variant
- **Outputs**: Multiple PDF plots in `analysis/02_QTL_calling/plots/`

#### `11.Compare_trans_results.R`
**Purpose**: Comprehensive comparison of trans-eQTLs and trans networks (744 lines)
- Loads trans-eQTL data from multiple sources:
  - MetaLCL bulk eQTLs
  - OneK1K single-cell eQTLs (8 cell types)
- Loads CRISPR trans network data:
  - Morris et al. (5 cis-genes: GFI1B, HHEX, IKZF1, NFE2, RUNX1)
  - Gasperini dataset (filtered to genes with trans-effects)
- Filters trans associations to >5Mbp from cis-gene TSS
- Calculates correlations between trans-cGene effects and trans-eGene effects
- Performs enrichment analysis (Fisher's test with resampling)
- Generates directionality tests (binomial test for concordant effect directions)
- Creates visualizations:
  - Barplots of number of trans-genes per cis-gene
  - Scatter plots with correlation statistics
  - Forest plots for correlation estimates and enrichment ORs
  - UpSet plots showing cell-type sharing
- Flips effect directions to harmonize with decreasing allele
- **Outputs**:
  - `Supplementary_Table9_trans_eGenes.txt` - All trans-eGene associations
  - `trans_comparison_results.txt` - Correlation and enrichment statistics
  - Multiple plots in `plots/supp_figs/` and `plots/figure_plots/`

#### `12.Compare_trans_results_enrichment.R`
**Purpose**: Test enrichment of trans-genes in gold standard disease genes (206 lines)
- Loads genome-wide gold standard genes (`gw_gold_genes.txt`)
- Queries BioMart for all protein-coding genes (background set)
- Defines enrichment function using resampling approach:
  - Simulates null distribution by random sampling
  - Calculates Fisher's exact test for observed vs. expected overlaps
- Compares trans-cGenes vs. trans-eGenes for enrichment in disease genes
- Analyzes both Morris and Gasperini trans networks
- Tests across multiple cis-genes and cell types
- Applies FDR correction for multiple testing
- **Output**: `supplementary_table_13.txt` - Enrichment statistics for trans-genes

## Key Data Sources

### CRISPR Functional Genomics:
- **STING-seq**: SNP-targeted CRISPRi screen
- **Gasperini et al.**: Enhancer targeting CRISPRi-seq
- **ENCODE**: Harmonized CRISPRi data

### eQTL Studies:
- **Interval**: Blood cell scRNA-seq eQTLs
- **OneK1K**: Single-cell immune eQTLs
- **GTEx**: Bulk tissue eQTLs
- **eQTL Catalogue**: Harmonized eQTL summary statistics
- **MAGE**: Multi-ancestry gene expression study

### GWAS Data:
- Fine-mapped GWAS variants with posterior inclusion probabilities (PIP)
- Multiple trait categories including autoimmune, complex diseases

## Output Files

### Main Outputs:
- `cres_with_grnas.txt` - All CRE-gene pairs from CRISPR experiments
- `cres_with_grna_eqtls_interval.txt` - CRE-eGene pairs with colocalization evidence
- `cres_with_grnas_power.txt` - CRE-gene pairs with statistical power estimates
- `finemap_snp_intersect_grna.txt` - Fine-mapped SNPs intersecting gRNA targets

### Intermediate Files:
- Gene distance rankings for each eQTL dataset
- BLAST alignment results
- Lifted-over genomic coordinates

## Software Requirements

### Modules:
- R (v4.3.1, v4.4.1)
- BLAST (v2.10.0)
- LiftOver

### R Packages:
- tidyverse
- data.table
- cowplot
- gridExtra
- fst

## Key Parameters

- **CRE merging window**: 4001 bp
- **TSS distance filter**: ±1 Mb
- **GWAS p-value threshold**: < 1×10⁻⁵
- **eQTL nominal p-value**: < 1×10⁻³
- **STING-seq FDR threshold**: Q < 0.1
- **Statistical power threshold**: ≥80%

## Notes

- All genomic coordinates are processed in both GRCh37/hg19 and GRCh38
- BLAST alignment parameters optimized for short gRNA sequences (word_size=20)
- Some scripts depend on completion of upstream steps
- Array jobs enable parallel processing of cell types for efficiency
- The pipeline includes automatic job monitoring for dependent steps

## Authors

Lab: Lappalainen Lab
Analysis: sghatan
Project: STING-seq eQTL Overlap Analysis
