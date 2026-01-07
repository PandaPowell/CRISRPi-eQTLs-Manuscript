# Post-Colocalization Analysis Scripts

This directory contains scripts for analyzing the overlap and comparison between CRISPRi target genes (cGenes) and eQTL target genes (eGenes) following colocalization analysis.

## R Scripts

### Gene Expression and Characterization

- **Gene_expr_all_cell_types_beta.R** - Compares gene expression levels between cGenes and eGenes across multiple cell types using OneK1K and Gasperini K562 single-cell RNA-seq data. Generates boxplots and density plots showing expression differences.

- **Heatmap_by_study_type.R** - Creates heatmaps displaying gene expression patterns across different eQTL studies for genes identified as cGenes versus eGenes.

- **Context_specific_exp.R** - Analyzes context-specific gene expression differences, determining what proportion of cGenes/eGenes are expressed in specific cell types or contexts.

### Distance and Effect Size Analysis

- **distance_vs_effect_size.R** - Examines the relationship between TSS distance and CRISPR effect size for both STING-seq and Gasperini datasets, testing whether distance influences regulatory effect magnitude.

- **distance_vs_logfc.R** - Performs regression analysis plotting log fold change versus TSS distance, adjusting for number of cells, to model the distance-effect size relationship.

### Chromatin Architecture

- **Compare_HiC_targets.R** - Compares CRE-gene pairs supported by Hi-C chromatin interaction data between cGenes and eGenes using K562 Hi-C data. Creates barplots showing proportions of genes with Hi-C support.

- **intersection_of_ATAC_peaks.R** - Analyzes ATAC-seq peak overlaps with CREs across different cell types. Compares chromatin accessibility patterns between CREs targeting cGenes versus eGenes.

### Gene Constraint and Function

- **cgene_egene_weissman.R** - Compares the proportion of cGenes and eGenes with high haploinsufficiency scores (pLI > 0.9) using Weissman gene constraint data.

### Gold Standard Gene Analysis

- **number_gold_gene_closest.R** - Determines how many gold standard causal genes from Mendelian disease and UK Biobank burden tests are identified as closest genes to CREs by distance rank.

- **gold_genes_and_power.R** - Comprehensive analysis of gold standard gene enrichment in cGenes vs eGenes. Calculates enrichment compared to random protein-coding genes, analyzes statistical power limitations, and creates Venn diagrams showing overlap between methods.

### Visualization

- **Locus_heatmap_plot.R** - Generates detailed locus-specific heatmaps for gold standard genes, showing CRISPRi results and eQTL associations across multiple studies and cell types. Highlights target genes and closest genes.

- **split_heatmaps.R** - Creates split heatmaps for gold genes with separate columns for CRISPRi and eQTL data when both methods identify the gene, allowing method comparison.

- **gwas_heatmaps.r** - Generates heatmaps showing which GWAS traits are associated with specific gold genes through CRISPRi and/or eQTL evidence.

### Statistical Power

- **power_by_coloc.R** - Analyzes statistical power of eQTL studies as a function of sample size and effect size. Compares power across different eQTL datasets and demonstrates the relationship between sample size and number of colocalizations detected.

- **Precision_Recall.R** - Calculates precision and recall metrics for different gene prediction methods (cGenes, eGenes, Hi-C, ABC-Max, closest gene, TWAS) using gold standard genes as the truth set. Generates precision-recall plots.

### Comprehensive Analysis

- **figure_plots.R** - Main comprehensive script generating multiple analysis figures including:
  - TSS distance distributions
  - Number of genes per CRE
  - Gene characteristics (expression, pLI, number of enhancers)
  - Closest gene rankings
  - MAF distributions
  - PP.H4 colocalization probabilities
  - Power comparisons
  - Upset plots showing eQTL dataset overlaps

- **GWAS_specific_mappings.R** - Maps GWAS trait-specific CRE-gene relationships, analyzing pleiotropy (number of traits per CRE), creating trait-specific barplots, and examining concordance between CRISPRi and eQTL trait associations using rank correlation.

## Shell Scripts

- **Run_Extract_MetaLCL.sh** - SLURM batch job script that converts a Jupyter notebook to Python and executes Extract_MetaLCL.py to extract trans-eQTL data from MetaLCL datasets. Requires 16 CPUs and 40GB memory.

- **intersection_of_ATAC_peaks.sh** - Uses bedtools to:
  1. Create BED file of unique CREs from processed data
  2. Intersect CREs with K562 ATAC-seq peaks
  3. Intersect CREs with ATAC-seq peaks from multiple cell types
  4. Calculate overlap statistics

## Python Scripts

- **Extract_MetaLCL.py** - Filters MetaLCL eQTL parquet files to extract trans-eQTL associations for cis-genes identified in Gasperini and Morris studies. Uses DuckDB for efficient parquet querying and parallel processing with joblib to handle large datasets.

## Key Concepts

- **cGenes**: Genes identified as CRISPRi targets (genes whose expression changes when nearby CREs are perturbed)
- **eGenes**: Genes identified through eQTL analysis (genes with expression associated with genetic variants)
- **Gold Standard Genes**: Curated set of high-confidence causal genes from Mendelian disease studies and UK Biobank burden tests
- **CREs**: Cis-regulatory elements (enhancers/regulatory regions)
- **Colocalization**: Statistical method to identify shared causal variants between GWAS and eQTL signals

## Dependencies

### R Packages
- data.table
- tidyverse (dplyr, tidyr, ggplot2, stringr)
- cowplot
- gridExtra
- ggrepel
- viridis
- Seurat
- scPower
- biomaRt
- readxl
- svglite
- eulerr
- UpSetR
- patchwork
- scales

### Python Packages
- duckdb
- pandas
- joblib

### System Tools
- bedtools (v2.31.0)
- SLURM (for job submission)

## Usage Notes

Most R scripts are designed to be run from the Post_coloc_analysis directory and expect specific input files in subdirectories:
- `processed_data/` - Processed CRE-gene mappings
- `CRISPR_data/` - CRISPRi experimental data
- `eQTL_catalogue/` - eQTL summary statistics
- `power_results/` - Statistical power calculations
- `plots/` - Output directory for figures

Scripts typically output plots to `plots/` or `plots/interval/` subdirectories.

## Analysis Workflow

1. Start with processed colocalization results in `processed_data/`
2. Run gene characterization scripts to compare properties
3. Calculate statistical power using `power_by_coloc.R`
4. Generate gold standard enrichment analyses
5. Create visualization plots for publication
6. Run GWAS-specific analyses for trait pleiotropy

## Contact

For questions about specific analyses or scripts, refer to the methods sections of the associated manuscript or contact the repository maintainer.
