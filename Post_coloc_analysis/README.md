# Post-Colocalization Analysis Pipeline

Comprehensive pipeline for analyzing the overlap between CRISPR-based functional genomics experiments and eQTL studies, with a focus on GWAS colocalization.

## Overview

This pipeline integrates:
- **CRISPR functional genomics data**: STING-seq, Gasperini et al., ENCODE CRISPRi
- **eQTL studies**: Interval, OneK1K, GTEx, eQTL Catalogue, MAGE
- **GWAS fine-mapping**: Statistical fine-mapping with posterior inclusion probabilities
- **Validation datasets**: Hi-C, ABC predictions, ATAC-seq, ChIP-seq

The analysis identifies cis-regulatory elements (CREs) that colocalize with GWAS signals and compares CRISPR-validated target genes (cGenes) with eQTL-identified genes (eGenes).

## Repository Structure

```
Post_coloc_analysis/
├── scripts/              # Main analysis scripts (see detailed README)
│   ├── gene_distances/   # Gene distance calculations
│   ├── power_scripts/    # Statistical power analysis
│   ├── plotting_scripts/ # Visualization scripts
│   └── util/             # Utility functions
├── 00.intersect_data/    # Intersection analysis scripts
│   └── scripts/          # BEDTools intersections
├── CRISPR_data/          # CRISPR data processing
│   └── gasperini_2019/   # Gasperini dataset scripts
├── ABC/                  # ABC model predictions
├── Hi_C/                 # Hi-C interaction data
│   └── K562/             # K562 cell line Hi-C
├── encode_spaceseq/      # ENCODE spacer sequence mapping
└── gasperini_spaceseq/   # Gasperini spacer sequence mapping
```

## Documentation

Detailed documentation is available:
- **Sequential workflow**: `scripts/READ_sequential_scripts.md` - Step-by-step pipeline execution
- **Utility scripts**: `scripts/README_non_sequential_scripts.md` - Analysis and visualization scripts

## Quick Start

### 1. CRISPR Data Processing
```bash
# Map ENCODE gRNA sequences to genome
bash scripts/01.determine_encode_targetseq.sh

# Map Gasperini gRNA sequences to genome
bash scripts/01.determine_gasperini_targetseq.sh
```

### 2. Run Intersection Analysis
```bash
# Execute all intersections (CRISPR × GWAS × eQTL)
bash scripts/02.run_intersections.sh
```

### 3. Calculate Gene Distances
```bash
# Calculate distance ranks for eQTL datasets
bash scripts/03.Run_gene_distances.sh
```

### 4. Format and Merge Data
```bash
# Integrate all datasets
Rscript scripts/04.format_merge_data_interval.R
```

### 5. Calculate Statistical Power
```bash
# Single-cell eQTL power (high memory cell types)
sbatch scripts/05.Run_sceQTL_power.high.mem.sh

# Single-cell eQTL power (low memory cell types)
sbatch scripts/05.Run_sceQTL_power.low.mem.sh

# Bulk eQTL power
sbatch scripts/05.Run_bulk_eQTL_power.sh
```

### 6. Process Power Results
```bash
# Format and integrate power estimates
Rscript scripts/06.process_format_power.R
```

### 7. Compare Targets
```bash
# Compare cGenes with eGenes
Rscript scripts/07.compare_targets.R
```

### 8. Trans-eQTL Analysis
```bash
# Compare CRISPR trans networks with trans-eQTLs
Rscript scripts/08.CompareTransNet_transeQTL.R
Rscript scripts/09.IdentifyGeneswTransNetw.R
Rscript scripts/10.Plot_trans_Res_allPairs.R
Rscript scripts/11.Compare_trans_results.R
Rscript scripts/12.Compare_trans_results_enrichment.R
```

## Key Data Sources

### CRISPR Functional Genomics
- **STING-seq**: SNP-targeted CRISPRi screen in K562 cells
- **Gasperini et al. (2019)**: Enhancer-targeting CRISPRi-seq
- **ENCODE**: Harmonized CRISPRi data across multiple labs

### eQTL Studies
- **Interval**: Blood cell scRNA-seq eQTLs
- **OneK1K**: Single-cell immune cell eQTLs (8 cell types)
- **GTEx**: Bulk tissue eQTLs (whole blood)
- **eQTL Catalogue**: Harmonized eQTL summary statistics
- **MAGE**: Multi-ancestry gene expression study

### Validation Data
- **Hi-C**: Chromatin interaction data (K562)
- **ABC Model**: Activity-by-Contact enhancer predictions
- **ATAC-seq**: Chromatin accessibility
- **ChIP-seq**: Transcription factor binding sites

## Requirements

### Software
- R (≥4.3.1)
- Python (≥3.8)
- BLAST (≥2.10.0)
- BEDTools (≥2.29.0)
- LiftOver
- SLURM job scheduler

### R Packages
```r
tidyverse
data.table
cowplot
gridExtra
fst
biomaRt
UpSetR
ggrepel
ComplexHeatmap
```

### Python Packages
```python
pandas
numpy
pybedtools
```

## Key Parameters

- **CRE merging window**: 4001 bp
- **TSS distance filter**: ±1 Mb
- **GWAS p-value threshold**: < 1×10⁻⁵
- **eQTL nominal p-value**: < 1×10⁻³
- **STING-seq FDR threshold**: Q < 0.1
- **Statistical power threshold**: ≥80%
- **Trans-eQTL distance**: >5 Mbp from cis-gene TSS

## Main Outputs

### CRE-Gene Associations
- `cres_with_grnas.txt` - All CRISPR-validated CRE-gene pairs
- `cres_with_grna_eqtls_interval.txt` - CRE-eGene pairs with colocalization
- `cres_with_grnas_power.txt` - CRE-gene pairs with power estimates

### Intersection Results
- `finemap_snp_intersect_grna.txt` - GWAS variants × CRISPR targets
- Dataset-specific intersections in `00.intersect_data/results/`

### Trans-eQTL Analysis
- `genes_cis_w_trans.txt` - Cis-genes with trans networks
- `Supplementary_Table9_trans_eGenes.txt` - Trans-eQTL associations
- `trans_comparison_results.txt` - Correlation and enrichment statistics

## Genome Builds

The pipeline handles both genome builds:
- **GRCh37/hg19**: GWAS fine-mapping, GTEx, OneK1K
- **GRCh38**: ENCODE data, newer datasets
- Automatic liftover between builds where needed

## Citation

If you use this pipeline or data, please cite:

- **STING-seq**: [Your publication]
- **Gasperini et al. (2019)**: Nature Genetics. doi:10.1038/s41588-018-0315-8
- **OneK1K**: Nature. doi:10.1038/s41586-023-06229-1
- **GTEx**: Science. doi:10.1126/science.aaz1776

## Contact

For questions or issues, please contact:
- Lab: Lappalainen Lab
- Maintainer: sghatan

## License

Please consult with the lab before using or distributing this code.
