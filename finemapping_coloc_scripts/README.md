# Fine-mapping and Colocalization Pipeline

A comprehensive pipeline for statistical fine-mapping of GWAS signals and colocalization analysis with eQTL datasets using SuSiE and coloc methods.

**Author:** Sam Ghatan
**Last Updated:** June 12, 2024

## Overview

This pipeline performs:
1. Statistical fine-mapping of GWAS summary statistics using SuSiE
2. Colocalization analysis between GWAS and eQTL signals
3. Generation of locus plots and credible set visualizations
4. Integration with multiple eQTL datasets (GTEx, OneK1K, MAGE, Interval)

## Table of Contents

- [Requirements](#requirements)
- [Pipeline Workflow](#pipeline-workflow)
- [Input Data](#input-data)
- [Usage](#usage)
- [Output](#output)
- [Scripts Description](#scripts-description)

## Requirements

### R Packages
```r
coloc
susieR
data.table
tidyverse
foreign
purrr
Rfast
geni.plots
bigsnpr
arrow
biomaRt
```

### System Requirements
- R version 4.4.1 or higher
- SLURM job scheduler (for batch submissions)
- Minimum 20-50GB RAM depending on analysis

## Pipeline Workflow

### Step 1: Fine-mapping with SuSiE

Perform statistical fine-mapping on GWAS summary statistics to identify credible sets of causal variants.

**Scripts:**
- `01.susie_finemap.R` - Core fine-mapping script
- `01.susie_finemap_full_stats.R` - Full statistics version
- `01.susie_finemap_Interval.R` - Interval-specific version
- `01.submit_susieR.sh` - SLURM submission script
- `01.submit_susieR_full_stat.sh` - Full stats submission script

**Input:**
- GWAS summary statistics (formatted)
- LD matrices (per chromosome, per region)
- Region definitions file

**Output:**
- Fine-mapped variants with posterior inclusion probabilities (PIPs)
- Credible sets (95% coverage by default)
- Locus plots
- BED files of fine-mapped variants

### Step 2: Colocalization with GTEx

Test for colocalization between GWAS signals and GTEx eQTLs.

**Scripts:**
- `02.susie_coloc_GTEx.R` - GTEx colocalization
- `02.submit_susie_coloc_gtex.sh` - SLURM submission

**Input:**
- GWAS summary statistics
- GTEx eQTL summary statistics
- LD matrices for both datasets

**Output:**
- Colocalization results with PP.H4 probabilities
- Stacked locus plots (GWAS + eQTL)

### Step 3: Colocalization with OneK1K

Test for colocalization with cell-type-specific eQTLs from OneK1K dataset.

**Scripts:**
- `03.susie_coloc_OneK1K.R` - OneK1K colocalization
- `03.submit_susie_coloc.sh` - SLURM submission
- `03.parallel_coloc.sh` - Parallel processing script

**Cell Types Analyzed:**
- NK cells
- B cells
- CD4+ T cells
- CD8+ T cells
- Dendritic cells
- Monocytes
- Other cell types

### Step 4: Colocalization with BCX/MAGE

Colocalization with Blood Cell Consortium (BCX) and MAGE datasets.

**Scripts:**
- `04.BCX_MAGE_coloc.R` - BCX/MAGE colocalization
- `04.BCX_eQTL_catalogue.R` - eQTL Catalogue integration
- `04.submit_BCX_coloc.sh` - SLURM submission

### Step 5: Merge Results

Combine colocalization results across chromosomes and datasets.

**Scripts:**
- `05.merge_coloc_results.sh` - Merge results across chromosomes and cell types

### Step 6: Interval Colocalization

**Scripts:**
- `06.Interval_coloc.R` - Interval dataset colocalization
- `06.submit_interval_coloc.sh` - SLURM submission
- `Interval_finemap.R` - Interval-specific fine-mapping

## Input Data

### GWAS Summary Statistics Format

Required columns:
```
variant          # Format: chr:pos:ref:alt
Chr              # Chromosome
Pos              # Position (hg19)
minor_allele     # Minor allele
minor_AF         # Minor allele frequency
n_complete_samples  # Sample size
beta             # Effect size
se               # Standard error
pval             # P-value
```

### eQTL Summary Statistics Format

Required columns (TensorQTL format):
```
phenotype_id     # Gene ID
variant_id       # Variant ID
af               # Allele frequency
ma_samples       # Minor allele samples
ma_count         # Minor allele count
pval_nominal     # Nominal p-value
beta             # Effect size
slope_se         # Standard error
```

### LD Matrices

- Pre-computed LD matrices in PLINK binary format (.bed/.bim/.fam)
- Organized by chromosome and genomic region
- Directory structure: `data/[DATASET]_LDmatrices/chr[CHR]/[LOWER].[UPPER]/`

### Region Definitions

Tab-separated file with columns:
```
chr    lower    upper
```

## Usage

### 1. Fine-mapping

```bash
# Submit fine-mapping jobs for all chromosomes
sbatch 01.submit_susieR.sh
```

### 2. Colocalization

```bash
# GTEx colocalization
sbatch 02.submit_susie_coloc_gtex.sh

# OneK1K colocalization
sbatch 03.submit_susie_coloc.sh

# BCX/MAGE colocalization
sbatch 04.submit_BCX_coloc.sh

# Interval colocalization
sbatch 06.submit_interval_coloc.sh
```

### 3. Merge Results

```bash
sbatch 05.merge_coloc_results.sh
```

### Running Individual Scripts

```bash
# Fine-mapping example
Rscript 01.susie_finemap.R \
    data/UKBB_sumstats/30000_formatted.tsv \
    results/gwas_regions/merged_blood_trait_regions.txt \
    data/UKBB_LDmatrices/ \
    22  # chromosome

# Colocalization example
Rscript 02.susie_coloc_GTEx.R \
    GTEx \
    22 \
    data/UKBB_sumstats/30000_formatted.tsv \
    data/eQTL_sumstats/GTEx_whlbld.sumstats.chr22.hg19.txt
```

## Output

### Fine-mapping Results

```
results/UKBB_SuSiE_finemap/[GWAS_ID]/
├── [GWAS_ID]_chr[CHR]_finemap_results.txt    # Fine-mapped variants
├── [GWAS_ID]_chr[CHR]_finemap_results.bed    # BED format
└── credible_set_plots/                        # Visualization
```

**Key output columns:**
- `pip` - Posterior inclusion probability
- `lbf_cs[1-10]` - Log Bayes factors for credible sets
- `cs` - Credible set assignment
- `cond_indep` - Conditionally independent signal ID

### Colocalization Results

```
results/Coloc_results_V2/[DATASET]/[GWAS_ID]/
├── chr[CHR]_coloc_results.txt                 # Per-chromosome results
└── Locus_plots/                               # Stacked locus plots
```

**Key output columns:**
- `PP.H0.abf` - PP for no association
- `PP.H1.abf` - PP for GWAS only
- `PP.H2.abf` - PP for eQTL only
- `PP.H3.abf` - PP for both traits, different causal variants
- `PP.H4.abf` - PP for colocalization (shared causal variant)
- `SNP.PP.H4` - Variant-level posterior probability
- `gwas_credible_set` - GWAS credible set ID
- `eqtl_credible_set` - eQTL credible set ID
- `method` - Colocalization method (susie_coloc or coloc)

### Interpretation

- **PP.H4 > 0.8**: Strong evidence for colocalization
- **PP.H4 > 0.5**: Moderate evidence for colocalization
- **PP.H3 > 0.8**: Both traits associated but likely different causal variants

## Scripts Description

### Main Analysis Scripts

| Script | Description |
|--------|-------------|
| `01.susie_finemap.R` | Fine-maps GWAS signals using SuSiE-RSS |
| `01.susie_finemap_full_stats.R` | Fine-mapping with full summary statistics |
| `01.susie_finemap_Interval.R` | Interval-specific fine-mapping |
| `02.susie_coloc_GTEx.R` | Colocalization with GTEx whole blood eQTLs |
| `03.susie_coloc_OneK1K.R` | Colocalization with OneK1K cell-type eQTLs |
| `04.BCX_MAGE_coloc.R` | Colocalization with BCX/MAGE datasets |
| `04.BCX_eQTL_catalogue.R` | Integration with eQTL Catalogue |
| `06.Interval_coloc.R` | Colocalization with Interval study |
| `Interval_finemap.R` | Interval dataset fine-mapping helper |

### Submission Scripts

| Script | Description |
|--------|-------------|
| `01.submit_susieR.sh` | Submit fine-mapping jobs |
| `01.submit_susieR_full_stat.sh` | Submit full stats fine-mapping |
| `02.submit_susie_coloc_gtex.sh` | Submit GTEx colocalization |
| `03.submit_susie_coloc.sh` | Submit OneK1K colocalization |
| `03.parallel_coloc.sh` | Parallel colocalization processing |
| `04.submit_BCX_coloc.sh` | Submit BCX colocalization |
| `06.submit_interval_coloc.sh` | Submit Interval colocalization |
| `05.merge_coloc_results.sh` | Merge results across datasets |

### Grid Computing Scripts

| Script | Description |
|--------|-------------|
| `launch_coloc_grid.sh` | Launch colocalization on grid |
| `run_coloc_grid.sh` | Execute grid jobs |
| `launch_coloc_grid_test.sh` | Test grid launch |
| `run_coloc_grid_test.sh` | Test grid execution |

### Legacy Scripts

The `old_scripts/` directory contains earlier versions and development scripts:
- `susie_coloc.R` / `susie_coloc_beta.R` / `susie_coloc_beta.2.0.R` - Previous colocalization versions
- `Calculate_LD.R` / `Check_LDblocks.R` - LD computation utilities
- `Find_LD_proxies.R` - Proxy SNP identification
- `locus_plot.R` / `stacked_locus_plots.R` - Visualization scripts
- `finemapping_examples.R` - Tutorial examples

## Methods

### Statistical Fine-mapping (SuSiE)

This pipeline uses **Sum of Single Effects (SuSiE)** regression for fine-mapping:
- Identifies multiple causal variants per locus
- Provides credible sets with specified coverage (default 95%)
- Returns posterior inclusion probabilities (PIPs) for each variant
- Uses summary statistics and LD matrices (SuSiE-RSS)

**Key parameters:**
- `L = 10`: Maximum number of causal variants per region
- `prior_variance = 50`: Prior effect size variance
- `estimate_residual_variance = FALSE/TRUE`: Depends on trait type
- `coverage = 0.95`: Credible set coverage

### Colocalization (coloc/SuSiE-coloc)

Two methods are used:
1. **SuSiE-coloc**: Performs fine-mapping and colocalization jointly, identifying which credible sets colocalize
2. **coloc.abf**: Approximate Bayes Factor approach when SuSiE fails to converge

Both methods test five hypotheses:
- H0: No association with either trait
- H1: Association with trait 1 only
- H2: Association with trait 2 only
- H3: Association with both traits, different causal variants
- H4: Association with both traits, shared causal variant

## Configuration

Key paths are configured within each script. Update these for your environment:

```r
wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap"
gwas_ld_path = "data/UKBB_LDmatrices/"
eqtl_ld_path = "data/OneK1K_LDmatrices/"
output_dir = "results/Coloc_results_V2/"
```

## Troubleshooting

### SuSiE Convergence Issues
If SuSiE doesn't converge:
- Script automatically reduces `L` from 10 to 1
- If still failing, falls back to standard coloc.abf

### Memory Requirements
- Fine-mapping: 50GB per chromosome
- Colocalization: 20GB per job
- Adjust `#SBATCH --mem` in submission scripts if needed

### LD Matrix Issues
Ensure LD matrices are:
- Properly formatted (PLINK binary)
- Matching the summary statistics genome build (hg19)
- Complete (no missing regions)

## Citation

If you use this pipeline, please cite:

- **SuSiE:** Wang et al. (2020) JRSS-B. doi:10.1111/rssb.12388
- **coloc:** Giambartolomei et al. (2014) PLoS Genet. doi:10.1371/journal.pgen.1004383
- **SuSiE-coloc:** Zou et al. (2022) bioRxiv. doi:10.1101/2022.10.14.512286

## Contact

For questions or issues, please contact the repository maintainer or open an issue on GitHub.

## License

Please consult with the lab before using or distributing this code.
