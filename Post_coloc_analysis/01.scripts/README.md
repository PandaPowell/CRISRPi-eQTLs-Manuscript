Variant-to-Gene Mapping: A Comparative Analysis of CRISPRi and eQTL Approaches

This repository contains code and analyses from our manuscript "CRISPRi perturbation screens and eQTLs provide complementary and distinct insights into GWAS target genes". The project systematically compares expression quantitative trait locus (eQTL) studies and pooled CRISPRi screens for mapping target genes of non-coding GWAS variants, integrating data from 44 blood trait GWAS, 50 eQTL datasets, and 14 CRISPRi studies.
Repository Structure
01–04: GWAS–CRE Target Definition

    01.determine_gasperini_targetseq.sh: Parses CRISPRi gRNA target sequences from Gasperini et al.

    02.determine_encode_targetseq.sh: Parses gRNA targets from ENCODE CRISPRi screens.

    03.run_intersections.sh: Intersects GWAS fine-mapped SNPs with CRISPRi gRNA regions.

    04.calculate_gwas_gene_distance.sh: Calculates distances between GWAS variants and gene TSSs.

05: Gene Distance and eQTL Data Extraction

    05.extract_eQTL_catalogue_sumstats.R: Extracts and formats summary statistics from the eQTL Catalogue.

    05.gtex_gene_distances.R: Calculates gene distances for GTEx data.

    05.MAGE_gene_distances.R: Extracts MAGE LCL gene–variant distances.

    05.Onek1k_gene_distances.R: Calculates gene distances for OneK1K sc-eQTL dataset.

    05.gwas_gene_distance.R: Generates GWAS variant–gene distance tables.

    05.Run_gene_distances.sh: Wrapper script to run all gene distance calculations.

06–07: Bulk eQTL Power Estimation

    06.format_merge_data.R: Formats and merges all data

    07.calculate_bulk_eqtl_power.R: Performs bulk eQTL power calculations.

    07.Run_bulk_eQTL_power.sh: Runs the bulk eQTL power estimation pipeline.

08: sc-eQTL Power Estimation

    08.calculate_sceqtl_power.R: Calculates sc-eQTL power using scPower.

    08.Run_sceQTL_power.high.mem.sh / low.mem.sh: Submit jobs with different memory footprints.

    08.Compare_targets.r: Compares eQTL and CRISPRi target mappings.

09–11: Target Gene Analysis and Power Integration

    09.Compare_HiC_targets.R: Compares enhancer–gene pairs with Hi-C data.

    10.process_format_power.R: Formats power estimates for downstream analysis.

    11.gold_genes_and_power.R: Assesses overlap with gold-standard gene list and integrates with power metrics.

12–14: Trans-regulatory Network Analysis

    12.CompareTransNet_transeQTL.R: Evaluates overlap between CRISPRi trans networks and trans-eQTLs.

    12.Mendelian_gene_enrichment.ipynb: Assesses enrichment of Mendelian genes.

    13.IdentifyGeneswTransNetw.R: Identifies cis-genes with known trans networks.

    14.Plot_trans_Res_allPairs.R: Generates plots comparing trans-cGenes and trans-eGenes.

15–17: MetaLCL and Trans-eQTL Post-processing

    15.Extract_MetaLCL.py: Extracts trans-eQTL hits from MetaLCL summary stats.

    15.Run_Extract_MetaLCL.sh: Shell script wrapper for above.

    16.Compare_trans_results.R: Merges and compares CRISPRi and eQTL trans results.

    17.Compare_trans_results_enrichment.R: Performs enrichment analysis on trans overlaps.

Visualization and Figures

    Figure_plots.r: Generates main figures for the manuscript.

    Extended_Data_Figures.R: Scripts for creating extended data figures.

    Heatmap_by_study_type.R: Plots heatmaps showing eQTL signal by study type.

    Locus_heatmap_plot.R: Generates locus-specific heatmaps of gene associations.

Utilities

    grna_liftover.sh: Lifts over gRNA coordinates to match reference genome builds.

Usage

All scripts are modular and organized by step. Please see individual scripts for specific dependencies and inputs. Most analyses were run on HPC clusters and rely on R (v4.1+), Python 3, and Bash. Power calculations use the scPower and DESeq2 frameworks.
