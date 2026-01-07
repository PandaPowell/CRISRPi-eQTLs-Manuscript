#!/bin/bash

# First reformat Hi-C file so that interacting region is in bed format

awk '{print $9"\t"$10"\t"$11"\t"$8}' K562.hg19.AllInteractions.SP4.FDR0.1.txt | tail -n +2 > K562.hg19.AllInteractions.SP4.FDR0.1.bed

sort -k1,1 -k2,2n K562.hg19.AllInteractions.SP4.FDR0.1.bed > K562.hg19.AllInteractions.SP4.FDR0.1.sorted.bed

tail -n +2 ../../cres_with_grnas.txt | awk -F, '{print "chr"$10"\t"$12-1"\t"$12"\t"$11}' | sort -u | sort -k1,1 -k2,2n > cres_w_grnas.bed

module load bedtools/2.25.0

echo -e "Interactor_Chr\tInteractor_Start\tInteractor_End\tInteractorID\tgrna_chr\tgrna_target_lower\tgrna_target_upper\tgrna_target\tdistance" > cres_w_grnas_HiC_interactions.bed

bedtools closest -d -wb -a K562.hg19.AllInteractions.SP4.FDR0.1.sorted.bed -b  cres_w_grnas.bed | \
awk '$9 != "-1"' | awk '$9 < 2000' | sort -k9,9 >> cres_w_grnas_HiC_interactions.bed
