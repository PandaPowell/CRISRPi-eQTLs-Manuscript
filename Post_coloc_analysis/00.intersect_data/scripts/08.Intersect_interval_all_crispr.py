#!/usr/bin/env python
# coding: utf-8

# In[23]:


#!/usr/bin/env python3
import os
import pandas as pd

# --- config ---
EQTL = "Interval"
BASE_DIR = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

RESULTS_ROOT = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/Coloc_results_V2"
CREDSET_ROOT = os.path.join(BASE_DIR, "sting_seq_credible_sets")

OUT_DIR = os.path.join(BASE_DIR, "stingseq_eqtl_overlap")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_FP = os.path.join(OUT_DIR, f"{EQTL}_coloc_stingseq.txt")
# -------------

merged_chunks = []

for chr_i in range(1, 23):  # chr1..chr22
    chr_name = f"chr{chr_i}"
    for gwas_id in range(30000, 30301, 10):  # 30000, 30010, ..., 30300
        coloc_fp = f"{RESULTS_ROOT}/{EQTL}/{gwas_id}/{chr_name}_coloc_results.txt"
        cred_fp  = f"{CREDSET_ROOT}/{gwas_id}/{gwas_id}_{chr_name}_stingseq_credset.txt"

        if not (os.path.exists(coloc_fp) and os.path.exists(cred_fp)):
            continue

        try:
            # coloc is CSV (commas)
            coloc = pd.read_csv(coloc_fp, sep=",", dtype=str, low_memory=False)

            # sting_seq file is typically tab/whitespace delimited; try tab first, then whitespace
            try:
                sting_seq = pd.read_csv(cred_fp, sep="\t", dtype=str, header=0, low_memory=False)
            except Exception:
                sting_seq = pd.read_csv(cred_fp, sep=r"\s+", dtype=str, header=0, engine="python", low_memory=False)

            if "SNP" not in coloc.columns or "finemap_snp" not in sting_seq.columns:
                continue  # skip if required columns are missing

            m = pd.merge(
                coloc,
                sting_seq,
                left_on="SNP",
                right_on="finemap_snp",
                how="inner"
            )
            if not m.empty:
                merged_chunks.append(m)

        except Exception as e:
            # Skip problematic pairs, keep going
            print(f"Skipping {gwas_id} {chr_name}: {e}")

# Write all results to a single file
if merged_chunks:
    final_df = pd.concat(merged_chunks, ignore_index=True)
    final_df.to_csv(OUT_FP, index=False)
    print(f"Wrote {len(final_df)} rows to {OUT_FP}")
else:
    # Still create an empty file if nothing merged
    pd.DataFrame().to_csv(OUT_FP, index=False)
    print(f"No merges produced rows. Created empty file at {OUT_FP}")

