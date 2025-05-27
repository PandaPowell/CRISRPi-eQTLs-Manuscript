#!/usr/bin/env python
# coding: utf-8

# In[1]:


import duckdb
import glob
import pandas as pd
from joblib import Parallel, delayed

# Load Gasperini cis-genes with trans hits
gas_df = pd.read_csv("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/genes_cis_w_trans.txt", sep='\t', header=None)
gas_cis_genes = set(gas_df[0].tolist())  # Convert to sorted list

# Load morris cis-genes with trans hits
trans_network_john = pd.read_excel("/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/science.adh7699_table_s4.xlsx", skiprows=2)

morrs_cis_genes = ["GFI1B", "NFE2", "IKZF1", "HHEX", "RUNX1"]

# Load MAGE fine-mapped cis-genes
file_path = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/MAGE/MAGE.v1.0.data/QTL_results/eQTL_results/eQTL_finemapping_results/eQTL_finemapping.significantAssociations.MAGE.v1.0.txt.gz"
mage = pd.read_csv(file_path, compression='gzip', sep='\t', low_memory=False)

# Filter for cis-genes
filtered_mage = mage[mage["geneSymbol"].isin(gas_cis_genes.union(morrs_cis_genes))]


# In[ ]:


# Directory containing Parquet files
parquet_dir = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/MetaLCL/"
parquet_files = glob.glob(parquet_dir + "*.parquet")

# Convert chromosome names safely
filtered_mage = filtered_mage.assign(variantChrom=filtered_mage["variantChrom"].str.replace("chr", "", regex=True))

# Generate WHERE clause for filtering
conditions = " OR ".join(
    [f"(CAST(chromosome AS STRING) = '{row.variantChrom}' AND position = {row.variantPosition})"
     for _, row in filtered_mage.iterrows()]
)

def process_parquet(file):
    """Process a single Parquet file and return filtered results."""
    print(f"Processing file {file}", flush=True)
    try:
        query = f"""
        SELECT * FROM parquet_scan('{file}')
        WHERE {conditions}
        """
        df = duckdb.query(query).to_df()

        if not df.empty:
            df["chromosome"] = df["chromosome"].astype(str)
            filtered_mage["variantChrom"] = filtered_mage["variantChrom"].astype(str)

            # Merge additional columns
            df_merged = df.merge(filtered_mage, left_on=["chromosome", "position"],
                                 right_on=["variantChrom", "variantPosition"], how="left")

            df_merged.drop(columns=["variantChrom", "variantPosition"], inplace=True)

            if "nlog10p" in df_merged.columns and "molecular_trait_id" in df_merged.columns and "geneSymbol" in df_merged.columns:
                df_unique = df_merged.loc[df_merged.groupby(["molecular_trait_id", "geneSymbol"])["nlog10p"].idxmax()]
                df_unique = df_unique.rename(columns={"geneSymbol": "cis_gene"})
                return df_unique

        return None

    except Exception as e:
        print(f"Error processing {file}: {e}")
        return None

# Parallel execution
num_cores = 16  # Adjust this based on available resources
filtered_data = Parallel(n_jobs=num_cores)(delayed(process_parquet)(file) for file in parquet_files)

# Filter out None values
filtered_data = [df for df in filtered_data if df is not None]

# Combine results
if filtered_data:
    final_df = pd.concat(filtered_data, ignore_index=True)
    final_df.to_parquet("MetaLCL_extracted_transQTLs.parquet", index=False)
    print(f"Extracted {len(final_df)} rows and saved to 'MetaLCL_extracted_transQTLs.parquet'.")
else:
    print("No matching rows found in any file.")


# In[2]:


# Query the first 5 rows directly
result = duckdb.sql("SELECT * FROM 'MetaLCL_extracted_transQTLs.parquet' LIMIT 5").df()

# Display result
print(result)


# In[ ]:




