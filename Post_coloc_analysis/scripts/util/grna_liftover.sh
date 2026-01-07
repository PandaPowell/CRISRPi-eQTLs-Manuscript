#!/bin/bash
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=liftover_rsids

module load bedtools

# -------- paths you may want to edit --------
CSV="processed_data/cres_with_grnas.txt"   # the file you showed (comma-separated)
DBSNP37_DIR="$home/genome/dbsnp151_GRCh37p13"
DBSNP38_DIR="$home/genome/dbsnp151_GRCh38p7"
OUTDIR="temp"
CHRS="$(seq 1 22)"          # add 'X Y' if needed: CHRS="$(seq 1 22) X Y"
# --------------------------------------------

mkdir -p "$OUTDIR"

# Combined output header
COMBINED="${OUTDIR}/hg19_hg38_rsids_all.csv"
echo "hg19_chr,hg19_start,hg19_end,ensembl_id,gene_name,grna_pos,rsid,hg38_chr,hg38_start,hg38_end" > "$COMBINED"

for i in $CHRS; do
  echo ">> chr${i}"

  # 1) Build hg19 BED from your CSV (skip header), using SNP in column 7: CHR:POS:REF:ALT
  #    Columns (1-based):
  #      4 = ensembl_id, 5 = gene_name, 7 = finemap_snp_intersect_grna ("CHR:POS:REF:ALT")
  #    Keep only rows for chr == $i. Output A (6 cols):
  #      chr  start  end  ensembl_id  gene_name  snp_id_string
  awk -F',' -v i="$i" '
    NR>1 {
      # parse CHR:POS:REF:ALT from col 7
      split($7, a, ":"); chr=a[1]; pos=a[2];
      # keep only this chromosome; POS must be numeric
      if (chr == i && pos ~ /^[0-9]+$/) {
        # BED is 0-based start, 1-based end
        printf("chr%s\t%d\t%d\t%s\t%s\t%s\n", chr, pos-1, pos, $4, $5, $7);
      }
    }
  ' "$CSV" \
    | sort -k1,1 -k2,2n -k3,3n -u > "${OUTDIR}/hg19_chr${i}.bed"

  # 2) Intersect with dbSNP GRCh37 to fetch rsIDs at those loci
  #    B is expected: col4 = rsID (standard dbSNP BED)
  bedtools intersect -wa -wb \
    -a "${OUTDIR}/hg19_chr${i}.bed" \
    -b "${DBSNP37_DIR}/sorted_bed_chr_${i}.bed.gz" \
    > "${OUTDIR}/hg19_rsids_chr${i}.txt"

  if [[ ! -s "${OUTDIR}/hg19_rsids_chr${i}.txt" ]]; then
    echo "   (no rsIDs found at these loci on chr${i})"
    continue
  fi

  # 3) Join by rsID to dbSNP GRCh38 to get hg38 coordinates
  #    In the -wa -wb output:
  #      A: $1..$6 = hg19_chr, hg19_start, hg19_end, ensembl_id, gene_name, grna_pos
  #      B: $7..   = chr, start, end, rsID, ...
  #      => rsID is overall $10  (6 A cols + B col4)
  awk -v FS='\t' -v OFS=',' '
    FNR==NR {
      # Build map from rsID -> A-fields (6 cols)
      rs = $10
      if (rs ~ /^rs[0-9]+$/) {
        key = rs
        a   = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6  # hg19_chr,hg19_start,hg19_end,ensembl_id,gene_name,grna_pos
        map[key] = a
      }
      next
    }
    {
      # dbSNP GRCh38: chr=$1, start=$2, end=$3, rsID=$4
      rs38 = $4
      if (rs38 in map) {
        split(map[rs38], x, OFS)
        # print: A6 + rsid + hg38 3
        print x[1], x[2], x[3], x[4], x[5], x[6], rs38, $1, $2, $3
      }
    }
  ' "${OUTDIR}/hg19_rsids_chr${i}.txt" <(zcat "${DBSNP38_DIR}/bed_chr_${i}.bed.gz") \
    > "${OUTDIR}/hg19_hg38_rsids_chr${i}.csv"

  if [[ -s "${OUTDIR}/hg19_hg38_rsids_chr${i}.csv" ]]; then
    cat "${OUTDIR}/hg19_hg38_rsids_chr${i}.csv" >> "$COMBINED"
  else
    echo "   (no rsID matches to hg38 on chr${i})"
  fi
done

echo "Done. Combined: $COMBINED"
