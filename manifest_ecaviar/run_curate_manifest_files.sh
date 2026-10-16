#!/bin/bash
#SBATCH --job-name=manifest
#SBATCH --partition=normal
#SBATCH --mem=5G

set -euo pipefail

module purge
if [[ -n "${VIRTUAL_ENV_ACTIVATE:-}" ]]; then
  source "$VIRTUAL_ENV_ACTIVATE"
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base_dir="$repo_dir"
trait="POAG_EUR"
sig_genes_directory="$base_dir/identify_sig_genes/$trait/"
manifest_df_dir="$base_dir/manifest_files/"
manifest_chunk_out_dir="$base_dir/manifest_files/chunked_manifest_files/"
plink_dir="$base_dir/data/1000Genomes/NYGC_HG38/Plink_Files/"
plink_subset="$base_dir/data/1000Genomes/NYGC_HG38_unrelated/sample_subset.txt"
gwas_file="$base_dir/data/gwas/gwas_aligned_UPD_1KG_NYGC_B38.txt"

python3 curate_manifest_files.py --sig_genes_directory "$sig_genes_directory" \
  --gwas_trait "$trait" \
  --gwas_file "$gwas_file" \
  --plink_directory "$plink_dir" \
  --plink_subset "$plink_subset" \
  --manifest_dir "$manifest_df_dir" \
  --manifest_chunks_dir "$manifest_chunk_out_dir"
