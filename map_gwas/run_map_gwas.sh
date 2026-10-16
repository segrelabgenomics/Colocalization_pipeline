#!/bin/bash
#SBATCH --partition=normal
#SBATCH --mem=0

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base_dir="$repo_dir/map_gwas/"
variant_file="1000G_NYGC_unique_variants_w_header.txt" 
#or 'GTEx_V8_Unique_Variants.txt' 
#or '1000G_unique_variants_w_header.txt'
reference_panel="1KG" # 1KG or GTEx or VCF
input_file="$repo_dir/data/gwas/gwas.txt"
output_file="${input_file%.txt}_aligned_UPD_${reference_panel}_NYGC_B38.txt"

if [[ -n "${VIRTUAL_ENV_ACTIVATE:-}" ]]; then
  source "$VIRTUAL_ENV_ACTIVATE"
fi

python3 map_gwas.py \
    --gwas_file "$input_file" \
    --variants_file "$base_dir$variant_file" \
    --ref "$reference_panel" \
    --output_name "$output_file" \
    --chr_name 'CHR' \
    --bp_name 'POS' \
    --build 'b38' \
    --a1 'EA' \
    --a2 'NEA' \
    --beta 'BETA' \
    --slope 'SE' \
    --variant_col_in_reference 'variant'
