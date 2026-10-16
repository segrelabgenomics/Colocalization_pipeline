#!/bin/bash
#SBATCH --job-name=ld_interval
#SBATCH --partition=normal
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4

set -euo pipefail

if [[ -n "${VIRTUAL_ENV_ACTIVATE:-}" ]]; then
  source "$VIRTUAL_ENV_ACTIVATE"
fi

module purge
module load PLINK

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base_dir="$repo_dir"
trait="POAG_EUR"
loci="$1"
ld_interval="$base_dir/ld_interval/ld_intervals/$trait"
plink_dir="$base_dir/data/1000Genomes/NYGC_HG38/Plink_Files/"
plink_subset="$base_dir/data/1000Genomes/NYGC_HG38_unrelated/sample_subset.txt"
second_dir="$base_dir/data/GTEx/V8_release/genotypes/PLINK_files_v8/"
second_subset="$base_dir/data/GTEx/V8_release/genotypes/samples_subset.txt"

mkdir -p "$ld_interval"
# consider arguments for additional window sizes

python3 main_ld_pipeline.py --gwas_variant "$loci" \
  --gwas_trait "$trait" \
  --ld_interval_directory "$ld_interval" \
  --plink_directory "$plink_dir" \
  --plink_subset "$plink_subset" \
  --second_plink_directory "$second_dir" \
  --second_plink_subset "$second_subset"


