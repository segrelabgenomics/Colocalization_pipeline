#!/bin/bash
#SBATCH --job-name=eCAV
#SBATCH --partition=normal
#SBATCH --mem=35G
#SBATCH --export=ALL

mkdir -p logs

exec > "logs/results_${11}.txt"
exec 2> "logs/error_${11}.txt"

echo "Running eCAVIAR job for CustomJobID: ${11}"
echo "GWAS Trait: $4, Tissue: $5, Gene: $6, QTL: $3"

module purge
module load caviar
module load PLINK #PLINK/1.90-beta-7.11-x86_64 
module load Python/3.12.12-GCCcore-13.3.0

if [[ -n "${VIRTUAL_ENV_ACTIVATE:-}" ]]; then
  source "$VIRTUAL_ENV_ACTIVATE"
fi

main="colocalization_pipeline_just_ecaviar.py"
echo $5

python3 "$main" --gwas_variant "$1" \
  --gwas_file "$2" \
  --qtl_type "$3" \
  --gwas_trait "$4" \
  --tissue "$5" \
  --gene "$6" \
  --plink_directory "$7" \
  --plink_subset "$8" \
  --output_directory "$9" \
  --ld_directory "${10}"

echo "Job completed successfully."
