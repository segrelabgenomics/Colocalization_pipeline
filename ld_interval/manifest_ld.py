import pandas as pd

# Read the input file (no header)
variants = pd.read_csv("inputs/loci.txt", sep="\t", header=None)

# Construct the sbatch commands
commands = ["sbatch ./run_ld_interval.sh " + str(variant) for variant in variants[0]]

# Write to output file
with open("manifest_variants_POAG_EUR.sh", "w") as f:
    for cmd in commands:
        f.write(cmd + "\n")

