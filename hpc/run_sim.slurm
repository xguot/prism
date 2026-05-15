#!/bin/bash
#SBATCH --job-name=smriti_sim
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --output=hpc/logs/sim_%j.out
#SBATCH --error=hpc/logs/sim_%j.err

# Load Rivanna environment modules
# Note: Specific versions may vary based on cluster updates
module purge
module load gcc/11.4.0
module load openblas/0.3.24
module load R/4.3.1

# Ensure the log directory exists
mkdir -p hpc/logs

# Execute the simulation scaling to 1,000 replications
# We use a fixed seed for reproducibility across HPC nodes
Rscript tests/simulation_study.R --reps 1000 --seed 42
