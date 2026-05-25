#!/bin/bash
#SBATCH -n 4
#SBATCH -t 7-00:00:00
#SBATCH --array=1-4
#SBATCH --partition=standard
#SBATCH -A tongx
#SBATCH -J smriti_fiml
#SBATCH -o hpc/logs/fiml_%A_%a.out
#SBATCH -e hpc/logs/fiml_%A_%a.err

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.3.1

# Ensure log directory exists
mkdir -p hpc/logs

# Execute the simulation
# SLURM_ARRAY_TASK_ID can be used if you have split the simulation into multiple files
# or if the R script handles the array ID to subset the simulation grid.
# For now, we point it to the full grid script, but note it will run 4 times.
Rscript hpc/simulation_full_grid.R
