#!/bin/bash
#SBATCH -n 4
#SBATCH -t 7-00:00:00
#SBATCH --partition=standard
#SBATCH -A tongx
#SBATCH -J smriti_sim
#SBATCH -o hpc/logs/sim_%j.out
#SBATCH -e hpc/logs/sim_%j.err
#SBATCH --mem=16G

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.3.1

# Ensure the log directory exists
mkdir -p hpc/logs

# Execute the simulation
Rscript hpc/simulation_full_grid.R
