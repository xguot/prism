#!/bin/bash
#SBATCH -p standard
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --array=1-90
#SBATCH -A tongx
#SBATCH -J smriti_beast
#SBATCH -o hpc/logs/beast_%A_%a.out
#SBATCH -e hpc/logs/beast_%A_%a.err
#SBATCH --mem=16G

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.3.1

# Ensure log directory exists
mkdir -p hpc/logs

# Neutralize multi-threaded BLAS to prevent contention during R-level parallelization
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Execute the simulation for the current array task (condition)
echo "Starting Smriti Beast Mode Task ${SLURM_ARRAY_TASK_ID} on $(hostname)..."
Rscript hpc/simulation_beast.R
