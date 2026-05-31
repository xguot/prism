#!/bin/bash
#SBATCH -p standard
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --array=1-30
#SBATCH -A tongx
#SBATCH -J smriti_tune
#SBATCH -o hpc/logs/tune_%A_%a.out
#SBATCH -e hpc/logs/tune_%A_%a.err
#SBATCH --mem=16G

# ── Tuning Array: rapid FIML-vs-Smriti hyperparameter sweep ─────────────────
# Splits the 30-condition grid (5 N × 3 miss × 1 dist × 2 mech) across
# 30 array tasks, one condition each. Each task loops over 5 lambda values
# internally, producing a single tune_results_<task>.rds file.
#
# After the array completes, aggregate with:
#   Rscript -e 'files <- list.files("sim_results", pattern="tune_results_.*\\\\.rds",
#               full.names=TRUE); all <- do.call(rbind, lapply(files, readRDS));
#               saveRDS(all, "sim_results/tune_results.rds")'

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.4.1

# Ensure log directory exists
mkdir -p hpc/logs
mkdir -p sim_results

# Neutralize multi-threaded BLAS to prevent contention
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "Starting Smriti Tuning Task ${SLURM_ARRAY_TASK_ID} on $(hostname)..."
Rscript hpc/simulation_tune.R
