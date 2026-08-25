#!/bin/bash
#SBATCH -p standard
#SBATCH -c 8
#SBATCH -t 18:00:00
#SBATCH --array=1-384
#SBATCH -A tongx
#SBATCH -J prism_prod
#SBATCH -o simulations/logs/prod_%A_%a.out
#SBATCH -e simulations/logs/prod_%A_%a.err
#SBATCH --mem=16G

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.4.1

# Local R library (no root access on Rivanna)
export R_LIBS_USER=~/R/rivanna-lib

# Log directory must exist BEFORE the job starts (sbatch resolves -o/-e at
# launch time); submit from the repository root
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Execute the simulation for the current array task (condition)
# Grid: 2 models x 6 N x 4 miss x 4 dist x 2 mech = 384 conditions
echo "Starting prism Production Task ${SLURM_ARRAY_TASK_ID} on $(hostname)..."
Rscript simulations/simulation.R
