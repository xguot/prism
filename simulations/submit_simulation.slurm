#!/bin/bash
# PRISM production array on Rivanna (pure runner; installs happen in
# scripts/rsync.sh, which also creates simulations/logs before submission).
#
# Workflow:
#   bash scripts/rsync.sh
#   ssh rivanna "cd ~/scratch/prism && sbatch simulations/submit_simulation.slurm"
#
# The --array size must equal the condition grid in simulation.R:
# 2 models x 6 N x 4 miss x 4 dist x 2 mech = 384 conditions.
# simulation.R stops with a clear error on any out-of-range task id instead
# of silently running the full grid in one job.
#SBATCH -p standard
#SBATCH -c 8
#SBATCH -t 18:00:00
#SBATCH --array=1-384
#SBATCH -A tongx
#SBATCH -J prism_prod
#SBATCH -o simulations/logs/prod_%A_%a.out
#SBATCH -e simulations/logs/prod_%A_%a.err
#SBATCH --mem=16G

set -euo pipefail

# Load Rivanna environment modules
module purge
module load goolf/11.4.0_4.1.4
module load R/4.4.1

# Do not set R_LIBS_USER: UVA's R profile rewrites library paths when it is
# present and can break R startup (missing 'utils'). simulation.R prepends
# ~/R/rivanna-lib to .libPaths itself, and the package was installed there
# by scripts/rsync.sh.
unset R_LIBS_USER R_LIBS_SITE R_LIBS || true

# Slurm starts the job in the directory sbatch was invoked from (the
# repository root in the rsync workflow).  SLURM_SUBMIT_DIR points there
# regardless of how the script path is resolved, so use it when present.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "${SLURM_SUBMIT_DIR}"
fi

# Fail fast with a clear message instead of an R fatal error if the
# repository layout is not where expected
if [ ! -f simulations/simulation.R ]; then
  echo "ERROR: simulations/simulation.R not found. Submit from the repository"
  echo "root: ssh rivanna \"cd ~/scratch/prism && sbatch simulations/submit_simulation.slurm\""
  exit 1
fi

# Neutralize multi-threaded BLAS to prevent contention during R-level
# parallelization across the 8 allocated cores
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "Starting prism Production Task ${SLURM_ARRAY_TASK_ID} on $(hostname)..."
Rscript simulations/simulation.R
