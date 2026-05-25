#!/bin/bash

# Rivanna Submission Wrapper for Smriti Simulation
# Adapts the Polar-Delta template for CPU-based longitudinal refinement.

PARTITION=${1:-standard}
CPUS=${2:-8}
TIME=${3:-12:00:00}
ACCOUNT=${4:-tongx}

SLURM_SCRIPT="hpc/smriti_job.slurm"

# Ensure log directory exists
mkdir -p hpc/logs

# Generate SLURM script
cat << EOB > $SLURM_SCRIPT
#!/bin/bash
#SBATCH -p $PARTITION
#SBATCH -c $CPUS
#SBATCH -t $TIME
#SBATCH -A $ACCOUNT
#SBATCH -J smriti-sim
#SBATCH -o hpc/logs/smriti_%j.out
#SBATCH -e hpc/logs/smriti_%j.err
#SBATCH --mem=16G

# Load Rivanna environment modules for R/C++ backend
module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.3.1

# Execution
mkdir -p ~/R/rivanna-lib
export R_LIBS_USER=~/R/rivanna-lib

echo "Starting Smriti Longitudinal Refinement Pipeline on $(hostname)..."
Rscript hpc/simulation_full_grid.R
EOB

chmod +x $SLURM_SCRIPT
sbatch $SLURM_SCRIPT
