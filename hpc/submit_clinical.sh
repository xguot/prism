#!/bin/bash

# Rivanna Submission Wrapper for ADNI Clinical Validation
# Adapts the Polar-Delta template for empirical longitudinal analysis.

PARTITION=${1:-standard}
CPUS=${2:-8}
TIME=${3:-12:00:00}
ACCOUNT=${4:-tongx}

SLURM_SCRIPT="hpc/clinical_job.slurm"

# Ensure log directory exists
mkdir -p hpc/logs

# Generate SLURM script
cat << EOB > $SLURM_SCRIPT
#!/bin/bash
#SBATCH -p $PARTITION
#SBATCH -c $CPUS
#SBATCH -t $TIME
#SBATCH -A $ACCOUNT
#SBATCH -J smriti-clinical
#SBATCH -o hpc/logs/clinical_%j.out
#SBATCH -e hpc/logs/clinical_%j.err
#SBATCH --mem=32G

# Load Rivanna environment modules
module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.3.1

# Execution
export R_LIBS_USER=~/R/rivanna-lib

echo "Starting ADNI Clinical Structural Validation on \$(hostname)..."
Rscript empirical/adni_structural_validation.R
EOB

chmod +x $SLURM_SCRIPT
sbatch $SLURM_SCRIPT
