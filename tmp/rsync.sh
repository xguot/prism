#!/bin/bash
# rsync smriti to Rivanna HPC scratch space
# Usage: bash tmp/rsync.sh [netbadge_id]
#
# If netbadge_id is omitted, defaults to 'xguo' (change as needed).
#
# Excludes:
#   - Compiled objects and shared libs (will rebuild on Rivanna)
#   - Previous sim_results (fresh run)
#   - .codewhale session data
#   - .git directory (optional — uncomment to include for reproducibility)
#   - manuscript_figures (will regenerate)
#   - tmp/

set -euo pipefail

NETBADGE="${1:-xguo}"
REMOTE="${NETBADGE}@rivanna.hpc.virginia.edu"
DEST="~/scratch/smriti"

# ── Source directory (this repo root) ───────────────────────────────────────
SRC="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== rsync smriti → Rivanna ==="
echo "  Source : ${SRC}"
echo "  Remote : ${REMOTE}:${DEST}"
echo ""

# Ensure remote scratch directory exists
ssh "${REMOTE}" "mkdir -p ${DEST}/hpc/logs ${DEST}/sim_results ${DEST}/sim_raw_data"

# Rsync with sensible exclusions
rsync -avz --progress \
  --exclude='.git/' \
  --exclude='.codewhale/' \
  --exclude='tmp/' \
  --exclude='src/*.o' \
  --exclude='src/*.so' \
  --exclude='src/smriti.so' \
  --exclude='sim_results/prod_results*.rds' \
  --exclude='sim_results/tune_results*.rds' \
  --exclude='manuscript_figures/' \
  --exclude='.Rproj.user/' \
  --exclude='.Rhistory' \
  --exclude='.RData' \
  --exclude='.DS_Store' \
  "${SRC}/" "${REMOTE}:${DEST}/"

echo ""
echo "=== rsync complete ==="
echo ""
echo "Next steps on Rivanna:"
echo "  1. ssh ${REMOTE}"
echo "  2. cd ${DEST}"
echo "  3. Rscript hpc/install_deps.R               # install R dependencies"
echo "  4. Rscript -e 'install.packages(\".\", repos=NULL, type=\"source\")'  # install smriti"
echo "  5. sbatch hpc/run_tune_array.slurm          # tuning study first"
echo "  6. sbatch hpc/run_production_array.slurm     # production run"
