#!/bin/bash
# rsync prism to Rivanna HPC scratch space
# Usage: bash scripts/rsync.sh
#
# Excludes:
#   - Compiled objects and shared libs (will rebuild on Rivanna)
#   - Previous sim_results (fresh run)
#   - .codewhale session data
#   - .git directory (optional — uncomment to include for reproducibility)
#   - docs/ (will regenerate)
#   - scripts/

set -euo pipefail

REMOTE="rivanna"
DEST="~/scratch/prism"

# ── Source directory (this repo root) ───────────────────────────────────────
SRC="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== rsync prism → Rivanna ==="
echo "  Source : ${SRC}"
echo "  Remote : ${REMOTE}:${DEST}"
echo ""

# Ensure remote scratch directory exists
ssh "${REMOTE}" "mkdir -p ${DEST}/simulations/logs ${DEST}/sim_results ${DEST}/sim_raw_data"

# Rsync with sensible exclusions
rsync -avz --progress \
  --exclude='.git/' \
  --exclude='.codewhale/' \
  --exclude='scripts/' \
  --exclude='src/*.o' \
  --exclude='src/*.so' \
  --exclude='src/prism.so' \
  --exclude='sim_results/prod_results*.rds' \
  --exclude='sim_results/tune_results*.rds' \
  --exclude='docs/' \
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
echo "  3. Rscript simulations/install_deps.R           # install R dependencies"
echo "  4. Rscript -e 'install.packages(\".\", repos=NULL, type=\"source\")'  # install prism"
echo "  5. sbatch simulations/run_tune_array.slurm      # tuning study first"
echo "  6. sbatch simulations/run_production_array.slurm # production run"
