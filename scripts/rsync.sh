#!/bin/bash
# Sync prism to Rivanna and install the synced source into ~/R/rivanna-lib.
#
# Usage: bash scripts/rsync.sh
#
# Afterwards, submit the production array from the cluster:
#   ssh rivanna "cd ~/scratch/prism && sbatch simulations/submit_simulation.slurm"
#
# Excludes:
#   - libs/ (local macOS build; the package is rebuilt from source remotely)
#   - src build artifacts
#   - sim_results/, sim_raw_data/, figs/ (cluster-side outputs stay remote)
#   - .codewhale session data, scripts/, docs/, tags, pdfs, local junk

set -euo pipefail

REMOTE="rivanna"
DEST="~/scratch/prism"
SRC="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== rsync prism -> Rivanna ==="
echo "  Source : ${SRC}"
echo "  Remote : ${REMOTE}:${DEST}"
echo ""

rsync -avz --progress \
  --exclude='.git/' \
  --exclude='.codewhale/' \
  --exclude='scripts/' \
  --exclude='libs/' \
  --exclude='src/*.o' \
  --exclude='src/*.so' \
  --exclude='sim_results/' \
  --exclude='sim_raw_data/' \
  --exclude='figs/' \
  --exclude='docs/' \
  --exclude='tags' \
  --exclude='*.pdf' \
  --exclude='*.aux' \
  --exclude='*.log' \
  --exclude='.Rproj.user/' \
  --exclude='.Rhistory' \
  --exclude='.RData' \
  --exclude='.DS_Store' \
  "${SRC}/" "${REMOTE}:${DEST}/"

echo ""
echo "=== install synced source on Rivanna ==="
ssh "${REMOTE}" "
  module purge >/dev/null 2>&1 || true
  module load goolf/11.4.0_4.1.4 >/dev/null 2>&1 || true
  module load R/4.4.1 >/dev/null 2>&1
  # Do not set R_LIBS_USER here: UVA's R profile rewrites library paths when
  # it is present and can break R startup (missing 'utils', so no
  # install.packages). install_deps.R and simulation.R manage the local
  # library themselves via lib paths.
  unset R_LIBS_USER R_LIBS_SITE R_LIBS
  cd ${DEST}
  mkdir -p simulations/logs sim_results sim_raw_data
  Rscript simulations/install_deps.R
  R CMD INSTALL -l ~/R/rivanna-lib .
"

echo ""
echo "=== done ==="
echo "Submit the production array (384 conditions):"
echo "  ssh ${REMOTE} \"cd ${DEST} && sbatch simulations/submit_simulation.slurm\""
