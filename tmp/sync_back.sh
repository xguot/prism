#!/bin/bash
# rsync simulation results back from Rivanna
# Usage: bash tmp/sync_back.sh

set -euo pipefail

REMOTE="rivanna"
SRC="~/scratch/smriti"

echo "=== rsync results ← Rivanna ==="
echo "  Remote : ${REMOTE}:${SRC}"
echo ""

# Pull aggregated production results (if already combined)
rsync -avz --progress \
  "${REMOTE}:${SRC}/sim_results/prod_results.rds" \
  "sim_results/"

# Pull per-task production files (if not yet aggregated)
rsync -avz --progress \
  --include='prod_results_*.rds' \
  --exclude='*' \
  "${REMOTE}:${SRC}/sim_results/" \
  "sim_results/"

# Pull tuning results
rsync -avz --progress \
  --include='tune_results*.rds' \
  --exclude='*' \
  "${REMOTE}:${SRC}/sim_results/" \
  "sim_results/"

echo ""
echo "=== rsync complete ==="
echo ""
echo "If per-task files were pulled, aggregate locally:"
echo "  Rscript -e 'files <- list.files(\"sim_results\", pattern=\"prod_results_.*[.]rds\", full.names=TRUE); all <- do.call(rbind, lapply(files, readRDS)); saveRDS(all, \"sim_results/prod_results.rds\")'"
echo ""
echo "Then run:"
echo "  Rscript manuscript_analysis.R"
echo "  Rscript manuscript_plots.R"
