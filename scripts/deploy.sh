#!/bin/bash
# One-step verified deploy script to Rivanna
# Syncs code, verifies source fix, re-installs cleanly, verifies installed binary, runs manual test, and submits.

set -euo pipefail

REMOTE="rivanna"
DEST="~/scratch/prism"
SRC="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== 1. Syncing codebase to Rivanna ==="
ssh "${REMOTE}" "mkdir -p ${DEST}/simulations/logs ${DEST}/sim_results ${DEST}/sim_raw_data"

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
echo "=== 2. Running Remote Verification & Deployment ==="
ssh "${REMOTE}" "bash -s" << 'EOF'
set -euo pipefail

DEST="${HOME}/scratch/prism"
cd $DEST

echo "-> [0/5] Loading cluster modules..."
module purge && module load goolf/11.4.0_4.1.4 && module load R/4.4.1 && module load java/21.0.11
export R_LIBS_USER=~/R/rivanna-lib


echo "-> [1/5] Verifying the source on Rivanna has the fix..."
if ! grep -q "as.matrix(implied" R/stochastic_fiml_impute.R; then
    echo "ERROR: Source does not contain the fix! Rsync must have failed."
    exit 1
fi
echo "✓ Source fix verified."

echo "-> [2/5] Reinstalling package cleanly..."
Rscript -e 'try(remove.packages("prism", lib="~/R/rivanna-lib"), silent=TRUE)'
Rscript -e 'install.packages(".", repos=NULL, type="source", lib="~/R/rivanna-lib")'

echo "-> [3/5] Verifying the INSTALLED binary has the fix..."
if ! Rscript -e 'library(prism, lib.loc="~/R/rivanna-lib"); cat(grep("fitted", deparse(body(prism_mi)), value=TRUE))' | grep -q "fitted"; then
    echo "ERROR: Installed binary does not contain the fix! The install silently failed or used a cached version."
    exit 1
fi
echo "✓ Installed binary fix verified."

echo "-> [4/5] Running manual test..."
TEST_OUT=$(Rscript -e '
suppressPackageStartupMessages({
  library(prism, lib.loc="~/R/rivanna-lib")
  library(lavaan)
  library(MASS)
})
set.seed(123)
n <- 100
latent <- MASS::mvrnorm(n, mu=c(6,2), Sigma=matrix(c(1,0,0,1),2,2))
df <- data.frame(T1=latent[,1]+rnorm(n), T2=latent[,1]+latent[,2]+rnorm(n), T3=latent[,1]+2*latent[,2]+rnorm(n), T4=latent[,1]+3*latent[,2]+rnorm(n))
cut <- sort(df$T1)[16]
df$T2[df$T1<cut] <- NA
df$T3[df$T1<cut] <- NA
df$T4[df$T1<cut] <- NA
fit <- growth("i=~1*T1+1*T2+1*T3+1*T4; s=~0*T1+1*T2+2*T3+3*T4; i~~s; i~~i; s~~s", data=df, missing="fiml")
imp <- stochastic_fiml_impute(df, fit)
cat("SUCCESS: rows =", nrow(imp), "\n")
')

echo "$TEST_OUT"

if echo "$TEST_OUT" | grep -q "SUCCESS"; then
    echo "✓ Test passed."

    echo "-> [5/5] Submitting simulation array..."
    # sbatch resolves the -o/-e log paths at launch time; the directory
    # must already exist when the job is submitted
    mkdir -p simulations/logs
    sbatch simulations/submit_simulation.slurm
else
    echo "ERROR: Manual test failed. Array submission aborted."
    exit 1
fi
EOF
