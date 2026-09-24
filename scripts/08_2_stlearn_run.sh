#!/usr/bin/env bash
# Run 08_2_stlearn_run.py (grid mode) for every sample, one process + log per sample.
# Usage: bash 08_2_stlearn_run.sh [sample ...]   (no args = all <sample>.h5ad in IN_DIR)
set -uo pipefail
SCRIPT_DIR=/mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium/scripts
BASE=/mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium
IN_DIR=$BASE/data/08_stlearn
OUT_DIR=$BASE/data/08_stlearn/cci
FIG_DIR=$BASE/figures/08_stlearn
LOG_DIR=$BASE/logs/08_stlearn
N_CPUS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$LOG_DIR"

if [ $# -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for f in "$IN_DIR"/*.h5ad; do SAMPLES+=("$(basename "$f" .h5ad)"); done
fi
echo "Samples: ${SAMPLES[*]}"

FAILED=()
for s in "${SAMPLES[@]}"; do
    log="$LOG_DIR/${s}_grid250_$(date +%Y%m%d_%H%M).log"
    echo "[$(date +%T)] $s -> $log"
    python -u "$SCRIPT_DIR/run_stlearn_cci.py" \
        --samples "$s" \
        --database connectomeDB2020_lit \
        --input-dir  "$IN_DIR" \
        --output-dir "$OUT_DIR" \
        --fig-dir    "$FIG_DIR" \
        --grid-n 250 \
        --grid-neighbours 1.5 \
        --min-spots 20 \
        --n-pairs 10000 \
        --n-cpus "$N_CPUS" \
        --annotation annotation \
        --cci-min-spots 5 \
        --n-perms 1000 \
        > "$log" 2>&1 || FAILED+=("$s")
done

echo "[$(date +%T)] Done."
[ ${#FAILED[@]} -eq 0 ] || { echo "Failed: ${FAILED[*]} (see $LOG_DIR)"; exit 1; }
