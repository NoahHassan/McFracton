#!/bin/bash
# Example Slurm array job: one temperature sweep per task, all from one config file.
#
# Each task overrides only what makes it different from its neighbours - here the coupling KU -
# and writes its own log, so nothing two tasks do can collide. Adapt the partition, the time limit
# and the module lines to your cluster before using it.
#
#   sbatch --array=0-7 scripts/slurm_array.sh config/spiderweb.cfg
#
#SBATCH --job-name=mcfracton
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

CONFIG="${1:?usage: sbatch --array=0-N scripts/slurm_array.sh <config file>}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="$REPO/build/linux-release/mcf_run"
RESULTS="$REPO/results"

if [ ! -x "$RUN" ]; then
    echo "mcf_run not built: run cmake --preset linux-release && cmake --build --preset linux-release" >&2
    exit 1
fi

mkdir -p "$RESULTS" logs

# The run records which commit produced it. mcf_run does not compile the hash in, so that neither
# build system needs a generated header - it is passed in here instead.
GIT_HASH="$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)"

# One value of the coupling per array index. Replace this with whatever your scan varies; any key
# that "mcf_run --list" shows, or any NumericalParams field, can be overridden the same way.
KU_VALUES=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
KU="${KU_VALUES[$TASK]}"

# A seed per task, so the tasks are independent but the whole array stays reproducible.
SEED=$((1000 + TASK))

OUT="$RESULTS/$(basename "$CONFIG" .cfg)_KU=${KU}_task=${TASK}.txt"

echo "task $TASK: KU=$KU seed=$SEED -> $OUT"

srun "$RUN" --config "$CONFIG" \
    --KU="$KU" \
    --seed="$SEED" \
    --out="$OUT" \
    --git_hash="$GIT_HASH"
