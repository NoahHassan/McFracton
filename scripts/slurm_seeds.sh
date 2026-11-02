#!/bin/bash
# Slurm array job: the same run once per seed, so the tasks are independent chains of one
# parameter set and the error bars come from the spread across them. Adapt the partition, the time
# limit and the module lines to your cluster before using it.
#
#   mkdir -p logs
#   sbatch --array=0-15 scripts/slurm_seeds.sh config/spiderweb_corrected.cfg
#
# Slurm does not create the directory of --output itself, hence the mkdir. Anything after the config file is passed on to mcf_run, so one setting can be changed for a whole
# array without editing the file:
#
#   sbatch --array=0-15 scripts/slurm_seeds.sh config/spiderweb_corrected.cfg --KU=0.5
#
#SBATCH --job-name=mcfracton
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

CONFIG="${1:?usage: sbatch --array=0-N scripts/slurm_seeds.sh <config file> [--key=value ...]}"
shift
TASK="${SLURM_ARRAY_TASK_ID:-0}"

# Slurm runs a copy of this script from its spool directory, so the repository is taken to be the
# directory the job was submitted from: submit from the repository root.
REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RUN="$REPO/build/linux-release/mcf_run"
RESULTS="$REPO/results"

if [ ! -x "$RUN" ]; then
    echo "mcf_run not built: run cmake --preset linux-release && cmake --build --preset linux-release" >&2
    exit 1
fi

mkdir -p "$RESULTS"

# The run records which commit produced it. mcf_run does not compile the hash in, so that neither
# build system needs a generated header - it is passed in here instead.
GIT_HASH="$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)"

# A seed per task, so the tasks are independent but the whole array stays reproducible.
SEED=$((1000 + TASK))

# The array job id keeps two submissions of the same config from overwriting each other.
OUT="$RESULTS/$(basename "$CONFIG" .cfg)_job=${SLURM_ARRAY_JOB_ID:-local}_seed=${SEED}.txt"

echo "task $TASK: seed=$SEED -> $OUT"

srun "$RUN" --config "$CONFIG" "$@" \
    --seed="$SEED" \
    --out="$OUT" \
    --git_hash="$GIT_HASH"
