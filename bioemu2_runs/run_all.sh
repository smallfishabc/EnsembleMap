#!/usr/bin/env bash
# Sample every system in systems/, or those matching the prefix given as $1.
#
#   ./run_all.sh            all of them
#   ./run_all.sh p53_       just the p53 ones
#
# Unsteered and one fixed seed throughout: the comparison is between sequences, so the
# sampling noise should be shared rather than varied.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${1:-}"
NUM_SAMPLES="${NUM_SAMPLES:-2000}"
BATCH_SIZE="${BATCH_SIZE:-64}"
SEED="${SEED:-0}"
OUT_ROOT="${OUT_ROOT:-$HERE/runs}"

mkdir -p "$OUT_ROOT"

for yaml in "$HERE"/systems/${PREFIX}*.yaml; do
    name="$(basename "$yaml" .yaml)"
    out="$OUT_ROOT/$name"
    if [ -d "$out" ]; then
        echo "skip  $name  (output directory exists)"
        continue
    fi
    echo "run   $name"
    bioemu2 sample-equilibrium \
        --system "$yaml" \
        --out "$out" \
        --num-samples "$NUM_SAMPLES" \
        --batch-size "$BATCH_SIZE" \
        --sampler edm \
        --seed "$SEED"
done

echo "done. output under $OUT_ROOT"
