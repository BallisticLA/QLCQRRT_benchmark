#!/bin/bash
# run_orth_gap_benchmark.sh
# Re-generates all CQRRT_diagnostic data for the orth-gap figure.
#
# Produces one diagnostic_*.csv per matrix in:
#   results/orth-gap/
#
# Usage: bash run_orth_gap_benchmark.sh [path-to-CQRRT_diagnostic-binary]
#
# Defaults to the standard build path if no argument given.

BENCH=${1:-/home/mymel/RandNLA/RandNLA-project/build/benchmark-build/CQRRT_diagnostic}
OUTDIR=$(dirname "$0")/../results/orth-gap
MTX=/home/mymel/matlab/QLCQRRT_benchmark/input_matrices/photogrammetry2/photogrammetry2.mtx

mkdir -p "$OUTDIR"

# Synthetic matrices: m=5000, n=1000, density=0.05, d=2n, nnz=4, runs=5
for kappa in 1e2 1e4 1e6 1e8; do
    echo "--- Synthetic kappa=$kappa ---"
    $BENCH double "$OUTDIR" gen 5000 1000 $kappa 0.05 2.0 5 4
done

# photogrammetry2 (4472 x 936)
echo "--- photogrammetry2 ---"
$BENCH double "$OUTDIR" "$MTX" 2.0 5 4

echo "=== Done. CSVs in $OUTDIR ==="
