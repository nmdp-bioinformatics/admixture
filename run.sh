#!/bin/sh
# Run NMF decomposition for K=2..12 and generate the error curve.
# Usage: ./run.sh [freqs_dir]
#   freqs_dir defaults to "freqs/"

FREQS_DIR="${1:-freqs/}"

for k in $(seq 2 12); do
    echo "running NMF for k = ${k}"
    python nmf.py "${FREQS_DIR}" -k ${k} --min-value 1e-7 --out-dir nmf_out_k${k}
done

echo "generating error curve..."
python error_curve.py "${FREQS_DIR}"
echo "done"
