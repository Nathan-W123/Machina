#!/bin/sh
# (Re)start combined_rim3 (v2 model + std weight 0.25 on the 8 test parts) in
# the background; resumable. Log: work/combined_rim3/run.log; marker: work/combined_rim3/done
cd "$(dirname "$0")/../../.." || exit 1
if pgrep -f "python/scripts/combined_rim[3]\.py all" >/dev/null; then
    echo "already running"; exit 0
fi
W=benchmarks/springback_fine_ml/work/combined_rim3
mkdir -p "$W"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH="$PWD/python" \
    nohup setsid python3 python/scripts/combined_rim3.py all \
    >> "$W/run.log" 2>&1 < /dev/null &
echo "started pid $!"
