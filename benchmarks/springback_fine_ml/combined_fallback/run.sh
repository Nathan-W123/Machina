#!/bin/sh
# (Re)start the ML-first-shot + engine-fallback run; resumable.
# Log: work/combined_fallback/run.log; marker: work/combined_fallback/done
cd "$(dirname "$0")/../../.." || exit 1
if pgrep -f "python/scripts/combined_fallbac[k]\.py all" >/dev/null; then
    echo "already running"; exit 0
fi
W=benchmarks/springback_fine_ml/work/combined_fallback
mkdir -p "$W"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH="$PWD/python" \
    nohup setsid python3 python/scripts/combined_fallback.py all \
    >> "$W/run.log" 2>&1 < /dev/null &
echo "started pid $!"
