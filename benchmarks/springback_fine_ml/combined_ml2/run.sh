#!/bin/sh
# (Re)start the ML first shot + ML second round run; resumable.
# Log: work/combined_ml2/run.log; marker: work/combined_ml2/done
cd "$(dirname "$0")/../../.." || exit 1
if pgrep -f "python/scripts/combined_ml[2]\.py all" >/dev/null; then
    echo "already running"; exit 0
fi
W=benchmarks/springback_fine_ml/work/combined_ml2
mkdir -p "$W"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH="$PWD/python" \
    nohup setsid python3 python/scripts/combined_ml2.py all \
    >> "$W/run.log" 2>&1 < /dev/null &
echo "started pid $!"
