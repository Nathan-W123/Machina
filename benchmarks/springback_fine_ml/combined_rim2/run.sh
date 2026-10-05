#!/bin/sh
# (Re)start combined_rim2 (fit-part compensations + 8 SparLab runs -> refit ->
# test compensations + 8 SparLab runs -> report) in the background, detached.
# Resumable: finished steps are skipped. 4 SparLab runs at a time, 1 thread each.
# Log: work/combined_rim2/run.log; completion marker: work/combined_rim2/done
cd "$(dirname "$0")/../../.." || exit 1
if pgrep -f "python/scripts/combined_rim[2]\.py all" >/dev/null; then
    echo "already running"; exit 0
fi
W=benchmarks/springback_fine_ml/work/combined_rim2
mkdir -p "$W"
rm -f "$W/done"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH="$PWD/python" \
    nohup setsid python3 python/scripts/combined_rim2.py all \
    >> "$W/run.log" 2>&1 < /dev/null &
echo "started pid $!"
