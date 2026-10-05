#!/bin/sh
# (Re)start the combined rim-pass + transfer-model run in the background. It
# resumes from the run cache (work/runs), the rim-pass data set (work/data_rim)
# and the per-step JSON in work/combined_rim. Log: work/combined_rim.log
cd "$(dirname "$0")/../../.." || exit 1
if pgrep -f "^python3 python/scripts/combined_rim_ml[.]py all" >/dev/null; then
    echo "already running"; exit 0
fi
mkdir -p benchmarks/springback_fine_ml/work
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH="$PWD/python" \
    nohup setsid python3 python/scripts/combined_rim_ml.py all "$@" \
    >> benchmarks/springback_fine_ml/work/combined_rim.log 2>&1 < /dev/null &
echo "started pid $!"
