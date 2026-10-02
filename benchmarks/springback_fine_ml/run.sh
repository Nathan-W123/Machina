#!/bin/sh
# (Re)start the ML stage in the background; resumes from the run cache, the
# data set, the model bundles and the saved compensations. Arguments are passed
# on (default: --stages 18 26). Log: benchmarks/springback_fine_ml/work/ml.log
cd "$(dirname "$0")/../.." || exit 1
mkdir -p benchmarks/springback_fine_ml/work
[ $# -eq 0 ] && set -- --stages 18 26
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$PWD/python" nohup setsid \
    python3 python/scripts/springback_fine_ml.py --out benchmarks/springback_fine_ml \
    --work benchmarks/springback_fine_ml/work --workers 4 --threads 4 "$@" \
    >> benchmarks/springback_fine_ml/work/ml.log 2>&1 < /dev/null &
echo "started pid $!"
