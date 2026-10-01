#!/bin/sh
# (Re)start the fine-mesh benchmark in the background; resumes from the run cache.
cd "$(dirname "$0")/../.." || exit 1
mkdir -p benchmarks/springback_fine/work
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$PWD/python" nohup setsid \
    python3 python/scripts/springback_fine_benchmark.py --out benchmarks/springback_fine \
    --work benchmarks/springback_fine/work --workers 4 "$@" \
    >> benchmarks/springback_fine/work/benchmark.log 2>&1 < /dev/null &
echo "started pid $!"
