#!/bin/sh
# (Re)start the combined compensation + verification on the 8 test parts,
# 4 at a time; each part resumes from its JSON in work/combined_rim.
# Logs: work/combined_rim/all8_<part>.log
cd "$(dirname "$0")/../../.." || exit 1
pgrep -f "combined_rim_all8[_]queue" >/dev/null && { echo "already running"; exit 0; }
W=benchmarks/springback_fine_ml/work/combined_rim
nohup setsid sh -c 'echo combined_rim_all8_queue >/dev/null; for p in dome-s2026-0000 dome-s2026-0001 elliptic_cone-s2026-0000 elliptic_cone-s2026-0001 pyramid-s2026-0000 pyramid-s2026-0001 truncated_cone-s2026-0000 truncated_cone-s2026-0001; do echo $p; done | xargs -P 4 -I{} sh -c "COMBINED_PART={} OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PRECOMP_ML_THREADS=1 PYTHONPATH=$PWD/python python3 python/scripts/combined_rim_ml.py one >> '"$W"'/all8_{}.log 2>&1"; echo ALL8_DONE >> '"$W"'/all8_done' >/dev/null 2>&1 < /dev/null &
echo started
