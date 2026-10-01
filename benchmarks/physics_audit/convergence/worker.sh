#!/bin/bash
# One worker of the study's run queue: takes the first case of queue.txt that
# is neither done, failed nor claimed (claim = atomic mkdir of cases/<name>.claim),
# runs it with one thread, repeats; exits when nothing is left.
# Start at most two workers (the machine is shared): nohup ./worker.sh A &
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE" || exit 1
export PYTHONPATH="$HERE/../../../python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while true; do
  picked=""
  while read -r name; do
    [ -z "$name" ] && continue
    case "$name" in \#*) continue;; esac
    [ -e "cases/$name/DONE" ] && continue
    [ -e "cases/$name/FAILED" ] && continue
    if mkdir "cases/$name.claim" 2>/dev/null; then picked=$name; break; fi
  done < queue.txt
  [ -z "$picked" ] && { echo "$(date -u +%FT%TZ) worker $1: queue empty"; exit 0; }
  echo "$(date -u +%FT%TZ) worker $1: start $picked (load $(cut -d' ' -f1-3 /proc/loadavg))"
  python3 run_case.py "$picked"
  echo "$(date -u +%FT%TZ) worker $1: end $picked rc=$?"
  rmdir "cases/$picked.claim" 2>/dev/null
done
