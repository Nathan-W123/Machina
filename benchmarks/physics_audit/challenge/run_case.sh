#!/bin/bash
# usage: run_case.sh <case>   (single thread, writes cases/<case>/{output,log,DONE|FAILED})
c=/home/user/wt/audit/benchmarks/physics_audit/challenge/cases/$1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
s=$(date +%s)
/home/user/wt/audit/build/bin/sparlab_form --config $c/deck/deck.json --output $c/output --threads 1 --strict-config > $c/sparlab_form.log 2>&1
rc=$?
echo "rc=$rc wall_s=$(( $(date +%s) - s ))" > $c/$([ $rc = 0 ] && echo DONE || echo FAILED)
