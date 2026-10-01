#!/bin/bash
# run one BM3 case: ./run_case.sh cases/<name>   (1 thread; log in <case>/log.txt)
c=$1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE
start=$(date +%s)
/home/user/wt/audit/build/bin/sparlab_form --config $c/deck.json --output $c/output --threads 1 --strict-config --no-vtk --verbosity info > $c/log.txt 2>&1
echo "exit $? wall $(( $(date +%s) - start )) s" >> $c/log.txt
