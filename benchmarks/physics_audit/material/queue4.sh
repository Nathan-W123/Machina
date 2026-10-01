#!/bin/bash
# waits for batches 2 and 3, then runs batch 4
cd /home/user/wt/audit/benchmarks/physics_audit/material
while pgrep -f "queue3.sh" > /dev/null || pgrep -f "run_batch.py truncated_cone-s2026-0000:E63 " > /dev/null; do sleep 15; done
python3 run_batch.py truncated_cone-s2026-0000:combo_lit pyramid-s2026-0001:combo_lit truncated_cone-s2026-0000:scaled090 truncated_cone-s2026-0000:t095 pyramid-s2026-0001:hill_sb110 > log_batch4.txt 2>&1
