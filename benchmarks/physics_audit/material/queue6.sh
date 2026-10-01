#!/bin/bash
# waits for batch 5, then runs batch 6 (Bauschinger and yield locus on the IM5 element)
cd /home/user/wt/audit/benchmarks/physics_audit/material
while pgrep -f "queue[345].sh" > /dev/null; do sleep 15; done
python3 run_batch.py truncated_cone-s2026-0000:hill_sb110_im5 truncated_cone-s2026-0000:kin_mild_im5 pyramid-s2026-0001:base_im5 pyramid-s2026-0001:kin_shutov_im5 truncated_cone-s2026-0000:combo_lit_im5 > log_batch6.txt 2>&1
