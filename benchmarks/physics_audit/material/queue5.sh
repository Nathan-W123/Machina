#!/bin/bash
# waits for batch 4, then runs batch 5 (a third test part)
cd /home/user/wt/audit/benchmarks/physics_audit/material
while pgrep -f "queue[34].sh" > /dev/null; do sleep 15; done
python3 run_batch.py elliptic_cone-s2026-0001:base elliptic_cone-s2026-0001:hill_sb110 elliptic_cone-s2026-0001:kin_shutov elliptic_cone-s2026-0001:hard_voce_coer elliptic_cone-s2026-0001:vm elliptic_cone-s2026-0001:E63 elliptic_cone-s2026-0001:combo_lit > log_batch5.txt 2>&1
