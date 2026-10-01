#!/bin/bash
# waits for batch 2, then runs batch 3 (at most 2 simulations at a time overall)
cd /home/user/wt/audit/benchmarks/physics_audit/material
while pgrep -f "run_batch.py truncated_cone-s2026-0000:E63" > /dev/null; do sleep 15; done
python3 run_batch.py truncated_cone-s2026-0000:base_im5 truncated_cone-s2026-0000:E63_im5 pyramid-s2026-0001:E63 pyramid-s2026-0001:kin_shutov pyramid-s2026-0001:vm pyramid-s2026-0001:hard_voce_coer pyramid-s2026-0001:hard_pow_iadicola pyramid-s2026-0001:ys115 pyramid-s2026-0001:ys85 pyramid-s2026-0001:mu02 truncated_cone-s2026-0000:kin_shutov_im5 truncated_cone-s2026-0000:vm_im5 truncated_cone-s2026-0000:hard_voce_coer_im5 pyramid-s2026-0001:Edeg > log_batch3.txt 2>&1
