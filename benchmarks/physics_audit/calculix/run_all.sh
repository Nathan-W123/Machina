#!/bin/bash
# Every run of this lens, in the order used (each line is one simulation pair;
# run at most two at a time, 1 thread each). Runs whose outputs exist are skipped
# by case_punch.py (SparLab) / rerun if CalculiX did not finish.
cd "$(dirname "$0")"
export PYTHONPATH=../../../python OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# test 3: homogeneous large strain (CCX_INC=400 TAG=_inc400 for the fine-increment CalculiX run)
python3 case_homogeneous.py
CCX_INC=400 TAG=_inc400 python3 case_homogeneous.py
# tests 1+2, strip, benchmark mesh (first auditor; CalculiX unload by retraction)
python3 case_punch.py strip_h2_std_mu0.1 kind=strip elem=std ccx_unload=retract
python3 case_punch.py strip_h2_nobbar_mu0.1 kind=strip elem=std_nobbar ccx_unload=retract
python3 case_punch.py strip_h2_im_mu0.1 kind=strip elem=im ccx_unload=retract
python3 case_punch.py strip_h2_im_mu0.1_remove kind=strip h=2 nz=2 elem=im
# strip mesh / element / increment study
python3 case_punch.py strip_h1_std_nz4 kind=strip h=1 nz=4 elem=std
python3 case_punch.py strip_h2_std_nz2_tp5 kind=strip h=2 nz=2 elem=std tp=5
python3 case_punch.py strip_h2_im_nz4 kind=strip h=2 nz=4 elem=im
python3 case_punch.py strip_h1_im_nz2 kind=strip h=1 nz=2 elem=im
python3 case_punch.py strip_h2_std_nz2_trav1 kind=strip h=2 nz=2 elem=std travel=1.0
python3 case_punch.py strip_h2_std_nz4_tp5 kind=strip h=2 nz=4 elem=std tp=5
python3 case_punch.py strip_h0.5_std_nz4 kind=strip h=0.5 nz=4 elem=std
python3 case_punch.py strip_h0.5_im_nz8 kind=strip h=0.5 nz=8 elem=im
python3 case_punch.py strip_h1_im_nz4 kind=strip h=1 nz=4 elem=im
python3 case_punch.py strip_h0.5_im_nz4 kind=strip h=0.5 nz=4 elem=im
# C3D20R references (the node-to-surface ones, ref_h0.5_nz4 / ref_h0.25_nz8, diverged)
python3 strip_reference.py ccx ref_h0.5_nz4_s2s h=0.5 nz=4 contact=s2s
python3 strip_reference.py ccx ref_h0.25_nz8_s2s_retract h=0.25 nz=8 contact=s2s unload=retract
python3 strip_reference.py compare ref_h0.5_nz4_s2s
python3 strip_reference.py compare ref_h0.25_nz8_s2s_retract
python3 strip_regions.py ref_h0.5_nz4_s2s
python3 strip_regions.py ref_h0.25_nz8_s2s_retract
# 3D plate, plunge 3 + drag 3 mm
python3 case_punch.py plate_h2_im_mu0.1_drag3 kind=plate elem=im depth=3 drag=3
python3 case_punch.py plate_h2_nobbar_mu0.1_drag3 kind=plate elem=std_nobbar depth=3 drag=3
python3 case_punch.py plate_h2_std_mu0.1_drag3 kind=plate elem=std depth=3 drag=3 only=sparlab
python3 case_punch.py plate_h2_std_mu0.1_drag3_trav1 kind=plate elem=std depth=3 drag=3 travel=1.0 only=sparlab
python3 case_punch.py plate_h2_im_mu0.1_drag3_trav1 kind=plate elem=im depth=3 drag=3 travel=1.0 only=sparlab
python3 case_punch.py plate_h2_std_mu0.1_drag3_tp5 kind=plate elem=std depth=3 drag=3 tp=5 only=sparlab
python3 case_punch.py plate_h1_std_mu0.1_drag3 kind=plate h=1 nz=2 elem=std depth=3 drag=3 only=sparlab
python3 case_punch.py plate_h1_im_mu0.1_drag3 kind=plate h=1 nz=2 elem=im depth=3 drag=3 only=sparlab
python3 case_punch.py plate_h1_std_nz4_mu0.1_drag3 kind=plate h=1 nz=4 elem=std depth=3 drag=3 only=sparlab
python3 plate_compare.py plate_h2_im_mu0.1_drag3
python3 plate_compare.py plate_h2_nobbar_mu0.1_drag3
python3 plate_mesh.py
