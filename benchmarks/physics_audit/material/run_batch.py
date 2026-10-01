#!/usr/bin/env python3
"""Run material-lens variants: `run_batch.py PART:VARIANT [PART:VARIANT ...] [--workers 2]`.

Variant specs come from variants.py (VARIANTS). Each run is single-threaded
(OMP/OpenBLAS 1 thread); at most `--workers` (<= 2) run at a time. Results
are appended to results.jsonl as each run finishes; a finished run
(runs/<part>/<variant>/COMPLETE) is measured again, not re-run.
"""
import argparse
import concurrent.futures as cf
import json
import os
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import matlens            # noqa: E402
from variants import VARIANTS   # noqa: E402


def job(pv):
    part, name = pv
    t0 = time.time()
    rec = matlens.run_variant(part, name, VARIANTS[name])
    rec["wall_s"] = time.time() - t0
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    jobs = [tuple(j.split(":")) for j in a.jobs]
    for p, v in jobs:
        if v not in VARIANTS or p not in matlens.PARTS:
            sys.exit(f"unknown {p}:{v}")
    # variants whose deck depends on another run (e_of_ep) run after it
    first = [j for j in jobs if not VARIANTS[j[1]].get("e_of_ep")]
    later = [j for j in jobs if VARIANTS[j[1]].get("e_of_ep")]
    for batch in (first, later):
        with cf.ProcessPoolExecutor(max_workers=min(2, a.workers)) as ex:
            for rec in ex.map(job, batch):
                print(time.strftime("%H:%M:%S"), json.dumps(
                    {k: rec.get(k) for k in ("part", "variant", "ok", "wall_s", "error")}),
                    (rec.get("dev_part") or {}).get("rms"), flush=True)


if __name__ == "__main__":
    main()
