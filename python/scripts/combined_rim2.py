#!/usr/bin/env python3
"""Combined compensation v2 (SparLab simulation): the DSIF rim pass + transfer
MLP + surrogate DA + optimiser of `combined_rim_ml.py`, with two fixes applied
TOGETHER (no ablation):

  fix 2 - a stronger std penalty in the optimiser: `surrogate_optimize(...,
     std_weight=STD_WEIGHT)` = 1.0 (v1: the default 0.25), chosen a priori
     (the std then counts as much as the residual), not tuned on any test part;
  fix 1 - the transfer correction also learns from the shapes the optimiser
     proposes: for the 8 role=fit TRAIN parts of combined_rim's selection.csv,
     the combined command (v1 transfer model + DA + optimiser with the new
     penalty, same masks / setup as `combined_rim_ml.combined`) is simulated
     with the DSIF + rim-pass setup and added to <work>/data_rim as a
     "compensated" sample (commanded = that command, formed = simulated).
     The correction is refitted on the 8 uncompensated + 8 compensated fit
     samples (the 4 holdout parts stay out); the correction form ("linear",
     as v1, or the existing regularised "gbm") is picked by leave-one-part-out
     error on those TRAIN samples only (both samples of a part leave together;
     "gbm" only if more than 5 % below "linear").

Then the combined compensation of the 8 TEST parts with the v2 model and the
new penalty, one verifying SparLab run each, and a report.

Steps (mode `all`; 4 SparLab runs at a time, 1 thread each; every step
persisted, a re-run skips finished work - a recorded failure is not retried):
  1. fit-part commands  <CW2>/fit/<part>-combined.{json,npz} (v1 model)
  2. fit-part runs      <CW2>/fit/<part>-sim.json (+ -formed.npz, + the data
                        set sample <part>-rimc1); the run cache is <work>/runs
  3. refit              <CW2>/refit.json, <CW2>/model/transfer_mlp_v2
  4. test commands      <CW2>/test/<part>-combined.{json,npz} (v2 model)
  5. test runs          <CW2>/test/<part>-verify.json (+ -formed.npz)
  6. report             benchmarks/springback_fine_ml/combined_rim2/
                        {results.csv, README.md, ...}; marker <CW2>/done
The ML computations (DA + optimiser, ~2-5 min each) run one at a time (a
lock), so the optimiser's wall-time cap means the same in every slot.

usage (worktree root; OPENBLAS_NUM_THREADS=1 PYTHONPATH=python):
    python3 python/scripts/combined_rim2.py all
    python3 python/scripts/combined_rim2.py status
    python3 python/scripts/combined_rim2.py report
    python3 python/scripts/combined_rim2.py dry --part dome-s2026-0010
        (compute one fit-part command into <CW2>/dry/, no SparLab run)
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional

import _bootstrap  # noqa: F401

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import springback_benchmark as sb  # noqa: E402
import smoke_optimizer_edge as se  # noqa: E402
import support_validation as sv  # noqa: E402
import combined_rim_ml as cr  # noqa: E402

MM = 1e-3
LABEL = "SparLab simulation"
OUT_ML = se.OUT_ML
WORK = se.WORK
OUT = OUT_ML / "combined_rim2"
CW2 = WORK / "combined_rim2"
CW1 = cr.CW
OUT1 = cr.OUT
COLS = se.COLS
#: SparLab runs at a time and threads per run. One run on 4 threads finishes
#: well inside the container's ~1 h restart interval (4 runs on 1 thread each
#: take ~85 min and were lost to restarts). Threads are an execution field:
#: the deck hash and the run cache do not change.
SLOTS = int(os.environ.get("COMBINED_SLOTS", "1"))
SIM_THREADS = int(os.environ.get("COMBINED_SIM_THREADS", "4"))
#: fix 2: the optimiser's std weight (v1: the default 0.25); chosen a priori
STD_WEIGHT = 1.0
STD_WEIGHT_V1 = 0.25
#: the correction forms tried for v2 (picked by LOPO on TRAIN data)
CANDIDATES = ("linear", "gbm")
#: a richer form than "linear" is chosen only if its LOPO error is lower by this
RICHER_MARGIN = 0.05
#: the 8 held-out TEST parts of stage_n26 (combined_rim/run_all8.sh)
TEST = ["dome-s2026-0000", "dome-s2026-0001", "elliptic_cone-s2026-0000",
        "elliptic_cone-s2026-0001", "pyramid-s2026-0000", "pyramid-s2026-0001",
        "truncated_cone-s2026-0000", "truncated_cone-s2026-0001"]

ML_LOCK = threading.Lock()


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def comp_sample_id(pid: str) -> str:
    return f"{pid}-rimc1"


def selection() -> pd.DataFrame:
    f = CW1 / "selection.csv"
    if not f.exists():
        raise SystemExit(f"{f} missing: run combined_rim_ml.py all first")
    sel = pd.read_csv(f)
    assert not set(sel.point_id) & set(TEST), "a TEST part in the training selection"
    assert (sel.split == "train").all()
    return sel


# ---------------------------------------------------------------------------
# the models
# ---------------------------------------------------------------------------
def model_v1(ctx: cr.Ctx, sel: pd.DataFrame):
    """combined_rim's transfer model (loaded from its bundle; refitted only
    if the bundle is missing - the same data, the same fit)."""
    return cr.transfer_model(ctx, sel)


def fit_samples(ctx: cr.Ctx, sel: pd.DataFrame) -> List[Any]:
    ds = ctx.dataset()
    out = []
    for p in sel[sel.role == "fit"].point_id:
        for sid in (cr.sample_id(p), comp_sample_id(p)):
            if ds.has(sid):
                out.append(ds.load(sid))
    return out


def lopo(ctx: cr.Ctx, samples: List[Any], correction: str,
         train_kinds=("uncompensated", "compensated")) -> Dict[str, Any]:
    """Leave-one-PART-out on TRAIN rim samples: for each part, a TransferModel
    (`correction`, with its own inner validation, as `transfer_surrogate`)
    fitted on the OTHER parts' samples of `train_kinds`, scored on every
    sample of the left-out part (dz RMS over the table's nodes). Returns
    per-sample rows and means by kind."""
    from precomp.ml.dataset import build_table
    from precomp.ml.models import TransferModel
    from precomp.ml.surrogate import model_prior

    base = ctx.mlp()
    tab = build_table(samples, base.features,
                      points_per_sample=ctx.cfg["surrogate"]["points_per_sample"],
                      mask=base.domain, rng=ctx.cfg["seed"], prior=model_prior(base.model))
    kind = {s.sample_id: s.kind for s in samples}
    rk = np.array([kind[g] for g in tab.groups])
    parts = sorted(set(tab.part_id))
    mu_b, _ = base.model.predict(tab.X)
    rows = []
    for g in parts:
        tr = (tab.part_id != g) & np.isin(rk, train_kinds)
        te = tab.part_id == g
        m = TransferModel(base.model, correction, seed=ctx.cfg["seed"])
        m.fit(tab.X[tr], tab.y[tr], tab.part_id[tr], feature_names=tab.feature_names)
        mu, _ = m.predict(tab.X[te])
        y, gs, b = tab.y[te], tab.groups[te], mu_b[te]
        for sid in sorted(set(gs)):
            k = gs == sid
            rows.append({"part": g, "sample_id": sid, "kind": kind[sid],
                         "accepted": bool(m.accepted),
                         "base_rms_mm": float(np.sqrt(np.mean((b[k] - y[k]) ** 2))) / MM,
                         "corrected_rms_mm": float(np.sqrt(np.mean((mu[k] - y[k]) ** 2))) / MM})
    df = pd.DataFrame(rows)
    out: Dict[str, Any] = {"correction": correction, "train_kinds": list(train_kinds),
                           "n_parts": len(parts), "per_sample": rows}
    for kd in ("all", "uncompensated", "compensated"):
        d = df if kd == "all" else df[df.kind == kd]
        if len(d):
            out[f"mean_{kd}"] = {"n": int(len(d)), "base_rms_mm": float(d.base_rms_mm.mean()),
                                 "corrected_rms_mm": float(d.corrected_rms_mm.mean())}
    return out


def refit(ctx: cr.Ctx, sel: pd.DataFrame) -> Dict[str, Any]:
    """Pick the correction form by LOPO on the 16 TRAIN fit samples, fit the
    v2 transfer model on them, and compare old vs new (TRAIN data only)."""
    from precomp.ml import evaluate_surrogate, load_model, save_model
    from precomp.ml.surrogate import transfer_surrogate

    f = CW2 / "refit.json"
    fit = fit_samples(ctx, sel)
    ids = sorted(s.sample_id for s in fit)
    bundle = CW2 / "model" / "transfer_mlp_v2"
    if f.exists() and (bundle / "manifest.json").exists():
        info = json.loads(f.read_text())
        if info.get("samples") == ids:
            return info
    n_comp = sum(s.kind == "compensated" for s in fit)
    log(f"refit: {len(fit)} fit samples ({n_comp} compensated)")
    if n_comp < 2:
        raise RuntimeError(f"only {n_comp} compensated fit samples: cannot refit")
    t0 = time.perf_counter()
    cands = {c: lopo(ctx, fit, c) for c in CANDIDATES}
    # the old protocol (v1: linear on the uncompensated samples only), scored
    # on the same left-out parts and samples
    old = lopo(ctx, fit, "linear", train_kinds=("uncompensated",))
    # parsimony, fixed a priori: the richer form only if it beats "linear" by
    # more than RICHER_MARGIN (relative)
    e = {c: cands[c]["mean_all"]["corrected_rms_mm"] for c in CANDIDATES}
    best = "linear"
    for c in CANDIDATES:
        if c != "linear" and e[c] < (1.0 - RICHER_MARGIN) * e[best]:
            best = c
    log("refit LOPO (mean over left-out samples, mm): old "
        f"{old['mean_all']['corrected_rms_mm']:.4f}; "
        + ", ".join(f"{c} {cands[c]['mean_all']['corrected_rms_mm']:.4f}" for c in CANDIDATES)
        + f" -> {best}")
    sur = transfer_surrogate(ctx.mlp(), fit, correction=best,
                             points_per_sample=ctx.cfg["surrogate"]["points_per_sample"],
                             seed=ctx.cfg["seed"])
    save_model(sur, bundle, {"created_at": now_iso(),
                             "notes": "stage_n26 MLP (backing plate, penalty 3) + "
                                      f"{best} TransferModel correction on 8 uncompensated + "
                                      f"{n_comp} compensated DSIF + rim-pass runs (penalty "
                                      "10) of the 8 fit parts - SparLab simulation"},
               overwrite=True)
    # the 4 holdout parts (uncompensated runs; never in either fit)
    hold = cr.rim_samples(ctx, list(sel[sel.role == "holdout"].point_id))
    held = {}
    if hold:
        for name, m in (("base", ctx.mlp()), ("v1", model_v1(ctx, sel)), ("v2", sur)):
            ev = evaluate_surrogate(m, hold, split="held_out", assess=False, check_unseen=False)
            s = ev.summary()["held_out"]
            held[name] = {"n_samples": s["n_samples"], "err_rms_mean_mm": 1e3 * s["err_rms_mean_m"]}
    # v1 (as used, fitted on the 8 uncompensated) on the 8 compensated runs:
    # not in its fit (their parts were)
    comp = [s for s in fit if s.kind == "compensated"]
    ev = evaluate_surrogate(model_v1(ctx, sel), comp, split="held_out", assess=False,
                            check_unseen=False)
    s = ev.summary()["held_out"]
    v1_on_comp = {"n_samples": s["n_samples"], "err_rms_mean_mm": 1e3 * s["err_rms_mean_m"],
                  "dz_rms_mean_mm": 1e3 * s["dz_rms_mean_m"]}
    v = sur.training["transfer"]
    info = {"samples": ids, "n_compensated": n_comp, "candidates": cands, "old_protocol": old,
            "chosen": best, "rule": "mean LOPO dz RMS over the left-out fit samples (both samples "
            "of a part leave together; TRAIN data only); 'gbm' only if more than "
            f"{100 * RICHER_MARGIN:.0f} % below 'linear'",
            "v2_accepted": v["accepted"], "v2_validation": v["validation"],
            "holdout_uncompensated": held, "v1_on_compensated_fit_runs": v1_on_comp,
            "time_s": time.perf_counter() - t0}
    sb.jdump(f, info)
    return info


def model_v2(ctx: cr.Ctx, sel: pd.DataFrame):
    from precomp.ml import load_model

    refit(ctx, sel)
    return load_model(CW2 / "model" / "transfer_mlp_v2")


# ---------------------------------------------------------------------------
# a combined command (as combined_rim_ml.combined, with the new std weight)
# ---------------------------------------------------------------------------
def command(ctx: cr.Ctx, sur, pid: str, d: Path, max_time: float, model_name: str):
    from precomp._util import PrecompError
    from precomp.fea.support import command_upper_bound
    from precomp.geometry.heightmap import HeightMap
    from precomp.ml import surrogate_compensate
    from precomp.ml.compensate import surrogate_optimize

    fpred, fcmd = d / f"{pid}-combined.json", d / f"{pid}-combined.npz"
    if fpred.exists() and fcmd.exists():
        return json.loads(fpred.read_text()), HeightMap.load(fcmd)
    d.mkdir(parents=True, exist_ok=True)
    target = ctx.points[pid].target()
    hold, adjust, band, strip = cr.masks(ctx, target)
    sc, da = ctx.cfg["surrogate"], ctx.cfg["da"]
    kw = dict(max_wall_angle_deg=da["max_wall_angle_deg"],
              upper_bound=command_upper_bound(ctx.setup, target),
              hold_mask=hold, adjust_mask=adjust, tool_radius=ctx.setup.tool_radius)
    dkw = dict(iterations=sc["iterations"], stagnation=sc["stagnation"], alpha=da["alpha"],
               direction=da["direction"], smoothing=da["smoothing"])
    with ML_LOCK:
        log(f"command {pid} ({model_name}): start")
        override = False
        t0 = time.perf_counter()
        try:
            dres = surrogate_compensate(target, ctx.setup, sur, **dkw, **kw)
        except PrecompError as exc:
            if "envelope" not in str(exc):
                raise
            override = True
            dres = surrogate_compensate(target, ctx.setup, sur, allow_out_of_envelope=True,
                                        **dkw, **kw)
        t_da = time.perf_counter() - t0
        dres.compensated.save(d / f"{pid}-combined-da.npz")
        t0 = time.perf_counter()
        ores = surrogate_optimize(target, ctx.setup, sur, start=dres, max_time_s=max_time,
                                  std_weight=STD_WEIGHT, allow_out_of_envelope=override, **kw)
        t_opt = time.perf_counter() - t0
        comp = ores.compensated
        ores.save(d / f"{pid}-combined-optimizer")
        mu, _ = sur.predict_deviation(target, ctx.setup, target)

    def pred(res):
        dv = sv.deviation(res.predicted, target, ctx.cfg["upper_band_depth"])
        return {k: dv[k] for k in COLS}

    pm = target.mask & (target.z < -1e-6)
    std = None if ores.std is None else float(np.sqrt(np.mean(np.asarray(ores.std)[pm] ** 2)))
    info = {"point_id": pid, "model": model_name, "std_weight": STD_WEIGHT,
            "data_source": f"transfer-model prediction ({model_name}); SparLab simulation data",
            "setup": "smoke_optimizer_edge.edge_setup(True) (dsif + rim pass, penalty 10)",
            "transfer_accepted": sur.training["transfer"]["accepted"],
            "transfer_correction": sur.training["transfer"]["correction"],
            "envelope_override": override,
            "masks": {"band_nodes": int(band.sum()), "strip_nodes": int(strip.sum()),
                      "adjust_nodes": int(adjust.sum()), "hold_nodes": int(hold.sum())},
            "target_predicted": {k: v for k, v in sv.deviation(
                target.with_z(target.z + mu), target, ctx.cfg["upper_band_depth"]).items()
                if k in COLS},
            "da": {**pred(dres), "stopped": dres.stopped, "best_iteration": dres.best_iteration,
                   "time_s": t_da, "commanded_depth_mm": dres.compensated.depth / MM},
            "optimizer": {**pred(ores), "stopped": ores.stopped, "time_s": t_opt,
                          "evaluations": len(ores.history), "in_envelope": ores.in_envelope,
                          "stop_reason": ores.history[-1].get("stop_reason")
                          if ores.history else None,
                          "std_rms_mm": None if std is None else std / MM,
                          "commanded_depth_mm": comp.depth / MM,
                          "max_abs_change_from_da_mm":
                              float(np.abs(comp.z - dres.compensated.z).max()) / MM}}
    comp.save(fcmd)
    sb.jdump(fpred, info)
    log(f"command {pid} ({model_name}): predicted DA {info['da']['vert_rms_mm']:.4f}, "
        f"optimiser {info['optimizer']['vert_rms_mm']:.4f} mm ({info['optimizer']['stop_reason']}"
        f", {t_opt:.0f} s, change from DA {info['optimizer']['max_abs_change_from_da_mm']:.3f} mm)")
    return info, comp


# ---------------------------------------------------------------------------
# a SparLab run of a command (DSIF + rim pass, penalty 10)
# ---------------------------------------------------------------------------
def add_comp_sample(ctx: cr.Ctx, pid: str, comp, formed, rec: Dict[str, Any]) -> None:
    from precomp.ml import Sample

    p = ctx.points[pid]
    target = p.target()
    s = Sample(comp_sample_id(pid), comp.copy(), formed, ctx.setup.to_dict(), "sim",
               cr.FIDELITY, "compensated", p.part.to_dict(), pid, target, None,
               {"created_at": now_iso(), "deck_hash": rec.get("deck_hash"),
                "sparlab_version": rec.get("sparlab_version"),
                "runtime_s": rec.get("runtime_s"),
                "note": "combined command (v1 transfer MLP + DA + optimiser, std_weight "
                        f"{STD_WEIGHT}), DSIF + rim pass, penalty 10 (combined_rim2.py)"})
    with ctx.lock:
        ds = ctx.dataset()
        if not ds.has(s.sample_id):
            ds.append(s)
            log(f"data set: added {s.sample_id}")


def sim(ctx: cr.Ctx, pid: str, comp, f: Path, as_sample: bool) -> Dict[str, Any]:
    from precomp.fea.runner import FormingError, simulate
    from precomp.geometry.heightmap import HeightMap

    ff = f.with_name(f.stem + "-formed.npz")
    if f.exists():
        rec = json.loads(f.read_text())
        if as_sample and rec.get("ok") and ff.exists():
            add_comp_sample(ctx, pid, comp, HeightMap.load(ff), rec)   # (if missing)
        return rec
    target = ctx.points[pid].target()
    rec: Dict[str, Any] = {"point_id": pid, "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    log(f"sim {pid}: start (dsif + rim pass, penalty 10) -> {f.name}")
    t0 = time.perf_counter()
    try:
        os.environ["OPENBLAS_NUM_THREADS"] = str(SIM_THREADS)  # the subprocess's BLAS
        res = simulate(dataclasses.replace(ctx.setup, threads=SIM_THREADS), comp,
                       WORK / "runs", target=target)
        formed = res.formed_surface(-1, grid=target.grid)
        prov = res.provenance
        rec.update(ok=True, deck_hash=prov.get("key"), cache_hit=prov.get("cache_hit"),
                   sparlab_version=prov.get("sparlab_version"))
        rec.update(sv.physics(res))
        rec.update(sv.deviation(formed, target, ctx.cfg["upper_band_depth"]))
        formed.save(ff)
        if as_sample:
            add_comp_sample(ctx, pid, comp, formed, rec)
    except FormingError as exc:
        r = {k: v for k, v in (exc.record or {}).items()
             if isinstance(v, (str, int, float, bool))}
        rec.update(ok=False, error=str(exc)[:2000], termination=r.get("termination")
                   or r.get("reason"), record=r)
    rec["wall_s"] = time.perf_counter() - t0
    sb.jdump(f, rec)
    log(f"sim {pid}: {'ok' if rec['ok'] else 'FAILED ' + str(rec.get('termination'))} "
        f"[{rec['wall_s']:.0f} s] vert RMS {rec.get('vert_rms_mm', float('nan')):.3f} mm")
    return rec


# ---------------------------------------------------------------------------
# the scheduler
# ---------------------------------------------------------------------------
def cmd_all(args) -> None:
    t_start = time.time()
    sel = selection()
    fit = list(sel[sel.role == "fit"].sort_values("data_runtime_p3_s",
                                                  ascending=False).point_id)
    log(f"combined_rim2: std_weight {STD_WEIGHT}; fit parts {fit}; test parts {TEST}")
    ctx = cr.Ctx()
    CW2.mkdir(parents=True, exist_ok=True)
    jobs = ([{"name": f"fit {p}", "kind": "fit", "pid": p} for p in fit]
            + [{"name": f"test {p}", "kind": "test", "pid": p} for p in TEST])
    fit_names = {f"fit {p}" for p in fit}
    done: set = set()
    running: set = set()
    cv = threading.Condition()
    models: Dict[str, Any] = {}
    mlock = threading.Lock()

    def get_model(name):
        with mlock:
            if name not in models:
                models[name] = model_v1(ctx, sel) if name == "v1" else model_v2(ctx, sel)
            return models[name]

    def runnable(j):
        return j["kind"] == "fit" or fit_names <= done

    def worker(i: int) -> None:
        while True:
            with cv:
                while True:
                    todo = [j for j in jobs if j["name"] not in done | running]
                    if not todo:
                        return
                    nxt = next((j for j in todo if runnable(j)), None)
                    if nxt is not None:
                        running.add(nxt["name"])
                        break
                    cv.wait(30)
            pid = nxt["pid"]
            try:
                if nxt["kind"] == "fit":
                    d = CW2 / "fit"
                    sim_f = d / f"{pid}-sim.json"
                    _, comp = command(ctx, None if (d / f"{pid}-combined.npz").exists()
                                      and (d / f"{pid}-combined.json").exists()
                                      else get_model("v1"), pid, d, args.max_time,
                                      "v1: transfer MLP, linear, 8 uncompensated fit runs")
                    sim(ctx, pid, comp, sim_f, as_sample=True)
                else:
                    d = CW2 / "test"
                    cached = (d / f"{pid}-combined.npz").exists() and \
                        (d / f"{pid}-combined.json").exists()
                    sur = None if cached else get_model("v2")
                    _, comp = command(ctx, sur, pid, d, args.max_time,
                                      "v2: transfer MLP refitted on 8 uncompensated + 8 "
                                      "compensated fit runs")
                    sim(ctx, pid, comp, d / f"{pid}-verify.json", as_sample=False)
            except Exception as exc:  # logged; the others go on
                log(f"{nxt['name']}: ERROR {type(exc).__name__}: {exc}")
                traceback.print_exc()
            with cv:
                running.discard(nxt["name"])
                done.add(nxt["name"])
                cv.notify_all()

    ths = [threading.Thread(target=worker, args=(i,), daemon=True) for i in range(SLOTS)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    try:
        refit(ctx, sel)     # (already done unless every test command was cached)
    except Exception as exc:
        log(f"refit: ERROR {exc}")
    cmd_report(args)
    n_ok = sum(json.loads(f.read_text()).get("ok", False)
               for f in (CW2 / "test").glob("*-verify.json"))
    (CW2 / "done").write_text(json.dumps({"finished": now_iso(), "test_runs_ok": n_ok,
                                          "this_invocation_s": time.time() - t_start}) + "\n")
    log(f"done ({n_ok}/8 test runs ok)")


def cmd_dry(args) -> None:
    """One fit-part command into <CW2>/dry (no SparLab run): a smoke test."""
    sel = selection()
    ctx = cr.Ctx()
    pid = args.part or sel[sel.role == "fit"].point_id.iloc[0]
    info, comp = command(ctx, model_v1(ctx, sel), pid, CW2 / "dry", args.max_time, "v1 (dry)")
    print(json.dumps({k: info[k] for k in ("da", "optimizer", "target_predicted")}, indent=1,
                     default=str))


def cmd_status(args) -> None:
    sel = selection()
    for kind, pids, nm in (("fit", list(sel[sel.role == "fit"].point_id), "sim"),
                           ("test", TEST, "verify")):
        for p in pids:
            c = (CW2 / kind / f"{p}-combined.json").exists()
            f = CW2 / kind / f"{p}-{nm}.json"
            x = json.loads(f.read_text()) if f.exists() else {}
            print(f"{kind:4s} {p:28s} command {'yes' if c else 'no '}  sim "
                  + ((("ok" if x.get("ok") else "FAILED")
                      + f" {x.get('wall_s', 0):.0f} s vert {x.get('vert_rms_mm', np.nan):.3f}")
                     if x else "pending"))
    print("refit.json:", "present" if (CW2 / "refit.json").exists() else "pending")
    print("done:", "yes" if (CW2 / "done").exists() else "no")


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _j(f: Path) -> Dict[str, Any]:
    return json.loads(f.read_text()) if f.exists() else {}


def cmd_report(args) -> None:
    sel = selection()
    a8 = pd.read_csv(OUT1 / "all8_results.csv")
    meth = {"uncompensated": "plate_unc", "ML-MLP": "mlp_plate",
            "combined (rim pass + transfer MLP + optimiser)": "v1"}
    rows = []
    for p in TEST:
        r: Dict[str, Any] = {"point_id": p}
        for m, key in meth.items():
            x = a8[(a8.point_id == p) & (a8.method == m)]
            for c in ("vert_rms_mm", "upper_band_rms_mm", "deep_rms_mm"):
                r[f"{key}_{c}"] = float(x[c].iloc[0]) if len(x) else np.nan
        p1 = _j(CW1 / f"{p}-combined.json").get("optimizer", {})
        p2i = _j(CW2 / "test" / f"{p}-combined.json")
        p2 = p2i.get("optimizer", {})
        v2 = _j(CW2 / "test" / f"{p}-verify.json")
        for c in ("vert_rms_mm", "upper_band_rms_mm", "deep_rms_mm"):
            r[f"v1_pred_{c}"] = p1.get(c, np.nan)
            r[f"v2_pred_{c}"] = p2.get(c, np.nan)
            r[f"v2_{c}"] = v2.get(c, np.nan) if v2.get("ok") else np.nan
        r["v2_status"] = ("ok" if v2.get("ok") else "FAILED: " + str(v2.get("termination")
                          or v2.get("error"))[:200]) if v2 else "pending"
        r["v2_deck_hash"] = v2.get("deck_hash")
        r["v2_envelope_override"] = p2i.get("envelope_override")
        r["v2_optimizer_stop"] = p2.get("stop_reason")
        r["v2_change_from_da_mm"] = p2.get("max_abs_change_from_da_mm")
        rows.append(r)
    df = pd.DataFrame(rows)
    df["data_source"] = LABEL
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "results.csv", index=False, float_format="%.4f")
    # fit-part runs
    fr = []
    for p in sel[sel.role == "fit"].point_id:
        ci = _j(CW2 / "fit" / f"{p}-combined.json")
        x = _j(CW2 / "fit" / f"{p}-sim.json")
        u = _j(CW1 / "runs" / f"{p}.json")
        fr.append({"point_id": p, "rim_unc_vert_rms_mm": u.get("vert_rms_mm"),
                   "pred_vert_rms_mm": ci.get("optimizer", {}).get("vert_rms_mm"),
                   "sim_vert_rms_mm": x.get("vert_rms_mm") if x.get("ok") else None,
                   "sim_upper_band_rms_mm": x.get("upper_band_rms_mm") if x.get("ok") else None,
                   "sim_deep_rms_mm": x.get("deep_rms_mm") if x.get("ok") else None,
                   "status": ("ok" if x.get("ok") else "FAILED") if x else "pending",
                   "runtime_s": x.get("runtime_s"), "deck_hash": x.get("deck_hash")})
    frd = pd.DataFrame(fr)
    frd["data_source"] = LABEL
    frd.to_csv(OUT / "fit_runs.csv", index=False, float_format="%.4f")
    rf = _j(CW2 / "refit.json")
    if rf:
        sb.jdump(OUT / "refit.json", rf)

    def f3(v):
        try:
            v = float(v)
        except (TypeError, ValueError):
            return "-"
        return "-" if not np.isfinite(v) else f"{v:.3f}"

    n2 = int(df.v2_vert_rms_mm.notna().sum())
    L = ["# Combined compensation v2 on the 8 test parts: retrained transfer correction + "
         "stronger uncertainty penalty", "",
         f"**{LABEL}** (sparlab_form), penalty 10, DSIF + rim pass "
         "(`smoke_optimizer_edge.edge_setup(True)`). Predictions are the surrogate's, not "
         "results. 8 test parts, one verifying run each: a small sample.", "",
         "**Fix 1 and fix 2 are applied together (no ablation)**: the difference v1 -> v2 "
         "cannot be split between them.", "",
         f"* fix 2: optimiser std weight {STD_WEIGHT} (v1: {STD_WEIGHT_V1}) on the mean squared "
         "model std, chosen a priori (the std then weighs as much as the residual); not tuned "
         "on the test parts. The transfer model has no conformal calibration (every rim-pass "
         "part is a TRAIN part of the base), so the raw ensemble std is used;",
         "* fix 1: the 8 role=fit TRAIN parts were compensated with the v1 model + the new "
         "penalty, simulated, and added to work/data_rim as `compensated` samples; the "
         "correction was refitted on 8 uncompensated + 8 compensated fit samples (the 4 "
         "holdout parts and all test parts kept out); its form picked by leave-one-part-out on "
         "TRAIN data only.", ""]
    L += ["## Test parts: vertical RMS of formed - target [mm]", "",
          "plate columns: backing plate (stage_n26 headline); v1 = combined_rim "
          f"(transfer v1, std weight {STD_WEIGHT_V1}); v2 = this run.", "",
          "| part | plate uncompensated | ML-MLP plate | combined v1 | combined v2 | "
          "v1 predicted | v2 predicted |", "|---|--:|--:|--:|--:|--:|--:|"]
    for _, r in df.iterrows():
        v2s = f3(r.v2_vert_rms_mm) if r.v2_status == "ok" else r.v2_status
        L.append(f"| {r.point_id} | {f3(r.plate_unc_vert_rms_mm)} | "
                 f"{f3(r.mlp_plate_vert_rms_mm)} | {f3(r.v1_vert_rms_mm)} | {v2s} | "
                 f"{f3(r.v1_pred_vert_rms_mm)} | {f3(r.v2_pred_vert_rms_mm)} |")
    m = df.mean(numeric_only=True)
    L.append(f"| **mean** | **{f3(m.plate_unc_vert_rms_mm)}** | **{f3(m.mlp_plate_vert_rms_mm)}**"
             f" | **{f3(m.v1_vert_rms_mm)}** | **{f3(m.v2_vert_rms_mm)}** | "
             f"**{f3(m.v1_pred_vert_rms_mm)}** | **{f3(m.v2_pred_vert_rms_mm)}** |")
    L += ["", f"v2 simulated on {n2} of 8 parts (means over the available ones).", "",
          "| mean over parts [mm] | plate uncompensated | ML-MLP plate | combined v1 | "
          "combined v2 |", "|---|--:|--:|--:|--:|"]
    for c, lab in (("upper_band_rms_mm", "upper-band RMS"), ("deep_rms_mm", "deep RMS")):
        L.append(f"| {lab} | {f3(m[f'plate_unc_{c}'])} | {f3(m[f'mlp_plate_{c}'])} | "
                 f"{f3(m[f'v1_{c}'])} | {f3(m[f'v2_{c}'])} |")
    L += ["", "Predicted vs simulated (mean vertical RMS): v1 "
          f"{f3(m.v1_pred_vert_rms_mm)} predicted / {f3(m.v1_vert_rms_mm)} simulated; v2 "
          f"{f3(m.v2_pred_vert_rms_mm)} predicted / {f3(m.v2_vert_rms_mm)} simulated. "
          "Predicted upper-band / deep RMS: v1 "
          f"{f3(m.v1_pred_upper_band_rms_mm)} / {f3(m.v1_pred_deep_rms_mm)}, v2 "
          f"{f3(m.v2_pred_upper_band_rms_mm)} / {f3(m.v2_pred_deep_rms_mm)}.", ""]
    if rf:
        L += ["## Transfer correction: old vs new (TRAIN rim-pass data only)", "",
              "Leave-one-part-out over the 8 fit parts: each part's samples (uncompensated and "
              "compensated) left out together, the correction refitted on the other 7 parts, "
              "dz prediction error RMS [mm] on the left-out samples (mean). 'old' = v1's "
              "protocol (linear, uncompensated samples only).", "",
              "| correction | trained on | all 16 | uncompensated 8 | compensated 8 |",
              "|---|---|--:|--:|--:|"]

        def lrow(lab, x, on):
            return (f"| {lab} | {on} | {f3(x['mean_all']['corrected_rms_mm'])} | "
                    f"{f3(x.get('mean_uncompensated', {}).get('corrected_rms_mm'))} | "
                    f"{f3(x.get('mean_compensated', {}).get('corrected_rms_mm'))} |")
        o = rf["old_protocol"]
        L.append(f"| base MLP (no correction) | - | {f3(o['mean_all']['base_rms_mm'])} | "
                 f"{f3(o['mean_uncompensated']['base_rms_mm'])} | "
                 f"{f3(o['mean_compensated']['base_rms_mm'])} |")
        L.append(lrow("old: linear", o, "uncompensated"))
        for c, x in rf["candidates"].items():
            L.append(lrow(f"new: {c}" + (" (chosen)" if c == rf["chosen"] else ""), x,
                          "uncompensated + compensated"))
        h = rf.get("holdout_uncompensated", {})
        L += ["", f"Chosen: **{rf['chosen']}** ({rf['rule']}); v2 accepted by its own "
              f"validation: {rf['v2_accepted']}. "
              f"v1 (as used) on the 8 compensated fit runs (not in its fit): dz error "
              f"{f3(rf['v1_on_compensated_fit_runs']['err_rms_mean_mm'])} mm."]
        if h:
            L += ["", "4 holdout TRAIN parts (uncompensated rim-pass runs, in neither fit), "
                  f"dz error [mm]: base {f3(h['base']['err_rms_mean_mm'])}, v1 "
                  f"{f3(h['v1']['err_rms_mean_mm'])}, v2 {f3(h['v2']['err_rms_mean_mm'])}."]
        L.append("")
    L += ["## The 8 compensated fit runs (TRAIN, data for fix 1)", "",
          "v1 model + std weight 1.0, DSIF + rim pass, penalty 10; vertical RMS [mm].", "",
          "| part | rim pass uncompensated | predicted | simulated | upper band | deep | "
          "status |", "|---|--:|--:|--:|--:|--:|---|"]
    for _, r in frd.iterrows():
        L.append(f"| {r.point_id} | {f3(r.rim_unc_vert_rms_mm)} | {f3(r.pred_vert_rms_mm)} | "
                 f"{f3(r.sim_vert_rms_mm)} | {f3(r.sim_upper_band_rms_mm)} | "
                 f"{f3(r.sim_deep_rms_mm)} | {r.status} |")
    L += ["", "Made by `python/scripts/combined_rim2.py` (run: `benchmarks/springback_fine_ml/"
          "combined_rim2/run.sh`). v1 results: `benchmarks/springback_fine_ml/combined_rim/`. "
          "Per-part rows: `results.csv`, fit runs: `fit_runs.csv`, correction fit: "
          "`refit.json`.", ""]
    (OUT / "README.md").write_text("\n".join(L))
    print("\n".join(L))


def main() -> int:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("PRECOMP_ML_THREADS", "1")
    try:
        import torch
        torch.set_num_threads(1)
    except ImportError:
        pass
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["all", "dry", "status", "report"])
    ap.add_argument("--part", default=None, help="dry: the fit part")
    ap.add_argument("--max-time", type=float, default=240.0,
                    help="optimiser wall-time cap [s] (v1's)")
    a = ap.parse_args()
    {"all": cmd_all, "dry": cmd_dry, "status": cmd_status, "report": cmd_report}[a.what](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
