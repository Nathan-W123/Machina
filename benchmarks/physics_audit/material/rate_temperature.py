#!/usr/bin/env python3
"""Order-of-magnitude strain rate and heating of the benchmark's forming (lens part c).

Inputs from the base cone run (plastic strain per element, tool path
length) and stated assumptions (feed speed, contact length, number of passes
that strain a point, Taylor-Quinney factor, AA5754 rho, c). The simulation
itself is rate- and temperature-independent; this estimates how far the real
process departs from that, in flow stress, and converts it to shape with the
measured sensitivities (sensitivity.csv: ys85 / ys115 / scaled090).
Writes rate_temperature.json.
"""
import json
from pathlib import Path

import numpy as np

import matlens

HERE = Path(__file__).resolve().parent
PART = "truncated_cone-s2026-0000"


def main():
    run = HERE / "runs" / PART / "base"
    summ = json.loads((run / "output" / "summary.json").read_text())
    el = np.loadtxt(run / "output" / "step_3_release_elements.csv", delimiter=",", skiprows=1)
    ep = el[:, 1]
    m = matlens.library_material()
    a = np.linspace(0, ep.max(), 2001)
    w_curve = np.concatenate([[0], np.cumsum(0.5 * (m.flow_stress(a[1:]) + m.flow_stress(a[:-1]))
                                             * np.diff(a))])
    w_max = float(np.interp(ep.max(), a, w_curve))           # J/m^3 at the most strained element
    rho, c, beta_tq = 2670.0, 900.0, 0.9                       # AA5754: kg/m^3, J/(kg K)
    k_th = 130.0                                               # W/(m K)
    alpha_th = k_th / (rho * c)
    # assumptions (not measured here): the tool strains a wall point over
    # n_pass passes (step-down 1 mm, wall 38.6 deg -> 1.25 mm pitch; contact
    # footprint 2.5-4 mm -> 2-3 passes); contact length along the path 2-3 mm;
    # feed 17-100 mm/s (1-6 m/min, a common ISF range; RoboForming's speeds
    # are not public)
    n_pass = (2.0, 3.0)
    L_c = (2e-3, 3e-3)
    v = (17e-3, 100e-3)
    d_eps = (ep.max() / n_pass[1], ep.max() / n_pass[0])
    rate_lo = d_eps[0] * v[0] / L_c[1]
    rate_hi = d_eps[1] * v[1] / L_c[0]
    t_contact = (L_c[0] / v[1], L_c[1] / v[0])
    diff_len = [float(np.sqrt(alpha_th * t)) for t in t_contact]
    dT_adiabatic_total = beta_tq * w_max / (rho * c)
    dT_adiabatic_pass = (dT_adiabatic_total / n_pass[1], dT_adiabatic_total / n_pass[0])
    # AA5754 at room temperature: negative strain-rate sensitivity (dynamic
    # strain ageing); |m| <= 0.005 assumed as a bound (Pandey et al. 2013 report
    # the sign, not a value in the text) -> flow stress change over the decades
    # between a quasi-static coupon (1e-3 /s) and the process
    decades = (np.log10(rate_lo / 1e-3), np.log10(rate_hi / 1e-3))
    m_srs = 0.005
    dsig_rate = [float(m_srs * np.log(10) * d) for d in decades]
    out = {
        "part": PART, "max_eq_plastic_strain": float(ep.max()),
        "tool_path_length_m": summ["tools"][0]["trajectory_length_m"],
        "assumptions": {"feed_m_per_s": v, "contact_length_m": L_c, "passes_straining_a_point": n_pass,
                        "rho_kg_m3": rho, "c_J_kgK": c, "taylor_quinney": beta_tq, "k_W_mK": k_th,
                        "srs_bound_m": m_srs},
        "strain_increment_per_pass": d_eps,
        "strain_rate_range_per_s": (float(rate_lo), float(rate_hi)),
        "decades_above_quasistatic_coupon_1e-3": [float(x) for x in decades],
        "flow_stress_change_from_rate_bound_rel": dsig_rate,
        "plastic_work_max_J_m3": w_max,
        "adiabatic_temperature_rise_total_K": float(dT_adiabatic_total),
        "adiabatic_temperature_rise_per_pass_K": [float(x) for x in dT_adiabatic_pass],
        "contact_time_s": t_contact, "thermal_diffusion_length_in_contact_time_m": diff_len,
        "sheet_thickness_m": 1e-3,
    }
    (HERE / "rate_temperature.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
