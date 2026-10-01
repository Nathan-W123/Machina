"""Variant specs of the material lens (see matlens.apply_variant / build_variant).

Every entry changes ONE thing relative to "base" (the benchmark's library
AA5754-O: E 70 GPa, yield 100 MPa, Hollomon K 420 MPa n 0.30 fitted by linear
+ Voce, Hill48 r 0.75/0.70/0.80, isotropic hardening only; friction 0.1;
standard Hex8 x 2 layers) unless its name says "combo" or ends in "_im5".
Literature sources of the numbers: README.md, section 'Literature'.
"""
import math

import numpy as np
from scipy.optimize import minimize_scalar, nnls

import matlens

MPA = 1e6


def voce_fit_with_backstresses(backstresses, lo=0.01, hi=0.6):
    """Isotropic (sy0, H, Q, delta) such that iso + sum of AF backstresses
    reproduces the library's monotonic uniaxial curve over plastic strain
    [lo, hi] (least squares; sy0, H, Q >= 0). backstresses: [(C [Pa], gamma)]."""
    m = matlens.library_material()
    a = np.linspace(0.0, hi, 601)
    base = m.flow_stress(a)
    X = sum(C / g * (1 - np.exp(-g * a)) for C, g in backstresses)
    w = a >= lo

    def solve(ld):
        d = math.exp(ld)
        A = np.column_stack([np.ones(w.sum()), a[w], 1 - np.exp(-d * a[w])])
        return nnls(A, (base - X)[w])

    best = minimize_scalar(lambda ld: solve(ld)[1], bounds=(math.log(0.1), math.log(500.0)),
                           method="bounded", options={"xatol": 1e-8})
    (sy0, H, Q), _ = solve(best.x)
    return float(sy0), float(H), float(Q), float(math.exp(best.x))


def _kin(backstresses):
    """A Chaboche variant: the AF backstresses as SparLab `backstresses` (deck
    patch: precomp's Material carries one AF term only) and the isotropic part
    refitted so that the monotonic curve stays the library's beyond eps_p 0.01."""
    sy0, H, Q, d = voce_fit_with_backstresses(backstresses)
    return {"material": {"yield_stress": sy0, "hardening_modulus": H, "saturation_stress": Q,
                         "saturation_rate": d},
            "deck": {"material": {"plasticity": {
                "backstresses": [{"modulus": float(C), "recovery": float(g)}
                                 for C, g in backstresses]}}}}


def _kin_on(backstresses, iso):
    """Backstresses with the isotropic part refitted so that the monotonic curve
    is the given linear + Voce law `iso` (Material fields) beyond eps_p 0.01."""
    a = np.linspace(0.0, 0.6, 601)
    target = (iso["yield_stress"] + iso["hardening_modulus"] * a
              + iso["saturation_stress"] * (1 - np.exp(-iso["saturation_rate"] * a)))
    X = sum(C / g * (1 - np.exp(-g * a)) for C, g in backstresses)
    w = a >= 0.01

    def solve(ld):
        d = math.exp(ld)
        A = np.column_stack([np.ones(w.sum()), a[w], 1 - np.exp(-d * a[w])])
        return nnls(A, (target - X)[w])

    best = minimize_scalar(lambda ld: solve(ld)[1], bounds=(math.log(0.1), math.log(500.0)),
                           method="bounded", options={"xatol": 1e-8})
    (sy0, H, Q), _ = solve(best.x)
    return {"material": {"yield_stress": float(sy0), "hardening_modulus": float(H),
                         "saturation_stress": float(Q), "saturation_rate": float(math.exp(best.x))},
            "deck": {"material": {"plasticity": {
                "backstresses": [{"modulus": float(C), "recovery": float(g)}
                                 for C, g in backstresses]}}}}


# Shutov-Kreissig two-backstress fit of Laurent et al. (2009) cyclic simple
# shear data on 1 mm AA5754-O (Shutov 2016, arXiv:1510.09163, Table 1:
# c1 = 115.5 MPa, kappa1 = 0.0885 1/MPa, c2 = 11500 MPa, kappa2 = 0.01676 1/MPa),
# mapped to Armstrong-Frederick in uniaxial terms: C = 1.5 c, gamma = c kappa sqrt(3/2)
# (saturation C/gamma = sqrt(3/2)/kappa: 13.8 and 73.1 MPa).
SHUTOV_AF = [(1.5 * 115.5 * MPA, 115.5 * 0.0885 * math.sqrt(1.5)),
             (1.5 * 11500.0 * MPA, 11500.0 * 0.01676 * math.sqrt(1.5))]
# A weak Bauschinger effect (Coer et al. 2018: "reversed shear results showed a
# very weak Bauschinger effect"): one backstress saturating at 20 MPa.
MILD_AF = [(2000.0 * MPA, 100.0)]


def _voce(Y0, Ysat, Cy):
    return {"yield_stress": Y0 * MPA, "hardening_modulus": 0.0,
            "saturation_stress": (Ysat - Y0) * MPA, "saturation_rate": Cy}


def _hollomon(K, n, sy):
    from precomp.materials import Material
    m = Material.from_hollomon("tmp", 70e9, 0.33, 2670.0, K=K * MPA, n=n, yield_stress=sy * MPA)
    return {"yield_stress": m.yield_stress, "hardening_modulus": m.hardening_modulus,
            "saturation_stress": m.saturation_stress, "saturation_rate": m.saturation_rate}


def _hill_stress_ratios(s45, s90, sb):
    """Hill48 by stress ratios (anisotropy.stress_ratios) instead of r-values."""
    return {"material": {"r0": None, "r45": None, "r90": None},
            "deck": {"material": {"plasticity": {
                "yield_criterion": "hill48",
                "anisotropy": {"stress_ratios": {"sigma_45": s45, "sigma_90": s90,
                                                 "sigma_biaxial": sb},
                               "rolling_direction": [1, 0, 0], "sheet_normal": [0, 0, 1]}}}}}


IM5 = {"deck": {"model": {"element_formulation": "incompatible_modes",
                          "integration": {"thickness_points": 5, "thickness_direction": "z"}}}}

# E(eps_p) = E0 - (E0 - Ea)(1 - exp(-xi eps_p)) (Yoshida-Uemori's chord-modulus
# law) with Ea = 0.85 E0, xi = 20: -15 % at saturation, -13 % at eps_p = 0.1;
# Al alloys: -10 % (Cleveland & Ghosh 2002, AA6022-T4) to -20 % (Nietsch et al.
# 2024, EN AW-6082-T6). Per element from the base run's final plastic strain.
EDEG = {"e_of_ep": {"Ea": 0.85 * 70e9, "xi": 20.0, "from": "base"}}


def _with(*specs):
    out = {}
    for s in specs:
        for k, v in s.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                matlens.deep_merge(out[k], v)
            else:
                out[k] = matlens.copy.deepcopy(v)
    return out


VM = {"material": {"r0": None, "r45": None, "r90": None}}
VOCE_COER = {"material": _voce(102.75, 292.14, 13.50)}          # Coer et al. 2018, RD

VARIANTS = {
    "base": {},
    # --- elasticity
    "E63": {"material": {"youngs_modulus": 63e9}},                  # -10 % (chord modulus)
    "Edeg": EDEG,                                                    # E(eps_p) per element
    # --- initial yield (whole curve shifted by +-15 MPa)
    "ys85": {"material": {"yield_stress": 85 * MPA}},
    "ys115": {"material": {"yield_stress": 115 * MPA}},
    # --- hardening curve (literature fits, full curve replaced)
    "hard_voce_coer": VOCE_COER,
    "hard_pow_iadicola": {"material": _hollomon(474.0, 0.317, 94.1)},  # Iadicola et al. 2008, RD
    # --- yield locus
    "vm": VM,                                                        # von Mises
    "hill_iadicola": {"material": {"r0": 0.69, "r45": 0.73, "r90": 0.87}},
    "hill_coer": {"material": {"r0": 0.663, "r45": 0.860, "r90": 0.717}},
    # stress-based Hill48 with the biaxial elongation Iadicola et al. measured
    # at 10-15 % strain (sigma_b / sigma_0 ~ 1.1; implies r0 ~ 1.5)
    "hill_sb110": _hill_stress_ratios(90.9 / 94.1, 92.1 / 94.1, 1.10),
    # --- kinematic hardening (monotonic curve kept beyond eps_p = 0.02)
    "kin_shutov": _kin(SHUTOV_AF),
    "kin_mild": _kin(MILD_AF),
    # --- friction
    "mu005": {"setup": {"friction": 0.05}},
    "mu02": {"setup": {"friction": 0.2}},
    # --- combination: every literature-supported change at once (von Mises,
    # Coer Voce curve with the Shutov backstresses refitted onto it, E(eps_p))
    "combo_lit": _with(VM, _kin_on(SHUTOV_AF, VOCE_COER["material"]), EDEG),
    # --- checks of the scale-invariance argument and of the sheet thickness
    # every stress parameter and E scaled by 0.9: shape should not change
    "scaled090": {"material": {"youngs_modulus": 63e9, "yield_stress": 90 * MPA,
                               "hardening_modulus": 0.9 * matlens.library_material().hardening_modulus,
                               "saturation_stress": 0.9 * matlens.library_material().saturation_stress}},
    "t095": {"setup": {"thickness": 0.95e-3}},                        # -5 % sheet thickness
    # --- element check: incompatible-mode Hex8, 5 points through each layer
    "base_im5": IM5,
    "E63_im5": _with({"material": {"youngs_modulus": 63e9}}, IM5),
    "kin_shutov_im5": _with(_kin(SHUTOV_AF), IM5),
    "hard_voce_coer_im5": _with(VOCE_COER, IM5),
    "vm_im5": _with(VM, IM5),
    "kin_mild_im5": _with(_kin(MILD_AF), IM5),
    "combo_lit_im5": _with(VM, _kin_on(SHUTOV_AF, VOCE_COER["material"]),
                           {"e_of_ep": dict(EDEG["e_of_ep"], **{"from": "base_im5"})}, IM5),
    "hill_sb110_im5": _with(_hill_stress_ratios(90.9 / 94.1, 92.1 / 94.1, 1.10), IM5),
}
