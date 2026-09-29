"""Materials: the SparLab block, the library, power-law conversion."""

import json

import numpy as np
import pytest

from precomp.materials import NOMINAL_SOURCE, Material, get_material, library

CURRENT_KEYS = {"yield_stress", "hardening_modulus", "saturation_stress", "saturation_rate",
                "kinematic_hardening_modulus"}


def test_to_sparlab_has_exactly_the_current_keys_unless_extras_are_set():
    m = Material("plain", 70e9, 0.33, 2700.0, 150e6, hardening_modulus=1e9,
                 saturation_stress=50e6, saturation_rate=20.0,
                 kinematic_hardening_modulus=2e8)
    doc = m.to_sparlab()
    assert set(doc) == {"name", "youngs_modulus", "poisson_ratio", "density", "plasticity"}
    assert set(doc["plasticity"]) == CURRENT_KEYS
    assert doc["plasticity"]["yield_stress"] == 150e6
    full = m.replace(r0=0.7, r45=0.8, r90=0.9, af_C=5e9, af_gamma=40.0).to_sparlab()
    p = full["plasticity"]
    # docs/configuration.md, material.plasticity: Hill48 by r-values, Chaboche backstresses
    assert p["yield_criterion"] == "hill48"
    assert p["anisotropy"] == {"r0": 0.7, "r45": 0.8, "r90": 0.9,
                               "rolling_direction": [1.0, 0.0, 0.0],
                               "sheet_normal": [0.0, 0.0, 1.0]}
    assert p["backstresses"] == [{"modulus": 5e9, "recovery": 40.0}]
    assert set(p) - CURRENT_KEYS == {"yield_criterion", "anisotropy", "backstresses"}


def test_invalid_materials_are_refused():
    with pytest.raises(ValueError, match="poisson"):
        Material("x", 70e9, 0.5, 2700.0, 1e8)
    with pytest.raises(ValueError, match="saturation_rate"):
        Material("x", 70e9, 0.3, 2700.0, 1e8, saturation_stress=1e7)
    with pytest.raises(ValueError, match="Armstrong"):
        Material("x", 70e9, 0.3, 2700.0, 1e8, af_C=1e9)
    with pytest.raises(ValueError, match="Hill48"):
        Material("x", 70e9, 0.3, 2700.0, 1e8, r0=1.0)
    with pytest.raises(ValueError, match="af_C"):              # SparLab: modulus > 0
        Material("x", 70e9, 0.3, 2700.0, 1e8, af_C=0.0, af_gamma=10.0)


def test_library_is_complete_marked_and_serialisable():
    lib = library()
    assert set(lib) == {"AA5754-O", "AA2024-O", "AA6061-T6", "AA7075-O", "DC04",
                        "Ti-6Al-4V-annealed"}
    for name, m in lib.items():
        assert m.source == NOMINAL_SOURCE
        assert m.hardening_source["rms_relative_error"] < 0.05
        json.dumps(m.to_sparlab())
        assert Material.from_dict(json.loads(json.dumps(m.to_dict()))) == m
        a = np.linspace(0, 1, 50)
        assert np.all(np.diff(m.uniaxial_stress(a)) >= 0)          # hardening, never softening
    assert get_material("dc04").name == "DC04"
    with pytest.raises(KeyError):
        get_material("unobtainium")


def test_swift_fit_starts_at_the_yield_point_and_follows_the_curve():
    K, eps0, n = 530e6, 0.012, 0.22
    m = Material.from_swift("steel", 210e9, 0.3, 7850.0, K, eps0, n, fit_range=(0.0, 0.5))
    assert m.yield_stress == pytest.approx(K * eps0 ** n, rel=1e-14)
    a = np.linspace(0.0, 0.5, 101)
    rel = np.abs(m.flow_stress(a) - K * (eps0 + a) ** n) / (K * (eps0 + a) ** n)
    assert rel.max() == pytest.approx(m.hardening_source["max_relative_error"], rel=0.05)
    assert rel.max() < 0.05
    assert m.hardening_source["law"] == "swift"


def test_hollomon_is_swift_through_the_yield_point():
    m = Material.from_hollomon("al", 70e9, 0.33, 2700.0, K=420e6, n=0.3, yield_stress=100e6)
    assert m.yield_stress == pytest.approx(100e6, rel=1e-12)
    assert m.hardening_source["eps0"] == pytest.approx((100 / 420) ** (1 / 0.3))
