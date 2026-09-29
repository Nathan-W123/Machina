/// \file test_plasticity.cpp
/// \brief J2 plasticity: the return map at a material point, the
///        elastoplastic element, and the solver's history.
///
/// Every check is exact or a derivative check. Uniaxial stress is a
/// proportional path, on which the backward-Euler radial return with linear
/// hardening is exact for any step, so stresses are compared with the
/// closed-form curves to round-off - monotonic, reversed (Bauschinger) and
/// with Voce saturation. Consistent tangents are compared with central
/// differences of the stress (3-D, plane strain, plane stress) and of the
/// element's internal force (with and without mean dilatation).
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Elastoplastic.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/material/Plasticity.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

IsotropicMaterial plastic_steel(Scalar h = 0.0, Scalar hk = 0.0, Scalar q = 0.0,
                                Scalar delta = 0.0) {
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  PlasticityParameters p;
  p.yield_stress = 250.0e6;
  p.hardening_modulus = h;
  p.kinematic_hardening_modulus = hk;
  p.saturation_stress = q;
  p.saturation_rate = delta;
  m.set_plasticity(p);
  return m;
}

/// Uniaxial stress sigma_11 at the axial strain e11: the lateral strains
/// found by Newton on the return so that sigma_22 = sigma_33 = 0.
PlasticResponse uniaxial(const IsotropicMaterial& m, Scalar e11, const PlasticState& committed,
                         Scalar delta_t = 0.0) {
  const Scalar th = m.thermal_expansion() * delta_t;
  Vector6 strain = Vector6::Zero();
  strain(0) = e11;
  for (int i = 1; i < 3; ++i) {
    strain(i) = committed.plastic_strain(i) + th -
                m.poisson_ratio() * (e11 - committed.plastic_strain(0) - th);
  }
  PlasticResponse r;
  for (int it = 0; it < 50; ++it) {
    r = j2_return(m, StressState::ThreeDimensional, strain, committed, delta_t);
    const Eigen::Vector2d residual(r.stress(1), r.stress(2));
    if (residual.cwiseAbs().maxCoeff() <=
        1.0e-13 * std::max(std::abs(r.stress(0)), m.plasticity().yield_stress)) {
      return r;
    }
    Eigen::Matrix2d j;
    j << r.tangent(1, 1), r.tangent(1, 2), r.tangent(2, 1), r.tangent(2, 2);
    const Eigen::Vector2d step = j.lu().solve(residual);
    strain(1) -= step(0);
    strain(2) -= step(1);
  }
  FAIL("the uniaxial-stress iteration did not converge");
  return r;
}

/// The uniaxial stress of a Voce material at total strain e (monotonic):
/// sigma = sigma_y(e - sigma / E), by bisection.
Scalar voce_uniaxial(const IsotropicMaterial& m, Scalar e) {
  const PlasticityParameters& p = m.plasticity();
  Scalar lo = p.yield_stress;
  Scalar hi = m.youngs_modulus() * e;
  for (int it = 0; it < 200; ++it) {
    const Scalar mid = 0.5 * (lo + hi);
    (mid - p.yield(e - mid / m.youngs_modulus()) > 0.0 ? hi : lo) = mid;
  }
  return 0.5 * (lo + hi);
}

Vector6 random_strain(Scalar amplitude, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  Vector6 v;
  for (int i = 0; i < 6; ++i) v(i) = amplitude * dist(gen);
  return v;
}

/// The consistent tangent against central differences of the stress, over
/// the strain components the idealisation takes.
Scalar tangent_error(const IsotropicMaterial& m, StressState state, const Vector6& strain,
                     const PlasticState& committed, Scalar delta_t = 0.0) {
  const PlasticResponse r = j2_return(m, state, strain, committed, delta_t);
  // (The plane-stress iteration leaves sigma_33 below 1e-12 of the stress,
  // which a step of 1e-8 keeps well below the tolerance.)
  const Scalar h = 1.0e-8;
  Scalar worst = 0.0;
  for (int j = 0; j < 6; ++j) {
    if (state != StressState::ThreeDimensional && (j == 4 || j == 5)) continue;
    if (state == StressState::PlaneStress && j == 2) continue;
    Vector6 up = strain;
    Vector6 um = strain;
    up(j) += h;
    um(j) -= h;
    const Vector6 fd = (j2_return(m, state, up, committed, delta_t).stress -
                        j2_return(m, state, um, committed, delta_t).stress) /
                       (2.0 * h);
    Vector6 column = r.tangent.col(j);
    worst = std::max(worst, (fd - column).cwiseAbs().maxCoeff());
  }
  return worst / r.tangent.cwiseAbs().maxCoeff();
}

}  // namespace

TEST_CASE("plasticity parameters are validated and the hardening laws are consistent",
          "[plasticity][material]") {
  IsotropicMaterial m(200.0e9, 0.3, 7800.0, "steel");
  PlasticityParameters p;
  p.yield_stress = -1.0;
  REQUIRE_THROWS_AS(m.set_plasticity(p), ConfigError);
  p.yield_stress = 250.0e6;
  p.saturation_stress = 100.0e6;  // no rate
  REQUIRE_THROWS_AS(m.set_plasticity(p), ConfigError);
  p = PlasticityParameters();
  p.hardening_modulus = 1.0e9;  // hardening of an elastic material
  REQUIRE_THROWS_AS(m.set_plasticity(p), ConfigError);
  p.yield_stress = std::numeric_limits<Scalar>::quiet_NaN();
  REQUIRE_THROWS_AS(m.set_plasticity(p), ConfigError);

  // sigma_y' is the derivative of sigma_y, and the stored isotropic
  // hardening energy the integral of sigma_y - sigma_y0.
  p = PlasticityParameters();
  p.yield_stress = 250.0e6;
  p.hardening_modulus = 1.5e9;
  p.saturation_stress = 120.0e6;
  p.saturation_rate = 25.0;
  REQUIRE_NOTHROW(m.set_plasticity(p));
  for (const Scalar a : {0.0, 0.003, 0.04, 0.2}) {
    const Scalar h = 1.0e-7;
    REQUIRE((p.yield(a + h) - p.yield(a - (a > 0 ? h : 0.0))) / (a > 0 ? 2.0 * h : h) ==
            Approx(p.yield_slope(a)).epsilon(1e-6));
    if (a > 0.0) {
      REQUIRE((p.isotropic_energy(a + h) - p.isotropic_energy(a - h)) / (2.0 * h) ==
              Approx(p.yield(a) - p.yield_stress).epsilon(1e-6));
    }
  }
  REQUIRE(p.isotropic_energy(0.0) == 0.0);
  // The copy with a new modulus keeps the plasticity.
  REQUIRE(m.with_youngs_modulus(100.0e9).plasticity().saturation_rate == 25.0);
}

TEST_CASE("uniaxial stress with linear hardening is exact for any step",
          "[plasticity][material]") {
  // sigma = sigma_y0 + (H + H_kin) eps_p on a monotonic path, so
  // sigma = E (sigma_y0 + (H + H_kin) eps) / (E + H + H_kin).
  struct Case {
    Scalar h;
    Scalar hk;
  };
  for (const Case c : {Case{0.0, 0.0}, Case{2.0e9, 0.0}, Case{0.0, 5.0e9}, Case{1.0e9, 3.0e9}}) {
    const IsotropicMaterial m = plastic_steel(c.h, c.hk);
    const Scalar e = m.youngs_modulus();
    const Scalar hh = c.h + c.hk;
    const Scalar target = 0.01;  // eight times the yield strain
    const Scalar exact = e * (250.0e6 + hh * target) / (e + hh);
    // One step.
    const PlasticResponse one = uniaxial(m, target, PlasticState());
    INFO("H = " << c.h << ", H_kin = " << c.hk);
    REQUIRE(one.yielding);
    REQUIRE(one.stress(0) == Approx(exact).epsilon(1e-12));
    // Forty steps, each from the state the last one committed.
    PlasticState state;
    PlasticResponse r;
    for (int k = 1; k <= 40; ++k) {
      r = uniaxial(m, target * k / 40.0, state);
      state = r.state;
    }
    REQUIRE(r.stress(0) == Approx(exact).epsilon(1e-12));
    const Scalar eps_p = target - exact / e;
    REQUIRE(state.plastic_strain(0) == Approx(eps_p).epsilon(1e-10));
    REQUIRE(state.plastic_strain(1) == Approx(-0.5 * eps_p).epsilon(1e-10));
    REQUIRE(state.equivalent_plastic_strain == Approx(eps_p).epsilon(1e-10));
    REQUIRE(von_mises_stress(r.stress) == Approx(exact).epsilon(1e-12));
    // Stored energy: elastic sigma^2 / 2E plus H alpha^2 / 2 plus the back
    // stress's H_kin eps_p^2 / 2.
    REQUIRE(r.energy == Approx(exact * exact / (2.0 * e) + 0.5 * hh * eps_p * eps_p)
                            .epsilon(1e-10));
    // The lateral contraction is elastic plus plastic, -nu sigma/E - eps_p/2.
    REQUIRE(uniaxial(m, target, PlasticState()).state.plastic_strain(2) ==
            Approx(-0.5 * eps_p).epsilon(1e-10));
  }
}

TEST_CASE("a reversed uniaxial path shows isotropic and kinematic hardening exactly",
          "[plasticity][material]") {
  const Scalar sy = 250.0e6;
  const Scalar e1 = 0.01;
  for (const bool kinematic : {false, true}) {
    const Scalar hmod = 10.0e9;
    const IsotropicMaterial m = kinematic ? plastic_steel(0.0, hmod) : plastic_steel(hmod, 0.0);
    const Scalar e = m.youngs_modulus();
    const Scalar s1 = e * (sy + hmod * e1) / (e + hmod);
    const Scalar ep1 = e1 - s1 / e;
    // Reverse yield: isotropic at -(sigma_y0 + H eps_p1); kinematic at
    // H_kin eps_p1 - sigma_y0 (the Bauschinger effect).
    const Scalar reverse_yield = kinematic ? hmod * ep1 - sy : -(sy + hmod * ep1);
    PlasticState state = uniaxial(m, e1, PlasticState()).state;
    // Unload halfway to reverse yield: elastic.
    const Scalar e_mid = e1 - 0.5 * (s1 - reverse_yield) / e;
    PlasticResponse r = uniaxial(m, e_mid, state);
    INFO((kinematic ? "kinematic" : "isotropic"));
    REQUIRE_FALSE(r.yielding);
    REQUIRE(r.stress(0) == Approx(0.5 * (s1 + reverse_yield)).epsilon(1e-12));
    state = r.state;
    // Past reverse yield to -e1, in one step.
    r = uniaxial(m, -e1, state);
    REQUIRE(r.yielding);
    // Isotropic: sigma = -(sigma_y0 + H alpha), alpha = eps_p1 + (eps_p1 -
    // eps_p), so sigma = E (H e - sigma_y0 - 2 H eps_p1) / (E + H).
    // Kinematic: sigma = H_kin eps_p - sigma_y0, so sigma = E (H_kin e -
    // sigma_y0) / (E + H_kin).
    const Scalar exact = kinematic ? e * (hmod * -e1 - sy) / (e + hmod)
                                   : e * (hmod * -e1 - sy - 2.0 * hmod * ep1) / (e + hmod);
    REQUIRE(r.stress(0) == Approx(exact).epsilon(1e-12));
    // The kinematic material's surface kept its size; the isotropic one's grew.
    const Scalar radius = von_mises_stress(r.stress - r.state.back_stress);
    REQUIRE(radius == Approx(kinematic ? sy : sy + hmod * r.state.equivalent_plastic_strain)
                          .epsilon(1e-10));
  }
}

TEST_CASE("Voce saturation is followed exactly on a uniaxial path", "[plasticity][material]") {
  const IsotropicMaterial m = plastic_steel(1.0e9, 0.0, 150.0e6, 30.0);
  for (const Scalar target : {0.002, 0.02, 0.2}) {
    const Scalar exact = voce_uniaxial(m, target);
    REQUIRE(uniaxial(m, target, PlasticState()).stress(0) == Approx(exact).epsilon(1e-11));
    PlasticState state;
    PlasticResponse r;
    for (int k = 1; k <= 25; ++k) {
      r = uniaxial(m, target * k / 25.0, state);
      state = r.state;
    }
    REQUIRE(r.stress(0) == Approx(exact).epsilon(1e-11));
  }
}

TEST_CASE("the consistent tangent is the derivative of the return", "[plasticity][material]") {
  // Every hardening mechanism at once, from a committed state with plastic
  // strain, back stress and accumulated strain, on a non-proportional
  // increment that yields again.
  const IsotropicMaterial m = plastic_steel(1.5e9, 4.0e9, 100.0e6, 20.0);
  for (const StressState state :
       {StressState::ThreeDimensional, StressState::PlaneStrain, StressState::PlaneStress}) {
    Vector6 first = random_strain(4.0e-3, 3u);
    Vector6 second = first + random_strain(3.0e-3, 5u);
    if (state != StressState::ThreeDimensional) {
      first(4) = first(5) = second(4) = second(5) = 0.0;
      if (state == StressState::PlaneStrain) first(2) = second(2) = 0.0;
    }
    const PlasticResponse r1 = j2_return(m, state, first, PlasticState(), 20.0);
    REQUIRE(r1.yielding);
    const PlasticResponse r2 = j2_return(m, state, second, r1.state, 35.0);
    INFO(to_string(state));
    REQUIRE(r2.yielding);
    REQUIRE(tangent_error(m, state, second, r1.state, 35.0) <= 1.0e-6);
    // The elastic branch too: a small step back inside the surface.
    const Vector6 back = second - 0.1 * (second - first);
    REQUIRE_FALSE(j2_return(m, state, back, r2.state).yielding);
    REQUIRE(tangent_error(m, state, back, r2.state) <= 1.0e-6);
    // The tangent is symmetric (associative flow, J2).
    REQUIRE((r2.tangent - r2.tangent.transpose()).cwiseAbs().maxCoeff() <=
            1.0e-12 * r2.tangent.cwiseAbs().maxCoeff());
  }
}

TEST_CASE("plane stress returns with sigma_33 = 0 and the 3-D state at its eps_33",
          "[plasticity][material]") {
  const IsotropicMaterial m = plastic_steel(2.0e9, 1.0e9);
  Vector6 strain = Vector6::Zero();
  strain(0) = 6.0e-3;
  strain(1) = -1.0e-3;
  strain(3) = 4.0e-3;
  const PlasticResponse r = j2_return(m, StressState::PlaneStress, strain, PlasticState());
  REQUIRE(r.yielding);
  REQUIRE(std::abs(r.stress(2)) <= 1.0e-11 * std::abs(r.stress(0)));
  REQUIRE(r.state.thickness_strain == r.strain_33);
  Vector6 full = strain;
  full(2) = r.strain_33;
  const PlasticResponse check = j2_return(m, StressState::ThreeDimensional, full, PlasticState());
  REQUIRE((check.stress - r.stress).cwiseAbs().maxCoeff() <= 1.0e-12 * check.stress.norm());
  // The condensed tangent has no thickness row or column.
  REQUIRE(r.tangent.row(2).cwiseAbs().maxCoeff() == 0.0);
  REQUIRE(r.tangent.col(2).cwiseAbs().maxCoeff() == 0.0);
  // Below yield it is the plane-stress elasticity matrix.
  const PlasticResponse elastic =
      j2_return(m, StressState::PlaneStress, 0.05 * strain, PlasticState());
  REQUIRE_FALSE(elastic.yielding);
  const Matrix3 d = m.plane_stress_matrix();
  const int rows[3] = {0, 1, 3};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      REQUIRE(elastic.tangent(rows[i], rows[j]) == Approx(d(i, j)).margin(1e-6 * d(0, 0)));
    }
  }
}

TEST_CASE("the loading flag gives the continuum tangent at zero increment",
          "[plasticity][material]") {
  // After a plastic step the point sits on its surface. At the same strain
  // the return is elastic (Delta gamma = 0) and, the point having yielded,
  // its tangent is the continuum elastoplastic one: the limit of the
  // consistent tangent as the increment vanishes.
  const Scalar h = 2.0e9;
  const Scalar hk = 1.0e9;
  const IsotropicMaterial m = plastic_steel(h, hk);
  Vector6 strain = Vector6::Zero();
  strain(0) = 5.0e-3;
  strain(3) = 2.0e-3;
  const PlasticResponse r1 = j2_return(m, StressState::ThreeDimensional, strain, PlasticState());
  REQUIRE(r1.state.loading);
  const PlasticResponse again = j2_return(m, StressState::ThreeDimensional, strain, r1.state);
  REQUIRE_FALSE(again.yielding);
  REQUIRE(again.state.loading);
  REQUIRE((again.stress - r1.stress).cwiseAbs().maxCoeff() <= 1.0e-6);
  const Scalar g = m.shear_modulus();
  Vector6 xi = r1.stress - r1.state.back_stress;
  const Scalar mean = (xi(0) + xi(1) + xi(2)) / 3.0;
  for (int i = 0; i < 3; ++i) xi(i) -= mean;
  const Scalar norm = std::sqrt(xi.head(3).squaredNorm() + 2.0 * xi.tail(3).squaredNorm());
  const Vector6 n = xi / norm;
  const Matrix6 continuum = m.three_dimensional_matrix() -
                            2.0 * g / (1.0 + (h + hk) / (3.0 * g)) * (n * n.transpose());
  REQUIRE((again.tangent - continuum).cwiseAbs().maxCoeff() <=
          1.0e-9 * continuum.cwiseAbs().maxCoeff());
  // A small further increment along the loading direction follows it.
  Vector6 d = n;
  d.tail(3) *= 2.0;  // engineering shear
  d *= 1.0e-9;
  const PlasticResponse further =
      j2_return(m, StressState::ThreeDimensional, strain + d, r1.state);
  REQUIRE(further.yielding);
  REQUIRE((further.stress - r1.stress - continuum * d).norm() <=
          1.0e-4 * (continuum * d).norm());
  // A step elsewhere returns elastically with the elastic tangent and
  // clears the flag.
  const PlasticResponse unload =
      j2_return(m, StressState::ThreeDimensional, 0.9 * strain, r1.state);
  REQUIRE_FALSE(unload.state.loading);
  REQUIRE((unload.tangent - m.three_dimensional_matrix()).cwiseAbs().maxCoeff() <= 1.0e-3);
}

TEST_CASE("a temperature change enters as a free strain", "[plasticity][material]") {
  IsotropicMaterial m = plastic_steel();
  m.set_thermal(1.2e-5, 20.0, 50.0);
  // Free expansion is stress-free however large.
  Vector6 free = Vector6::Zero();
  free.head(3).setConstant(1.2e-5 * 400.0);
  const PlasticResponse r = j2_return(m, StressState::ThreeDimensional, free, PlasticState(), 400.0);
  REQUIRE(r.stress.cwiseAbs().maxCoeff() <= 1.0e-6);
  REQUIRE_FALSE(r.yielding);
  // A bar held axially: sigma = -E alpha dT until it yields, then -sigma_y
  // (no hardening).
  const Scalar dt_elastic = 50.0;
  REQUIRE(uniaxial(m, 0.0, PlasticState(), dt_elastic).stress(0) ==
          Approx(-200.0e9 * 1.2e-5 * dt_elastic).epsilon(1e-12));
  const PlasticResponse hot = uniaxial(m, 0.0, PlasticState(), 150.0);
  REQUIRE(hot.yielding);
  REQUIRE(hot.stress(0) == Approx(-250.0e6).epsilon(1e-12));
  // Cooling back to the reference leaves a tensile residual stress,
  // E alpha dT - sigma_y, below yield (E alpha dT < 2 sigma_y).
  const PlasticResponse cold = uniaxial(m, 0.0, hot.state, 0.0);
  REQUIRE_FALSE(cold.yielding);
  REQUIRE(cold.stress(0) == Approx(200.0e9 * 1.2e-5 * 150.0 - 250.0e6).epsilon(1e-10));
  // Heated further instead (E alpha dT > 2 sigma_y), it yields back to +sigma_y.
  const PlasticResponse hotter = uniaxial(m, 0.0, PlasticState(), 300.0);
  const PlasticResponse cooled = uniaxial(m, 0.0, hotter.state, 0.0);
  REQUIRE(cooled.yielding);
  REQUIRE(cooled.stress(0) == Approx(250.0e6).epsilon(1e-12));
}

namespace {

/// A single-element model of each type, distorted, with a plastic material.
FemModel one_element(ElementType type, StressState state, const IsotropicMaterial& m) {
  StructuredMeshSpec spec;
  spec.nx = spec.ny = spec.nz = 2;
  spec.lx = 1.0;
  spec.ly = 0.8;
  spec.lz = 0.6;
  Mesh mesh = type == ElementType::Quad4   ? make_perturbed_quad_mesh(spec, 0.2)
              : type == ElementType::Tri3  ? make_perturbed_tri_mesh(spec, 0.2)
              : type == ElementType::Hex8  ? make_perturbed_hex_mesh(spec, 0.2)
              : type == ElementType::Tet4  ? make_perturbed_tet_mesh(spec, 0.2)
                                           : make_perturbed_tet10_mesh(spec, 0.2);
  const Scalar thickness = mesh.dim() == 2 ? 0.05 : 1.0;
  return FemModel(std::move(mesh), m, thickness, state, IntegrationOptions());
}

Vector element_displacement(const FemModel& model, Index e, Scalar amplitude, unsigned seed) {
  const int dim = model.dim();
  const Matrix x0 = model.mesh().element_coordinates(e);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  // A homogeneous stretch and shear plus a random part.
  Vector ue(dim * x0.cols());
  for (Eigen::Index a = 0; a < x0.cols(); ++a) {
    for (int k = 0; k < dim; ++k) {
      const Scalar homogeneous = k == 0 ? 0.8 * x0(0, a) + 0.3 * x0(1, a) : -0.2 * x0(k, a);
      ue(dim * a + k) = amplitude * (homogeneous + 0.5 * dist(gen));
    }
  }
  return ue;
}

}  // namespace

TEST_CASE("the elastoplastic element tangent is the derivative of its internal force",
          "[plasticity][element]") {
  const IsotropicMaterial m = plastic_steel(1.0e9, 2.0e9, 50.0e6, 40.0);
  struct Case {
    ElementType type;
    StressState state;
  };
  const std::vector<Case> cases = {
      {ElementType::Quad4, StressState::PlaneStrain}, {ElementType::Quad4, StressState::PlaneStress},
      {ElementType::Tri3, StressState::PlaneStrain},  {ElementType::Hex8, StressState::ThreeDimensional},
      {ElementType::Tet4, StressState::ThreeDimensional},
      {ElementType::Tet10, StressState::ThreeDimensional}};
  for (const Case& c : cases) {
    FemModel model = one_element(c.type, c.state, m);
    for (const bool bbar : {false, true}) {
      const Index e = 0;
      const std::vector<PlasticState> virgin(
          static_cast<std::size_t>(elastoplastic_points(model)));
      // Load into the plastic range, commit, then load on differently.
      const Vector u1 = element_displacement(model, e, 6.0e-3, 7u);
      const ElastoplasticElement first =
          elastoplastic_element(model, e, u1, virgin, bbar, nullptr, 0.0, false);
      REQUIRE(first.yielding_points > 0);
      const Vector u2 = u1 + element_displacement(model, e, 3.0e-3, 9u);
      const ElastoplasticElement base =
          elastoplastic_element(model, e, u2, first.states, bbar, nullptr, 0.0, true);
      const Scalar h = 1.0e-8;
      Scalar worst = 0.0;
      for (Eigen::Index j = 0; j < u2.size(); ++j) {
        Vector up = u2;
        Vector um = u2;
        up(j) += h;
        um(j) -= h;
        const Vector fd =
            (elastoplastic_element(model, e, up, first.states, bbar, nullptr, 0.0, false)
                 .internal_force -
             elastoplastic_element(model, e, um, first.states, bbar, nullptr, 0.0, false)
                 .internal_force) /
            (2.0 * h);
        worst = std::max(worst, (fd - base.tangent.col(j)).cwiseAbs().maxCoeff());
      }
      INFO(to_string(c.type) << " " << to_string(c.state) << (bbar ? " B-bar" : ""));
      REQUIRE(base.yielding_points > 0);
      REQUIRE(worst <= 1.0e-6 * base.tangent.cwiseAbs().maxCoeff());
    }
  }
}

TEST_CASE("the element's thermal load rate is the derivative of its internal force",
          "[plasticity][element]") {
  IsotropicMaterial m = plastic_steel(2.0e9);
  m.set_thermal(1.2e-5, 20.0, 50.0);
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    for (const StressState state :
         {StressState::PlaneStress, StressState::PlaneStrain, StressState::ThreeDimensional}) {
      if ((type == ElementType::Quad4) != (state != StressState::ThreeDimensional)) continue;
      FemModel model = one_element(type, state, m);
      Vector temperature(model.mesh().num_nodes());
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        temperature(node) = 20.0 + 300.0 * model.mesh().node(node).x() + 100.0;
      }
      const std::vector<PlasticState> virgin(
          static_cast<std::size_t>(elastoplastic_points(model)));
      const Vector ue = element_displacement(model, 0, 2.0e-3, 3u);
      const Scalar lambda = 0.8;
      // With finite kinematics the change of the Green thermal strain,
      // alpha dT (1 + alpha dT / 2), is the rate.
      for (const Kinematics kinematics : {Kinematics::SmallStrain, Kinematics::Finite}) {
        for (const bool bbar : {false, true}) {
          const ElastoplasticElement base = elastoplastic_element(
              model, 0, ue, virgin, bbar, &temperature, lambda, false, kinematics);
          REQUIRE(base.yielding_points > 0);
          const Scalar h = 1.0e-7;
          const Vector fd = (elastoplastic_element(model, 0, ue, virgin, bbar, &temperature,
                                                   lambda + h, false, kinematics)
                                 .internal_force -
                             elastoplastic_element(model, 0, ue, virgin, bbar, &temperature,
                                                   lambda - h, false, kinematics)
                                 .internal_force) /
                            (2.0 * h);
          INFO(to_string(type) << " " << to_string(state) << " " << to_string(kinematics)
                               << (bbar ? " mean dilatation" : ""));
          REQUIRE((fd - base.thermal_force_rate).cwiseAbs().maxCoeff() <=
                  1.0e-6 * base.thermal_force_rate.cwiseAbs().maxCoeff());
        }
      }
    }
  }
}

TEST_CASE("an elastic element is the linear element, and mean dilatation leaves a constant "
          "dilatation alone",
          "[plasticity][element]") {
  const IsotropicMaterial elastic(200.0e9, 0.3, 7800.0, "elastic");
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    const StressState state =
        type == ElementType::Quad4 ? StressState::PlaneStrain : StressState::ThreeDimensional;
    FemModel model = one_element(type, state, elastic);
    const std::vector<PlasticState> virgin(static_cast<std::size_t>(elastoplastic_points(model)));
    const Vector ue = element_displacement(model, 0, 1.0e-3, 5u);
    const ElastoplasticElement r =
        elastoplastic_element(model, 0, ue, virgin, false, nullptr, 0.0, true);
    const Matrix k = model.element().stiffness(model.mesh().element_coordinates(0),
                                               model.constitutive(),
                                               model.dim() == 2 ? model.thickness() : 1.0,
                                               model.integration());
    INFO(to_string(type));
    REQUIRE((r.tangent - k).cwiseAbs().maxCoeff() <= 1.0e-12 * k.cwiseAbs().maxCoeff());
    REQUIRE((r.internal_force - k * ue).cwiseAbs().maxCoeff() <=
            1.0e-12 * (k.cwiseAbs() * ue.cwiseAbs()).maxCoeff());
    REQUIRE(r.energy == Approx(0.5 * ue.dot(k * ue)).epsilon(1e-12));
    // A homogeneous strain has a constant dilatation, which B-bar keeps.
    const Matrix x0 = model.mesh().element_coordinates(0);
    Vector homogeneous(ue.size());
    for (Eigen::Index a = 0; a < x0.cols(); ++a) {
      for (int kk = 0; kk < model.dim(); ++kk) {
        homogeneous(model.dim() * a + kk) = 1.0e-3 * (kk == 0 ? x0(0, a) - 0.4 * x0(1, a)
                                                              : 0.7 * x0(kk, a));
      }
    }
    const Vector f_std =
        elastoplastic_element(model, 0, homogeneous, virgin, false, nullptr, 0.0, false)
            .internal_force;
    const Vector f_bar =
        elastoplastic_element(model, 0, homogeneous, virgin, true, nullptr, 0.0, false)
            .internal_force;
    REQUIRE((f_std - f_bar).cwiseAbs().maxCoeff() <= 1.0e-10 * f_std.cwiseAbs().maxCoeff());
  }
}

TEST_CASE("the finite-kinematics elastoplastic tangent is the derivative of its internal force",
          "[plasticity][element]") {
  // Green-Lagrange strain and second Piola-Kirchhoff stress, with the
  // geometric stiffness and (averaged) the Green-strain mean dilatation's
  // extra term, at a large displacement with a rigid turn.
  const IsotropicMaterial m = plastic_steel(1.0e9, 2.0e9, 50.0e6, 40.0);
  struct Case {
    ElementType type;
    StressState state;
  };
  const std::vector<Case> cases = {
      {ElementType::Quad4, StressState::PlaneStrain}, {ElementType::Quad4, StressState::PlaneStress},
      {ElementType::Tri3, StressState::PlaneStrain},  {ElementType::Hex8, StressState::ThreeDimensional},
      {ElementType::Tet10, StressState::ThreeDimensional}};
  for (const Case& c : cases) {
    FemModel model = one_element(c.type, c.state, m);
    const int dim = model.dim();
    const Matrix x0 = model.mesh().element_coordinates(0);
    Matrix rot = Matrix::Identity(dim, dim);
    rot(0, 0) = rot(1, 1) = std::cos(0.5);
    rot(0, 1) = -std::sin(0.5);
    rot(1, 0) = std::sin(0.5);
    Vector turn(dim * x0.cols());
    const Matrix rigid = (rot - Matrix::Identity(dim, dim)) * x0;
    for (Eigen::Index a = 0; a < x0.cols(); ++a) {
      for (int k = 0; k < dim; ++k) turn(dim * a + k) = rigid(k, a);
    }
    for (const bool averaged : {false, true}) {
      const std::vector<PlasticState> virgin(
          static_cast<std::size_t>(elastoplastic_points(model)));
      const Vector u1 = turn + element_displacement(model, 0, 6.0e-3, 7u);
      const ElastoplasticElement first = elastoplastic_element(
          model, 0, u1, virgin, averaged, nullptr, 0.0, false, Kinematics::Finite);
      REQUIRE(first.yielding_points > 0);
      const Vector u2 = u1 + element_displacement(model, 0, 3.0e-3, 9u);
      const ElastoplasticElement base = elastoplastic_element(
          model, 0, u2, first.states, averaged, nullptr, 0.0, true, Kinematics::Finite);
      const Scalar h = 1.0e-8;
      Scalar worst = 0.0;
      for (Eigen::Index j = 0; j < u2.size(); ++j) {
        Vector up = u2;
        Vector um = u2;
        up(j) += h;
        um(j) -= h;
        const Vector fd = (elastoplastic_element(model, 0, up, first.states, averaged, nullptr,
                                                 0.0, false, Kinematics::Finite)
                               .internal_force -
                           elastoplastic_element(model, 0, um, first.states, averaged, nullptr,
                                                 0.0, false, Kinematics::Finite)
                               .internal_force) /
                          (2.0 * h);
        worst = std::max(worst, (fd - base.tangent.col(j)).cwiseAbs().maxCoeff());
      }
      INFO(to_string(c.type) << " " << to_string(c.state) << (averaged ? " averaged" : ""));
      REQUIRE(base.yielding_points > 0);
      REQUIRE(worst <= 1.0e-6 * base.tangent.cwiseAbs().maxCoeff());
    }
  }
}

TEST_CASE("finite-kinematics plasticity is objective and tends to the linear stiffness",
          "[plasticity][element]") {
  const IsotropicMaterial m = plastic_steel(2.0e9, 1.0e9);
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    const StressState state =
        type == ElementType::Quad4 ? StressState::PlaneStrain : StressState::ThreeDimensional;
    FemModel model = one_element(type, state, m);
    const int dim = model.dim();
    const Matrix x0 = model.mesh().element_coordinates(0);
    const std::vector<PlasticState> virgin(static_cast<std::size_t>(elastoplastic_points(model)));
    const Vector u = element_displacement(model, 0, 5.0e-3, 3u);
    // x -> R x: the same Green strain, so the same state and energy.
    Matrix rot = Matrix::Identity(dim, dim);
    rot(0, 0) = rot(1, 1) = 0.0;
    rot(0, 1) = -1.0;
    rot(1, 0) = 1.0;
    Vector rotated(u.size());
    for (Eigen::Index a = 0; a < x0.cols(); ++a) {
      const Vector x = x0.col(a) + u.segment(dim * a, dim);
      rotated.segment(dim * a, dim) = rot * x - x0.col(a);
    }
    for (const bool averaged : {false, true}) {
      const ElastoplasticElement r1 =
          elastoplastic_element(model, 0, u, virgin, averaged, nullptr, 0.0, false,
                                Kinematics::Finite);
      const ElastoplasticElement r2 =
          elastoplastic_element(model, 0, rotated, virgin, averaged, nullptr, 0.0, false,
                                Kinematics::Finite);
      INFO(to_string(type) << (averaged ? " averaged" : ""));
      REQUIRE(r1.yielding_points > 0);
      REQUIRE(r2.yielding_points == r1.yielding_points);
      REQUIRE(r2.energy == Approx(r1.energy).epsilon(1e-10));
      for (std::size_t q = 0; q < r1.states.size(); ++q) {
        REQUIRE(r2.states[q].equivalent_plastic_strain ==
                Approx(r1.states[q].equivalent_plastic_strain).margin(1e-14));
        REQUIRE((r2.states[q].back_stress - r1.states[q].back_stress).norm() <=
                1.0e-8 * std::max(1.0, r1.states[q].back_stress.norm()));
      }
      // The internal force turns with the body: f_rotated = R f per node.
      for (Eigen::Index a = 0; a < x0.cols(); ++a) {
        const Vector expected = rot * r1.internal_force.segment(dim * a, dim);
        REQUIRE((r2.internal_force.segment(dim * a, dim) - expected).norm() <=
                1.0e-9 * r1.internal_force.cwiseAbs().maxCoeff());
      }
    }
    // At zero displacement, with the element's own operator, the tangent is
    // the linear stiffness.
    const ElastoplasticElement zero =
        elastoplastic_element(model, 0, Vector::Zero(u.size()), virgin, false, nullptr, 0.0,
                              true, Kinematics::Finite);
    const Matrix k = model.element().stiffness(x0, model.constitutive(),
                                               dim == 2 ? model.thickness() : 1.0,
                                               model.integration());
    REQUIRE((zero.tangent - k).cwiseAbs().maxCoeff() <= 1.0e-12 * k.cwiseAbs().maxCoeff());
  }
}

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

SelectorGroup box(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
                  Scalar zmin = -kInf, Scalar zmax = kInf) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  s.zmin = zmin;
  s.zmax = zmax;
  g.members.push_back(s);
  g.name = "box";
  return g;
}

/// A bar along x on a distorted mesh: x held on the left face, rigid
/// motion removed at points of it that the uniform lateral contraction
/// leaves in place.
FemModel bar(ElementType type, const IsotropicMaterial& m, Scalar length, Scalar height,
             Scalar width) {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 2;
  spec.nz = 2;
  spec.lx = length;
  spec.ly = height;
  spec.lz = width;
  Mesh mesh = type == ElementType::Quad4  ? make_perturbed_quad_mesh(spec, 0.25)
              : type == ElementType::Hex8 ? make_perturbed_hex_mesh(spec, 0.25)
                                          : make_perturbed_tet10_mesh(spec, 0.25);
  const int dim = mesh.dim();
  FemModel model(std::move(mesh), m, dim == 2 ? width : 1.0,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint left;
  left.region = box(-kInf, 0.0);
  left.fix_x = true;
  model.constraints().push_back(left);
  DisplacementConstraint origin;
  origin.region = box(-kInf, 0.0, -kInf, 0.0, -kInf, 0.0);
  origin.fix_y = true;
  origin.fix_z = dim == 3;
  model.constraints().push_back(origin);
  if (dim == 3) {
    DisplacementConstraint turn;
    turn.region = box(-kInf, 0.0, height, kInf, -kInf, 0.0);
    turn.fix_z = true;
    model.constraints().push_back(turn);
  }
  return model;
}

NonlinearMonitor right_reaction(Scalar length) {
  NonlinearMonitor mon;
  mon.name = "force";
  mon.region = box(length, kInf);
  mon.component = 0;
  mon.quantity = NonlinearMonitor::Quantity::Reaction;
  return mon;
}

}  // namespace

TEST_CASE("a bar pulled and pushed back follows the exact cyclic curve on distorted meshes",
          "[plasticity][solver]") {
  // Uniaxial stress under a prescribed end displacement along the load path
  // 0 -> 1 -> -1 with kinematic hardening: elastic, plastic, elastic
  // unloading to the Bauschinger reverse yield H_kin eps_p1 - sigma_y, and
  // reverse plastic flow. A homogeneous state is exact on any mesh, so
  // every step's end force and every element's stress must match.
  const Scalar hk = 10.0e9;
  const IsotropicMaterial m = plastic_steel(0.0, hk);
  const Scalar e = m.youngs_modulus();
  const Scalar sy = 250.0e6;
  const Scalar length = 1.0;
  const Scalar height = 0.2;
  const Scalar width = 0.1;
  const Scalar e1 = 0.01;
  const Scalar s1 = e * (sy + hk * e1) / (e + hk);
  const Scalar ep1 = e1 - s1 / e;
  const auto exact = [&](Scalar strain, bool reversed) {
    if (!reversed) return strain * e <= sy ? e * strain : e * (sy + hk * strain) / (e + hk);
    const Scalar elastic = s1 - e * (e1 - strain);
    return elastic >= hk * ep1 - sy ? elastic : e * (hk * strain - sy) / (e + hk);
  };
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    {
      const Kinematics kin = Kinematics::SmallStrain;
      FemModel model = bar(type, m, length, height, width);
      DisplacementConstraint pull;
      pull.region = box(length, kInf);
      pull.set(0, true, e1 * length);
      model.constraints().push_back(pull);
      LoadCaseSpec lc;
      lc.name = "cycle";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      NonlinearOptions options;
      options.kinematics = kin;
      options.steps = 8;
      options.load_path = {1.0, -1.0};
      options.residual_tolerance = 1.0e-11;
      options.displacement_tolerance = 1.0e-11;
      options.monitors.push_back(right_reaction(length));
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      INFO(to_string(type));
      REQUIRE(r.completed);
      REQUIRE(r.load_factor == -1.0);
      REQUIRE(r.plastic);
      const Scalar area = height * width;
      bool reversed = false;
      int reverse_yield_steps = 0;
      for (const NonlinearStep& s : r.steps) {
        if (s.load_factor == 1.0) reversed = true;
        const Scalar expected = exact(s.load_factor * e1, reversed && s.load_factor < 1.0);
        REQUIRE(s.monitors[0] == Approx(expected * area).epsilon(1e-9));
        if (reversed && s.load_factor < 1.0 && s.yielding_points > 0) ++reverse_yield_steps;
      }
      REQUIRE(reverse_yield_steps > 0);
      const Scalar final_stress = exact(-e1, true);
      for (Index el = 0; el < model.mesh().num_elements(); ++el) {
        REQUIRE(r.element_cauchy(0, el) == Approx(final_stress).epsilon(1e-9));
        REQUIRE(std::abs(r.element_cauchy(1, el)) <= 1.0e-8 * std::abs(final_stress));
      }
      REQUIRE(r.plastic_points == r.total_points);
      // eps_p at the end: the strain less the elastic part.
      REQUIRE(r.max_plastic_strain ==
              Approx(ep1 + (ep1 - (-e1 - final_stress / e))).epsilon(1e-8));
    }
  }
}

TEST_CASE("a homogeneous finite elastoplastic deformation is exact on distorted meshes",
          "[plasticity][solver]") {
  // Every boundary node carries u = lambda (F - I) X with a 3 % stretch, a
  // shear and a rigid turn of 0.1 rad. The interior must follow, and every
  // point must hold the return of the same Green strains in the solver's
  // own steps. (A larger turn is no test: lambda (F - I) runs along the
  // chord of the rotation, which compresses the body - by 3 % at mid-path
  // for 0.5 rad, a pressure of 5 GPa whose geometric stiffness outweighs the
  // plastic tangent, and the homogeneous state is rightly found unstable.)
  const IsotropicMaterial m = plastic_steel(1.0e9, 2.0e9, 60.0e6, 30.0);
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 3;
    Mesh mesh = type == ElementType::Quad4  ? make_perturbed_quad_mesh(spec, 0.25)
                : type == ElementType::Hex8 ? make_perturbed_hex_mesh(spec, 0.25)
                                            : make_perturbed_tet10_mesh(spec, 0.25);
    const int dim = mesh.dim();
    const StressState state = dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional;
    Matrix stretch = Matrix::Identity(dim, dim);
    stretch(0, 0) = 1.03;
    stretch(1, 1) = 0.985;
    stretch(0, 1) = 0.01;
    if (dim == 3) stretch(2, 2) = 1.005;
    Matrix rot = Matrix::Identity(dim, dim);
    rot(0, 0) = rot(1, 1) = std::cos(0.1);
    rot(0, 1) = -std::sin(0.1);
    rot(1, 0) = std::sin(0.1);
    const Matrix f = rot * stretch;
    FemModel model(mesh, m, 1.0, state, IntegrationOptions());
    std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
    for (const Mesh::BoundaryFace& face : model.mesh().boundary_faces()) {
      for (Index node : face.nodes) on_boundary[static_cast<std::size_t>(node)] = 1;
    }
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      if (!on_boundary[static_cast<std::size_t>(node)]) continue;
      const Vector x = model.mesh().coordinates().col(node);
      const Vector target = (f - Matrix::Identity(dim, dim)) * x;
      DisplacementConstraint bc;
      Selector s;
      s.kind = SelectorKind::NodeIds;
      s.ids.assign(1, node);
      bc.region.members.push_back(s);
      for (int k = 0; k < dim; ++k) bc.set(k, true, target(k));
      model.constraints().push_back(bc);
    }
    LoadCaseSpec lc;
    lc.name = "patch";
    lc.prescribed_displacement_only = true;
    model.load_case_specs().push_back(lc);
    model.finalize();
    Assembler assembler(model);
    NonlinearOptions options;
    options.steps = 4;
    options.residual_tolerance = 1.0e-11;
    options.displacement_tolerance = 1.0e-11;
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    INFO(to_string(type));
    REQUIRE(r.completed);
    // By default Q4 (plane strain) and Hex8 average the dilatation, Tet10 not.
    REQUIRE(r.mean_dilatation == (type != ElementType::Tet10));
    Scalar err = 0.0;
    for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
      const Vector x = model.mesh().coordinates().col(node);
      err = std::max(err, (r.displacement.segment(node * dim, dim) -
                           (f - Matrix::Identity(dim, dim)) * x)
                              .norm());
    }
    REQUIRE(err <= 1.0e-10);
    // The point history along the solver's steps.
    PlasticState point;
    PlasticResponse response;
    Matrix h;
    for (const NonlinearStep& s : r.steps) {
      h = s.load_factor * (f - Matrix::Identity(dim, dim));
      Matrix e_t = 0.5 * (h + h.transpose() + h.transpose() * h);
      Vector6 e6 = Vector6::Zero();
      e6(0) = e_t(0, 0);
      e6(1) = e_t(1, 1);
      e6(3) = 2.0 * e_t(0, 1);
      if (dim == 3) {
        e6(2) = e_t(2, 2);
        e6(4) = 2.0 * e_t(1, 2);
        e6(5) = 2.0 * e_t(2, 0);
      }
      response = j2_return(m, state, e6, point);
      point = response.state;
    }
    REQUIRE(response.yielding);
    Matrix3 f3 = Matrix3::Identity();
    f3.topLeftCorner(dim, dim) += h;
    Matrix3 s3;
    const Vector6& sv = response.stress;
    s3 << sv(0), sv(3), sv(5), sv(3), sv(1), sv(4), sv(5), sv(4), sv(2);
    const Matrix3 cauchy = f3 * s3 * f3.transpose() / f3.determinant();
    for (Index el = 0; el < model.mesh().num_elements(); ++el) {
      REQUIRE(r.element_cauchy(0, el) == Approx(cauchy(0, 0)).epsilon(1e-8));
      REQUIRE(r.element_cauchy(1, el) == Approx(cauchy(1, 1)).epsilon(1e-8));
      REQUIRE(r.element_plastic_strain(el) ==
              Approx(point.equivalent_plastic_strain).epsilon(1e-8));
    }
    REQUIRE(r.equilibrium.relative_force_error <= 1.0e-10);
  }
}

TEST_CASE("heated with finite kinematics a free body expands stress-free and a homogeneous "
          "state holds the exact return",
          "[plasticity][solver][thermal]") {
  // With finite kinematics the return takes away the Green strain of the free
  // thermal stretch 1 + alpha dT, alpha dT (1 + alpha dT / 2), so a body
  // heated on statically determinate supports takes u = alpha dT x without
  // stress; a linear thermal strain would leave its stretch 1.8e-5 short
  // here. Moved homogeneously as well, every point holds the return of the
  // same Green strain less that thermal strain - which, volumetric, sets the
  // pressure.
  const Scalar alpha = 1.2e-5;
  const Scalar t_ref = 20.0;
  const Scalar dt = 500.0;
  IsotropicMaterial m = plastic_steel(1.0e9, 2.0e9, 60.0e6, 30.0);
  m.set_thermal(alpha, t_ref, 50.0);
  const auto heat = [&]() {
    LoadCaseSpec lc;
    lc.name = "heat";
    lc.temperature.source = TemperatureSpec::Source::Uniform;
    lc.temperature.uniform = t_ref + dt;
    return lc;
  };
  const auto hold = [](const Vector& x, int dim, const Vector& value, const std::array<bool, 3>& on) {
    DisplacementConstraint bc;
    Selector s;
    s.kind = SelectorKind::NearestNode;
    s.point = Vector3::Zero();
    s.point.head(dim) = x;
    bc.region.members.push_back(s);
    for (int k = 0; k < dim; ++k) {
      if (on[static_cast<std::size_t>(k)]) bc.set(k, true, value(k));
    }
    return bc;
  };
  NonlinearOptions options;
  options.kinematics = Kinematics::Finite;
  options.steps = 4;
  options.residual_tolerance = 1.0e-11;
  options.displacement_tolerance = 1.0e-11;
  for (const ElementType type : {ElementType::Quad4, ElementType::Hex8, ElementType::Tet10}) {
    INFO(to_string(type));
    StructuredMeshSpec spec;
    spec.nx = spec.ny = spec.nz = 3;
    const auto mesh = [&]() {
      return type == ElementType::Quad4  ? make_perturbed_quad_mesh(spec, 0.25)
             : type == ElementType::Hex8 ? make_perturbed_hex_mesh(spec, 0.25)
                                         : make_perturbed_tet10_mesh(spec, 0.25);
    };
    {
      // Free expansion (plane stress for the Q4: its thickness is free too).
      FemModel model(mesh(), m, type == ElementType::Quad4 ? 0.01 : 1.0,
                     type == ElementType::Quad4 ? StressState::PlaneStress
                                                : StressState::ThreeDimensional,
                     IntegrationOptions());
      const int dim = model.dim();
      const Vector zero = Vector::Zero(dim);
      model.constraints().push_back(hold(Vector::Zero(dim), dim, zero, {true, true, true}));
      Vector corner = Vector::Zero(dim);
      corner(0) = 1.0;
      model.constraints().push_back(hold(corner, dim, zero, {false, true, true}));
      if (dim == 3) {
        corner.setZero();
        corner(1) = 1.0;
        model.constraints().push_back(hold(corner, dim, zero, {false, false, true}));
      }
      model.load_case_specs().push_back(heat());
      model.finalize();
      Assembler assembler(model);
      const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
      REQUIRE(r.completed);
      REQUIRE(r.plastic);
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        const Vector3 x = model.mesh().node(node);
        for (int k = 0; k < dim; ++k) {
          REQUIRE(r.displacement(dim * node + k) == Approx(alpha * dt * x(k)).margin(1.0e-12));
        }
      }
      REQUIRE(r.plastic_points == 0);
      REQUIRE(r.element_von_mises.cwiseAbs().maxCoeff() <= 1.0);
    }
    {
      // The heated homogeneous deformation: a 3 % stretch, a shear and a turn
      // of 0.1 rad on the boundary, the temperature rising with the load.
      FemModel model(mesh(), m, 1.0,
                     type == ElementType::Quad4 ? StressState::PlaneStrain
                                                : StressState::ThreeDimensional,
                     IntegrationOptions());
      const int dim = model.dim();
      const StressState state = model.stress_state();
      Matrix stretch = Matrix::Identity(dim, dim);
      stretch(0, 0) = 1.03;
      stretch(1, 1) = 0.985;
      stretch(0, 1) = 0.01;
      if (dim == 3) stretch(2, 2) = 1.005;
      Matrix rot = Matrix::Identity(dim, dim);
      rot(0, 0) = rot(1, 1) = std::cos(0.1);
      rot(0, 1) = -std::sin(0.1);
      rot(1, 0) = std::sin(0.1);
      const Matrix f = rot * stretch;
      std::vector<char> on_boundary(static_cast<std::size_t>(model.mesh().num_nodes()), 0);
      for (const Mesh::BoundaryFace& face : model.mesh().boundary_faces()) {
        for (Index node : face.nodes) on_boundary[static_cast<std::size_t>(node)] = 1;
      }
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        if (!on_boundary[static_cast<std::size_t>(node)]) continue;
        const Vector x = model.mesh().coordinates().col(node);
        DisplacementConstraint bc;
        Selector s;
        s.kind = SelectorKind::NodeIds;
        s.ids.assign(1, node);
        bc.region.members.push_back(s);
        const Vector target = (f - Matrix::Identity(dim, dim)) * x;
        for (int k = 0; k < dim; ++k) bc.set(k, true, target(k));
        model.constraints().push_back(bc);
      }
      model.load_case_specs().push_back(heat());
      model.finalize();
      Assembler assembler(model);
      const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
      REQUIRE(r.completed);
      Scalar err = 0.0;
      for (Index node = 0; node < model.mesh().num_nodes(); ++node) {
        const Vector x = model.mesh().coordinates().col(node);
        err = std::max(err, (r.displacement.segment(node * dim, dim) -
                             (f - Matrix::Identity(dim, dim)) * x)
                                .norm());
      }
      REQUIRE(err <= 1.0e-10);
      // The point history along the solver's steps, the temperature change
      // lambda dT entering as the Green strain of its stretch.
      PlasticState point;
      PlasticResponse response;
      Matrix h;
      for (const NonlinearStep& s : r.steps) {
        h = s.load_factor * (f - Matrix::Identity(dim, dim));
        const Matrix e_t = 0.5 * (h + h.transpose() + h.transpose() * h);
        Vector6 e6 = Vector6::Zero();
        e6(0) = e_t(0, 0);
        e6(1) = e_t(1, 1);
        e6(3) = 2.0 * e_t(0, 1);
        if (dim == 3) {
          e6(2) = e_t(2, 2);
          e6(4) = 2.0 * e_t(1, 2);
          e6(5) = 2.0 * e_t(2, 0);
        }
        const Scalar change = s.load_factor * dt;
        const Scalar green_change = ((1.0 + alpha * change) * (1.0 + alpha * change) - 1.0) /
                                    (2.0 * alpha);
        response = j2_return(m, state, e6, point, green_change);
        point = response.state;
      }
      REQUIRE(response.yielding);
      Matrix3 f3 = Matrix3::Identity();
      f3.topLeftCorner(dim, dim) += h;
      Matrix3 s3;
      const Vector6& sv = response.stress;
      s3 << sv(0), sv(3), sv(5), sv(3), sv(1), sv(4), sv(5), sv(4), sv(2);
      const Matrix3 cauchy = f3 * s3 * f3.transpose() / f3.determinant();
      for (Index el = 0; el < model.mesh().num_elements(); ++el) {
        REQUIRE(r.element_cauchy(0, el) == Approx(cauchy(0, 0)).epsilon(1e-8));
        REQUIRE(r.element_cauchy(1, el) == Approx(cauchy(1, 1)).epsilon(1e-8));
        REQUIRE(r.element_plastic_strain(el) ==
                Approx(point.equivalent_plastic_strain).epsilon(1e-8));
      }
    }
  }
}

TEST_CASE("load control stops at a plastic collapse and arc length agrees with it below",
          "[plasticity][solver]") {
  // A perfectly plastic bar under an end traction of 1.2 sigma_y collapses
  // at lambda = 1 / 1.2: beyond it there is no equilibrium.
  const IsotropicMaterial perfect = plastic_steel();
  const Scalar length = 1.0;
  const auto pulled = [&](const IsotropicMaterial& m, Scalar traction) {
    FemModel model = bar(ElementType::Quad4, m, length, 0.2, 0.1);
    LoadCaseSpec lc;
    lc.name = "pull";
    TractionLoadSpec t;
    t.region = box(length, kInf);
    t.traction = Vector3(traction, 0.0, 0.0);
    lc.tractions.push_back(t);
    model.load_case_specs().push_back(lc);
    model.finalize();
    return model;
  };
  {
    FemModel model = pulled(perfect, 1.2 * 250.0e6);
    Assembler assembler(model);
    NonlinearOptions options;
    options.kinematics = Kinematics::SmallStrain;
    NonlinearStaticAnalysis analysis(model, assembler, options);
    const NonlinearResult r = analysis.solve(0);
    REQUIRE_FALSE(r.completed);
    const Scalar collapse = 1.0 / 1.2;
    REQUIRE(r.load_factor <= collapse * (1.0 + 1e-12));
    // The collapse shows as the tangent losing its definiteness (singular
    // at collapse, round-off tips it) or as steps that cannot converge;
    // either way the run brackets it.
    const Scalar bound =
        std::isnan(r.critical_bound) ? r.unreached_load_factor : r.critical_bound;
    REQUIRE(bound >= collapse);
    REQUIRE(bound - r.load_factor <= 0.01);
    REQUIRE(r.termination.find("plastic collapse") != std::string::npos);
    REQUIRE(r.steps.size() < 100);
  }
  // With hardening the same bar carries the load; arc length and load
  // control end in the same (path-independent, uniaxial) state.
  const IsotropicMaterial hardening = plastic_steel(3.0e9);
  FemModel model = pulled(hardening, 1.2 * 250.0e6);
  Assembler assembler(model);
  NonlinearOptions options;
  options.kinematics = Kinematics::SmallStrain;
  options.residual_tolerance = 1.0e-11;
  options.displacement_tolerance = 1.0e-11;
  const NonlinearResult lc = NonlinearStaticAnalysis(model, assembler, options).solve(0);
  options.method = NonlinearOptions::Method::ArcLength;
  const NonlinearResult al = NonlinearStaticAnalysis(model, assembler, options).solve(0);
  REQUIRE(lc.completed);
  REQUIRE(al.completed);
  REQUIRE((lc.displacement - al.displacement).norm() <= 1.0e-9 * lc.displacement.norm());
  // sigma = 1.2 sigma_y, eps_p = (sigma - sigma_y) / H.
  REQUIRE(lc.max_plastic_strain == Approx(0.2 * 250.0e6 / 3.0e9).epsilon(1e-9));
}

TEST_CASE("small-strain kinematics with elastic materials is the linear analysis",
          "[plasticity][solver]") {
  CantileverCase c;
  c.poisson = 0.3;
  c.tip_load = -3.0e5;  // a large deflection, which small strain must ignore
  FemModel model = make_cantilever(c, 20, 4);
  Assembler assembler(model);
  StaticAnalysis linear(model, assembler);
  const Vector ul = linear.solve_all().front().displacement;
  NonlinearOptions options;
  options.kinematics = Kinematics::SmallStrain;
  options.steps = 3;
  options.residual_tolerance = 1.0e-12;
  options.displacement_tolerance = 1.0e-12;
  const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, options).solve(0);
  REQUIRE(r.completed);
  REQUIRE_FALSE(r.plastic);
  REQUIRE((r.displacement - ul).norm() <= 1.0e-10 * ul.norm());
  for (const NonlinearStep& s : r.steps) REQUIRE(s.iterations == 1);
  // It reports what it neglects: a rotation large enough to warn about.
  REQUIRE(r.max_rotation > 0.05);
  REQUIRE_FALSE(r.warnings.empty());
  // The neo-Hookean law needs finite kinematics.
  options.law = HyperelasticModel::NeoHookean;
  REQUIRE_THROWS_AS(NonlinearStaticAnalysis(model, assembler, options), ConfigError);
}

TEST_CASE("the plasticity block and the non-linear keys of a deck parse and validate",
          "[plasticity][config]") {
  // A plane-strain strip 1 m x 0.1 m of plastic steel under a tip force.
  const auto deck = [](const std::string& plasticity, const std::string& nonlinear) {
    return json::parse(R"({
      "mesh": { "type": "structured_quad", "nx": 20, "ny": 2, "lx": 1.0, "ly": 0.1 },
      "material": { "youngs_modulus": 200e9, "poisson_ratio": 0.3,
                    "plasticity": )" + plasticity + R"( },
      "model": { "thickness": 0.01, "stress_state": "plane_strain" },
      "boundary_conditions": [ { "fix": ["x", "y"], "region": { "box": { "xmax": 0.0 } } } ],
      "load_cases": [ { "name": "tip", "point_loads": [ { "force": [0.0, -8000.0],
            "region": { "box": { "xmin": 1.0 } } } ] } ],
      "nonlinear": )" + nonlinear + "}");
  };
  const std::string steel = R"({ "yield_stress": 250e6, "hardening_modulus": 2e9,
      "kinematic_hardening_modulus": 1e9, "saturation_stress": 80e6, "saturation_rate": 20 })";
  const Configuration config = parse_configuration(
      deck(steel, R"({ "enabled": true, "kinematics": "small_strain", "mean_dilatation": "all",
                       "steps": 4, "load_path": [1.0, 0.0] })"),
      "inline", true);
  const PlasticityParameters& p = config.material().plasticity();
  REQUIRE(p.yield_stress == 250.0e6);
  REQUIRE(p.hardening_modulus == 2.0e9);
  REQUIRE(p.kinematic_hardening_modulus == 1.0e9);
  REQUIRE(p.saturation_stress == 80.0e6);
  REQUIRE(p.saturation_rate == 20.0);
  const NonlinearOptions& o = config.nonlinear.options;
  REQUIRE(o.kinematics == Kinematics::SmallStrain);
  REQUIRE(o.mean_dilatation == MeanDilatation::All);
  REQUIRE(o.load_path == std::vector<Scalar>{1.0, 0.0});
  // Booleans stand for "all" and "none"; the default is "auto".
  REQUIRE(parse_configuration(deck(steel, R"({ "enabled": true, "mean_dilatation": false })"),
                              "inline", true)
              .nonlinear.options.mean_dilatation == MeanDilatation::None);
  REQUIRE(parse_configuration(deck(steel, R"({ "enabled": true })"), "inline", true)
              .nonlinear.options.mean_dilatation == MeanDilatation::Auto);

  // The run loads the strip past yield and unloads it: a permanent set,
  // plastic points, and the load path's end.
  FemModel model = build_model(config);
  Assembler assembler(model);
  const NonlinearResult r = NonlinearStaticAnalysis(model, assembler, o).solve(0);
  REQUIRE(r.completed);
  REQUIRE(r.load_factor == 0.0);
  REQUIRE(r.plastic);
  REQUIRE(r.plastic_points > 0);
  REQUIRE(r.max_plastic_strain > 0.0);
  REQUIRE(r.displacement.cwiseAbs().maxCoeff() > 0.0);  // the permanent set
  REQUIRE(r.steps.size() == 8);

  const auto refuses = [&](const std::string& plasticity, const std::string& nonlinear) {
    REQUIRE_THROWS_AS(parse_configuration(deck(plasticity, nonlinear), "inline", true),
                      ConfigError);
  };
  const std::string on = R"({ "enabled": true })";
  refuses(R"({ "yield_stress": -1.0 })", on);
  refuses(R"({ "hardening_modulus": 1e9 })", on);  // no yield stress
  refuses(R"({ "yield_stress": 250e6, "saturation_stress": 50e6 })", on);  // no rate
  refuses(R"({ "yield_stress": 250e6, "hardening": 1e9 })", on);  // unknown key
  refuses(steel, R"({ "enabled": true, "kinematics": "moderate" })");
  refuses(steel, R"({ "enabled": true, "mean_dilatation": "some" })");
  refuses(steel, R"({ "enabled": true, "load_path": [0.0, 1.0] })");   // starts where it is
  refuses(steel, R"({ "enabled": true, "load_path": [1.0, 1.0] })");   // no leg
  refuses(steel, R"({ "enabled": true, "method": "arc_length", "load_path": [1.0, 0.0] })");
  refuses(steel, R"({ "enabled": true, "load_path": [1.0, 0.0], "load_factors": [0.5] })");
  refuses(steel, R"({ "enabled": true, "kinematics": "small_strain",
                      "follower_pressure": true })");
  refuses(steel, R"({ "enabled": true, "kinematics": "small_strain",
                      "material_model": "neo_hookean" })");
  refuses(steel, R"({ "enabled": true, "material_model": "neo_hookean" })");  // plastic
}

TEST_CASE("the CalculiX export writes the hardening curve and one step per leg",
          "[plasticity][io]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 2;
  spec.nz = 2;
  spec.lx = 0.4;
  spec.ly = 0.1;
  spec.lz = 0.1;
  IsotropicMaterial steel(200.0e9, 0.3, 7800.0, "steel");
  PlasticityParameters p;
  p.yield_stress = 250.0e6;
  p.hardening_modulus = 2.0e9;
  p.saturation_stress = 100.0e6;
  p.saturation_rate = 25.0;
  steel.set_plasticity(p);
  FemModel model(make_structured_hex_mesh(spec), steel, 1.0, StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint root;
  root.region = box(-kInf, 0.0);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);
  DisplacementConstraint tip;
  tip.region = box(0.4, kInf);
  tip.set(1, true, -0.004);
  model.constraints().push_back(tip);
  LoadCaseSpec lc;
  lc.name = "bend";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();

  CalculixNonlinearExport nl;
  nl.load_cases = {0};
  nl.increments = 6;
  nl.nlgeom = false;
  nl.load_path = {1.0, -0.5};
  ensure_directory("results/_test_tmp");
  const std::vector<std::string> decks =
      write_calculix_decks(model, "results/_test_tmp/plastic", "unit", &nl);
  REQUIRE(decks.size() == 2);
  REQUIRE(decks[1].find("_small_strain.inp") != std::string::npos);
  std::ifstream in(decks[1]);
  std::stringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();
  // Two steps without NLGEOM, each of 6 fixed increments.
  REQUIRE(text.find("NLGEOM") == std::string::npos);
  std::size_t steps = 0;
  for (std::size_t at = text.find("*STATIC, DIRECT"); at != std::string::npos;
       at = text.find("*STATIC, DIRECT", at + 1)) {
    ++steps;
  }
  REQUIRE(steps == 2);
  // The second leg drives the tip to -0.5 of the prescribed value.
  REQUIRE(text.find(", 2, 2, 0.002") != std::string::npos);
  // The hardening table: exact nodes of sigma_y(alpha), from alpha = 0, up
  // to 10, with chords within 1e-4 Q of the curve.
  const std::size_t start = text.find("*PLASTIC\n");
  REQUIRE(start != std::string::npos);
  std::istringstream table(text.substr(start + 9));
  std::string line;
  std::vector<std::pair<Scalar, Scalar>> points;
  while (std::getline(table, line) && line[0] != '*') {
    const std::size_t comma = line.find(',');
    points.emplace_back(std::stod(line.substr(0, comma)), std::stod(line.substr(comma + 1)));
  }
  REQUIRE(points.size() > 10);
  REQUIRE(points.front().second == 0.0);
  REQUIRE(points.back().second == Approx(10.0));
  Scalar worst = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    REQUIRE(points[i].first == Approx(steel.plasticity().yield(points[i].second)).epsilon(1e-9));
    if (i == 0) continue;
    // The chord's largest departure, sampled.
    for (int k = 1; k < 10; ++k) {
      const Scalar s = k / 10.0;
      const Scalar a = points[i - 1].second + s * (points[i].second - points[i - 1].second);
      const Scalar chord = points[i - 1].first + s * (points[i].first - points[i - 1].first);
      worst = std::max(worst, std::abs(chord - steel.plasticity().yield(a)));
    }
  }
  REQUIRE(worst <= 1.0e-4 * p.saturation_stress);
}
