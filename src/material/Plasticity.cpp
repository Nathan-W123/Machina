#include "sparlab/material/Plasticity.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace sparlab {
namespace {

const Scalar kSqrt23 = std::sqrt(2.0 / 3.0);

/// The norm of a symmetric tensor given by its tensorial Voigt components.
Scalar tensor_norm(const Vector6& v) {
  return std::sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2) +
                   2.0 * (v(3) * v(3) + v(4) * v(4) + v(5) * v(5)));
}

/// The deviatoric stress 2G dev(eps) of an engineering-shear strain, as
/// tensorial Voigt components.
Vector6 deviatoric_stress(Scalar g, const Vector6& strain) {
  const Scalar mean = (strain(0) + strain(1) + strain(2)) / 3.0;
  Vector6 s;
  for (int i = 0; i < 3; ++i) s(i) = 2.0 * g * (strain(i) - mean);
  for (int i = 3; i < 6; ++i) s(i) = g * strain(i);  // 2G (gamma / 2)
  return s;
}

/// The squared norm of the deviator of an engineering-shear strain.
Scalar deviatoric_norm2(const Vector6& strain) {
  const Scalar mean = (strain(0) + strain(1) + strain(2)) / 3.0;
  Scalar out = 0.0;
  for (int i = 0; i < 3; ++i) out += (strain(i) - mean) * (strain(i) - mean);
  for (int i = 3; i < 6; ++i) out += 0.5 * strain(i) * strain(i);  // 2 (gamma / 2)^2
  return out;
}

/// The squared tensor norm of an engineering-shear strain.
Scalar strain_norm2(const Vector6& strain) {
  Scalar out = 0.0;
  for (int i = 0; i < 3; ++i) out += strain(i) * strain(i);
  for (int i = 3; i < 6; ++i) out += 0.5 * strain(i) * strain(i);
  return out;
}

/// The isotropic elasticity tensor in Voigt form: K m m^T + 2G (I - m m^T/3)
/// with G on the engineering-shear diagonal.
Matrix6 elastic_tangent(Scalar k, Scalar g) {
  Matrix6 c = Matrix6::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) c(i, j) = k - 2.0 * g / 3.0;
    c(i, i) = k + 4.0 * g / 3.0;
    c(3 + i, 3 + i) = g;
  }
  return c;
}

/// The 3-D return at a fully given strain.
PlasticResponse return_3d(const IsotropicMaterial& m, const Vector6& strain,
                          const PlasticState& committed, Scalar delta_t) {
  const PlasticityParameters& p = m.plasticity();
  const Scalar g = m.shear_modulus();
  const Scalar k = m.lame_lambda() + 2.0 * g / 3.0;
  const Scalar hk = p.kinematic_hardening_modulus;

  Vector6 elastic = strain - committed.plastic_strain;
  const Scalar thermal = m.thermal_expansion() * delta_t;
  if (thermal != 0.0) elastic.head(3).array() -= thermal;
  const Scalar volumetric = elastic(0) + elastic(1) + elastic(2);
  const Vector6 s_trial = deviatoric_stress(g, elastic);
  const Vector6 xi_trial = s_trial - committed.back_stress;
  const Scalar xi_norm = tensor_norm(xi_trial);
  const Scalar alpha_n = committed.equivalent_plastic_strain;
  const Scalar radius_n = kSqrt23 * p.yield(alpha_n);
  const Scalar f_trial = xi_norm - radius_n;

  PlasticResponse out;
  out.state = committed;
  out.state.loading = false;
  out.strain_33 = strain(2);
  // Round-off in the trial stress must not start plastic flow on the yield
  // surface itself.
  if (!p.enabled() || f_trial <= 1.0e-12 * radius_n) {
    out.stress = s_trial;
    out.stress.head(3).array() += k * volumetric;
    out.tangent = elastic_tangent(k, g);
    out.energy = 0.5 * k * volumetric * volumetric + g * deviatoric_norm2(elastic) +
                 p.isotropic_energy(alpha_n) + hk / 3.0 * strain_norm2(committed.plastic_strain);
    // Still on the surface after a plastic step: assume further loading
    // (the continuum tangent, Delta gamma = 0) - see PlasticState::loading.
    if (p.enabled() && committed.loading && f_trial > -1.0e-10 * radius_n && xi_norm > 0.0) {
      const Vector6 n = xi_trial / xi_norm;
      out.tangent -= 2.0 * g / (1.0 + (p.yield_slope(alpha_n) + hk) / (3.0 * g)) *
                     (n * n.transpose());
      out.state.loading = true;
    }
    return out;
  }

  // The consistency condition for the plastic multiplier.
  const Scalar stiffness = 2.0 * g + 2.0 / 3.0 * hk;
  Scalar dgamma = f_trial / (stiffness + 2.0 / 3.0 * p.yield_slope(alpha_n));
  if (p.saturation_stress != 0.0) {
    // Voce: g(dgamma) is convex and decreasing, so Newton from the left of the
    // root rises monotonically onto it.
    dgamma = 0.0;
    bool converged = false;
    for (int it = 0; it < 60; ++it) {
      const Scalar alpha = alpha_n + kSqrt23 * dgamma;
      const Scalar residual = xi_norm - stiffness * dgamma - kSqrt23 * p.yield(alpha);
      if (std::abs(residual) <= 1.0e-13 * radius_n) {
        converged = true;
        break;
      }
      dgamma += residual / (stiffness + 2.0 / 3.0 * p.yield_slope(alpha));
    }
    if (!converged) {
      std::ostringstream os;
      os << "the J2 return of material '" << m.name()
         << "' did not converge for the plastic multiplier";
      throw SolverError(os.str());
    }
  }
  const Scalar alpha = alpha_n + kSqrt23 * dgamma;
  const Vector6 n = xi_trial / xi_norm;

  out.yielding = true;
  out.stress = s_trial - 2.0 * g * dgamma * n;
  out.stress.head(3).array() += k * volumetric;
  PlasticState& st = out.state;
  st.loading = true;
  st.equivalent_plastic_strain = alpha;
  st.back_stress += 2.0 / 3.0 * hk * dgamma * n;
  for (int i = 0; i < 3; ++i) st.plastic_strain(i) += dgamma * n(i);
  for (int i = 3; i < 6; ++i) st.plastic_strain(i) += 2.0 * dgamma * n(i);  // engineering shear

  const Scalar theta = 1.0 - 2.0 * g * dgamma / xi_norm;
  const Scalar theta_bar = 1.0 / (1.0 + (p.yield_slope(alpha) + hk) / (3.0 * g)) - (1.0 - theta);
  Matrix6 c = Matrix6::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) c(i, j) = k - 2.0 * g * theta / 3.0;
    c(i, i) = k + 4.0 * g * theta / 3.0;
    c(3 + i, 3 + i) = g * theta;
  }
  c -= 2.0 * g * theta_bar * (n * n.transpose());
  out.tangent = c;

  Vector6 elastic_new = elastic;
  for (int i = 0; i < 3; ++i) elastic_new(i) -= dgamma * n(i);
  for (int i = 3; i < 6; ++i) elastic_new(i) -= 2.0 * dgamma * n(i);
  out.energy = 0.5 * k * volumetric * volumetric + g * deviatoric_norm2(elastic_new) +
               p.isotropic_energy(alpha) + hk / 3.0 * strain_norm2(st.plastic_strain);
  return out;
}

}  // namespace

PlasticResponse j2_return(const IsotropicMaterial& material, StressState state,
                          const Vector6& strain, const PlasticState& committed, Scalar delta_t) {
  if (state != StressState::PlaneStress) return return_3d(material, strain, committed, delta_t);

  // Plane stress: eps_33 such that sigma_33 = 0. The elastic predictor of
  // eps_33 is exact for an elastic step; a plastic one converges
  // quadratically with the consistent C_33,33 > 0.
  const Scalar g = material.shear_modulus();
  const Scalar lambda = material.lame_lambda();
  const Scalar thermal = material.thermal_expansion() * delta_t;
  const Vector6& ep = committed.plastic_strain;
  Vector6 trial = strain;
  trial(2) = ep(2) + thermal -
             lambda / (lambda + 2.0 * g) *
                 ((strain(0) - ep(0) - thermal) + (strain(1) - ep(1) - thermal));
  PlasticResponse out = return_3d(material, trial, committed, delta_t);
  const Scalar scale =
      std::max(out.stress.head(3).cwiseAbs().maxCoeff(), material.plasticity().yield_stress);
  int it = 0;
  while (std::abs(out.stress(2)) > 1.0e-12 * std::max(scale, std::numeric_limits<Scalar>::min())) {
    if (++it > 30) {
      std::ostringstream os;
      os << "the plane-stress J2 return of material '" << material.name()
         << "' did not reach sigma_33 = 0 (left " << out.stress(2) << " Pa)";
      throw SolverError(os.str());
    }
    trial(2) -= out.stress(2) / out.tangent(2, 2);
    out = return_3d(material, trial, committed, delta_t);
  }
  // Condense eps_33 out: d sigma_ab = (C_ab - C_a3 C_3b / C_33) d eps_b.
  const Matrix6 c = out.tangent;
  out.tangent = c - c.col(2) * c.row(2) / c(2, 2);
  out.tangent.row(2).setZero();
  out.tangent.col(2).setZero();
  out.state.thickness_strain = trial(2);
  out.strain_33 = trial(2);
  return out;
}

Scalar von_mises_stress(const Vector6& s) {
  const Scalar a = s(0) - s(1);
  const Scalar b = s(1) - s(2);
  const Scalar c = s(2) - s(0);
  return std::sqrt(0.5 * (a * a + b * b + c * c) +
                   3.0 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5)));
}

}  // namespace sparlab
