#include "sparlab/material/IsotropicMaterial.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <cmath>
#include <sstream>
#include <utility>

namespace sparlab {

IsotropicMaterial::IsotropicMaterial(Scalar youngs_modulus, Scalar poisson_ratio,
                                     Scalar density, std::string name)
    : e_(youngs_modulus), nu_(poisson_ratio), rho_(density), name_(std::move(name)) {
  if (!(e_ > 0.0)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a positive Young's modulus (got " << e_
       << " Pa)";
    throw ConfigError(os.str());
  }
  if (!(nu_ > -1.0 && nu_ < 0.5)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a Poisson ratio in (-1, 0.5) for a "
       << "positive-definite constitutive matrix (got " << nu_ << ")";
    throw ConfigError(os.str());
  }
  if (!(rho_ >= 0.0)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a non-negative density (got " << rho_
       << " kg/m^3)";
    throw ConfigError(os.str());
  }
}

Matrix3 IsotropicMaterial::plane_stress_matrix() const {
  Matrix3 d = Matrix3::Zero();
  const Scalar f = e_ / (1.0 - nu_ * nu_);
  d(0, 0) = f;
  d(0, 1) = f * nu_;
  d(1, 0) = f * nu_;
  d(1, 1) = f;
  d(2, 2) = f * (1.0 - nu_) / 2.0;
  return d;
}

Matrix3 IsotropicMaterial::plane_strain_matrix() const {
  Matrix3 d = Matrix3::Zero();
  const Scalar f = e_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
  d(0, 0) = f * (1.0 - nu_);
  d(0, 1) = f * nu_;
  d(1, 0) = f * nu_;
  d(1, 1) = f * (1.0 - nu_);
  d(2, 2) = f * (1.0 - 2.0 * nu_) / 2.0;
  return d;
}

Matrix6 IsotropicMaterial::three_dimensional_matrix() const {
  Matrix6 d = Matrix6::Zero();
  const Scalar lambda = lame_lambda();
  const Scalar g = shear_modulus();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) d(i, j) = lambda;
    d(i, i) = lambda + 2.0 * g;
    d(3 + i, 3 + i) = g;
  }
  return d;
}

void IsotropicMaterial::set_thermal(Scalar expansion, Scalar reference_temperature,
                                    Scalar conductivity) {
  if (!std::isfinite(expansion) || !std::isfinite(reference_temperature) ||
      !std::isfinite(conductivity)) {
    throw ConfigError("material '" + name_ + "' has a non-finite thermal property");
  }
  if (conductivity < 0.0) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a non-negative conductivity (got "
       << conductivity << " W/(m K))";
    throw ConfigError(os.str());
  }
  alpha_ = expansion;
  t_ref_ = reference_temperature;
  k_ = conductivity;
}

Vector IsotropicMaterial::thermal_strain(StressState state, Scalar delta_t) const {
  const Scalar e0 = alpha_ * delta_t;
  switch (state) {
    case StressState::PlaneStress: {
      Vector v(3);
      v << e0, e0, 0.0;
      return v;
    }
    case StressState::PlaneStrain: {
      // eps_zz = 0 restrains the out-of-plane expansion, which reappears as
      // nu alpha dT in each in-plane direction.
      Vector v(3);
      v << (1.0 + nu_) * e0, (1.0 + nu_) * e0, 0.0;
      return v;
    }
    case StressState::ThreeDimensional: {
      Vector v = Vector::Zero(6);
      v.head(3).setConstant(e0);
      return v;
    }
  }
  throw ConfigError("unhandled stress state");
}

Matrix IsotropicMaterial::constitutive(StressState state) const {
  switch (state) {
    case StressState::PlaneStress: return plane_stress_matrix();
    case StressState::PlaneStrain: return plane_strain_matrix();
    case StressState::ThreeDimensional: return three_dimensional_matrix();
  }
  throw ConfigError("unhandled stress state");
}

}  // namespace sparlab
