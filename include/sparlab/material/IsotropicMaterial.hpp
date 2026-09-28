/// \file IsotropicMaterial.hpp
/// \brief Linear-elastic isotropic material and its constitutive matrices.
///
/// Voigt ordering used throughout SparLab is
/// \f$ \{\sigma_{xx}, \sigma_{yy}, \sigma_{xy}\} \f$ in 2-D and
/// \f$ \{\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz},
/// \sigma_{zx}\} \f$ in 3-D, paired with the *engineering* strain vector
/// (\f$\gamma_{xy} = 2\varepsilon_{xy}\f$ and so on).
///
/// Plane stress (\f$\sigma_{zz}=0\f$):
/// \f[
///   \mathbf{D} = \frac{E}{1-\nu^2}
///   \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \tfrac{1-\nu}{2}
///   \end{bmatrix}
/// \f]
/// Plane strain (\f$\varepsilon_{zz}=0\f$):
/// \f[
///   \mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)}
///   \begin{bmatrix} 1-\nu & \nu & 0 \\ \nu & 1-\nu & 0 \\ 0 & 0 &
///   \tfrac{1-2\nu}{2}\end{bmatrix}
/// \f]
/// In plane strain the out-of-plane stress is
/// \f$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy})\f$, which StressRecovery
/// accounts for when forming the von Mises stress.
///
/// **Thermal strain.** A temperature change \f$\Delta T = T - T_{ref}\f$ gives
/// the isotropic free strain \f$\alpha\,\Delta T\f$. Its Voigt form depends on
/// the idealisation: \f$\alpha\Delta T\{1, 1, 0\}\f$ in plane stress, where
/// the free out-of-plane expansion leaves the in-plane law unchanged;
/// \f$(1+\nu)\,\alpha\Delta T\{1, 1, 0\}\f$ in plane strain, where the
/// restrained out-of-plane expansion \f$\varepsilon_{zz} = 0\f$ adds
/// \f$\nu\alpha\Delta T\f$ in each in-plane direction and
/// \f$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy}) - E\alpha\Delta T\f$; and
/// \f$\alpha\Delta T\{1,1,1,0,0,0\}\f$ in 3-D. The stress is
/// \f$\sigma = D(\varepsilon - \varepsilon_0)\f$.
///
/// Three-dimensional elasticity, with Lame constants
/// \f$\lambda = E\nu/((1+\nu)(1-2\nu))\f$ and \f$G = E/(2(1+\nu))\f$:
/// \f[
///   \mathbf{D} = \begin{bmatrix}
///     \lambda+2G & \lambda & \lambda & 0 & 0 & 0
///     \\ \lambda & \lambda+2G & \lambda & 0 & 0 & 0
///     \\ \lambda & \lambda & \lambda+2G & 0 & 0 & 0
///     \\ 0 & 0 & 0 & G & 0 & 0 \\ 0 & 0 & 0 & 0 & G & 0 \\ 0 & 0 & 0 & 0 & 0 & G
///   \end{bmatrix}
/// \f]
#pragma once

#include "sparlab/core/Types.hpp"

#include <string>

namespace sparlab {

/// Isotropic linear-elastic material described by \f$(E, \nu, \rho)\f$.
class IsotropicMaterial {
 public:
  /// \param youngs_modulus E [Pa], must be > 0.
  /// \param poisson_ratio nu [-], must lie in (-1, 0.5) for plane strain and
  ///        3-D elasticity and (-1, 1) excluding +-1 for plane stress. SparLab
  ///        restricts the input to (-1, 0.5) so every idealisation remains
  ///        well posed.
  /// \param density rho [kg/m^3], must be >= 0 (0 disables modal analysis).
  /// \param name optional label carried into result files.
  IsotropicMaterial(Scalar youngs_modulus, Scalar poisson_ratio, Scalar density,
                    std::string name = "material");

  Scalar youngs_modulus() const { return e_; }
  Scalar poisson_ratio() const { return nu_; }
  Scalar density() const { return rho_; }
  const std::string& name() const { return name_; }

  /// Shear modulus \f$G = E / (2(1+\nu))\f$ [Pa].
  Scalar shear_modulus() const { return e_ / (2.0 * (1.0 + nu_)); }

  /// First Lame constant \f$\lambda = E\nu/((1+\nu)(1-2\nu))\f$ [Pa].
  Scalar lame_lambda() const { return e_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_)); }

  /// Constitutive matrix for the requested idealisation [Pa]: 3 x 3 for the
  /// plane states, 6 x 6 for three-dimensional elasticity.
  Matrix constitutive(StressState state) const;

  /// Plane-stress constitutive matrix [Pa].
  Matrix3 plane_stress_matrix() const;

  /// Plane-strain constitutive matrix [Pa].
  Matrix3 plane_strain_matrix() const;

  /// Three-dimensional constitutive matrix [Pa].
  Matrix6 three_dimensional_matrix() const;

  /// Return a copy of this material with a scaled Young's modulus. Used by the
  /// aerospace material-stiffness sweep.
  IsotropicMaterial with_youngs_modulus(Scalar e) const {
    IsotropicMaterial copy(e, nu_, rho_, name_);  // validates e
    copy.alpha_ = alpha_;
    copy.t_ref_ = t_ref_;
    copy.k_ = k_;
    return copy;
  }

  /// Linear thermal expansion coefficient alpha [1/K]; 0 means no thermal
  /// strain.
  Scalar thermal_expansion() const { return alpha_; }
  /// Temperature at which the material carries no thermal strain [K, or
  /// deg C when every temperature of the deck is in deg C].
  Scalar reference_temperature() const { return t_ref_; }
  /// Thermal conductivity k [W/(m K)], used by the conduction solve.
  Scalar conductivity() const { return k_; }
  /// Set the thermal properties.
  /// 	hrows ConfigError for a negative conductivity or non-finite input.
  void set_thermal(Scalar expansion, Scalar reference_temperature, Scalar conductivity);

  /// Voigt thermal strain \f$\varepsilon_0\f$ of a temperature change for
  /// an idealisation (see the file comment): 3 components for the plane
  /// states, 6 in 3-D.
  Vector thermal_strain(StressState state, Scalar delta_t) const;

  /// Out-of-plane stress of a plane-strain state with in-plane stresses
  /// sx, sy at temperature change delta_t:
  /// \f$\nu(\sigma_{xx}+\sigma_{yy}) - E\alpha\Delta T\f$.
  Scalar plane_strain_sigma_zz(Scalar sx, Scalar sy, Scalar delta_t) const {
    return nu_ * (sx + sy) - e_ * alpha_ * delta_t;
  }

 private:
  Scalar e_;
  Scalar nu_;
  Scalar rho_;
  std::string name_;
  Scalar alpha_ = 0.0;
  Scalar t_ref_ = 0.0;
  Scalar k_ = 0.0;
};

}  // namespace sparlab
