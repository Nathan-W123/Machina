#include "sparlab/fem/Elastoplastic.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/TotalLagrangian.hpp"
#include "sparlab/material/Hyperelastic.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace sparlab {
namespace {

/// An operator or strain in the 3-D Voigt order {11, 22, 33, 12, 23, 31}: a
/// plane element's rows {11, 22, 12} go to rows 0, 1 and 3.
Matrix voigt6(const Matrix& b, int dim) {
  if (dim == 3) return b;
  Matrix out = Matrix::Zero(6, b.cols());
  out.row(0) = b.row(0);
  out.row(1) = b.row(1);
  out.row(3) = b.row(2);
  return out;
}

/// The dilatation row of a strain operator: the sum of its normal rows.
Vector dilatation_row(const Matrix& b6) {
  return (b6.row(0) + b6.row(1) + b6.row(2)).transpose();
}

/// Everything a point needs from the kinematics.
struct Point {
  Matrix b;              ///< 6 x nd: B, B_NL, or their mean-dilatation forms
  Vector6 strain;        ///< epsilon or E, likewise
  Matrix g;              ///< reference gradients dN/dX, dim x n
  Matrix h;              ///< displacement gradient, dim x dim
  Scalar weight = 0.0;   ///< w detJ t
  Scalar delta_t = 0.0;  ///< temperature change per unit load factor
};

struct Kinematic {
  std::vector<Point> points;
  bool averaged = false;  ///< mean dilatation applied
  Matrix gbar;            ///< finite + averaged: the volume average of G^T G, n x n
};

Kinematic kinematics_of(const FemModel& model, Index e, const Vector& ue, Kinematics kinematics,
                        bool mean_dilatation, const Vector* temperature) {
  const Mesh& mesh = model.mesh();
  const Element& element = model.element();
  const int dim = mesh.dim();
  const int n = mesh.nodes_per_elem();
  const Scalar t = dim == 2 ? model.thickness() : 1.0;
  const Matrix x0 = mesh.element_coordinates(e);
  const IsotropicMaterial& mat = model.material_of(e);
  const bool thermal = temperature != nullptr && temperature->size() > 0 &&
                       mat.thermal_expansion() != 0.0;
  Vector te;
  if (thermal) {
    te.resize(n);
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < n; ++a) te(a) = (*temperature)(nodes[a]);
  }
  Matrix u(dim, n);
  for (int a = 0; a < n; ++a) {
    for (int k = 0; k < dim; ++k) u(k, a) = ue(dim * a + k);
  }
  const bool finite = kinematics == Kinematics::Finite;

  Kinematic out;
  for (const IntegrationPoint& ip : element.integration_rule(model.integration())) {
    const StrainOperator op = element.strain_operator(x0, ip.point);
    Point p;
    p.g = reference_gradients(op, dim, n);
    p.h = u * p.g.transpose();
    if (finite) {
      const Matrix f = Matrix::Identity(dim, dim) + p.h;
      p.b = voigt6(green_lagrange_operator(f, p.g), dim);
      p.strain = voigt6(green_lagrange_voigt(p.h), dim);
    } else {
      p.b = voigt6(op.b, dim);
      p.strain = p.b * ue;
    }
    p.weight = ip.weight * op.detJ * t;
    if (thermal) {
      p.delta_t = element.shape_functions(ip.point).dot(te) - mat.reference_temperature();
    }
    out.points.push_back(std::move(p));
  }
  // Mean dilatation: every point's dilatation replaced by the element's
  // volume average, in the strain and in its variation.
  if (mean_dilatation && model.stress_state() != StressState::PlaneStress &&
      out.points.size() > 1) {
    out.averaged = true;
    Vector mean_row = Vector::Zero(dim * n);
    Scalar mean_trace = 0.0;
    Scalar volume = 0.0;
    if (finite) out.gbar = Matrix::Zero(n, n);
    for (const Point& p : out.points) {
      mean_row += p.weight * dilatation_row(p.b);
      mean_trace += p.weight * (p.strain(0) + p.strain(1) + p.strain(2));
      volume += p.weight;
      if (finite) out.gbar += p.weight * (p.g.transpose() * p.g);
    }
    mean_row /= volume;
    mean_trace /= volume;
    if (finite) out.gbar /= volume;
    for (Point& p : out.points) {
      const Vector correction = (mean_row - dilatation_row(p.b)) / 3.0;
      const Scalar shift = (mean_trace - (p.strain(0) + p.strain(1) + p.strain(2))) / 3.0;
      for (int r = 0; r < 3; ++r) {
        p.b.row(r) += correction.transpose();
        p.strain(r) += shift;
      }
    }
  }
  return out;
}

/// The temperature change the return sees at load factor lambda, and its
/// derivative with respect to lambda: lambda dT0 with small strain; with
/// finite kinematics the change whose linear thermal strain is the Green
/// strain of the free thermal stretch, dT (1 + alpha dT / 2).
void thermal_change(const IsotropicMaterial& mat, Kinematics kinematics, Scalar dt0,
                    Scalar lambda, Scalar& change, Scalar& rate) {
  const Scalar dt = lambda * dt0;
  if (kinematics == Kinematics::Finite) {
    change = thermal_green_lagrange_change(mat, dt);
    rate = dt0 * (1.0 + mat.thermal_expansion() * dt);
  } else {
    change = dt;
    rate = dt0;
  }
}

Matrix3 tensor3(const Vector6& v) {
  Matrix3 t;
  t << v(0), v(3), v(5), v(3), v(1), v(4), v(5), v(4), v(2);
  return t;
}

Vector6 voigt_of(const Matrix3& t) {
  Vector6 v;
  v << t(0, 0), t(1, 1), t(2, 2), t(0, 1), t(1, 2), t(2, 0);
  return v;
}

}  // namespace

std::string to_string(Kinematics kinematics) {
  return kinematics == Kinematics::SmallStrain ? "small_strain" : "finite";
}

Kinematics parse_kinematics(const std::string& text) {
  if (text == "finite") return Kinematics::Finite;
  if (text == "small_strain") return Kinematics::SmallStrain;
  throw ConfigError("unknown kinematics '" + text + "'; expected \"finite\" or \"small_strain\"");
}

int elastoplastic_points(const FemModel& model) {
  return static_cast<int>(model.element().integration_rule(model.integration()).size());
}

ElastoplasticElement elastoplastic_element(const FemModel& model, Index e, const Vector& ue,
                                           const std::vector<PlasticState>& committed,
                                           bool mean_dilatation, const Vector* temperature,
                                           Scalar temperature_scale, bool want_tangent,
                                           Kinematics kinematics) {
  if (model.dofs_per_node() != model.dim()) {
    throw ConfigError("the elastoplastic analysis is written for continuum elements");
  }
  const Kinematic kin = kinematics_of(model, e, ue, kinematics, mean_dilatation, temperature);
  if (committed.size() != kin.points.size()) {
    throw SolverError("elastoplastic element: the stored internal variables do not match "
                      "the integration rule");
  }
  const IsotropicMaterial& mat = model.material_of(e);
  const StressState state = model.stress_state();
  const int dim = model.dim();
  const int n = model.mesh().nodes_per_elem();
  const int nd = static_cast<int>(ue.size());
  const bool finite = kinematics == Kinematics::Finite;
  const bool thermal = std::any_of(kin.points.begin(), kin.points.end(),
                                   [](const Point& p) { return p.delta_t != 0.0; });

  ElastoplasticElement out;
  out.internal_force = Vector::Zero(nd);
  if (want_tangent) out.tangent = Matrix::Zero(nd, nd);
  if (thermal) out.thermal_force_rate = Vector::Zero(nd);
  out.states.reserve(kin.points.size());
  for (std::size_t q = 0; q < kin.points.size(); ++q) {
    const Point& p = kin.points[q];
    Scalar change = 0.0;
    Scalar rate = 0.0;
    thermal_change(mat, kinematics, p.delta_t, temperature_scale, change, rate);
    const PlasticResponse r = j2_return(mat, state, p.strain, committed[q], change);
    const Scalar w = p.weight;
    out.internal_force.noalias() += w * (p.b.transpose() * r.stress);
    out.energy += w * r.energy;
    if (want_tangent) {
      out.tangent.noalias() += w * (p.b.transpose() * r.tangent * p.b);
      if (finite) {
        // Geometric stiffness: S : Delta delta E, and with the averaged
        // dilatation (tr S / 3) (Delta delta tr-bar E - Delta delta tr E).
        Matrix s(dim, dim);
        if (dim == 3) {
          s = tensor3(r.stress);
        } else {
          s << r.stress(0), r.stress(3), r.stress(3), r.stress(1);
        }
        Matrix geo = p.g.transpose() * s * p.g;
        if (kin.averaged) {
          geo += (r.stress(0) + r.stress(1) + r.stress(2)) / 3.0 *
                 (kin.gbar - p.g.transpose() * p.g);
        }
        for (int a = 0; a < n; ++a) {
          for (int b = 0; b < n; ++b) {
            const Scalar v = w * geo(a, b);
            for (int k = 0; k < dim; ++k) out.tangent(dim * a + k, dim * b + k) += v;
          }
        }
      }
    }
    if (thermal && rate != 0.0) {
      // At fixed displacement the thermal strain alpha dT_e(lambda) m moves
      // the stress by -C alpha (d dT_e/d lambda) m, with the consistent (and,
      // in plane stress, condensed) tangent.
      Vector6 m = Vector6::Zero();
      m.head(3).setConstant(mat.thermal_expansion() * rate);
      out.thermal_force_rate.noalias() -= w * (p.b.transpose() * (r.tangent * m));
    }
    if (r.yielding) ++out.yielding_points;
    out.states.push_back(r.state);
  }
  if (want_tangent) out.tangent = 0.5 * (out.tangent + out.tangent.transpose());
  return out;
}

ElastoplasticStress elastoplastic_stress(const FemModel& model, Index e, const Vector& ue,
                                         const std::vector<PlasticState>& committed,
                                         bool mean_dilatation, const Vector* temperature,
                                         Scalar temperature_scale, Kinematics kinematics) {
  const Kinematic kin = kinematics_of(model, e, ue, kinematics, mean_dilatation, temperature);
  if (committed.size() != kin.points.size()) {
    throw SolverError("elastoplastic element: the stored internal variables do not match "
                      "the integration rule");
  }
  const IsotropicMaterial& mat = model.material_of(e);
  const StressState state = model.stress_state();
  const int dim = model.dim();
  const bool finite = kinematics == Kinematics::Finite;
  ElastoplasticStress out;
  out.min_jacobian = std::numeric_limits<Scalar>::infinity();
  for (std::size_t q = 0; q < kin.points.size(); ++q) {
    const Point& p = kin.points[q];
    Scalar change = 0.0;
    Scalar rate = 0.0;
    thermal_change(mat, kinematics, p.delta_t, temperature_scale, change, rate);
    // The committed state is the converged one at this displacement, so the
    // return reproduces it (an elastic check from it).
    const PlasticResponse r = j2_return(mat, state, p.strain, committed[q], change);
    Vector6 cauchy = r.stress;
    // The deformation measures use the element's own strain, not the
    // averaged one.
    const Matrix sym = 0.5 * (p.h + p.h.transpose());
    if (finite) {
      Matrix3 f = Matrix3::Identity();
      f.topLeftCorner(dim, dim) += p.h;
      // Plane stress: the thickness stretch of the Green strain the return
      // found, sqrt(1 + 2 E_33).
      if (state == StressState::PlaneStress) f(2, 2) = std::sqrt(1.0 + 2.0 * r.strain_33);
      const Scalar j = f.determinant();
      cauchy = voigt_of(f * tensor3(r.stress) * f.transpose() / j);
      const Matrix green = sym + 0.5 * p.h.transpose() * p.h;
      out.max_strain = std::max(out.max_strain, green.cwiseAbs().maxCoeff());
      if (state == StressState::PlaneStress) {
        out.max_strain = std::max(out.max_strain, std::abs(r.strain_33));
      }
      out.min_jacobian = std::min(out.min_jacobian, j);
    } else {
      out.max_strain = std::max(out.max_strain, sym.cwiseAbs().maxCoeff());
      if (state == StressState::PlaneStress) {
        out.max_strain = std::max(out.max_strain, std::abs(r.strain_33));
      }
      const Matrix skew = 0.5 * (p.h - p.h.transpose());
      out.max_rotation =
          std::max(out.max_rotation, std::sqrt(0.5 * skew.squaredNorm()));
      out.max_quadratic_strain = std::max(out.max_quadratic_strain,
                                          (0.5 * p.h.transpose() * p.h).cwiseAbs().maxCoeff());
    }
    out.cauchy += cauchy;
    out.piola_kirchhoff += r.stress;
    out.max_point_von_mises = std::max(out.max_point_von_mises, von_mises_stress(cauchy));
    out.max_equivalent_plastic_strain =
        std::max(out.max_equivalent_plastic_strain, r.state.equivalent_plastic_strain);
    if (r.state.equivalent_plastic_strain > 0.0) ++out.plastic_points;
  }
  const Scalar count = static_cast<Scalar>(kin.points.size());
  out.cauchy /= count;
  out.piola_kirchhoff /= count;
  out.von_mises = von_mises_stress(out.cauchy);
  if (!finite) out.min_jacobian = 1.0;
  return out;
}

}  // namespace sparlab
