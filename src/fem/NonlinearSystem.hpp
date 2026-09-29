/// \file NonlinearSystem.hpp
/// \brief The non-linear system of a load case - internal forces, loads that
///        follow the deformation, tangent, the committed history of the
///        elastoplastic points - and its factorisation, shared by the static
///        path-following (NonlinearStatic.cpp) and the non-linear transient
///        (Dynamics.cpp). Internal to the library.
#pragma once

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/FaceGeometry.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Elastoplastic.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/Loads.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/TotalLagrangian.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace sparlab {
namespace detail {

/// A boundary face under a pressure that follows it.
struct PressureFace {
  Index element = 0;
  std::vector<int> local;  ///< the face's nodes as local indices of the element
  Scalar pressure = 0.0;   ///< [Pa] at lambda = 1
};

/// True when a sparse matrix equals its transpose to round-off: no entry of
/// A - A^T exceeds 1e-12 times the largest entry of A.
inline bool numerically_symmetric(const SparseMatrix& a) {
  const SparseMatrix difference = SparseMatrix(a.transpose()) - a;
  Scalar top = 0.0;
  for (Eigen::Index k = 0; k < a.nonZeros(); ++k) top = std::max(top, std::abs(a.valuePtr()[k]));
  Scalar asymmetry = 0.0;
  for (Eigen::Index k = 0; k < difference.nonZeros(); ++k) {
    asymmetry = std::max(asymmetry, std::abs(difference.valuePtr()[k]));
  }
  return asymmetry <= 1.0e-12 * top;
}

/// Factorisation of the free-free tangent: LDL^T for a symmetric tangent
/// (with the inertia it reveals), LU for a non-symmetric one or when LDL^T
/// meets a vanishing pivot. A follower pressure makes the element tangents
/// non-symmetric, but where the rim of the loaded surface is held - a closed
/// surface, or one whose edges lie on supports or symmetry planes - the
/// non-symmetric parts cancel in the assembled free-free block; a tangent
/// declared non-symmetric is therefore tested, and factorised as symmetric
/// when it is.
class TangentFactor {
 public:
  bool factorize(const SparseMatrix& a, bool symmetric) {
    negative_ = -1;
    if (!symmetric) symmetric = numerically_symmetric(a);
    use_lu_ = !symmetric;
    if (symmetric) {
      ldlt_.compute(a);
      if (ldlt_.info() == Eigen::Success) {
        const Vector d = ldlt_.vectorD();
        const Scalar top = d.cwiseAbs().maxCoeff();
        if (top > 0.0 && d.cwiseAbs().minCoeff() > 1.0e-13 * top && d.allFinite()) {
          negative_ = static_cast<int>((d.array() < 0.0).count());
          return true;
        }
      }
      use_lu_ = true;
    }
    lu_.analyzePattern(a);
    lu_.factorize(a);
    return lu_.info() == Eigen::Success;
  }
  Vector solve(const Vector& b) {
    Vector x = use_lu_ ? Vector(lu_.solve(b)) : Vector(ldlt_.solve(b));
    return x;
  }
  int negative_pivots() const { return negative_; }
  /// The last tangent factorised was symmetric (LDL^T).
  bool symmetric() const { return !use_lu_; }
  std::string name() const { return use_lu_ ? "SparseLU" : "SimplicialLDLT"; }

 private:
  Eigen::SimplicialLDLT<SparseMatrix> ldlt_;
  Eigen::SparseLU<SparseMatrix> lu_;
  bool use_lu_ = false;
  int negative_ = -1;
};

/// Everything the path-following needs from the model at a state (u, lambda).
struct Evaluation {
  Vector internal;   ///< f_int [N]
  Vector external;   ///< f_ext at (u, lambda) [N]
  Vector residual;   ///< f_int - f_ext
  Vector load_rate;  ///< q = -dR/dlambda
  SparseMatrix tangent;
  Scalar energy = 0.0;
  /// Norm of the gross assembly, sum_e |f_e| per DOF plus |f_ext| (and the
  /// thermal forces below): the scale of the round-off in the residual,
  /// below which no tolerance can be met.
  Scalar gross = 0.0;
  /// Norm of the gross thermal forces, sum_e |lambda df_e/dlambda| per DOF:
  /// the load a temperature change exerts, which a free expansion balances
  /// with no external force or reaction at all.
  Scalar thermal = 0.0;
  /// Elastoplastic elements: the internal variables of every point after the
  /// return from the committed ones (empty for an elastic element), and the
  /// number of points whose return was plastic.
  std::vector<std::vector<PlasticState>> states;
  int yielding_points = 0;
};

/// Norm over the free DOFs of |K| |u|, the tangent's absolute entries times
/// the absolute displacement: each stored displacement is exact only to its
/// relative rounding eps, which the stiffness turns into forces of
/// eps |K| |u|. In a slender structure under large rotation this dominates:
/// the Tet10 elastica's residual stagnates at 0.3 eps times it (1.3e-5 N),
/// 180 times above the gross element forces' share.
inline Scalar stiffness_gross(const SparseMatrix& k, const Vector& u, const std::vector<Index>& free) {
  Vector product = Vector::Zero(k.rows());
  for (Eigen::Index col = 0; col < k.outerSize(); ++col) {
    const Scalar uc = std::abs(u(col));
    if (uc == 0.0) continue;
    for (SparseMatrix::InnerIterator it(k, col); it; ++it) {
      product(it.row()) += std::abs(it.value()) * uc;
    }
  }
  Scalar sum = 0.0;
  for (Index d : free) sum += product(d) * product(d);
  return std::sqrt(sum);
}

/// Round-off floor of a residual: a multiple of the machine epsilon times the
/// gross size of the sums that form it, and of the forces the rounding of
/// the displacement itself causes. Newton's residual stagnates at about 150
/// epsilon times the gross element forces on the plane cantilever (the
/// element forces are themselves sums of stress times area), so 1024 epsilon
/// sits safely above that; and at 0.3 epsilon |K||u| on the elastica, so 64
/// epsilon sits above that - both far below any useful tolerance.
///
/// A floor belongs to a state that is resolved to working precision. One
/// above a millionth of the residual's scale (`scale`) means the state is
/// not: past a limit or a plastic collapse load the displacement runs away
/// and inflates |K||u| without bound, until the "floor" exceeds the load
/// itself and any residual would pass. Such a floor is not accepted in place
/// of the tolerance (0 is returned).
inline constexpr Scalar kFloorLimit = 1.0e-6;
inline Scalar residual_floor(const Evaluation& ev, Scalar k_gross, Scalar scale) {
  const Scalar eps = std::numeric_limits<Scalar>::epsilon();
  const Scalar floor = std::max(1024.0 * eps * ev.gross, 64.0 * eps * k_gross);
  return floor <= kFloorLimit * scale ? floor : 0.0;
}

class NonlinearSystem {
 public:
  NonlinearSystem(const FemModel& model, const Assembler& assembler, std::size_t lc,
                  const NonlinearOptions& options)
      : model_(model), assembler_(assembler), options_(options),
        small_(options.kinematics == Kinematics::SmallStrain) {
    const LoadCaseSpec& spec = model.load_case_specs()[lc];
    const LoadCaseData& data = model.load_case_data(lc);
    const Index n = model.dofs().num_dofs();
    dead_ = data.mechanical.size() > 0 ? data.mechanical : Vector::Zero(n);
    if (data.body.size() > 0) dead_ += data.body;
    if (data.temperature.size() > 0) temperature_ = data.temperature;

    // Elastoplastic elements keep their points' internal variables, starting
    // virgin; mean dilatation where it relaxes a constraint.
    const Index ne = model.mesh().num_elements();
    const int points = elastoplastic_points(model);
    committed_.resize(static_cast<std::size_t>(ne));
    averaged_.assign(static_cast<std::size_t>(ne), 0);
    virgin_.assign(static_cast<std::size_t>(points), PlasticState());
    for (Index e = 0; e < ne; ++e) {
      if (!model.material_of(e).plasticity().enabled()) continue;
      if (model.dim() == 2 && !model.material_of(e).plasticity().plane_compatible()) {
        throw ConfigError("material '" + model.material_of(e).name() +
                          "': in a plane model the out-of-plane axis z must be the rolling, "
                          "transverse or normal direction of its Hill48 frame");
      }
      plastic_ = true;
      // Armstrong-Frederick recovery makes the consistent tangent
      // non-symmetric (factorised by LU, without the inertia test).
      if (!model.material_of(e).plasticity().symmetric_tangent()) symmetric_ = false;
      committed_[static_cast<std::size_t>(e)].assign(static_cast<std::size_t>(points),
                                                     PlasticState());
      const ElementType type = model.mesh().element_type();
      const bool wanted =
          options.mean_dilatation == MeanDilatation::All ||
          (options.mean_dilatation == MeanDilatation::Auto &&
           (type == ElementType::Quad4 || type == ElementType::Hex8));
      averaged_[static_cast<std::size_t>(e)] =
          wanted && points > 1 && model.stress_state() != StressState::PlaneStress;
    }

    // Small strain: every load acts on the undeformed geometry, as in the
    // linear analysis. Finite kinematics: a follower pressure leaves the dead
    // load and is integrated over the deformed faces.
    if (!small_ && options.follower_pressure && !spec.pressures.empty()) {
      LoadCaseSpec only;
      only.name = spec.name;
      only.pressures = spec.pressures;
      dead_ -= assemble_load_vector(model.mesh(), model.element(), only, model.thickness(),
                                    model.integration());
      const std::vector<Mesh::BoundaryFace> boundary = model.mesh().boundary_faces();
      const std::vector<std::vector<int>>& table = element_local_faces(model.mesh().element_type());
      for (const PressureLoadSpec& p : spec.pressures) {
        for (const Mesh::BoundaryFace& face : faces_in_region(model.mesh(), boundary, p.region)) {
          PressureFace pf;
          pf.element = face.element;
          pf.local = table[static_cast<std::size_t>(face.local_face)];
          pf.pressure = p.pressure;
          faces_[face.element].push_back(pf);
        }
      }
      symmetric_ = false;
    }
    // With finite kinematics a rotation acts at the deformed position: the
    // centrifugal load of the reference geometry leaves the dead load and
    // returns as lambda (f_c0 + omega^2 M_perp u).
    if (!small_ && spec.centrifugal.enabled) {
      LoadCaseSpec only;
      only.name = spec.name;
      only.centrifugal = spec.centrifugal;
      centrifugal_ = assemble_body_load_vector(model, only);
      dead_ -= centrifugal_;
      omega2_ = spec.centrifugal.angular_velocity * spec.centrifugal.angular_velocity;
      const Vector3 e = spec.centrifugal.axis.normalized();
      const int dim = model.dim();
      Matrix perp = Matrix::Identity(dim, dim);
      if (dim == 3) perp -= e * e.transpose();
      const int npe = model.mesh().nodes_per_elem();
      // M_perp = int rho N^T N dV (I - e e^T): the element mass matrix
      // (which carries the element's density) with each node pair's
      // isotropic block turned into the projection normal to the axis.
      spin_ = assembler.assemble_elementwise([&](Index el) {
        const Matrix& m = assembler.element_mass(el);  // with the density, dim x dim blocks
        Matrix me = Matrix::Zero(m.rows(), m.cols());
        for (int a = 0; a < npe; ++a) {
          for (int b = 0; b < npe; ++b) {
            me.block(dim * a, dim * b, dim, dim) = m(dim * a, dim * b) * perp;
          }
        }
        return me;
      });
    }
  }

  bool symmetric() const { return symmetric_; }
  const Vector& dead() const { return dead_; }
  /// Some material is elastoplastic.
  bool plastic() const { return plastic_; }
  /// Element e averages its dilatation.
  bool averaged(Index e) const { return averaged_[static_cast<std::size_t>(e)] != 0; }
  /// Element e is elastoplastic.
  bool elastoplastic(Index e) const { return !committed_[static_cast<std::size_t>(e)].empty(); }
  /// The internal variables of element e's points at the last converged
  /// step (the virgin state for an elastic element).
  const std::vector<PlasticState>& committed(Index e) const {
    const std::vector<PlasticState>& c = committed_[static_cast<std::size_t>(e)];
    return c.empty() ? virgin_ : c;
  }
  /// Make the internal variables of a converged evaluation the committed
  /// ones.
  void commit(Evaluation& ev) {
    if (!plastic_) return;
    for (std::size_t e = 0; e < committed_.size(); ++e) {
      if (!committed_[e].empty()) committed_[e] = std::move(ev.states[e]);
    }
    ev.states.clear();
  }
  /// The largest accumulated plastic strain of the committed state.
  Scalar max_plastic_strain() const {
    Scalar top = 0.0;
    for (const std::vector<PlasticState>& element : committed_) {
      for (const PlasticState& p : element) top = std::max(top, p.equivalent_plastic_strain);
    }
    return top;
  }

  Evaluation evaluate(const Vector& u, Scalar lambda, bool want_tangent) const {
    const Mesh& mesh = model_.mesh();
    const int dim = mesh.dim();
    const int npe = mesh.nodes_per_elem();
    const int nd = dim * npe;
    const Index ne = mesh.num_elements();
    const Index n = model_.dofs().num_dofs();
    const Scalar thickness = dim == 2 ? model_.thickness() : 1.0;
    const FaceShape shape = face_shape_of(mesh.element_type());
    const int points = shape == FaceShape::Tri6
                           ? std::max(model_.integration().edge_points, 3)
                           : model_.integration().edge_points;
    const Vector* temperature = temperature_.size() > 0 ? &temperature_ : nullptr;

    std::vector<Vector> f_el(static_cast<std::size_t>(ne));
    std::vector<Vector> rate_el(static_cast<std::size_t>(ne));
    std::vector<Vector> pressure_el(static_cast<std::size_t>(ne));
    std::vector<Matrix> k_el(want_tangent ? static_cast<std::size_t>(ne) : 0);
    std::vector<Scalar> energy_el(static_cast<std::size_t>(ne), 0.0);
    std::vector<std::vector<PlasticState>> states(plastic_ ? static_cast<std::size_t>(ne) : 0);
    std::vector<int> yielding(plastic_ ? static_cast<std::size_t>(ne) : 0, 0);
    std::string failure;
#ifdef SPARLAB_HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (Index e = 0; e < ne; ++e) {
      try {
        const Index* nodes = mesh.element_nodes(e);
        Vector ue(nd);
        for (int a = 0; a < npe; ++a) {
          for (int k = 0; k < dim; ++k) ue(dim * a + k) = u(nodes[a] * dim + k);
        }
        TotalLagrangianElement tl;
        if (small_ || elastoplastic(e)) {
          ElastoplasticElement ep = elastoplastic_element(
              model_, e, ue, committed(e), averaged(e), temperature, lambda, want_tangent,
              options_.kinematics);
          tl.internal_force = std::move(ep.internal_force);
          tl.tangent = std::move(ep.tangent);
          tl.thermal_force_rate = std::move(ep.thermal_force_rate);
          tl.energy = ep.energy;
          if (elastoplastic(e)) {
            states[static_cast<std::size_t>(e)] = std::move(ep.states);
            yielding[static_cast<std::size_t>(e)] = ep.yielding_points;
          }
        } else {
          tl = total_lagrangian_element(model_, e, ue, options_.law, temperature, lambda,
                                        want_tangent);
        }
        const auto it = faces_.find(e);
        if (it != faces_.end()) {
          const Matrix x0 = mesh.element_coordinates(e);
          Vector fp = Vector::Zero(nd);
          for (const PressureFace& face : it->second) {
            const int nf = static_cast<int>(face.local.size());
            Matrix xf(dim, nf);
            for (int a = 0; a < nf; ++a) {
              const int l = face.local[static_cast<std::size_t>(a)];
              xf.col(a) = x0.col(l) + ue.segment(dim * l, dim);
            }
            // The face load at unit load factor, and its stiffness at lambda.
            const Matrix f1 = face_pressure_forces(shape, xf, face.pressure, thickness, points);
            for (int a = 0; a < nf; ++a) {
              const int l = face.local[static_cast<std::size_t>(a)];
              fp.segment(dim * l, dim) += f1.col(a).head(dim);
            }
            if (want_tangent && lambda != 0.0) {
              const Matrix kp = face_pressure_stiffness(shape, xf, lambda * face.pressure,
                                                        thickness, points);
              for (int a = 0; a < nf; ++a) {
                const int la = face.local[static_cast<std::size_t>(a)];
                for (int b = 0; b < nf; ++b) {
                  const int lb = face.local[static_cast<std::size_t>(b)];
                  tl.tangent.block(dim * la, dim * lb, dim, dim) -=
                      kp.block(3 * a, 3 * b, dim, dim);
                }
              }
            }
          }
          pressure_el[static_cast<std::size_t>(e)] = std::move(fp);
        }
        f_el[static_cast<std::size_t>(e)] = std::move(tl.internal_force);
        rate_el[static_cast<std::size_t>(e)] = std::move(tl.thermal_force_rate);
        energy_el[static_cast<std::size_t>(e)] = tl.energy;
        if (want_tangent) k_el[static_cast<std::size_t>(e)] = std::move(tl.tangent);
      } catch (const std::exception& ex) {
#ifdef SPARLAB_HAVE_OPENMP
#pragma omp critical(sparlab_nonlinear_failure)
#endif
        if (failure.empty()) failure = ex.what();
      }
    }
    if (!failure.empty()) throw SolverError(failure);

    Evaluation out;
    out.internal = Vector::Zero(n);
    out.load_rate = dead_;
    Vector pressure = Vector::Zero(n);
    Vector gross = Vector::Zero(n);
    Vector thermal = Vector::Zero(n);  // sum_e |lambda thermal force rate|
    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = mesh.element_nodes(e);
      const Vector& fe = f_el[static_cast<std::size_t>(e)];
      const Vector& re = rate_el[static_cast<std::size_t>(e)];
      const Vector& pe = pressure_el[static_cast<std::size_t>(e)];
      for (int a = 0; a < npe; ++a) {
        for (int k = 0; k < dim; ++k) {
          const Index dof = nodes[a] * dim + k;
          out.internal(dof) += fe(dim * a + k);
          gross(dof) += std::abs(fe(dim * a + k));
          if (re.size() > 0) {
            out.load_rate(dof) -= re(dim * a + k);
            thermal(dof) += std::abs(lambda * re(dim * a + k));
          }
          if (pe.size() > 0) pressure(dof) += pe(dim * a + k);
        }
      }
      out.energy += energy_el[static_cast<std::size_t>(e)];
      if (plastic_) out.yielding_points += yielding[static_cast<std::size_t>(e)];
    }
    out.states = std::move(states);
    out.external = lambda * (dead_ + pressure);
    out.load_rate += pressure;
    if (centrifugal_.size() > 0) {
      const Vector rotation = centrifugal_ + omega2_ * (spin_ * u);
      out.external += lambda * rotation;
      out.load_rate += rotation;
    }
    out.residual = out.internal - out.external;
    out.gross = (gross + thermal + out.external.cwiseAbs()).norm();
    out.thermal = thermal.norm();
    if (want_tangent) {
      out.tangent = assembler_.assemble_elementwise(
          [&](Index e) -> Matrix { return k_el[static_cast<std::size_t>(e)]; });
      if (centrifugal_.size() > 0 && lambda != 0.0) out.tangent -= (lambda * omega2_) * spin_;
    }
    return out;
  }

 private:
  const FemModel& model_;
  const Assembler& assembler_;
  const NonlinearOptions& options_;
  bool small_ = false;
  bool plastic_ = false;
  std::vector<std::vector<PlasticState>> committed_;  ///< per element; empty if elastic
  std::vector<char> averaged_;                        ///< per element: mean dilatation
  std::vector<PlasticState> virgin_;                  ///< one per point, for elastic elements
  Vector dead_;
  Vector temperature_;
  std::map<Index, std::vector<PressureFace>> faces_;
  Vector centrifugal_;
  SparseMatrix spin_;
  Scalar omega2_ = 0.0;
  bool symmetric_ = true;
};

inline Vector restrict(const Vector& full, const std::vector<Index>& dofs) {
  Vector out(static_cast<Eigen::Index>(dofs.size()));
  for (std::size_t i = 0; i < dofs.size(); ++i) out(static_cast<Eigen::Index>(i)) = full(dofs[i]);
  return out;
}

inline void add_to(Vector& full, const std::vector<Index>& dofs, const Vector& part, Scalar factor) {
  for (std::size_t i = 0; i < dofs.size(); ++i) {
    full(dofs[i]) += factor * part(static_cast<Eigen::Index>(i));
  }
}

/// The free-free block of a full tangent (entries in the free DOFs' order).
inline SparseMatrix free_block(const SparseMatrix& full, const std::vector<Index>& free_dofs) {
  std::vector<Index> map(static_cast<std::size_t>(full.rows()), -1);
  for (std::size_t i = 0; i < free_dofs.size(); ++i) {
    map[static_cast<std::size_t>(free_dofs[i])] = static_cast<Index>(i);
  }
  TripletList triplets;
  triplets.reserve(static_cast<std::size_t>(full.nonZeros()));
  for (Eigen::Index col = 0; col < full.outerSize(); ++col) {
    const Index c = map[static_cast<std::size_t>(col)];
    if (c < 0) continue;
    for (SparseMatrix::InnerIterator it(full, col); it; ++it) {
      const Index r = map[static_cast<std::size_t>(it.row())];
      if (r >= 0) triplets.emplace_back(r, c, it.value());
    }
  }
  const auto size = static_cast<Eigen::Index>(free_dofs.size());
  SparseMatrix out(size, size);
  out.setFromTriplets(triplets.begin(), triplets.end());
  out.makeCompressed();
  return out;
}

/// The free-prescribed block applied to a prescribed-DOF vector.
inline Vector free_prescribed_times(const SparseMatrix& full, const std::vector<Index>& free_dofs,
                             const std::vector<Index>& fixed, const Vector& values) {
  Vector full_values = Vector::Zero(full.cols());
  for (std::size_t i = 0; i < fixed.size(); ++i) {
    full_values(fixed[i]) = values(static_cast<Eigen::Index>(i));
  }
  return restrict(full * full_values, free_dofs);
}


}  // namespace detail
}  // namespace sparlab
