#include "sparlab/fem/StaticAnalysis.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// Nodal vector (force, reaction, ...) at node `n` of a full-length vector
/// with `ndpn` DOFs per node: its translational part, padded with a zero z
/// component on a 2-D model.
Vector3 nodal_vector(const Vector& full, Index n, int dim, int ndpn) {
  Vector3 v = Vector3::Zero();
  for (int k = 0; k < dim; ++k) v(k) = full(n * ndpn + k);
  return v;
}

/// Rotational part (a nodal moment, a rotation) of node `n` of a full-length
/// vector; zero on a model whose nodes carry translations only.
Vector3 nodal_rotation_part(const Vector& full, Index n, int ndpn) {
  Vector3 v = Vector3::Zero();
  if (ndpn == kMaxDofsPerNode) {
    for (int k = 0; k < 3; ++k) v(k) = full(n * ndpn + 3 + k);
  }
  return v;
}

/// Moment of `f` acting at `x` about the origin. The 2-D branch keeps the
/// scalar expression of the plane formulation, which is the z component of
/// the cross product written out.
Vector3 moment_about_origin(const Vector3& x, const Vector3& f, int dim) {
  if (dim == 2) return Vector3(0.0, 0.0, x.x() * f.y() - x.y() * f.x());
  return x.cross(f);
}

Scalar magnitude(const Vector3& v, int dim) {
  return dim == 2 ? std::hypot(v.x(), v.y()) : std::hypot(v.x(), v.y(), v.z());
}

std::string vector_text(const Vector3& v, int dim) {
  std::ostringstream os;
  os << "(" << v.x() << ", " << v.y();
  if (dim == 3) os << ", " << v.z();
  os << ")";
  return os.str();
}

}  // namespace

StaticAnalysis::StaticAnalysis(const FemModel& model, const Assembler& assembler,
                               StaticAnalysisOptions options)
    : model_(model), assembler_(assembler), options_(options) {
  if (!model_.finalized()) {
    throw ModelError(
        "StaticAnalysis requires a finalised model; call FemModel::finalize() first");
  }
  if (options_.check_model) require_well_posed(model_);
  prescribed_ = model_.dofs().prescribed_vector();
  layout_.dim = model_.dim();
  layout_.dofs_per_node = model_.dofs_per_node();
  layout_.coordinates = &model_.mesh().coordinates();
  layout_.unknowns = &model_.dofs().free_dofs();
}

void StaticAnalysis::use_external_solver(LinearSolver* solver) {
  external_solver_ = solver;
  prepared_ = false;
}

LinearSolver& StaticAnalysis::active_solver() {
  if (external_solver_ != nullptr) return *external_solver_;
  if (!owned_solver_) owned_solver_ = make_linear_solver(options_.linear);
  return *owned_solver_;
}

const LinearSolver& StaticAnalysis::solver() const {
  if (external_solver_ != nullptr) return *external_solver_;
  if (!owned_solver_) throw ModelError("StaticAnalysis::solver() called before prepare()");
  return *owned_solver_;
}

int StaticAnalysis::last_iterations() const {
  return prepared_ ? solver().last_iterations() : 0;
}

Vector StaticAnalysis::free_guess(const Vector* initial_guess) const {
  if (initial_guess == nullptr || initial_guess->size() != model_.dofs().num_dofs()) {
    return Vector();
  }
  return model_.dofs().restrict_to_free(*initial_guess);
}

void StaticAnalysis::prepare(const Vector* stiffness_scale) {
  k_full_ = assembler_.assemble_stiffness(stiffness_scale);
  k_ff_ = assembler_.reduce_free_free(k_full_);
  if (model_.dofs().has_nonzero_prescribed()) {
    k_fp_ = assembler_.reduce_free_prescribed(k_full_);
  } else {
    k_fp_ = SparseMatrix(model_.dofs().num_free(), model_.dofs().num_constrained());
  }
  LinearSolver& solver = active_solver();
  solver.set_layout(layout_);
  solver.factorize(k_ff_);
  prepared_ = true;
}

Vector StaticAnalysis::solve_load_vector(const Vector& applied_force,
                                         const Vector* initial_guess) {
  if (!prepared_) prepare();
  if (applied_force.size() != model_.dofs().num_dofs()) {
    std::ostringstream os;
    os << "load vector has length " << applied_force.size() << " but the model has "
       << model_.dofs().num_dofs() << " DOFs";
    throw ModelError(os.str());
  }

  Vector rhs = model_.dofs().restrict_to_free(applied_force);
  if (model_.dofs().has_nonzero_prescribed()) {
    const auto& fixed = model_.dofs().constrained_dofs();
    Vector up(static_cast<Eigen::Index>(fixed.size()));
    for (std::size_t k = 0; k < fixed.size(); ++k) {
      up(static_cast<Eigen::Index>(k)) = prescribed_(fixed[k]);
    }
    rhs -= k_fp_ * up;
  }

  LinearSolver& solver = active_solver();
  const Vector guess = solver.iterative() ? free_guess(initial_guess) : Vector();
  const Vector uf = guess.size() > 0 ? solver.solve_from(rhs, guess) : solver.solve(rhs);
  if (!uf.allFinite()) {
    throw SolverError(
        "the linear solver returned a non-finite displacement field; the reduced "
        "stiffness matrix is singular or severely ill-conditioned");
  }
  return model_.dofs().expand(uf);
}

Vector StaticAnalysis::solve_homogeneous(const Vector& rhs, const Vector* initial_guess) {
  if (!prepared_) prepare();
  if (rhs.size() != model_.dofs().num_dofs()) {
    std::ostringstream os;
    os << "adjoint right-hand side has length " << rhs.size() << " but the model has "
       << model_.dofs().num_dofs() << " DOFs";
    throw ModelError(os.str());
  }
  LinearSolver& solver = active_solver();
  const Vector guess = solver.iterative() ? free_guess(initial_guess) : Vector();
  const Vector free_rhs = model_.dofs().restrict_to_free(rhs);
  const Vector reduced =
      guess.size() > 0 ? solver.solve_from(free_rhs, guess) : solver.solve(free_rhs);
  if (!reduced.allFinite()) {
    throw SolverError(
        "the adjoint solve returned a non-finite field; the reduced stiffness matrix is "
        "singular or severely ill-conditioned");
  }
  Vector full = Vector::Zero(model_.dofs().num_dofs());
  const auto& free = model_.dofs().free_dofs();
  for (std::size_t k = 0; k < free.size(); ++k) full(free[k]) = reduced(static_cast<Eigen::Index>(k));
  return full;
}

StaticSolution StaticAnalysis::build_solution(const std::string& name, Scalar weight,
                                              const Vector& applied_force,
                                              const LoadCaseData* data) {
  const int dim = model_.dim();
  const int ndpn = model_.dofs_per_node();
  StaticSolution sol;
  sol.load_case_name = name;
  sol.weight = weight;
  sol.displacement = solve_load_vector(applied_force);

  // Residual of the reduced system, computed from the assembled blocks.
  {
    Vector rhs = model_.dofs().restrict_to_free(applied_force);
    if (model_.dofs().has_nonzero_prescribed()) {
      const auto& fixed = model_.dofs().constrained_dofs();
      Vector up(static_cast<Eigen::Index>(fixed.size()));
      for (std::size_t k = 0; k < fixed.size(); ++k) {
        up(static_cast<Eigen::Index>(k)) = prescribed_(fixed[k]);
      }
      rhs -= k_fp_ * up;
    }
    const Vector uf = model_.dofs().restrict_to_free(sol.displacement);
    sol.scaled_residual = scaled_residual(k_ff_, uf, rhs);
    sol.backward_error = backward_error(k_ff_, uf, rhs);
    if (!residual_accepted(sol.scaled_residual, sol.backward_error,
                           options_.linear.residual_tolerance)) {
      std::ostringstream os;
      os << "load case '" << name << "': the linear solve left a scaled residual of "
         << sol.scaled_residual << ", above the recorded tolerance "
         << options_.linear.residual_tolerance << ", and a backward error of "
         << sol.backward_error << ", above round-off (" << kRoundoffBackwardError << ")";
      throw SolverError(os.str());
    }
  }
  sol.solver_iterations = active_solver().last_iterations();
  sol.solver_name = active_solver().name();

  // Reactions from the full residual r = K u - f.
  const Vector residual = k_full_ * sol.displacement - applied_force;
  sol.reactions = Vector::Zero(model_.dofs().num_dofs());
  for (Index d : model_.dofs().constrained_dofs()) sol.reactions(d) = residual(d);

  sol.compliance = applied_force.dot(sol.displacement);
  sol.strain_energy = 0.5 * sol.displacement.dot(k_full_ * sol.displacement);
  if (data != nullptr && data->thermal.size() > 0) {
    // Elastic energy of sigma = D (B u - eps0): the thermal part of the load
    // does work against the free expansion, not against the stiffness.
    sol.strain_energy += data->thermal_self_energy - data->thermal.dot(sol.displacement);
  }

  // Peak displacement magnitude and its node.
  const Index nn = model_.mesh().num_nodes();
  for (Index n = 0; n < nn; ++n) {
    const Scalar mag = magnitude(nodal_vector(sol.displacement, n, dim, ndpn), dim);
    if (mag > sol.max_displacement_magnitude) {
      sol.max_displacement_magnitude = mag;
      sol.max_displacement_node = n;
    }
  }

  // Global force and moment balance about the origin. On a shell or beam
  // model the nodal moments (applied and reacted) add to the moments of the
  // forces. Both balances are measured against the gross size of the applied
  // loads, sum |f_n| and sum |x_n| |f_n|, the scale of the round-off in their
  // sums: the resultant of a self-equilibrated load - a thermal strain, or a
  // self-weight carried by a traction - is itself round-off. A case driven by
  // prescribed displacements alone applies no load: its reactions balance
  // among themselves, measured against the gross force |K| |u| they are
  // formed from at their DOFs, which bounds their rounding - a prescribed
  // rigid motion (a shaken base) strains nothing, and its reactions are that
  // rounding alone.
  EquilibriumCheck& eq = sol.equilibrium;
  Vector gross = Vector::Zero(k_full_.rows());
  for (Eigen::Index col = 0; col < k_full_.outerSize(); ++col) {
    const Scalar uc = std::abs(sol.displacement(col));
    for (SparseMatrix::InnerIterator it(k_full_, col); it; ++it) {
      gross(it.row()) += std::abs(it.value()) * uc;
    }
  }
  Vector gross_reactions = Vector::Zero(k_full_.rows());
  for (Index d : model_.dofs().constrained_dofs()) gross_reactions(d) = gross(d);
  Scalar force_scale = 0.0;
  Scalar applied_moment_scale = 0.0;
  Scalar reaction_scale = 0.0;
  Scalar reaction_moment_scale = 0.0;
  for (Index n = 0; n < nn; ++n) {
    const Vector3 x = model_.mesh().node(n);
    const Vector3 fa = nodal_vector(applied_force, n, dim, ndpn);
    const Vector3 fr = nodal_vector(sol.reactions, n, dim, ndpn);
    const Scalar gr = std::max(fr.norm(), nodal_vector(gross_reactions, n, dim, ndpn).norm());
    eq.applied_force += fa;
    eq.reaction_force += fr;
    eq.applied_moment += moment_about_origin(x, fa, dim);
    eq.reaction_moment += moment_about_origin(x, fr, dim);
    force_scale += fa.norm();
    applied_moment_scale += x.norm() * fa.norm();
    reaction_scale += gr;
    reaction_moment_scale += x.norm() * gr;
    if (ndpn == kMaxDofsPerNode) {
      const Vector3 ma = nodal_rotation_part(applied_force, n, ndpn);
      const Vector3 mr = nodal_rotation_part(sol.reactions, n, ndpn);
      eq.applied_moment += ma;
      eq.reaction_moment += mr;
      applied_moment_scale += ma.norm();
      reaction_moment_scale += mr.norm();
    }
  }
  if (force_scale == 0.0) force_scale = reaction_scale;
  if (applied_moment_scale == 0.0) applied_moment_scale = reaction_moment_scale;
  eq.force_residual = eq.applied_force + eq.reaction_force;
  eq.moment_residual = eq.applied_moment + eq.reaction_moment;
  eq.relative_force_error = eq.force_residual.norm() / std::max(force_scale, 1.0e-30);
  eq.relative_moment_error =
      eq.moment_residual.norm() / std::max(applied_moment_scale, 1.0e-30);

  if (force_scale > 0.0 && eq.relative_force_error > options_.equilibrium_tolerance) {
    std::ostringstream os;
    os << "load case '" << name << "': global force balance is violated. Applied "
       << vector_text(eq.applied_force, dim) << " N, reactions "
       << vector_text(eq.reaction_force, dim) << " N, relative error "
       << eq.relative_force_error << " exceeds the tolerance "
       << options_.equilibrium_tolerance;
    throw SolverError(os.str());
  }
  if (applied_moment_scale > 0.0 &&
      eq.relative_moment_error > options_.equilibrium_tolerance) {
    log::warn("load case '", name, "': moment balance residual ",
              (dim == 2 ? eq.moment_residual.z() : eq.moment_residual.norm()),
              " N m (relative ", eq.relative_moment_error,
              ") exceeds the equilibrium tolerance ", options_.equilibrium_tolerance);
  }

  return sol;
}

std::vector<StaticSolution> StaticAnalysis::solve_all(const Vector* stiffness_scale) {
  prepare(stiffness_scale);
  const std::vector<Vector>& loads = model_.load_vectors();
  const std::vector<LoadCaseSpec>& specs = model_.load_case_specs();

  std::vector<StaticSolution> solutions;
  solutions.reserve(loads.size());
  for (std::size_t l = 0; l < loads.size(); ++l) {
    solutions.push_back(
        build_solution(specs[l].name, specs[l].weight, loads[l], &model_.load_case_data(l)));
    log::debug("load case '", specs[l].name, "': compliance ",
               solutions.back().compliance, " J, max |u| ",
               solutions.back().max_displacement_magnitude, " m");
  }
  return solutions;
}

Scalar StaticAnalysis::weighted_compliance(const std::vector<StaticSolution>& solutions,
                                           const std::vector<Scalar>& weights) {
  if (solutions.size() != weights.size()) {
    throw ModelError("weighted_compliance: solution and weight counts differ");
  }
  Scalar c = 0.0;
  for (std::size_t l = 0; l < solutions.size(); ++l) {
    c += weights[l] * solutions[l].compliance;
  }
  return c;
}

}  // namespace sparlab
