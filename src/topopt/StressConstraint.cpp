#include "sparlab/topopt/StressConstraint.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {

void StressConstraintOptions::validate(const SimpOptions& simp) const {
  if (!enabled) return;
  if (!(limit > 0.0)) {
    std::ostringstream os;
    os << "stress.limit must be a positive allowable von Mises stress in Pa (got " << limit
       << ")";
    throw ConfigError(os.str());
  }
  if (!(p_norm >= 2.0 && p_norm <= 64.0)) {
    throw ConfigError("stress.p_norm must lie in [2, 64]");
  }
  if (!(relaxation > 0.0 && relaxation < simp.penalty)) {
    std::ostringstream os;
    os << "stress.relaxation (q = " << relaxation
       << ") must be positive and below the SIMP penalty (p = " << simp.penalty
       << "), otherwise the relaxed stress does not vanish with the density";
    throw ConfigError(os.str());
  }
  if (!(scaling_blend > 0.0 && scaling_blend <= 1.0)) {
    throw ConfigError("stress.scaling_blend must lie in (0, 1]");
  }
  if (!(feasibility_tolerance > 0.0)) {
    throw ConfigError("stress.feasibility_tolerance must be positive");
  }
}

Matrix von_mises_matrix(StressState state, Scalar poisson) {
  // Full six-component von Mises matrix and the map from the model's Voigt
  // vector onto (xx, yy, zz, xy, yz, zx).
  Matrix6 v6 = Matrix6::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) v6(i, j) = i == j ? 1.0 : -0.5;
    v6(3 + i, 3 + i) = 3.0;
  }
  const int nv = voigt_components(stress_state_dimension(state));
  Matrix map = Matrix::Zero(6, nv);
  if (state == StressState::ThreeDimensional) {
    map.setIdentity();
  } else {
    map(0, 0) = 1.0;
    map(1, 1) = 1.0;
    map(3, 2) = 1.0;
    if (state == StressState::PlaneStrain) {
      map(2, 0) = poisson;
      map(2, 1) = poisson;
    }
  }
  return map.transpose() * v6 * map;
}

StressConstraint::StressConstraint(const FemModel& model, const Assembler& assembler,
                                   const DensityFilter& filter,
                                   StressConstraintOptions options)
    : model_(model), assembler_(assembler), filter_(filter), options_(options) {
  if (!options_.enabled) throw ConfigError("StressConstraint built with stress.enabled = false");
  if (filter_.type() == FilterType::Sensitivity) {
    throw ConfigError(
        "stress constraints need an exact design-to-density map; use filter.type = "
        "density (or none), not the heuristic sensitivity filter");
  }
  if (model_.dofs().has_nonzero_prescribed()) {
    throw ConfigError(
        "stress constraints assume design-independent loads; non-zero prescribed "
        "displacements are not supported in a stress-constrained run");
  }
  if (model_.dofs_per_node() != model_.dim()) {
    throw ConfigError(
        "the aggregated stress constraint is formulated on continuum elements; a shell "
        "or beam model's stresses vary through the thickness or over the section, and "
        "the constraint does not support them");
  }
  v_ = von_mises_matrix(model_.stress_state(), model_.material().poisson_ratio());
  uniform_ = assembler_.uses_element_cache();
  if (uniform_) db_uniform_ = db_at_centre(0);
}

Matrix StressConstraint::db_at_centre(Index e) const {
  const NaturalPoint centre = model_.element().reference_centroid();
  const StrainOperator op =
      model_.element().strain_operator(model_.mesh().element_coordinates(e), centre);
  return model_.constitutive() * op.b;  // nv x edofs
}

StressEvaluation StressConstraint::evaluate(ComplianceObjective& objective,
                                            const ObjectiveEvaluation& eval,
                                            std::size_t load_case, Scalar scale,
                                            bool need_gradients) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * model_.dofs_per_node();
  const DofManager& dofs = model_.dofs();
  if (load_case >= eval.displacements.size()) {
    throw ConfigError("stress constraint asked for a load case the evaluation lacks");
  }
  if (eval.physical_density.size() != ne) {
    throw ConfigError("stress constraint received a density of the wrong length");
  }
  const Vector& u = eval.displacements[load_case];
  const Vector& rho = eval.physical_density;
  const Scalar q = options_.relaxation;
  const Scalar big_p = options_.p_norm;
  const Scalar limit = options_.limit;

  StressEvaluation out;
  out.scale = scale;
  out.solid_von_mises.setZero(ne);
  out.relaxed_von_mises.setZero(ne);

  // Solid-material stress at each element centre and its relaxed ratio.
  Vector ue(edofs);
  std::vector<Vector> sigma(static_cast<std::size_t>(ne));
  Vector s(ne);
  for (Index e = 0; e < ne; ++e) {
    dofs.gather(mesh.element_nodes(e), npe, u, ue);
    const Matrix& db = uniform_ ? db_uniform_ : db_at_centre(e);
    sigma[static_cast<std::size_t>(e)] = db * ue;
    const Vector& sig = sigma[static_cast<std::size_t>(e)];
    const Scalar vm = std::sqrt(std::max(sig.dot(v_ * sig), 0.0));
    out.solid_von_mises(e) = vm;
    out.relaxed_von_mises(e) = std::pow(std::max(rho(e), 0.0), q) * vm;
    s(e) = out.relaxed_von_mises(e) / limit;
  }
  out.max_relaxed_ratio = s.maxCoeff(&out.max_element);

  // p-norm, scaled by the maximum so large P cannot overflow.
  const Scalar smax = std::max(out.max_relaxed_ratio, 1.0e-300);
  Scalar sum = 0.0;
  for (Index e = 0; e < ne; ++e) sum += std::pow(s(e) / smax, big_p);
  out.p_norm_ratio = smax * std::pow(sum, 1.0 / big_p);
  out.constraint = scale * out.p_norm_ratio - 1.0;
  if (!std::isfinite(out.constraint)) {
    throw SolverError("the aggregated stress constraint evaluated to a non-finite value");
  }
  if (!need_gradients) return out;

  // dg_PN/ds_e and the adjoint right-hand side.
  const Scalar gpn = std::max(out.p_norm_ratio, 1.0e-300);
  Vector dg_ds(ne);
  for (Index e = 0; e < ne; ++e) dg_ds(e) = std::pow(s(e) / gpn, big_p - 1.0);

  Vector psi = Vector::Zero(model_.dofs().num_dofs());
  Vector explicit_term = Vector::Zero(ne);
  for (Index e = 0; e < ne; ++e) {
    const Scalar vm = out.solid_von_mises(e);
    const Scalar r = std::max(rho(e), 0.0);
    const Scalar rq = std::pow(r, q);
    // Explicit density dependence through rho^q (zero where rho = 0, since
    // q < 1 would otherwise blow up; the term is multiplied by rho^q anyway).
    explicit_term(e) = r > 0.0 ? dg_ds(e) * q * std::pow(r, q - 1.0) * vm / limit : 0.0;
    if (vm <= 0.0 || dg_ds(e) == 0.0) continue;
    const Matrix& db = uniform_ ? db_uniform_ : db_at_centre(e);
    const Vector& sig = sigma[static_cast<std::size_t>(e)];
    // d sigma_vm / d u_e = (V sigma)^T D B / sigma_vm, weighted.
    const Vector row = (dg_ds(e) * rq / (limit * vm)) * (db.transpose() * (v_ * sig));
    dofs.scatter_add(mesh.element_nodes(e), npe, row, psi);
  }

  // Adjoint solve with the factorisation of the current SIMP stiffness.
  const Vector lambda = objective.solve_adjoint(psi, static_cast<int>(load_case));
  const Vector dfactor = simp_stiffness_derivatives(rho, objective.simp());

  Vector dgpn(ne);
  Vector le(edofs);
  for (Index e = 0; e < ne; ++e) {
    dofs.gather(mesh.element_nodes(e), npe, u, ue);
    dofs.gather(mesh.element_nodes(e), npe, lambda, le);
    const Scalar implicit = -dfactor(e) * le.dot(assembler_.element_stiffness(e) * ue);
    dgpn(e) = explicit_term(e) + implicit;
  }
  out.dg_dphysical = scale * dgpn;
  // Through the projection (when on) and the filter.
  out.dg_dx = objective.chain_to_design(eval, out.dg_dphysical);
  return out;
}

Scalar StressConstraint::next_scale(Scalar previous_scale,
                                    const StressEvaluation& current) const {
  if (!(current.p_norm_ratio > 0.0)) return previous_scale;
  const Scalar target = current.max_relaxed_ratio / current.p_norm_ratio;
  return options_.scaling_blend * target + (1.0 - options_.scaling_blend) * previous_scale;
}

}  // namespace sparlab
