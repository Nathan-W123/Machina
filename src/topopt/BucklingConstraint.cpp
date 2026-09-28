#include "sparlab/topopt/BucklingConstraint.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// Adjoint slots of the buckling modes, clear of the per-load-case slots of
/// the stress constraint.
constexpr int kAdjointSlotBase = 1024;

}  // namespace

void BucklingConstraintOptions::validate() const {
  if (!enabled) return;
  if (!(min_load_factor > 0.0)) {
    std::ostringstream os;
    os << "buckling_constraint.min_load_factor must be positive (got " << min_load_factor
       << "); it is the required multiple of the design load";
    throw ConfigError(os.str());
  }
  if (num_modes < 1 || num_modes > 50) {
    throw ConfigError("buckling_constraint.num_modes must lie in [1, 50]");
  }
  if (!(ks_parameter >= 1.0 && ks_parameter <= 500.0)) {
    throw ConfigError("buckling_constraint.ks_parameter must lie in [1, 500]");
  }
  if (!(solid_threshold > 0.0 && solid_threshold < 1.0)) {
    throw ConfigError("buckling_constraint.solid_threshold must lie in (0, 1)");
  }
}

Vector buckling_stress_factors(const Vector& rho, Scalar penalty) {
  Vector out(rho.size());
  for (Eigen::Index e = 0; e < rho.size(); ++e) out(e) = std::pow(std::max(rho(e), 0.0), penalty);
  return out;
}

Vector buckling_stress_derivatives(const Vector& rho, Scalar penalty) {
  Vector out(rho.size());
  for (Eigen::Index e = 0; e < rho.size(); ++e) {
    const Scalar r = std::max(rho(e), 0.0);
    out(e) = r > 0.0 ? penalty * std::pow(r, penalty - 1.0) : 0.0;
  }
  return out;
}

BucklingConstraint::BucklingConstraint(const FemModel& model, const Assembler& assembler,
                                       BucklingConstraintOptions options)
    : model_(model), assembler_(assembler), options_(std::move(options)) {
  if (!options_.enabled) {
    throw ConfigError("BucklingConstraint built with buckling_constraint.enabled = false");
  }
  options_.validate();
  options_.eigen.num_modes = options_.num_modes;
  const std::vector<LoadCaseSpec>& specs = model_.load_case_specs();
  if (options_.load_cases.empty()) {
    for (std::size_t l = 0; l < specs.size(); ++l) cases_.push_back(l);
  } else {
    for (const std::string& name : options_.load_cases) {
      const auto it = std::find_if(specs.begin(), specs.end(),
                                   [&](const LoadCaseSpec& s) { return s.name == name; });
      if (it == specs.end()) {
        std::ostringstream os;
        os << "buckling_constraint.load_cases names '" << name
           << "', which is not a load case of the deck (";
        for (std::size_t l = 0; l < specs.size(); ++l) os << (l ? ", " : "") << specs[l].name;
        os << ")";
        throw ConfigError(os.str());
      }
      cases_.push_back(static_cast<std::size_t>(it - specs.begin()));
    }
  }
  subspace_.resize(specs.size());
}

BucklingEvaluation BucklingConstraint::evaluate(ComplianceObjective& objective,
                                                const ObjectiveEvaluation& eval,
                                                std::size_t load_case, bool need_gradients) {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * model_.dofs_per_node();
  if (load_case >= eval.displacements.size()) {
    throw ConfigError("buckling constraint asked for a load case the evaluation lacks");
  }
  const Vector& u = eval.displacements[load_case];
  const Vector& rho = eval.physical_density;
  const Scalar penalty = objective.simp().penalty;
  const Vector stress_scale = buckling_stress_factors(rho, penalty);
  const DofManager& dofs = model_.dofs();

  // Buckling eigenproblem of the current design, warm started.
  const SparseMatrix k_g = assemble_geometric_stiffness(model_, assembler_, u, &stress_scale);
  const FreeSolve solve = [&](const Vector& b) {
    return dofs.restrict_to_free(objective.solve_adjoint(dofs.expand(b)));
  };
  std::vector<char> solid(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) {
    solid[static_cast<std::size_t>(e)] = rho(e) >= options_.solid_threshold ? 1 : 0;
  }
  Matrix& warm = subspace_[load_case];
  const BucklingResult buckling =
      solve_buckling(model_, assembler_, objective.stiffness(), k_g, options_.eigen, solve,
                     warm.size() > 0 ? &warm : nullptr, &eval.stiffness_factors, &solid);
  warm = buckling.subspace;

  BucklingEvaluation out;
  out.load_factors = buckling.load_factors;
  out.solid_energy_fraction = buckling.solid_energy_fraction;
  out.iterations = buckling.iterations;
  out.transformed = buckling.transformed;
  out.sigma = buckling.sigma;
  out.no_positive_load_factor = buckling.no_positive_load_factor;
  const Eigen::Index m = out.load_factors.size();
  if (m == 0) {
    // Nothing buckles under this load: the constraint is inactive.
    out.ks = 0.0;
    out.constraint = -1.0;
    if (need_gradients) {
      out.dg_dphysical.setZero(ne);
      out.dg_dx = objective.chain_to_design(eval, out.dg_dphysical);
    }
    return out;
  }

  // KS aggregate of r_i = lambda_req / lambda_i.
  const Scalar p = options_.ks_parameter;
  const Vector r = options_.min_load_factor * out.load_factors.cwiseInverse();
  const Scalar r_max = r.maxCoeff();
  Vector w(m);
  for (Eigen::Index i = 0; i < m; ++i) w(i) = std::exp(p * (r(i) - r_max));
  const Scalar sum = w.sum();
  w /= sum;
  out.ks = r_max + std::log(sum) / p;
  out.constraint = out.ks - 1.0;
  if (!std::isfinite(out.constraint)) {
    throw SolverError("the aggregated buckling constraint evaluated to a non-finite value");
  }
  if (!need_gradients) return out;

  // Sensitivities of every aggregated load factor.
  const Vector dk = simp_stiffness_derivatives(rho, objective.simp());
  const Vector dg = buckling_stress_derivatives(rho, penalty);
  const Element& element = model_.element();
  const Matrix& d = model_.constitutive();
  out.dg_dphysical.setZero(ne);
  Vector ue(edofs);
  Vector pe(edofs);
  Vector ae(edofs);
  std::vector<Vector> ge(static_cast<std::size_t>(ne));
  for (Eigen::Index i = 0; i < m; ++i) {
    const Scalar lambda = out.load_factors(i);
    const Vector phi = buckling.mode_shapes.col(i);
    // Adjoint load g = sum_e (E_G / E_0) A_e^T g_e(phi_e).
    Vector rhs = Vector::Zero(dofs.num_dofs());
    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = mesh.element_nodes(e);
      dofs.gather(nodes, npe, phi, pe);
      Vector& g = ge[static_cast<std::size_t>(e)];
      g = element.geometric_stiffness_derivative(mesh.element_coordinates(e), d, pe, 1.0,
                                                 model_.thickness(), model_.integration());
      dofs.scatter_add(nodes, npe, g, rhs, stress_scale(e));
    }
    const Vector adjoint = objective.solve_adjoint(
        rhs, kAdjointSlotBase + static_cast<int>(64 * load_case) + static_cast<int>(i));
    Vector dlambda(ne);
    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = mesh.element_nodes(e);
      dofs.gather(nodes, npe, u, ue);
      dofs.gather(nodes, npe, phi, pe);
      dofs.gather(nodes, npe, adjoint, ae);
      const Matrix& ke = assembler_.element_stiffness(e);
      const Scalar elastic = pe.dot(ke * pe);
      const Scalar geometric = ge[static_cast<std::size_t>(e)].dot(ue);
      const Scalar implicit = ae.dot(ke * ue);
      dlambda(e) = lambda * dk(e) * elastic + lambda * lambda * dg(e) * geometric -
                   lambda * lambda * dk(e) * implicit;
    }
    // d r_i / d rho = -lambda_req / lambda_i^2 d lambda_i / d rho.
    out.dg_dphysical -= (w(i) * options_.min_load_factor / (lambda * lambda)) * dlambda;
    out.dlambda_dphysical.push_back(std::move(dlambda));
  }
  out.dg_dx = objective.chain_to_design(eval, out.dg_dphysical);
  return out;
}

}  // namespace sparlab
