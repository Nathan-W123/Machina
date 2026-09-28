#include "sparlab/topopt/Sensitivity.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Timer.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {

ComplianceObjective::ComplianceObjective(const FemModel& model,
                                         const Assembler& assembler,
                                         const DensityFilter& filter,
                                         const DesignDomain& domain, SimpOptions simp,
                                         StaticAnalysisOptions analysis_options)
    : model_(model),
      assembler_(assembler),
      filter_(filter),
      domain_(domain),
      simp_(simp),
      analysis_options_(analysis_options),
      weights_(model.normalised_weights()) {
  simp_.validate();
  if (filter_.num_elements() != model_.mesh().num_elements()) {
    throw ConfigError("the filter was built for a different mesh than the model");
  }
  // The SIMP stiffness floor keeps K_ff positive definite everywhere, so the
  // rigid-body diagnostics are run once by the caller on the solid model
  // instead of on every SIMP iterate.
  analysis_options_.check_model = false;
}

void ComplianceObjective::set_projection(Scalar beta, Scalar eta) {
  if (filter_.type() == FilterType::Sensitivity) {
    throw ConfigError(
        "the Heaviside projection needs the density filter (or none): the sensitivity "
        "filter modifies the gradient heuristically, so there is no chain rule to "
        "extend through the projection");
  }
  if (!(beta > 0.0) || !(eta > 0.0 && eta < 1.0)) {
    throw ConfigError("projection needs beta > 0 and eta in (0, 1)");
  }
  projection_ = true;
  beta_ = beta;
  eta_ = eta;
}

void ComplianceObjective::set_overhang(const OverhangFilter* overhang) {
  if (overhang != nullptr && filter_.type() == FilterType::Sensitivity) {
    throw ConfigError(
        "the overhang filter needs the density filter (or none): the sensitivity filter "
        "modifies the gradient heuristically, so there is no chain rule through it");
  }
  if (overhang != nullptr && overhang->num_elements() != model_.mesh().num_elements()) {
    throw ConfigError("the overhang filter was built for a different mesh than the model");
  }
  overhang_ = overhang;
}

Vector ComplianceObjective::physical_density(const Vector& x) const {
  Vector rho = filter_.to_physical(x).cwiseMax(0.0).cwiseMin(1.0);
  if (overhang_ != nullptr) rho = overhang_->apply(rho);
  return projection_ ? heaviside_project(rho, beta_, eta_) : rho;
}

Vector ComplianceObjective::physical_density_at(const Vector& x, Scalar eta) const {
  if (!projection_) {
    throw ConfigError("a density at another projection threshold needs the projection on");
  }
  Vector rho = filter_.to_physical(x).cwiseMax(0.0).cwiseMin(1.0);
  if (overhang_ != nullptr) rho = overhang_->apply(rho);
  return heaviside_project(rho, beta_, eta);
}

Vector ComplianceObjective::chain_at(const ObjectiveEvaluation& eval, Scalar eta,
                                     const Vector& d_dphysical) const {
  if (!projection_) {
    throw ConfigError("a gradient at another projection threshold needs the projection on");
  }
  Vector d = d_dphysical.cwiseProduct(heaviside_derivative(eval.printable_density, beta_, eta));
  if (overhang_ != nullptr) d = overhang_->pull_back(eval.filtered_density, d);
  return filter_.pull_back(d);
}

const SparseMatrix& ComplianceObjective::stiffness() const {
  if (!analysis_) {
    throw ModelError("the stiffness matrix was requested before the first evaluation");
  }
  return analysis_->stiffness();
}

Vector ComplianceObjective::chain_to_design(const ObjectiveEvaluation& eval,
                                            const Vector& d_dphysical) const {
  Vector d = eval.projection_derivative.size() == d_dphysical.size()
                 ? Vector(d_dphysical.cwiseProduct(eval.projection_derivative))
                 : d_dphysical;
  if (overhang_ != nullptr) d = overhang_->pull_back(eval.filtered_density, d);
  return filter_.pull_back(d);
}

ObjectiveEvaluation ComplianceObjective::evaluate(const Vector& x, bool need_gradients) {
  const Index ne = model_.mesh().num_elements();
  if (x.size() != ne) {
    std::ostringstream os;
    os << "design vector has length " << x.size() << " but the mesh has " << ne
       << " elements";
    throw ConfigError(os.str());
  }
  if (!x.allFinite()) throw ConfigError("design vector contains NaN or Inf");

  ObjectiveEvaluation out;
  // Filtering is a convex combination, so the result stays in [0,1] up to
  // round-off; clamp to keep pow() well defined.
  out.filtered_density = filter_.to_physical(x).cwiseMax(0.0).cwiseMin(1.0);
  out.printable_density =
      overhang_ != nullptr ? overhang_->apply(out.filtered_density) : out.filtered_density;
  if (projection_) {
    out.physical_density = heaviside_project(out.printable_density, beta_, eta_);
    out.projection_derivative = heaviside_derivative(out.printable_density, beta_, eta_);
  } else {
    out.physical_density = out.printable_density;
  }
  out.stiffness_factors = simp_stiffness_factors(out.physical_density, simp_);

  Timer timer;
  if (!solver_) solver_ = make_linear_solver(analysis_options_.linear);
  analysis_ = std::make_unique<StaticAnalysis>(model_, assembler_, analysis_options_);
  StaticAnalysis& analysis = *analysis_;
  analysis.use_external_solver(solver_.get());
  analysis.prepare(&out.stiffness_factors);
  out.assemble_seconds = timer.elapsed_seconds();
  adjoint_iterations_ = 0;
  const bool warm = analysis_options_.linear.warm_start && solver_->iterative();
  previous_displacements_.resize(model_.load_vectors().size());

  timer.reset();
  const std::vector<Vector>& loads = model_.load_vectors();
  out.displacements.reserve(loads.size());
  out.load_case_compliance.reserve(loads.size());

  const int npe = model_.mesh().nodes_per_elem();
  const int edofs = npe * model_.dofs_per_node();
  out.element_strain_energy.setZero(ne);
  Vector dc_dphys = Vector::Zero(ne);
  const Vector dfactor = need_gradients
                             ? simp_stiffness_derivatives(out.physical_density, simp_)
                             : Vector();

  Vector ue(edofs);
  for (std::size_t l = 0; l < loads.size(); ++l) {
    Vector& previous = previous_displacements_[l];
    const Vector u = analysis.solve_load_vector(
        loads[l], warm && previous.size() == loads[l].size() ? &previous : nullptr);
    out.linear_iterations += analysis.last_iterations();
    total_linear_iterations_ += analysis.last_iterations();
    if (warm) previous = u;
    ++num_solves_;
    const Scalar c = loads[l].dot(u);
    out.load_case_compliance.push_back(c);
    out.compliance += weights_[l] * c;
    out.displacements.push_back(u);

    for (Index e = 0; e < ne; ++e) {
      model_.dofs().gather(model_.mesh().element_nodes(e), npe, u, ue);
      // u_e^T K_e^0 u_e: the unit-density element strain energy times two.
      const Scalar quad = ue.dot(assembler_.element_stiffness(e) * ue);
      out.element_strain_energy(e) +=
          weights_[l] * 0.5 * out.stiffness_factors(e) * quad;
      if (need_gradients) dc_dphys(e) -= weights_[l] * dfactor(e) * quad;
    }
  }
  out.solve_seconds = timer.elapsed_seconds();

  out.volume = domain_.volume_of(out.physical_density);
  out.volume_fraction = out.volume / domain_.domain_volume();

  if (need_gradients) {
    out.dc_dphysical = dc_dphys;
    if (projection_ || overhang_ != nullptr) {
      // Through the projection, the overhang filter and the (linear) density
      // filter, in that order.
      out.dc_dx = chain_to_design(out, dc_dphys);
      out.dv_dx = chain_to_design(out, domain_.element_volumes());
    } else {
      out.dc_dx = filter_.transform_gradient(x, dc_dphys);
      out.dv_dx = filter_.pull_back(domain_.element_volumes());
    }

    if (out.dc_dx.maxCoeff() > 0.0) {
      // For minimum compliance the gradient must be non-positive everywhere.
      // A positive entry indicates a bug or a pathological filter; report it
      // rather than clipping silently.
      log::warn("compliance gradient has a positive entry (max = ",
                out.dc_dx.maxCoeff(),
                "), which is not expected for minimum-compliance SIMP; the optimality "
                "criteria update will clamp it");
    }
  }

  if (!std::isfinite(out.compliance)) {
    throw SolverError("the compliance objective evaluated to a non-finite value");
  }
  return out;
}

Vector ComplianceObjective::solve_adjoint(const Vector& rhs, int slot) {
  if (!analysis_) {
    throw ModelError("solve_adjoint called before the objective was evaluated");
  }
  ++num_solves_;
  const bool warm = slot >= 0 && analysis_options_.linear.warm_start && solver_->iterative();
  if (warm && static_cast<std::size_t>(slot) >= previous_adjoints_.size()) {
    previous_adjoints_.resize(static_cast<std::size_t>(slot) + 1);
  }
  const Vector* guess = nullptr;
  if (warm && previous_adjoints_[static_cast<std::size_t>(slot)].size() == rhs.size()) {
    guess = &previous_adjoints_[static_cast<std::size_t>(slot)];
  }
  Vector lambda = analysis_->solve_homogeneous(rhs, guess);
  adjoint_iterations_ += analysis_->last_iterations();
  total_linear_iterations_ += analysis_->last_iterations();
  if (warm) previous_adjoints_[static_cast<std::size_t>(slot)] = lambda;
  return lambda;
}

Scalar ComplianceObjective::compliance_at(const Vector& x) {
  return evaluate(x, /*need_gradients=*/false).compliance;
}

SensitivityCheckResult verify_sensitivities(ComplianceObjective& objective,
                                            const DesignDomain& domain, const Vector& x,
                                            const std::vector<Index>& elements,
                                            Scalar step, Scalar tolerance) {
  if (!(step > 0.0)) {
    throw ConfigError("finite-difference step must be positive");
  }
  SensitivityCheckResult result;
  result.step = step;
  result.tolerance = tolerance;

  const ObjectiveEvaluation base = objective.evaluate(x, /*need_gradients=*/true);
  result.gradient_infinity_norm = base.dc_dx.cwiseAbs().maxCoeff();

  std::vector<Index> targets = elements;
  if (targets.empty()) {
    targets.resize(static_cast<std::size_t>(domain.num_elements()));
    for (Index e = 0; e < domain.num_elements(); ++e) {
      targets[static_cast<std::size_t>(e)] = e;
    }
  }

  Scalar sum_sq_rel = 0.0;
  Vector xp = x;
  Vector xm = x;
  for (Index e : targets) {
    SensitivityCheckEntry entry;
    entry.element = e;
    entry.analytical = base.dc_dx(e);

    const Scalar lo = domain.lower_bounds()(e);
    const Scalar hi = domain.upper_bounds()(e);
    if (!domain.is_free(e)) {
      entry.excluded = true;
      entry.exclusion_reason =
          "passive variable (lower bound equals upper bound), so no admissible "
          "perturbation exists";
    } else if (x(e) - step < lo || x(e) + step > hi) {
      std::ostringstream os;
      os << "design variable " << x(e) << " is within the step " << step
         << " of its bounds [" << lo << ", " << hi
         << "], so a central difference would leave the feasible box";
      entry.excluded = true;
      entry.exclusion_reason = os.str();
    }

    if (entry.excluded) {
      ++result.num_excluded;
      result.entries.push_back(std::move(entry));
      continue;
    }

    xp(e) = x(e) + step;
    xm(e) = x(e) - step;
    const Scalar cp = objective.compliance_at(xp);
    const Scalar cm = objective.compliance_at(xm);
    xp(e) = x(e);
    xm(e) = x(e);

    entry.finite_difference = (cp - cm) / (2.0 * step);
    entry.absolute_error = std::abs(entry.analytical - entry.finite_difference);
    const Scalar scale =
        std::max({std::abs(entry.analytical), std::abs(entry.finite_difference),
                  1.0e-30});
    entry.relative_error = entry.absolute_error / scale;
    result.max_scaled_error = std::max(
        result.max_scaled_error,
        entry.absolute_error / std::max(scale, 1.0e-3 * result.gradient_infinity_norm));

    result.max_absolute_error = std::max(result.max_absolute_error, entry.absolute_error);
    result.max_relative_error = std::max(result.max_relative_error, entry.relative_error);
    sum_sq_rel += entry.relative_error * entry.relative_error;
    ++result.num_tested;
    result.entries.push_back(std::move(entry));
  }

  if (result.num_tested == 0) {
    throw ConfigError(
        "every candidate element was excluded from the finite-difference check; move "
        "the design point away from its bounds or reduce the step size");
  }
  result.rms_relative_error =
      std::sqrt(sum_sq_rel / static_cast<Scalar>(result.num_tested));

  // Directional derivative along the gradient restricted to tested elements.
  {
    Vector direction = Vector::Zero(x.size());
    for (const SensitivityCheckEntry& entry : result.entries) {
      if (!entry.excluded) direction(entry.element) = entry.analytical;
    }
    const Scalar dnorm = direction.norm();
    if (dnorm > 0.0) {
      direction /= dnorm;
      // Shrink the step so no component leaves the box.
      Scalar h = step;
      for (Index e = 0; e < domain.num_elements(); ++e) {
        if (direction(e) == 0.0) continue;
        const Scalar room = direction(e) > 0.0
                                ? (domain.upper_bounds()(e) - x(e)) / direction(e)
                                : (domain.lower_bounds()(e) - x(e)) / direction(e);
        h = std::min(h, std::abs(room));
      }
      if (h > 0.0) {
        const Scalar cp = objective.compliance_at((x + h * direction).eval());
        const Scalar cm = objective.compliance_at((x - h * direction).eval());
        const Scalar fd = (cp - cm) / (2.0 * h);
        const Scalar an = base.dc_dx.dot(direction);
        result.directional_relative_error =
            std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), 1.0e-30});
      }
    }
  }

  result.passed = result.max_relative_error <= tolerance;
  log::info("sensitivity check: step ", step, ", ", result.num_tested, " elements (",
            result.num_excluded, " excluded), max relative error ",
            result.max_relative_error, ", RMS ", result.rms_relative_error,
            ", directional ", result.directional_relative_error, " -> ",
            (result.passed ? "PASS" : "FAIL"), " at tolerance ", tolerance);
  return result;
}

}  // namespace sparlab
