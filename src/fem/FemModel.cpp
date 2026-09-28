#include "sparlab/fem/FemModel.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/HeatConduction.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/Loads.hpp"

#include <numeric>
#include <sstream>
#include <utility>

namespace sparlab {

std::string to_string(MassType type) {
  switch (type) {
    case MassType::Consistent: return "consistent";
    case MassType::Lumped: return "lumped";
  }
  return "unknown";
}

FemModel::FemModel(Mesh mesh, IsotropicMaterial material, Scalar thickness,
                   StressState stress_state, IntegrationOptions integration)
    : mesh_(std::move(mesh)),
      materials_{std::move(material)},
      thickness_(thickness),
      stress_state_(stress_state),
      integration_(integration),
      element_(make_element(mesh_.element_type())),
      d_{materials_.front().constitutive(stress_state)},
      dofs_(mesh_.num_nodes(), element_->dofs_per_node()) {
  if (stress_state_dimension(stress_state_) != mesh_.dim()) {
    std::ostringstream os;
    os << "stress state '" << to_string(stress_state_) << "' belongs to a "
       << stress_state_dimension(stress_state_) << "-D model but the mesh is "
       << mesh_.dim() << "-D (" << to_string(mesh_.element_type()) << " elements); use "
       << (mesh_.dim() == 3 ? "'three_dimensional'" : "'plane_stress' or 'plane_strain'");
    throw ConfigError(os.str());
  }
  if (!(thickness_ > 0.0)) {
    std::ostringstream os;
    os << "model thickness must be positive (got " << thickness_ << " m)";
    throw ConfigError(os.str());
  }
  if (mesh_.dim() == 3 && thickness_ != 1.0) {
    std::ostringstream os;
    os << "a 3-D model has no thickness (got " << thickness_
       << " m); leave model.thickness at its default of 1 for a solid mesh";
    throw ConfigError(os.str());
  }
  mesh_.validate();
}

void FemModel::set_material(const IsotropicMaterial& material) {
  materials_.front() = material;
  d_.front() = material.constitutive(stress_state_);
}

void FemModel::assign_material(const IsotropicMaterial& material,
                               const std::vector<Index>& elements) {
  if (finalized_) {
    throw ModelError("materials must be assigned before FemModel::finalize()");
  }
  const Index ne = mesh_.num_elements();
  for (Index e : elements) {
    if (e < 0 || e >= ne) {
      std::ostringstream os;
      os << "material '" << material.name() << "' assigned to element " << e
         << ", outside [0, " << ne - 1 << "]";
      throw ModelError(os.str());
    }
  }
  materials_.push_back(material);
  d_.push_back(material.constitutive(stress_state_));
  if (element_material_.empty()) element_material_.assign(static_cast<std::size_t>(ne), 0);
  const int id = static_cast<int>(materials_.size()) - 1;
  for (Index e : elements) element_material_[static_cast<std::size_t>(e)] = id;
}

void FemModel::set_conduction_solver(const LinearSolverOptions& options) {
  conduction_solver_ = std::make_shared<LinearSolverOptions>(options);
}

void FemModel::finalize(bool require_load_cases) {
  if (finalized_) return;
  if (load_case_specs_.empty() && require_load_cases) {
    throw ConfigError("model has no load cases; define at least one");
  }
  const Index constrained = apply_constraints(mesh_, constraints_, dofs_);
  log::info("applied ", constraints_.size(), " boundary condition group(s): ",
            constrained, " of ", dofs_.num_dofs(), " DOFs prescribed, ",
            dofs_.num_free(), " free");

  load_vectors_.clear();
  load_data_.clear();
  load_vectors_.reserve(load_case_specs_.size());
  load_data_.reserve(load_case_specs_.size());
  for (const LoadCaseSpec& spec : load_case_specs_) {
    if (!(spec.weight >= 0.0)) {
      std::ostringstream os;
      os << "load case '" << spec.name << "' has a negative weight (" << spec.weight
         << "); weights must be non-negative";
      throw ConfigError(os.str());
    }
    LoadCaseData data;
    data.mechanical = assemble_load_vector(mesh_, *element_, spec, thickness_, integration_);
    Vector total = data.mechanical;
    if (spec.has_body_loads()) {
      data.body = assemble_body_load_vector(*this, spec);
      for (Index n = 0; n < mesh_.num_nodes(); ++n) {
        for (int k = 0; k < mesh_.dim(); ++k) {
          data.body_resultant(k) += data.body(n * dofs_.dofs_per_node() + k);
        }
      }
      total += data.body;
    }
    if (spec.has_temperature()) {
      if (spec.temperature.source == TemperatureSpec::Source::Conduction) {
        LinearSolverOptions linear;
        if (conduction_solver_) linear = *conduction_solver_;
        const ConductionResult solved = solve_conduction(*this, spec.temperature.conduction, linear);
        data.temperature = solved.temperature;
        data.conduction = solved.summary;
        data.conduction_solved = true;
        log::info("load case '", spec.name, "': steady conduction, temperature ",
                  solved.summary.min_temperature, " to ", solved.summary.max_temperature,
                  " K, heat in ", solved.summary.applied_heat, " W, out through prescribed "
                  "temperatures ", solved.summary.prescribed_heat, " W (relative balance "
                  "error ", solved.summary.relative_balance_error, ")");
      } else {
        data.temperature = resolve_region_temperatures(mesh_, spec.temperature);
      }
      const ThermalLoad thermal = assemble_thermal_load(*this, data.temperature);
      data.thermal = thermal.force;
      data.thermal_self_energy = thermal.self_energy;
      total += data.thermal;
    }
    // A zero load vector is a mistake *unless* the case is driven by prescribed
    // displacements, which is how the patch test and any enforced-deflection
    // study work.
    if (total.norm() == 0.0 && !dofs_.has_nonzero_prescribed() &&
        !spec.prescribed_displacement_only) {
      log::warn("load case '", spec.name,
                "' has a zero resultant force vector and no non-zero prescribed "
                "displacement, so its solution is identically zero");
    }
    load_vectors_.push_back(std::move(total));
    load_data_.push_back(std::move(data));
  }

  if (!load_case_specs_.empty()) {
    Scalar weight_sum = 0.0;
    for (const LoadCaseSpec& spec : load_case_specs_) weight_sum += spec.weight;
    if (!(weight_sum > 0.0)) {
      throw ConfigError(
          "the sum of load-case weights is zero; the compliance objective would be "
          "identically zero");
    }
  }
  finalized_ = true;
}

const std::vector<Vector>& FemModel::load_vectors() const {
  if (!finalized_) {
    throw ModelError("load vectors requested before FemModel::finalize() was called");
  }
  return load_vectors_;
}

const LoadCaseData& FemModel::load_case_data(std::size_t l) const {
  if (!finalized_) {
    throw ModelError("load-case data requested before FemModel::finalize() was called");
  }
  if (l >= load_data_.size()) {
    std::ostringstream os;
    os << "load case " << l << " requested, the model has " << load_data_.size();
    throw ModelError(os.str());
  }
  return load_data_[l];
}

std::vector<Scalar> FemModel::normalised_weights() const {
  Scalar sum = 0.0;
  for (const LoadCaseSpec& spec : load_case_specs_) sum += spec.weight;
  if (!(sum > 0.0)) return {};
  std::vector<Scalar> w;
  w.reserve(load_case_specs_.size());
  for (const LoadCaseSpec& spec : load_case_specs_) w.push_back(spec.weight / sum);
  return w;
}

Vector FemModel::element_volumes() const {
  const Index ne = mesh_.num_elements();
  Vector v(ne);
  if (mesh_.dim() == 2) {
    for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_measure(e) * thickness_;
  } else {
    for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_measure(e);
  }
  return v;
}

Scalar FemModel::domain_volume() const { return element_volumes().sum(); }

}  // namespace sparlab
