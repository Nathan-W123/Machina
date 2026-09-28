#include "sparlab/fem/HeatConduction.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/FaceGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// Shape-function gradients (dim x num_nodes) read off a strain operator: the
/// normal-strain row i of B holds dN_a/dx_i in column dim * a + i.
Matrix gradients_from_b(const Matrix& b, int dim, int nodes) {
  Matrix g(dim, nodes);
  for (int a = 0; a < nodes; ++a) {
    for (int i = 0; i < dim; ++i) g(i, a) = b(i, dim * a + i);
  }
  return g;
}

int face_points(FaceShape shape, const IntegrationOptions& integration) {
  return shape == FaceShape::Tri6 ? std::max(integration.edge_points, 3)
                                  : integration.edge_points;
}

}  // namespace

Matrix element_conductivity(const FemModel& model, Index e) {
  const Mesh& mesh = model.mesh();
  const Element& element = model.element();
  const int nn = mesh.nodes_per_elem();
  const int dim = mesh.dim();
  const Scalar t = dim == 2 ? model.thickness() : 1.0;
  const Scalar k = model.material_of(e).conductivity();
  const Matrix coords = mesh.element_coordinates(e);
  Matrix ke = Matrix::Zero(nn, nn);
  for (const IntegrationPoint& ip : element.integration_rule(model.integration())) {
    const StrainOperator op = element.strain_operator(coords, ip.point);
    const Matrix g = gradients_from_b(op.b, dim, nn);
    ke.noalias() += (k * t * ip.weight * op.detJ) * (g.transpose() * g);
  }
  return 0.5 * (ke + ke.transpose());
}

ConductionResult solve_conduction(const FemModel& model, const ConductionSpec& spec,
                                  const LinearSolverOptions& linear_in) {
  if (model.dofs_per_node() != model.dim()) {
    throw ConfigError("the conduction solve is formulated for continuum meshes; a shell "
                      "or beam model takes its temperature from regions or a uniform value");
  }
  const Mesh& mesh = model.mesh();
  const Index nn = mesh.num_nodes();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int dim = mesh.dim();
  const Scalar t = dim == 2 ? model.thickness() : 1.0;
  for (const IsotropicMaterial& m : model.materials()) {
    if (!(m.conductivity() > 0.0)) {
      throw ConfigError("material '" + m.name() +
                        "' has no thermal conductivity; set material.conductivity [W/(m "
                        "K)] to solve for the temperature field");
    }
  }
  if (spec.prescribed.empty() && spec.convection.empty()) {
    throw ConfigError(
        "the conduction problem has neither a prescribed temperature nor a convection "
        "boundary: with fluxes and sources alone the temperature is determined only up to "
        "a constant (and has no steady state unless they balance). Prescribe a "
        "temperature somewhere or add convection");
  }

  ConductionResult result;
  ConductionSummary& sum = result.summary;
  TripletList triplets;
  triplets.reserve(static_cast<std::size_t>(ne) * npe * npe);
  Vector rhs = Vector::Zero(nn);
  Scalar source_heat = 0.0;
  Scalar flux_heat = 0.0;

  // Conduction.
  for (Index e = 0; e < ne; ++e) {
    const Matrix ke = element_conductivity(model, e);
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < npe; ++a) {
      for (int b = 0; b < npe; ++b) triplets.emplace_back(nodes[a], nodes[b], ke(a, b));
    }
  }

  // Volumetric sources: F_a = int N_a Q t dV with the stiffness rule.
  const Element& element = model.element();
  for (const RegionValue& source : spec.sources) {
    std::vector<Index> elements;
    if (source.whole_model) {
      for (Index e = 0; e < ne; ++e) elements.push_back(e);
    } else {
      elements = source.region.select_elements(mesh);
      if (elements.empty()) {
        throw ConfigError("heat source region '" + source.region.name +
                          "' selected no element");
      }
    }
    sum.source_elements += static_cast<int>(elements.size());
    for (Index e : elements) {
      const Matrix coords = mesh.element_coordinates(e);
      const Index* nodes = mesh.element_nodes(e);
      for (const IntegrationPoint& ip : element.integration_rule(model.integration())) {
        const StrainOperator op = element.strain_operator(coords, ip.point);
        const Vector n = element.shape_functions(ip.point);
        const Scalar w = source.value * t * ip.weight * op.detJ;
        for (int a = 0; a < npe; ++a) rhs(nodes[a]) += w * n(a);
        source_heat += w * n.sum();
      }
    }
  }

  // Surface flux and convection on boundary faces.
  const bool need_faces = !spec.fluxes.empty() || !spec.convection.empty();
  const std::vector<Mesh::BoundaryFace> faces =
      need_faces ? mesh.boundary_faces() : std::vector<Mesh::BoundaryFace>();
  const FaceShape shape = face_shape_of(mesh.element_type());
  const int points = face_points(shape, model.integration());
  const std::vector<std::vector<int>>& table = element_local_faces(mesh.element_type());
  const auto face_nodes = [&](const Mesh::BoundaryFace& face) {
    const Index* enodes = mesh.element_nodes(face.element);
    const std::vector<int>& fn = table[static_cast<std::size_t>(face.local_face)];
    std::vector<Index> out(fn.size());
    for (std::size_t a = 0; a < fn.size(); ++a) out[a] = enodes[fn[a]];
    return out;
  };
  for (const RegionValue& flux : spec.fluxes) {
    const std::vector<Mesh::BoundaryFace> hit = faces_in_region(mesh, faces, flux.region);
    if (hit.empty()) {
      throw ConfigError("heat flux region '" + flux.region.name +
                        "' matched no boundary face; it must contain every node of at "
                        "least one face on the mesh boundary");
    }
    sum.flux_faces += static_cast<int>(hit.size());
    for (const Mesh::BoundaryFace& face : hit) {
      const Matrix xf = element_face_coordinates(mesh, face.element, face.local_face);
      const Vector w = face_shape_integrals(shape, xf, t, points);
      const std::vector<Index> fnodes = face_nodes(face);
      for (std::size_t a = 0; a < fnodes.size(); ++a) {
        rhs(fnodes[a]) += flux.value * w(static_cast<Eigen::Index>(a));
      }
      flux_heat += flux.value * w.sum();
    }
  }
  struct ConvectionFace {
    std::vector<Index> nodes;
    Vector integrals;
    Scalar h = 0.0;
    Scalar ambient = 0.0;
  };
  std::vector<ConvectionFace> convection_faces;
  for (const ConvectionSpec& conv : spec.convection) {
    if (!(conv.film_coefficient > 0.0)) {
      throw ConfigError("convection region '" + conv.region.name +
                        "' needs a positive film coefficient h [W/(m^2 K)]");
    }
    const std::vector<Mesh::BoundaryFace> hit = faces_in_region(mesh, faces, conv.region);
    if (hit.empty()) {
      throw ConfigError("convection region '" + conv.region.name +
                        "' matched no boundary face");
    }
    sum.convection_faces += static_cast<int>(hit.size());
    for (const Mesh::BoundaryFace& face : hit) {
      const Matrix xf = element_face_coordinates(mesh, face.element, face.local_face);
      const Matrix m = face_shape_products(shape, xf, t, points);
      const Vector w = face_shape_integrals(shape, xf, t, points);
      const std::vector<Index> fnodes = face_nodes(face);
      for (std::size_t a = 0; a < fnodes.size(); ++a) {
        rhs(fnodes[a]) += conv.film_coefficient * conv.ambient * w(static_cast<Eigen::Index>(a));
        for (std::size_t b = 0; b < fnodes.size(); ++b) {
          triplets.emplace_back(fnodes[a], fnodes[b],
                                conv.film_coefficient *
                                    m(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)));
        }
      }
      convection_faces.push_back({fnodes, w, conv.film_coefficient, conv.ambient});
    }
  }

  SparseMatrix k(nn, nn);
  k.setFromTriplets(triplets.begin(), triplets.end());
  k.makeCompressed();

  // Prescribed temperatures, later regions winning.
  std::vector<char> fixed(static_cast<std::size_t>(nn), 0);
  Vector temperature = Vector::Zero(nn);
  for (const RegionValue& p : spec.prescribed) {
    const std::vector<Index> nodes = p.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("prescribed-temperature region '" + p.region.name +
                        "' selected no nodes");
    }
    for (Index n : nodes) {
      fixed[static_cast<std::size_t>(n)] = 1;
      temperature(n) = p.value;
    }
  }
  std::vector<Index> free_nodes;
  std::vector<Index> reduced(static_cast<std::size_t>(nn), -1);
  for (Index n = 0; n < nn; ++n) {
    if (fixed[static_cast<std::size_t>(n)]) {
      ++sum.prescribed_nodes;
    } else {
      reduced[static_cast<std::size_t>(n)] = static_cast<Index>(free_nodes.size());
      free_nodes.push_back(n);
    }
  }
  const Index nf = static_cast<Index>(free_nodes.size());
  if (nf > 0) {
    // K_ff T_f = F_f - K_fp T_p.
    TripletList kff;
    Vector b(nf);
    for (Index i = 0; i < nf; ++i) b(i) = rhs(free_nodes[static_cast<std::size_t>(i)]);
    for (Eigen::Index col = 0; col < k.outerSize(); ++col) {
      const Index rc = reduced[static_cast<std::size_t>(col)];
      for (SparseMatrix::InnerIterator it(k, col); it; ++it) {
        const Index rr = reduced[static_cast<std::size_t>(it.row())];
        if (rr < 0) continue;
        if (rc >= 0) {
          kff.emplace_back(rr, rc, it.value());
        } else {
          b(rr) -= it.value() * temperature(col);
        }
      }
    }
    SparseMatrix a(nf, nf);
    a.setFromTriplets(kff.begin(), kff.end());
    a.makeCompressed();
    LinearSolverOptions linear = linear_in;
    // The scalar problem has no rigid-body near-null space for multigrid;
    // it is a fraction of the structural size and is factorised directly.
    linear.type = LinearSolverType::SimplicialLdlt;
    const std::unique_ptr<LinearSolver> solver = make_linear_solver(linear);
    solver->factorize(a);
    const Vector tf = solver->solve(b);
    if (!tf.allFinite()) {
      throw SolverError("the conduction solve returned a non-finite temperature field");
    }
    sum.scaled_residual = scaled_residual(a, tf, b);
    if (sum.scaled_residual > linear.residual_tolerance) {
      std::ostringstream os;
      os << "the conduction solve left a scaled residual of " << sum.scaled_residual
         << ", above the tolerance " << linear.residual_tolerance;
      throw SolverError(os.str());
    }
    for (Index i = 0; i < nf; ++i) temperature(free_nodes[static_cast<std::size_t>(i)]) = tf(i);
    sum.solver = solver->name();
  } else {
    sum.solver = "none (every node prescribed)";
  }

  // Heat balance: r = K T - F is the heat entering at the prescribed nodes.
  const Vector r = k * temperature - rhs;
  Scalar prescribed_in = 0.0;
  for (Index n = 0; n < nn; ++n) {
    if (!fixed[static_cast<std::size_t>(n)]) continue;
    prescribed_in += r(n);
    if (r(n) > 0.0) {
      sum.prescribed_inflow += r(n);
    } else {
      sum.prescribed_outflow -= r(n);
    }
  }
  Scalar convection_in = 0.0;
  for (const ConvectionFace& f : convection_faces) {
    for (std::size_t a = 0; a < f.nodes.size(); ++a) {
      convection_in += f.h * (f.ambient - temperature(f.nodes[a])) *
                       f.integrals(static_cast<Eigen::Index>(a));
    }
  }
  // Convection uses the consistent face matrix; its row sums equal the
  // integrals above, so this counts the same heat the solve balanced.
  sum.applied_heat = source_heat + flux_heat + convection_in;
  sum.prescribed_heat = -prescribed_in;
  const Scalar scale =
      std::max({std::abs(sum.applied_heat), std::abs(sum.prescribed_heat), std::abs(source_heat),
                std::abs(flux_heat), sum.prescribed_inflow, sum.prescribed_outflow, 1.0e-300});
  sum.relative_balance_error = std::abs(sum.applied_heat - sum.prescribed_heat) / scale;
  sum.min_temperature = temperature.minCoeff();
  sum.max_temperature = temperature.maxCoeff();
  if (sum.relative_balance_error > 1.0e-6) {
    std::ostringstream os;
    os << "the conduction solution violates the heat balance: " << sum.applied_heat
       << " W enter through sources, fluxes and convection, " << sum.prescribed_heat
       << " W leave through prescribed temperatures (relative error "
       << sum.relative_balance_error << ")";
    throw SolverError(os.str());
  }
  result.temperature = std::move(temperature);
  return result;
}

}  // namespace sparlab
