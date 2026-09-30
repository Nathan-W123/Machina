#include "sparlab/elements/Element.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Hex8.hpp"
#include "sparlab/elements/Hex8Incompatible.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Tet10.hpp"
#include "sparlab/elements/Tet4.hpp"
#include "sparlab/elements/Tri3.hpp"

#include <Eigen/Cholesky>

#include <sstream>

namespace sparlab {
namespace {

/// Shape-function gradients (dim x num_nodes) read off a strain operator: the
/// normal-strain row i of B holds dN_a/dx_i in column dim * a + i for every
/// element in the library.
Matrix gradients_from_b(const Matrix& b, int dim, int nodes) {
  Matrix g(dim, nodes);
  for (int a = 0; a < nodes; ++a) {
    for (int i = 0; i < dim; ++i) g(i, a) = b(i, dim * a + i);
  }
  return g;
}

/// Stress tensor (dim x dim) from a Voigt vector (xx, yy, xy) or
/// (xx, yy, zz, xy, yz, zx).
Matrix stress_tensor(const Vector& s, int dim) {
  Matrix t(dim, dim);
  if (dim == 2) {
    t << s(0), s(2),
         s(2), s(1);
  } else {
    t << s(0), s(3), s(5),
         s(3), s(1), s(4),
         s(5), s(4), s(2);
  }
  return t;
}

void check_geometric_inputs(const Element& element, const Matrix& d, const Vector& v,
                            const char* what) {
  if (d.rows() != element.num_voigt() || d.cols() != element.num_voigt()) {
    std::ostringstream os;
    os << to_string(element.type()) << " " << what << " expects a " << element.num_voigt()
       << " x " << element.num_voigt() << " constitutive matrix, received " << d.rows()
       << " x " << d.cols();
    throw ModelError(os.str());
  }
  if (v.size() != element.num_dofs()) {
    std::ostringstream os;
    os << to_string(element.type()) << " " << what << " expects an element vector of "
       << element.num_dofs() << " entries, received " << v.size();
    throw ModelError(os.str());
  }
}

}  // namespace

std::string to_string(ElementFormulation formulation) {
  switch (formulation) {
    case ElementFormulation::Standard: return "standard";
    case ElementFormulation::IncompatibleModes: return "incompatible_modes";
  }
  return "standard";
}

ElementFormulation parse_element_formulation(const std::string& text) {
  if (text == "standard") return ElementFormulation::Standard;
  if (text == "incompatible_modes") return ElementFormulation::IncompatibleModes;
  throw ConfigError("unknown element formulation '" + text +
                    "'; expected \"standard\" or \"incompatible_modes\"");
}

Matrix Element::internal_mode_gradients(const Matrix& /*coords*/,
                                        const NaturalPoint& /*point*/) const {
  return Matrix(dim(), 0);
}

Matrix Element::internal_strain_operator(const Matrix& coords, const NaturalPoint& point) const {
  const Matrix g = internal_mode_gradients(coords, point);
  const int nd = dim();
  const auto m = static_cast<int>(g.cols());
  Matrix b = Matrix::Zero(num_voigt(), nd * m);
  for (int a = 0; a < m; ++a) {
    const int c = nd * a;
    if (nd == 2) {
      b(0, c + 0) = g(0, a);
      b(1, c + 1) = g(1, a);
      b(2, c + 0) = g(1, a);
      b(2, c + 1) = g(0, a);
    } else {
      b(0, c + 0) = g(0, a);
      b(1, c + 1) = g(1, a);
      b(2, c + 2) = g(2, a);
      b(3, c + 0) = g(1, a);
      b(3, c + 1) = g(0, a);
      b(4, c + 1) = g(2, a);
      b(4, c + 2) = g(1, a);
      b(5, c + 0) = g(2, a);
      b(5, c + 2) = g(0, a);
    }
  }
  return b;
}

InternalCondensation Element::condense_internal(const Matrix& coords, const Matrix& d,
                                                Scalar thickness,
                                                const IntegrationOptions& opts) const {
  InternalCondensation out;
  const int ni = num_internal_dofs();
  if (ni == 0) {
    out.coupling = Matrix(0, num_dofs());
    out.inverse = Matrix(0, 0);
    return out;
  }
  const Scalar t = dim() == 2 ? thickness : 1.0;
  Matrix kaa = Matrix::Zero(ni, ni);
  Matrix kau = Matrix::Zero(ni, num_dofs());
  for (const IntegrationPoint& ip : integration_rule(opts)) {
    const StrainOperator op = strain_operator(coords, ip.point);
    const Matrix bt = internal_strain_operator(coords, ip.point);
    const Scalar w = t * ip.weight * op.detJ;
    const Matrix dbt = d * bt;
    kaa.noalias() += w * (bt.transpose() * dbt);
    kau.noalias() += w * (dbt.transpose() * op.b);
  }
  kaa = 0.5 * (kaa + kaa.transpose());
  const Eigen::LLT<Matrix> llt(kaa);
  if (llt.info() != Eigen::Success) {
    throw SolverError(to_string(type()) +
                      ": the stiffness of the internal modes is not positive definite; the "
                      "constitutive matrix is not positive definite or the element is "
                      "degenerate");
  }
  out.coupling = llt.solve(kau);
  out.inverse = llt.solve(Matrix::Identity(ni, ni));
  return out;
}

Matrix Element::condensed_strain_operator(const Matrix& coords, const NaturalPoint& point,
                                          const InternalCondensation& condensation) const {
  const StrainOperator op = strain_operator(coords, point);
  if (condensation.coupling.rows() == 0) return op.b;
  return op.b - internal_strain_operator(coords, point) * condensation.coupling;
}

const std::vector<int>& Element::face_nodes(int local_face) const {
  const std::vector<std::vector<int>>& faces = element_local_faces(type());
  if (local_face < 0 || local_face >= static_cast<int>(faces.size())) {
    std::ostringstream os;
    os << to_string(type()) << " local face index " << local_face << " is outside [0, "
       << faces.size() - 1 << "]";
    throw MeshError(os.str());
  }
  return faces[static_cast<std::size_t>(local_face)];
}

Matrix Element::geometric_stiffness(const Matrix& coords, const Matrix& d, const Vector& ue,
                                    Scalar stress_scale, Scalar thickness,
                                    const IntegrationOptions& opts) const {
  check_geometric_inputs(*this, d, ue, "geometric stiffness");
  const int nd = dim();
  const int nn = num_nodes();
  const Scalar t = nd == 2 ? thickness : 1.0;
  // The scalar block G^T sigma G is shared by the dim displacement components.
  Matrix block = Matrix::Zero(nn, nn);
  for (const IntegrationPoint& ip : integration_rule(opts)) {
    const StrainOperator op = strain_operator(coords, ip.point);
    const Vector sigma = stress_scale * (d * (op.b * ue));
    const Matrix g = gradients_from_b(op.b, nd, nn);
    block.noalias() += (t * ip.weight * op.detJ) * (g.transpose() * stress_tensor(sigma, nd) * g);
  }
  Matrix kg = Matrix::Zero(num_dofs(), num_dofs());
  for (int a = 0; a < nn; ++a) {
    for (int b = 0; b < nn; ++b) {
      const Scalar value = 0.5 * (block(a, b) + block(b, a));
      for (int k = 0; k < nd; ++k) kg(nd * a + k, nd * b + k) = value;
    }
  }
  return kg;
}

Matrix Element::geometric_stiffness_of_stress(const Matrix& coords,
                                              const std::vector<Vector>& stresses,
                                              Scalar thickness,
                                              const IntegrationOptions& opts) const {
  if (dofs_per_node() != dim()) {
    throw ModelError(to_string(type()) +
                     " supplies its own geometric stiffness; the stress-field form is "
                     "for continuum elements");
  }
  const std::vector<IntegrationPoint> rule = integration_rule(opts);
  if (stresses.size() != rule.size()) {
    std::ostringstream os;
    os << to_string(type()) << " geometric stiffness received " << stresses.size()
       << " point stresses for a rule of " << rule.size() << " points";
    throw ModelError(os.str());
  }
  const int nd = dim();
  const int nn = num_nodes();
  const Scalar t = nd == 2 ? thickness : 1.0;
  Matrix block = Matrix::Zero(nn, nn);
  for (std::size_t q = 0; q < rule.size(); ++q) {
    const IntegrationPoint& ip = rule[q];
    const StrainOperator op = strain_operator(coords, ip.point);
    const Matrix g = gradients_from_b(op.b, nd, nn);
    block.noalias() +=
        (t * ip.weight * op.detJ) * (g.transpose() * stress_tensor(stresses[q], nd) * g);
  }
  Matrix kg = Matrix::Zero(num_dofs(), num_dofs());
  for (int a = 0; a < nn; ++a) {
    for (int b = 0; b < nn; ++b) {
      const Scalar value = 0.5 * (block(a, b) + block(b, a));
      for (int k = 0; k < nd; ++k) kg(nd * a + k, nd * b + k) = value;
    }
  }
  return kg;
}

Vector Element::geometric_stiffness_derivative(const Matrix& coords, const Matrix& d,
                                               const Vector& phi, Scalar stress_scale,
                                               Scalar thickness,
                                               const IntegrationOptions& opts) const {
  check_geometric_inputs(*this, d, phi, "geometric-stiffness derivative");
  const int nd = dim();
  const int nn = num_nodes();
  const Scalar t = nd == 2 ? thickness : 1.0;
  // Displacement-gradient matrix H(i, k) = d phi_k / d x_i of the mode.
  Matrix modes(nn, nd);
  for (int a = 0; a < nn; ++a) {
    for (int k = 0; k < nd; ++k) modes(a, k) = phi(nd * a + k);
  }
  Vector out = Vector::Zero(num_dofs());
  for (const IntegrationPoint& ip : integration_rule(opts)) {
    const StrainOperator op = strain_operator(coords, ip.point);
    const Matrix h = gradients_from_b(op.b, nd, nn) * modes;
    const Matrix phi_tensor = h * h.transpose();  // Phi_ij = sum_k H_ik H_jk
    Vector hat(num_voigt());
    if (nd == 2) {
      hat << phi_tensor(0, 0), phi_tensor(1, 1), 2.0 * phi_tensor(0, 1);
    } else {
      hat << phi_tensor(0, 0), phi_tensor(1, 1), phi_tensor(2, 2), 2.0 * phi_tensor(0, 1),
          2.0 * phi_tensor(1, 2), 2.0 * phi_tensor(2, 0);
    }
    out.noalias() += (stress_scale * t * ip.weight * op.detJ) * (op.b.transpose() * (d * hat));
  }
  return out;
}

std::unique_ptr<Element> make_element(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return std::make_unique<Quad4Element>();
    case ElementType::Hex8: return std::make_unique<Hex8Element>();
    case ElementType::Tri3: return std::make_unique<Tri3Element>();
    case ElementType::Tet4: return std::make_unique<Tet4Element>();
    case ElementType::Tet10: return std::make_unique<Tet10Element>();
  }
  throw ConfigError("no element implementation registered for the requested type");
}

std::unique_ptr<Element> make_element(ElementType type, const IntegrationOptions& opts) {
  if (opts.formulation == ElementFormulation::Standard) return make_element(type);
  if (type != ElementType::Hex8) {
    throw ConfigError("the \"" + to_string(opts.formulation) +
                      "\" element formulation is available for Hex8 meshes only; the mesh "
                      "is " + to_string(type));
  }
  return std::make_unique<Hex8IncompatibleElement>();
}

}  // namespace sparlab
