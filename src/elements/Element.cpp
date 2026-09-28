#include "sparlab/elements/Element.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Hex8.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Tet10.hpp"
#include "sparlab/elements/Tet4.hpp"
#include "sparlab/elements/Tri3.hpp"

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

}  // namespace sparlab
