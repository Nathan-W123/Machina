/// \file Contact.hpp
/// \brief Unilateral contact in the non-linear static analysis: a surface of
///        the model against a rigid obstacle given analytically, or against
///        another surface of the model; frictionless or with Coulomb friction.
///
/// **Assumptions.** Small displacements and small sliding: the contact
/// geometry - the normals, the pairing of the two surfaces and the weights
/// below - is that of the reference configuration, and the gap is linear in
/// the displacement, which is the classical Signorini problem. It holds with
/// the `small_strain` kinematics of the non-linear analysis (linear
/// elasticity or J2 plasticity), and the analysis refuses `finite`
/// kinematics with contact. Linear elements only (Q4, Tri3, Hex8, Tet4): the
/// face of a Tet10 has zero corner weights `int N dA`, for which the dual
/// basis below does not exist.
///
/// **Discretisation (dual mortar).** On the slave surface the contact
/// pressure is interpolated with the dual basis \f$\psi_j\f$ of each face
/// (Wohlmuth), biorthogonal to its shape functions:
/// \f$\int \psi_j N_k\,dA = \delta_{jk} D_j\f$, \f$D_j = \int N_j\,dA\f$.
/// The contact force on slave node j is then \f$D_j\,p_j\,\nu_j\f$ alone
/// (\f$\nu_j\f$ the direction of the pressure on the slave body, into it),
/// and on master node l \f$-\sum_j M_{jl}\,p_j\,\nu_j\f$ with
/// \f$M_{jl} = \int \psi_j N_l^{m}\,dA\f$ over the slave surface, the master
/// shape function taken at the point the slave normal field projects to
/// (segment-by-segment in 2-D; on an auxiliary plane per slave face with
/// polygon clipping in 3-D, Puso and Laursen). The weighted gap
/// \f[
///   \tilde g_j = \tilde g_{0j} + D_j\,\nu_j\cdot u_j - \sum_l M_{jl}\,\nu_j\cdot u_l
/// \f]
/// must stay non-negative, \f$\tilde g_{0j} = \int\psi_j g_0\,dA\f$ the
/// weighted initial gap. Against a rigid obstacle \f$\nu_j\f$ is the
/// obstacle's normal at the node, \f$\tilde g_j = D_j\,(g(X_j) + \nu_j \cdot
/// (u_j - \lambda d))\f$ with the obstacle moved by \f$\lambda d\f$ - exact on
/// a flat obstacle, and a consistent nodal rule on a curved one.
///
/// **Solution.** The pressure is condensed: slave node j's equilibrium gives
/// \f$p_j = \nu_{F}\cdot R_{F}/(D_j |\nu_F|^2)\f$ from the residual
/// \f$R = f_{int} - f_{ext}\f$ of its free components F, so Newton's method
/// works on the displacements alone. A semismooth Newton method (primal-dual
/// active set, Hueber and Wohlmuth) decides the contact status from the
/// complementarity function: node j is in contact where
/// \f$p_j - c\,\tilde g_j / D_j > 0\f$, \f$c = E / h\f$. In contact its
/// normal equation becomes \f$\tilde g_j = 0\f$; open, \f$p_j = 0\f$. With
/// Coulomb friction (coefficient \f$\mu\f$) and the tangential traction
/// \f$\lambda_\tau\f$ taken in the direction of the slip, a contact node
/// sticks where \f$|\lambda_\tau + c\,\tilde u_\tau / D_j| < \mu p_j\f$ -
/// its weighted slip increment over the step \f$\tilde u_\tau\f$ is then zero
/// - and slips otherwise, with \f$\lambda_\tau = \mu p_j\f$ along
/// \f$\lambda_\tau + c\,\tilde u_\tau / D_j\f$. The master equations take the
/// condensed contact forces of their slave nodes. The system is not
/// symmetric; it is solved by sparse LU. For linear elasticity without
/// friction the active set is exact after finitely many iterations and the
/// solution then exact to round-off.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Selector.hpp"

#include <string>
#include <vector>

namespace sparlab {

class FemModel;

/// A rigid obstacle given analytically. Its signed distance function is
/// positive on the body's side of its surface.
struct RigidObstacle {
  enum class Kind { Plane, Cylinder, Sphere };
  Kind kind = Kind::Plane;
  /// Plane: a point of it. Cylinder: a point of its axis. Sphere: its
  /// centre [m].
  Vector3 point = Vector3::Zero();
  /// Plane: its normal, pointing out of the obstacle towards the body.
  /// Cylinder: its axis (3-D; in 2-D the cylinder is a circle in the model
  /// plane and the axis is z). Normalised on use.
  Vector3 direction = Vector3::UnitY();
  Scalar radius = 0.0;  ///< cylinder, sphere [m]
  /// The body lies inside the cylinder or sphere - a rigid cavity - rather
  /// than outside it.
  bool inside = false;
  /// Rigid translation of the obstacle at load factor 1 [m]; along the path
  /// it has moved by lambda times this.
  Vector3 motion = Vector3::Zero();

  /// Signed distance of `x` from the surface (positive on the body's side)
  /// and, in `normal`, its gradient: the unit normal pointing to the body's
  /// side. `dim` 2 measures a cylinder in the x-y plane.
  /// \throws ModelError at the axis of a cylinder or the centre of a sphere,
  ///         where the normal is undefined.
  Scalar gap(const Vector3& x, Vector3& normal, int dim) const;
};

std::string to_string(RigidObstacle::Kind kind);
/// "plane", "cylinder" or "sphere".
RigidObstacle::Kind parse_obstacle_kind(const std::string& text);

/// One contact pair: a slave surface against a rigid obstacle or a master
/// surface.
struct ContactPairSpec {
  std::string name;
  /// The slave surface: every boundary face whose nodes all lie in the
  /// region (as a pressure selects its faces).
  SelectorGroup slave;
  /// Against the rigid obstacle below, or (false) the master surface.
  bool rigid = true;
  RigidObstacle obstacle;
  SelectorGroup master;
  /// Coulomb friction coefficient; 0 is frictionless.
  Scalar friction = 0.0;
};

struct ContactOptions {
  bool enabled = false;
  std::vector<ContactPairSpec> pairs;
  /// Scales the complementarity parameter c = E / h of the active-set test
  /// (E the largest Young's modulus next to a slave node, h its face size).
  /// It weighs gap against pressure while the set is still changing; the
  /// converged solution does not depend on it.
  Scalar complementarity = 1.0;
  /// Master faces are searched for within this multiple of a slave face's
  /// size (plus its initial gap to them).
  Scalar search_factor = 2.0;
};

enum class ContactStatus { Open, Stick, Slip };
std::string to_string(ContactStatus status);

/// One slave node at a converged state.
struct ContactNodeResult {
  Index node = 0;
  std::size_t pair = 0;
  Scalar weight = 0.0;                ///< D_j = int N_j dA [m^2] (times the thickness in 2-D)
  Vector3 normal = Vector3::Zero();   ///< nu_j: direction of the pressure on the slave body
  Scalar gap = 0.0;                   ///< weighted gap / D_j [m]; 0 in contact
  Scalar pressure = 0.0;              ///< [Pa], >= 0
  Vector3 traction = Vector3::Zero(); ///< tangential traction on the slave body [Pa]
  Vector3 slip = Vector3::Zero();     ///< accumulated slip of slave relative to master [m]
  ContactStatus status = ContactStatus::Open;
};

/// A contact pair at a converged state.
struct ContactPairResult {
  std::string name;
  bool rigid = true;
  Scalar friction = 0.0;
  int nodes = 0;          ///< slave nodes that take part
  int excluded = 0;       ///< slave nodes left out (see `warnings` of the analysis)
  int active = 0;         ///< nodes in contact (sticking or slipping)
  int sticking = 0;       ///< with friction
  int slipping = 0;
  Scalar area = 0.0;      ///< sum of D_j over the nodes in contact [m^2]
  /// Resultant contact force on the slave body [N]: sum of D_j (p_j nu_j + t_j).
  Vector3 force = Vector3::Zero();
  Scalar max_pressure = 0.0;       ///< [Pa]
  Scalar min_gap = 0.0;            ///< smallest gap over the slave nodes [m] (negative: penetration)
  Scalar max_slip = 0.0;           ///< largest accumulated slip [m]
};

/// The contact conditions of a model: the geometry of every pair, built once
/// in the reference configuration, and the linearised equations of one
/// Newton iteration.
class ContactProblem {
 public:
  /// \throws ConfigError for an empty surface, a quadratic (Tet10) face,
  ///         a node on two slave surfaces, a node that is both slave and
  ///         master, a non-positive radius or a negative friction coefficient.
  ContactProblem(const FemModel& model, const ContactOptions& options);

  /// One slave node that takes part in contact.
  struct Node {
    Index node = 0;
    std::size_t pair = 0;
    Scalar weight = 0.0;            ///< D_j
    Vector3 normal = Vector3::Zero();  ///< nu_j
    Scalar initial_gap = 0.0;       ///< weighted: g~0_j [m^3 in 3-D]
    /// Mortar pair: master nodes and M_jl, with sum_l M_jl = D_j.
    std::vector<std::pair<Index, Scalar>> masters;
    /// Free components of the node, the normal restricted to them and an
    /// orthonormal basis of the free directions normal to it (the
    /// tangential directions; one in 2-D, two in 3-D with every component
    /// free).
    std::vector<int> free;
    Vector3 free_normal = Vector3::Zero();
    std::vector<Vector3> tangents;
    Scalar complementarity = 0.0;   ///< c [Pa/m]
    Scalar friction = 0.0;          ///< mu
    Vector3 motion = Vector3::Zero();  ///< obstacle motion at lambda = 1 (rigid pair)
  };

  const std::vector<Node>& nodes() const { return nodes_; }
  const ContactOptions& options() const { return options_; }
  /// Per pair: slave nodes left out, with the reason (for the warnings).
  const std::vector<std::string>& exclusions() const { return exclusions_; }
  bool frictional() const { return frictional_; }

  /// The weighted gap of node i at the full displacement `u` and load
  /// factor `lambda`.
  Scalar weighted_gap(std::size_t i, const Vector& u, Scalar lambda) const;

  /// The condensed Newton system of one iteration over the free DOFs.
  struct Linearization {
    SparseMatrix matrix;        ///< free x free, not symmetric
    Vector rhs;                 ///< -(residual of the condensed equations)
    std::vector<ContactStatus> status;  ///< per node of nodes()
    /// Norm over the nodes in contact of their contact forces
    /// D_j (p_j nu_j + t_j) and of the forces c g~_j that would close their
    /// gaps: part of the scale the residual is judged against.
    Scalar force_scale = 0.0;
  };

  /// Linearise at the state `u` (full, prescribed values included), load
  /// factor `lambda`, with the full residual `residual` = f_int - f_ext and
  /// the full tangent; `u_start` and `lambda_start` are the last converged
  /// state, from which the slip of the step is measured.
  Linearization linearize(const Vector& u, Scalar lambda, const Vector& residual,
                          const SparseMatrix& tangent, const Vector& u_start,
                          Scalar lambda_start) const;

  /// The contact state of every node at a converged state, with the slip
  /// accumulated up to the last commit plus that of the current step.
  std::vector<ContactNodeResult> node_results(const Vector& u, Scalar lambda,
                                              const Vector& residual,
                                              const std::vector<ContactStatus>& status,
                                              const Vector& u_start,
                                              Scalar lambda_start) const;
  /// Pair totals of node_results.
  std::vector<ContactPairResult> pair_results(const std::vector<ContactNodeResult>& nodes) const;

  /// Accept a converged step: accumulate the slip of the nodes in contact.
  void commit(const Vector& u, Scalar lambda, const std::vector<ContactStatus>& status,
              const Vector& u_start, Scalar lambda_start);

 private:
  /// Contact pressure and tangential traction (in the tangent basis) of
  /// node i from the residual.
  void tractions(const Node& n, const Vector& residual, Scalar& pressure,
                 Eigen::Vector2d& tangential) const;
  /// Weighted tangential slip of node i over the step, in its tangent basis.
  Eigen::Vector2d weighted_slip(const Node& n, const Vector& u, Scalar lambda,
                                const Vector& u_start, Scalar lambda_start) const;

  const FemModel& model_;
  ContactOptions options_;
  std::vector<Node> nodes_;
  std::vector<std::string> exclusions_;
  std::vector<int> excluded_;  ///< per pair: slave nodes left out
  std::vector<Vector3> slip_;  ///< accumulated slip per node, global [m]
  bool frictional_ = false;
};

}  // namespace sparlab
