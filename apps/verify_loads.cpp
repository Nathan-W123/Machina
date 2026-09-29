/// \file verify_loads.cpp
/// \brief Verification of the pressure, volume and thermal loads against exact
///        solutions of elasticity and heat conduction, under mesh refinement.
///
/// Studies:
///   * `lame-cylinder`    a thick cylinder under internal pressure in plane
///                        strain (Lame): Q4 and Tri3 meshes of a quarter
///                        section, and one layer of Hex8 and curved Tet10 cells
///                        held at u_z = 0;
///   * `rotating-disk`    an annular disk spinning about its axis in plane
///                        stress (Q4, Tri3) and a long spinning cylinder in
///                        plane strain (Hex8, Tet10): the centrifugal load;
///   * `thermal-cylinder` steady conduction in a thick cylinder with uniform
///                        heat generation, a fixed bore temperature and
///                        convection from its outer surface, then the thermal
///                        stress of that temperature in plane strain;
///   * `bimetal-strip`    a free two-material strip heated uniformly, whose
///                        curvature is Timoshenko's (1925);
///   * `self-weight`      a bar hanging under its own weight, whose quadratic
///                        exact field the Tet10 space contains.
///
/// Every reference is an exact solution of the continuum problem that the
/// model discretises - no beam or plate idealisation stands between them - so
/// the error must vanish under refinement at the rate of the element: O(h^2)
/// in the displacements and temperatures of the linear elements, O(h^3) for
/// the Tet10, whose edge nodes lie on the curved surfaces. The measured error
/// is the RMS nodal error relative to the RMS exact field, a discrete L2 norm,
/// whose rate the theory gives without the log factor and corner effects of
/// the maximum norm; the largest nodal error is reported beside it. Each study
/// measures the rate between the two finest meshes and checks it against the
/// element's order, less a margin of 0.3 for the pre-asymptotic range: a load
/// integrated wrongly leaves an error that does not vanish, and the measured
/// rate falls towards zero.
///
/// The curved models are quarter sections with symmetry supports (u_y = 0 on
/// theta = 0, u_x = 0 on theta = 90 deg). A 3-D section is one layer of cells
/// with u_z = 0 at every node, which is exactly plane strain.

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Exact axisymmetric solutions
// ---------------------------------------------------------------------------

/// Radial temperature change dT(r) = t0 + t_log ln(r/a) + t_r2 r^2 [K]: the
/// general steady solution of -k (T'' + T'/r) = Q with uniform generation Q
/// (t_r2 = -Q / 4k), which includes the logarithmic and the uniform profile.
struct RadialTemperature {
  Scalar a = 1.0;      ///< radius in the logarithm [m]
  Scalar t0 = 0.0;     ///< [K]
  Scalar t_log = 0.0;  ///< [K]
  Scalar t_r2 = 0.0;   ///< [K/m^2]

  Scalar operator()(Scalar r) const { return t0 + t_log * std::log(r / a) + t_r2 * r * r; }
  Scalar derivative(Scalar r) const { return t_log / r + 2.0 * t_r2 * r; }
  /// \f$I(r) = \int_a^r \Delta T(s)\,s\,ds\f$, in closed form.
  Scalar weighted_integral(Scalar r) const { return antiderivative(r) - antiderivative(a); }

 private:
  Scalar antiderivative(Scalar s) const {
    const Scalar s2 = s * s;
    return 0.5 * t0 * s2 + t_log * (0.5 * s2 * std::log(s / a) - 0.25 * s2) +
           0.25 * t_r2 * s2 * s2;
  }
};

/// Exact small-strain solution of a hollow disk in plane stress, or of a long
/// hollow cylinder in plane strain (eps_zz = 0), with inner radius a and outer
/// radius b, under pressures on its two surfaces, a steady rotation about its
/// axis (the body force rho omega^2 r) and a radial temperature change.
/// Equilibrium in terms of the radial displacement,
/// \f[
///   \frac{d}{dr}\Big[\frac{1}{r}\frac{d(r u)}{dr}\Big]
///     = c_T\,\frac{d\Delta T}{dr} - c_b\,\rho\omega^2 r ,
/// \f]
/// with \f$c_T = (1+\nu)\alpha\f$, \f$c_b = (1-\nu^2)/E\f$ in plane stress and
/// \f$c_T = (1+\nu)\alpha/(1-\nu)\f$, \f$c_b = (1+\nu)(1-2\nu)/(E(1-\nu))\f$
/// in plane strain, integrates to
/// \f[
///   u = \frac{c_T I(r) - c_b\rho\omega^2 (r^4 - a^4)/8}{r} + C_1 r + \frac{C_2}{r},
///   \qquad I(r) = \int_a^r \Delta T\,s\,ds ,
/// \f]
/// and \f$C_1, C_2\f$ follow from \f$\sigma_r(a) = -p_a\f$,
/// \f$\sigma_r(b) = -p_b\f$. Lame's cylinder, the rotating disk and the thermal
/// stresses of a long cylinder (Timoshenko and Goodier, Theory of Elasticity,
/// sections 33, 32 and 150) are special cases; the studies check the first
/// two against the textbook closed forms and every case against equilibrium
/// and its boundary conditions.
class AxisymmetricSolution {
 public:
  struct Data {
    bool plane_strain = false;
    Scalar youngs = 1.0;
    Scalar poisson = 0.0;
    Scalar alpha = 0.0;
    Scalar a = 1.0;
    Scalar b = 2.0;
    Scalar p_inner = 0.0;     ///< [Pa]
    Scalar p_outer = 0.0;     ///< [Pa]
    Scalar rho_omega2 = 0.0;  ///< rho omega^2 [N/m^4]
    RadialTemperature delta_t;
  };

  explicit AxisymmetricSolution(const Data& d) : d_(d) {
    const Scalar e = d.youngs;
    const Scalar nu = d.poisson;
    if (d.plane_strain) {
      c_t_ = (1.0 + nu) * d.alpha / (1.0 - nu);
      c_b_ = (1.0 + nu) * (1.0 - 2.0 * nu) / (e * (1.0 - nu));
    } else {
      c_t_ = (1.0 + nu) * d.alpha;
      c_b_ = (1.0 - nu * nu) / e;
    }
    // sigma_r is affine in (C1, C2): evaluate it for three sets of constants.
    Eigen::Matrix2d m;
    Eigen::Vector2d rhs;
    const Scalar radii[2] = {d.a, d.b};
    const Scalar pressures[2] = {d.p_inner, d.p_outer};
    for (int i = 0; i < 2; ++i) {
      const Scalar s0 = radial_stress(radii[i], 0.0, 0.0);
      m(i, 0) = radial_stress(radii[i], 1.0, 0.0) - s0;
      m(i, 1) = radial_stress(radii[i], 0.0, 1.0) - s0;
      rhs(i) = -pressures[i] - s0;
    }
    const Eigen::Vector2d c = m.fullPivLu().solve(rhs);
    c1_ = c(0);
    c2_ = c(1);
  }

  const Data& data() const { return d_; }
  Scalar u(Scalar r) const { return displacement(r, c1_, c2_); }
  Scalar sigma_r(Scalar r) const { return radial_stress(r, c1_, c2_); }
  Scalar sigma_theta(Scalar r) const {
    Scalar sr = 0.0;
    Scalar st = 0.0;
    stresses(r, c1_, c2_, sr, st);
    return st;
  }
  /// Out-of-plane stress: \f$\nu(\sigma_r+\sigma_\theta) - E\alpha\Delta T\f$ in
  /// plane strain, zero in plane stress.
  Scalar sigma_z(Scalar r) const {
    if (!d_.plane_strain) return 0.0;
    Scalar sr = 0.0;
    Scalar st = 0.0;
    stresses(r, c1_, c2_, sr, st);
    return d_.poisson * (sr + st) - d_.youngs * d_.alpha * d_.delta_t(r);
  }

  /// Largest residual of the radial equilibrium equation
  /// \f$\sigma_r' + (\sigma_r - \sigma_\theta)/r + \rho\omega^2 r = 0\f$ over
  /// the wall, by central differences, times (b - a) over the largest stress,
  /// and the boundary-condition mismatch over the same stress: a check of the
  /// reference itself.
  Scalar self_check() const {
    const Scalar span = d_.b - d_.a;
    const Scalar step = 1.0e-4 * span;
    Scalar scale = 0.0;
    for (int i = 0; i <= 50; ++i) {
      const Scalar r = d_.a + span * i / 50.0;
      scale = std::max({scale, std::abs(sigma_r(r)), std::abs(sigma_theta(r)),
                        std::abs(sigma_z(r))});
    }
    Scalar worst = std::max(std::abs(sigma_r(d_.a) + d_.p_inner),
                            std::abs(sigma_r(d_.b) + d_.p_outer)) / scale;
    for (int i = 1; i < 50; ++i) {
      const Scalar r = d_.a + span * i / 50.0;
      const Scalar dsr = (sigma_r(r + step) - sigma_r(r - step)) / (2.0 * step);
      const Scalar residual = dsr + (sigma_r(r) - sigma_theta(r)) / r + d_.rho_omega2 * r;
      worst = std::max(worst, std::abs(residual) * span / scale);
    }
    return worst;
  }

 private:
  Scalar displacement(Scalar r, Scalar c1, Scalar c2) const {
    const Scalar a4 = d_.a * d_.a * d_.a * d_.a;
    return (c_t_ * d_.delta_t.weighted_integral(r) -
            c_b_ * d_.rho_omega2 * (r * r * r * r - a4) / 8.0) / r +
           c1 * r + c2 / r;
  }
  Scalar displacement_derivative(Scalar r, Scalar c1, Scalar c2) const {
    const Scalar a4 = d_.a * d_.a * d_.a * d_.a;
    return c_t_ * (d_.delta_t(r) - d_.delta_t.weighted_integral(r) / (r * r)) -
           c_b_ * d_.rho_omega2 * (3.0 * r * r * r * r + a4) / (8.0 * r * r) + c1 -
           c2 / (r * r);
  }
  void stresses(Scalar r, Scalar c1, Scalar c2, Scalar& sr, Scalar& st) const {
    const Scalar e = d_.youngs;
    const Scalar nu = d_.poisson;
    const Scalar er = displacement_derivative(r, c1, c2);
    const Scalar et = displacement(r, c1, c2) / r;
    const Scalar dt = d_.delta_t(r);
    if (d_.plane_strain) {
      const Scalar f = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
      const Scalar thermal = e * d_.alpha * dt / (1.0 - 2.0 * nu);
      sr = f * ((1.0 - nu) * er + nu * et) - thermal;
      st = f * (nu * er + (1.0 - nu) * et) - thermal;
    } else {
      const Scalar f = e / (1.0 - nu * nu);
      const Scalar thermal = (1.0 + nu) * d_.alpha * dt;
      sr = f * (er + nu * et - thermal);
      st = f * (nu * er + et - thermal);
    }
  }
  Scalar radial_stress(Scalar r, Scalar c1, Scalar c2) const {
    Scalar sr = 0.0;
    Scalar st = 0.0;
    stresses(r, c1, c2, sr, st);
    return sr;
  }

  Data d_;
  Scalar c_t_ = 0.0;
  Scalar c_b_ = 0.0;
  Scalar c1_ = 0.0;
  Scalar c2_ = 0.0;
};

/// Largest difference between two radial profiles over the wall, relative to
/// the largest magnitude of the second.
Scalar profile_difference(Scalar a, Scalar b, const std::function<Scalar(Scalar)>& f,
                          const std::function<Scalar(Scalar)>& g) {
  Scalar diff = 0.0;
  Scalar scale = 0.0;
  for (int i = 0; i <= 40; ++i) {
    const Scalar r = a + (b - a) * i / 40.0;
    diff = std::max(diff, std::abs(f(r) - g(r)));
    scale = std::max(scale, std::abs(g(r)));
  }
  return diff / scale;
}

Selector box_selector() {
  Selector s;
  s.kind = SelectorKind::Box;
  return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Quarter sections of a cylinder (shared with verify_nonlinear.cpp)
// ---------------------------------------------------------------------------

Mesh quarter_annulus(const Mesh& unit, Scalar a, Scalar b) {
  Matrix x = unit.coordinates();
  for (Index n = 0; n < unit.num_nodes(); ++n) {
    const Scalar r = a + (b - a) * x(0, n);
    const Scalar eta = x(1, n);
    // Exact zeros on the two symmetry planes, so their supports and the
    // selectors see every node there.
    const Scalar c = std::abs(eta - 1.0) < 1.0e-12 ? 0.0 : std::cos(0.5 * kPi * eta);
    const Scalar s = std::abs(eta) < 1.0e-12 ? 0.0 : std::sin(0.5 * kPi * eta);
    x(0, n) = r * c;
    x(1, n) = r * s;
  }
  Mesh mapped(std::move(x), unit.connectivity(), unit.element_type());
  mapped.validate();
  return mapped;
}

Mesh sector_mesh(ElementType type, Index nr, Index nt, Scalar a, Scalar b) {
  StructuredMeshSpec ms;
  ms.nx = nr;
  ms.ny = nt;
  ms.nz = 1;
  ms.lx = 1.0;
  ms.ly = 1.0;
  ms.lz = (b - a) / static_cast<Scalar>(nr);
  switch (type) {
    case ElementType::Quad4: return quarter_annulus(make_structured_quad_mesh(ms), a, b);
    case ElementType::Tri3: return quarter_annulus(make_structured_tri_mesh(ms), a, b);
    case ElementType::Hex8: return quarter_annulus(make_structured_hex_mesh(ms), a, b);
    case ElementType::Tet4: return quarter_annulus(make_structured_tet_mesh(ms), a, b);
    case ElementType::Tet10: return quarter_annulus(make_structured_tet10_mesh(ms), a, b);
  }
  throw ConfigError("unhandled element type");
}

void add_quarter_supports(FemModel& model) {
  DisplacementConstraint theta0;
  theta0.region.name = "theta_0";
  Selector y0 = box_selector();
  y0.ymax = 0.0;
  theta0.region.members.push_back(y0);
  theta0.fix_y = true;
  model.constraints().push_back(theta0);

  DisplacementConstraint theta90;
  theta90.region.name = "theta_90";
  Selector x0 = box_selector();
  x0.xmax = 0.0;
  theta90.region.members.push_back(x0);
  theta90.fix_x = true;
  model.constraints().push_back(theta90);

  if (model.dim() == 3) {
    DisplacementConstraint plane;
    plane.region.name = "plane_strain";
    Selector all;
    all.kind = SelectorKind::All;
    plane.region.members.push_back(all);
    plane.fix_z = true;
    model.constraints().push_back(plane);
  }
}

SelectorGroup cylinder_surface(const std::string& name, Scalar radius, bool bore) {
  SelectorGroup g;
  g.name = name;
  Selector s;
  if (bore) {
    s.kind = SelectorKind::Circle;
    s.radius = radius;
  } else {
    s.kind = SelectorKind::Annulus;
    s.inner_radius = radius;
    s.radius = 2.0 * radius;
  }
  g.members.push_back(s);
  return g;
}

namespace {

SelectorGroup whole_model_region() {
  SelectorGroup g;
  g.name = "all";
  Selector s;
  s.kind = SelectorKind::All;
  g.members.push_back(s);
  return g;
}

/// A finalised model, its solution and its stresses.
struct SolvedModel {
  StaticSolution solution;
  StressField field;
  Vector temperature;  ///< empty without a temperature field
  ConductionSummary conduction;
  Index num_dofs = 0;
};

/// Solve with the sparse Cholesky factorisation, or with multigrid CG (to a
/// relative residual of 1e-12) for a solid model above 50 000 unknowns, whose
/// factorisation would dominate the run time.
SolvedModel solve_model(FemModel& model) {
  model.finalize();
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.type = model.dim() == 3 && model.dofs().num_free() > 50000
                            ? LinearSolverType::AmgCg
                            : LinearSolverType::SimplicialLdlt;
  StaticAnalysis analysis(model, assembler, options);
  SolvedModel out;
  out.solution = analysis.solve_all().front();
  const LoadCaseData& data = model.load_case_data(0);
  out.temperature = data.temperature;
  out.conduction = data.conduction;
  out.field = recover_stresses(model, assembler, out.solution.displacement, nullptr,
                               data.temperature.size() > 0 ? &data.temperature : nullptr);
  out.num_dofs = model.dofs().num_dofs();
  return out;
}

/// Errors of a quarter-section solution against an axisymmetric reference.
/// The RMS errors are relative discrete L2 norms over the nodes,
/// \f$\sqrt{\sum_n |u_h - u|^2 / \sum_n |u|^2}\f$; the max errors divide the
/// largest nodal error by the largest exact value.
struct SectionErrors {
  Scalar displacement = 0.0;      ///< RMS nodal displacement error
  Scalar displacement_max = 0.0;  ///< max nodal |u_h - u| / max |u|
  Scalar stress = 0.0;            ///< max element error of sigma_r, sigma_theta,
                                  ///< sigma_z at the centroid / max exact stress
  Scalar hoop = 0.0;              ///< the sigma_theta part of `stress`
  Scalar bore_displacement = 0.0;  ///< mean u_r of the bore nodes [m]
  Scalar temperature = 0.0;       ///< RMS nodal error of T - T_ref
  Scalar temperature_max = 0.0;   ///< max nodal |T_h - T| / (max T - min T)
};

SectionErrors section_errors(const FemModel& model, const SolvedModel& solved,
                             const AxisymmetricSolution& exact,
                             const std::function<Scalar(Scalar)>* temperature = nullptr) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  const Vector& u = solved.solution.displacement;
  const Scalar a = exact.data().a;
  SectionErrors out;
  Scalar u_scale = 0.0;
  Scalar u_error = 0.0;
  Scalar u_sq = 0.0;
  Scalar u_error_sq = 0.0;
  Scalar bore_sum = 0.0;
  Index bore_nodes = 0;
  Scalar t_min = 1.0e300;
  Scalar t_max = -1.0e300;
  Scalar t_error = 0.0;
  Scalar t_sq = 0.0;
  Scalar t_error_sq = 0.0;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    const Scalar r = std::hypot(x.x(), x.y());
    const Scalar c = x.x() / r;
    const Scalar s = x.y() / r;
    const Scalar ur = exact.u(r);
    const Scalar ux = u(dim * n);
    const Scalar uy = u(dim * n + 1);
    const Scalar err = std::hypot(ux - ur * c, uy - ur * s);
    u_scale = std::max(u_scale, std::abs(ur));
    u_error = std::max(u_error, err);
    u_sq += ur * ur;
    u_error_sq += err * err;
    if (r <= a * (1.0 + 1.0e-9)) {
      bore_sum += c * ux + s * uy;
      ++bore_nodes;
    }
    if (temperature != nullptr) {
      // The RMS error is taken relative to the temperature change T - T_ref,
      // which drives the thermal strain, not to the absolute temperature.
      const Scalar t = (*temperature)(r);
      const Scalar change = exact.data().delta_t(r);
      const Scalar t_err = solved.temperature(n) - t;
      t_min = std::min(t_min, t);
      t_max = std::max(t_max, t);
      t_error = std::max(t_error, std::abs(t_err));
      t_error_sq += t_err * t_err;
      t_sq += change * change;
    }
  }
  out.displacement = std::sqrt(u_error_sq / u_sq);
  out.displacement_max = u_error / u_scale;
  out.bore_displacement = bore_nodes > 0 ? bore_sum / static_cast<Scalar>(bore_nodes) : 0.0;
  if (temperature != nullptr) {
    out.temperature = std::sqrt(t_error_sq / t_sq);
    out.temperature_max = t_error / (t_max - t_min);
  }

  Scalar s_scale = 0.0;
  Scalar s_error = 0.0;
  Scalar hoop_error = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector3 xc = mesh.element_centroid(e);
    const Scalar r = std::hypot(xc.x(), xc.y());
    const Scalar c = xc.x() / r;
    const Scalar s = xc.y() / r;
    const Scalar sxx = solved.field.element_stress(0, e);
    const Scalar syy = solved.field.element_stress(1, e);
    const Scalar sxy = solved.field.element_stress(dim == 2 ? 2 : 3, e);
    const Scalar szz = dim == 3 ? solved.field.element_stress(2, e)
                       : solved.field.element_sigma_zz.size() > 0
                           ? solved.field.element_sigma_zz(e)
                           : 0.0;
    const Scalar srr_h = c * c * sxx + s * s * syy + 2.0 * c * s * sxy;
    const Scalar stt_h = s * s * sxx + c * c * syy - 2.0 * c * s * sxy;
    const Scalar srr = exact.sigma_r(r);
    const Scalar stt = exact.sigma_theta(r);
    const Scalar szz_exact = exact.sigma_z(r);
    s_scale = std::max({s_scale, std::abs(srr), std::abs(stt), std::abs(szz_exact)});
    hoop_error = std::max(hoop_error, std::abs(stt_h - stt));
    s_error = std::max({s_error, std::abs(srr_h - srr), std::abs(stt_h - stt),
                        std::abs(szz - szz_exact)});
  }
  out.stress = s_error / s_scale;
  out.hoop = hoop_error / s_scale;
  return out;
}

/// Order of accuracy of an element's nodal displacements (and temperatures).
int element_order(ElementType type) { return type == ElementType::Tet10 ? 3 : 2; }

/// Refinement ladders of the section studies: n_r cells through the wall,
/// twice as many round the quarter.
std::vector<Index> section_ladder(ElementType type) {
  switch (type) {
    case ElementType::Quad4:
    case ElementType::Tri3: return {4, 8, 16, 32, 64};
    case ElementType::Hex8: return {4, 8, 16, 32};
    case ElementType::Tet4: return {4, 8, 16, 32};
    case ElementType::Tet10: return {2, 4, 8, 16, 32};
  }
  return {};
}

/// Convergence bookkeeping of one element type.
struct Ladder {
  std::vector<Scalar> h;
  std::vector<Scalar> displacement;      ///< RMS
  std::vector<Scalar> displacement_max;
  std::vector<Scalar> stress;
  std::vector<Scalar> temperature;       ///< RMS
  std::vector<Scalar> temperature_max;
  std::vector<Scalar> extra;  ///< a study-specific scalar error

  static Scalar last_order(const std::vector<Scalar>& hs, const std::vector<Scalar>& e) {
    const std::size_t n = e.size();
    if (n < 2) return 0.0;
    return observed_order(hs[n - 2], e[n - 2], hs[n - 1], e[n - 1]);
  }
  Scalar order(const std::vector<Scalar>& e) const { return last_order(h, e); }
  /// Order between the latest two entries, or NaN for the first one.
  Scalar running_order(const std::vector<Scalar>& e) const {
    return e.size() < 2 ? std::nan("") : order(e);
  }
};

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

/// Run a quarter-section study over the element types and their ladders,
/// writing one CSV row per mesh. `build` makes the loaded model on a mesh.
/// Returns the ladders by element name.
std::map<std::string, Ladder> run_sections(
    const std::string& csv_path, const std::vector<ElementType>& types, Scalar a, Scalar b,
    const std::function<FemModel(Mesh)>& build,
    const std::function<const AxisymmetricSolution&(int dim)>& exact_for,
    const std::function<Scalar(Scalar)>* temperature, json::Value& block,
    const std::function<Scalar(const FemModel&, const SolvedModel&)>& extra = nullptr) {
  std::vector<std::string> header = {
      "element", "n_r", "n_theta", "h[m]", "num_dofs", "u_rms_error[-]", "u_rms_order[-]",
      "u_max_error[-]", "u_max_order[-]", "stress_error[-]", "stress_order[-]",
      "hoop_error[-]", "u_r_bore[m]", "u_r_bore_exact[m]", "solver"};
  if (temperature != nullptr) {
    header.push_back("T_rms_error[-]");
    header.push_back("T_rms_order[-]");
    header.push_back("T_max_error[-]");
    header.push_back("T_max_order[-]");
  }
  if (extra) {
    header.push_back("extra_error[-]");
    header.push_back("extra_order[-]");
  }
  CsvWriter csv(csv_path, header);
  std::map<std::string, Ladder> ladders;
  for (const ElementType type : types) {
    const std::string name = to_string(type);
    Ladder& ladder = ladders[name];
    json::Value records = json::Value::make_array();
    for (const Index nr : section_ladder(type)) {
      const Index nt = 2 * nr;
      Mesh mesh = sector_mesh(type, nr, nt, a, b);
      const int dim = mesh.dim();
      FemModel model = build(std::move(mesh));
      const SolvedModel solved = solve_model(model);
      const AxisymmetricSolution& exact = exact_for(dim);
      const SectionErrors err = section_errors(model, solved, exact, temperature);
      ladder.h.push_back((b - a) / static_cast<Scalar>(nr));
      ladder.displacement.push_back(err.displacement);
      ladder.displacement_max.push_back(err.displacement_max);
      ladder.stress.push_back(err.stress);
      if (temperature != nullptr) {
        ladder.temperature.push_back(err.temperature);
        ladder.temperature_max.push_back(err.temperature_max);
      }
      Scalar extra_error = 0.0;
      if (extra) {
        extra_error = extra(model, solved);
        ladder.extra.push_back(extra_error);
      }
      std::vector<std::string> row = {
          name,
          fmt(static_cast<Scalar>(nr)),
          fmt(static_cast<Scalar>(nt)),
          fmt(ladder.h.back(), 8),
          fmt(static_cast<Scalar>(solved.num_dofs), 9),
          fmt(err.displacement, 8),
          fmt(ladder.running_order(ladder.displacement), 4),
          fmt(err.displacement_max, 8),
          fmt(ladder.running_order(ladder.displacement_max), 4),
          fmt(err.stress, 8),
          fmt(ladder.running_order(ladder.stress), 4),
          fmt(err.hoop, 8),
          fmt(err.bore_displacement, 10),
          fmt(exact.u(a), 10),
          solved.solution.solver_name};
      if (temperature != nullptr) {
        row.push_back(fmt(err.temperature, 8));
        row.push_back(fmt(ladder.running_order(ladder.temperature), 4));
        row.push_back(fmt(err.temperature_max, 8));
        row.push_back(fmt(ladder.running_order(ladder.temperature_max), 4));
      }
      if (extra) {
        row.push_back(fmt(extra_error, 8));
        row.push_back(fmt(ladder.running_order(ladder.extra), 4));
      }
      csv.raw_row(row);
      json::Value rec = json::Value::make_object();
      rec.set("n_r", json::Value::make_number(static_cast<Scalar>(nr)));
      rec.set("n_theta", json::Value::make_number(static_cast<Scalar>(nt)));
      rec.set("num_dofs", json::Value::make_number(static_cast<Scalar>(solved.num_dofs)));
      rec.set("displacement_rms_error", json::Value::make_number(err.displacement));
      rec.set("displacement_max_error", json::Value::make_number(err.displacement_max));
      rec.set("stress_error", json::Value::make_number(err.stress));
      rec.set("hoop_stress_error", json::Value::make_number(err.hoop));
      rec.set("bore_radial_displacement_m", json::Value::make_number(err.bore_displacement));
      if (temperature != nullptr) {
        rec.set("temperature_rms_error", json::Value::make_number(err.temperature));
        rec.set("temperature_max_error", json::Value::make_number(err.temperature_max));
      }
      if (extra) rec.set("extra_error", json::Value::make_number(extra_error));
      rec.set("force_balance_error",
              json::Value::make_number(solved.solution.equilibrium.relative_force_error));
      records.push_back(rec);
    }
    json::Value entry = json::Value::make_object();
    entry.set("meshes", records);
    entry.set("displacement_rms_order",
              json::Value::make_number(ladder.order(ladder.displacement)));
    entry.set("displacement_max_order",
              json::Value::make_number(ladder.order(ladder.displacement_max)));
    entry.set("stress_order", json::Value::make_number(ladder.order(ladder.stress)));
    if (temperature != nullptr) {
      entry.set("temperature_rms_order",
                json::Value::make_number(ladder.order(ladder.temperature)));
      entry.set("temperature_max_order",
                json::Value::make_number(ladder.order(ladder.temperature_max)));
    }
    if (extra) entry.set("extra_order", json::Value::make_number(ladder.order(ladder.extra)));
    entry.set("expected_order", json::Value::make_number(element_order(type)));
    block.set(name, entry);
  }
  csv.close();
  return ladders;
}

/// Largest shortfall of the measured order below the element's order over
/// the element types, for the chosen error series; the note lists them all.
Scalar order_shortfall(const std::map<std::string, Ladder>& ladders,
                       const std::function<const std::vector<Scalar>&(const Ladder&)>& series,
                       const std::string& label, std::ostringstream& note) {
  Scalar worst = -1.0e300;
  for (const auto& kv : ladders) {
    const ElementType type = kv.first == "Quad4"  ? ElementType::Quad4
                             : kv.first == "Tri3" ? ElementType::Tri3
                             : kv.first == "Hex8" ? ElementType::Hex8
                             : kv.first == "Tet4" ? ElementType::Tet4
                                                  : ElementType::Tet10;
    const std::vector<Scalar>& e = series(kv.second);
    const Scalar p = kv.second.order(e);
    worst = std::max(worst, static_cast<Scalar>(element_order(type)) - p);
    note << kv.first << " " << label << " " << app::format(e.back(), 3) << " (order "
         << app::format(p, 3) << "); ";
  }
  return worst;
}

constexpr Scalar kOrderMargin = 0.3;

}  // namespace

// ---------------------------------------------------------------------------
// Lame's thick cylinder
// ---------------------------------------------------------------------------

StudyOutcome study_lame_cylinder(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.2;
  const Scalar p = 100.0e6;
  const Scalar e = 200.0e9;
  const Scalar nu = 0.3;
  AxisymmetricSolution::Data d;
  d.plane_strain = true;
  d.youngs = e;
  d.poisson = nu;
  d.a = a;
  d.b = b;
  d.p_inner = p;
  const AxisymmetricSolution exact(d);

  // Lame's closed form (plane strain): sigma_r = A - B/r^2,
  // sigma_theta = A + B/r^2, u = (1 + nu)/E [(1 - 2 nu) A r + B / r].
  const Scalar big_a = p * a * a / (b * b - a * a);
  const Scalar big_b = p * a * a * b * b / (b * b - a * a);
  const Scalar reference_check = std::max(
      {profile_difference(a, b, [&](Scalar r) { return exact.u(r); },
                          [&](Scalar r) {
                            return (1.0 + nu) / e * ((1.0 - 2.0 * nu) * big_a * r + big_b / r);
                          }),
       profile_difference(a, b, [&](Scalar r) { return exact.sigma_theta(r); },
                          [&](Scalar r) { return big_a + big_b / (r * r); }),
       profile_difference(a, b, [&](Scalar r) { return exact.sigma_r(r); },
                          [&](Scalar r) { return big_a - big_b / (r * r); })});
  const Scalar self_check = exact.self_check();

  const auto build = [&](Mesh mesh) {
    const int dim = mesh.dim();
    FemModel model(std::move(mesh), IsotropicMaterial(e, nu, 7850.0, "steel"), 1.0,
                   dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional,
                   IntegrationOptions());
    add_quarter_supports(model);
    LoadCaseSpec load;
    load.name = "internal_pressure";
    PressureLoadSpec pressure;
    pressure.region = cylinder_surface("bore", a, true);
    pressure.pressure = p;
    load.pressures.push_back(pressure);
    model.load_case_specs().push_back(load);
    return model;
  };
  json::Value block = json::Value::make_object();
  const std::map<std::string, Ladder> ladders = run_sections(
      path_join(out_dir, "lame_cylinder.csv"),
      {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet10}, a, b,
      build, [&](int) -> const AxisymmetricSolution& { return exact; }, nullptr, block);
  block.set("reference_vs_closed_form", json::Value::make_number(reference_check));
  block.set("reference_self_check", json::Value::make_number(self_check));
  block.set("kind", json::Value::make_string("verification (exact elasticity solution)"));
  block.set("note",
            json::Value::make_string(
                "Quarter section of a cylinder a = 0.1 m, b = 0.2 m under an internal "
                "pressure of 100 MPa in plane strain (E = 200 GPa, nu = 0.3), with symmetry "
                "supports; the 3-D sections are one cell deep with u_z = 0. Errors: the "
                "RMS nodal displacement error relative to the RMS exact displacement "
                "(checked), the largest nodal error over the largest exact displacement, "
                "and the largest error of sigma_r, sigma_theta, sigma_z at element "
                "centroids (the element's quadrature-point average) over the largest exact "
                "stress. Reference: Lame's solution."));
  summary.set("lame_cylinder", block);

  std::ostringstream note;
  const Scalar shortfall = order_shortfall(
      ladders, [](const Ladder& l) -> const std::vector<Scalar>& { return l.displacement; },
      "u error", note);
  note << "reference vs closed form " << app::format(reference_check, 2);
  StudyOutcome outcome;
  outcome.name = "thick cylinder under internal pressure vs Lame (plane strain)";
  outcome.kind = "verification";
  outcome.metric = "largest shortfall of the displacement convergence order";
  outcome.value = shortfall;
  outcome.tolerance = kOrderMargin;
  outcome.passed = shortfall <= kOrderMargin && reference_check < 1.0e-12 && self_check < 1.0e-6;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Rotating disk and cylinder
// ---------------------------------------------------------------------------

StudyOutcome study_rotating_disk(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.3;
  const Scalar rho = 7800.0;
  const Scalar omega = 600.0;
  const Scalar e = 200.0e9;
  const Scalar nu = 0.3;
  AxisymmetricSolution::Data d;
  d.youngs = e;
  d.poisson = nu;
  d.a = a;
  d.b = b;
  d.rho_omega2 = rho * omega * omega;
  d.plane_strain = false;
  const AxisymmetricSolution disk(d);
  d.plane_strain = true;
  const AxisymmetricSolution cylinder(d);

  // Plane-stress closed form (Timoshenko and Goodier, section 32):
  // u = rho w^2 r / (8 E) [(3 + nu)(1 - nu)(a^2 + b^2) + (3 + nu)(1 + nu) a^2 b^2 / r^2
  //                        - (1 - nu^2) r^2];
  // the plane-strain one is the same with E -> E / (1 - nu^2), nu -> nu / (1 - nu).
  const auto disk_u = [&](Scalar r, Scalar ee, Scalar n) {
    return rho * omega * omega * r / (8.0 * ee) *
           ((3.0 + n) * (1.0 - n) * (a * a + b * b) + (3.0 + n) * (1.0 + n) * a * a * b * b / (r * r) -
            (1.0 - n * n) * r * r);
  };
  const Scalar reference_check = std::max(
      profile_difference(a, b, [&](Scalar r) { return disk.u(r); },
                         [&](Scalar r) { return disk_u(r, e, nu); }),
      profile_difference(a, b, [&](Scalar r) { return cylinder.u(r); },
                         [&](Scalar r) {
                           return disk_u(r, e / (1.0 - nu * nu), nu / (1.0 - nu));
                         }));
  const Scalar self_check = std::max(disk.self_check(), cylinder.self_check());

  const Scalar disk_thickness = 0.01;
  const auto build = [&](Mesh mesh) {
    const int dim = mesh.dim();
    FemModel model(std::move(mesh), IsotropicMaterial(e, nu, rho, "steel"),
                   dim == 2 ? disk_thickness : 1.0,
                   dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                   IntegrationOptions());
    add_quarter_supports(model);
    LoadCaseSpec load;
    load.name = "rotation";
    load.centrifugal.enabled = true;
    load.centrifugal.angular_velocity = omega;
    load.centrifugal.axis = Vector3::UnitZ();
    load.centrifugal.point = Vector3::Zero();
    model.load_case_specs().push_back(load);
    return model;
  };
  json::Value block = json::Value::make_object();
  const std::map<std::string, Ladder> ladders = run_sections(
      path_join(out_dir, "rotating_disk.csv"),
      {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet10}, a, b,
      build,
      [&](int dim) -> const AxisymmetricSolution& { return dim == 2 ? disk : cylinder; },
      nullptr, block);
  block.set("reference_vs_closed_form", json::Value::make_number(reference_check));
  block.set("reference_self_check", json::Value::make_number(self_check));
  block.set("kind", json::Value::make_string("verification (exact elasticity solution)"));
  block.set("note",
            json::Value::make_string(
                "Quarter section, a = 0.1 m, b = 0.3 m, steel (rho = 7800 kg/m^3) spinning "
                "at 600 rad/s about z with free surfaces: a 10 mm disk in plane stress (Q4, "
                "Tri3) and a long cylinder in plane strain (Hex8, Tet10, one cell deep with "
                "u_z = 0). The centrifugal body force rho omega^2 r is assembled from the "
                "consistent mass. Reference: Timoshenko and Goodier, section 32."));
  summary.set("rotating_disk", block);

  std::ostringstream note;
  const Scalar shortfall = order_shortfall(
      ladders, [](const Ladder& l) -> const std::vector<Scalar>& { return l.displacement; },
      "u error", note);
  note << "reference vs closed form " << app::format(reference_check, 2);
  StudyOutcome outcome;
  outcome.name = "rotating disk (plane stress) and cylinder (plane strain) vs exact";
  outcome.kind = "verification";
  outcome.metric = "largest shortfall of the displacement convergence order";
  outcome.value = shortfall;
  outcome.tolerance = kOrderMargin;
  outcome.passed = shortfall <= kOrderMargin && reference_check < 1.0e-12 && self_check < 1.0e-6;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Conduction and thermal stress in a thick cylinder
// ---------------------------------------------------------------------------

StudyOutcome study_thermal_cylinder(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.2;
  const Scalar k = 45.0;          // W/(m K)
  const Scalar q = 1.0e6;         // W/m^3
  const Scalar t_bore = 400.0;    // K
  const Scalar h_film = 500.0;    // W/(m^2 K)
  const Scalar t_ambient = 300.0; // K
  const Scalar t_ref = 300.0;     // K
  const Scalar alpha = 1.2e-5;
  const Scalar e = 200.0e9;
  const Scalar nu = 0.3;

  // T(r) = -q r^2 / (4 k) + A ln(r / a) + B with T(a) = t_bore and
  // -k T'(b) = h (T(b) - T_inf).
  const Scalar big_b = t_bore + q * a * a / (4.0 * k);
  const Scalar big_a = (0.5 * q * b + h_film * (q * b * b / (4.0 * k) - big_b + t_ambient)) /
                       (h_film * std::log(b / a) + k / b);
  const std::function<Scalar(Scalar)> temperature = [=](Scalar r) {
    return -q * r * r / (4.0 * k) + big_a * std::log(r / a) + big_b;
  };
  const auto temperature_slope = [=](Scalar r) { return -q * r / (2.0 * k) + big_a / r; };
  // Check of the temperature: the convection condition and the heat balance
  // of the quarter per unit depth.
  const Scalar convection_mismatch =
      std::abs(-k * temperature_slope(b) - h_film * (temperature(b) - t_ambient)) /
      (h_film * std::abs(temperature(b) - t_ambient));
  const Scalar heat_in_bore = -k * temperature_slope(a) * 0.5 * kPi * a;  // W/m
  const Scalar heat_generated = q * 0.25 * kPi * (b * b - a * a);         // W/m
  const Scalar heat_out = h_film * (temperature(b) - t_ambient) * 0.5 * kPi * b;
  const Scalar balance_mismatch =
      std::abs(heat_in_bore + heat_generated - heat_out) / std::abs(heat_out);

  AxisymmetricSolution::Data d;
  d.plane_strain = true;
  d.youngs = e;
  d.poisson = nu;
  d.alpha = alpha;
  d.a = a;
  d.b = b;
  d.delta_t.a = a;
  d.delta_t.t0 = big_b - t_ref;
  d.delta_t.t_log = big_a;
  d.delta_t.t_r2 = -q / (4.0 * k);
  const AxisymmetricSolution exact(d);
  const Scalar self_check = exact.self_check();
  const Scalar temperature_check = profile_difference(
      a, b, [&](Scalar r) { return d.delta_t(r) + t_ref; }, temperature);

  const auto build = [&](Mesh mesh) {
    const int dim = mesh.dim();
    IsotropicMaterial steel(e, nu, 7850.0, "steel");
    steel.set_thermal(alpha, t_ref, k);
    FemModel model(std::move(mesh), steel, 1.0,
                   dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional,
                   IntegrationOptions());
    add_quarter_supports(model);
    LoadCaseSpec load;
    load.name = "heated";
    load.temperature.source = TemperatureSpec::Source::Conduction;
    RegionValue bore;
    bore.region = cylinder_surface("bore", a, true);
    bore.value = t_bore;
    load.temperature.conduction.prescribed.push_back(bore);
    ConvectionSpec film;
    film.region = cylinder_surface("outer", b, false);
    film.film_coefficient = h_film;
    film.ambient = t_ambient;
    load.temperature.conduction.convection.push_back(film);
    RegionValue source;
    source.region = whole_model_region();
    source.whole_model = true;
    source.value = q;
    load.temperature.conduction.sources.push_back(source);
    model.load_case_specs().push_back(load);
    LinearSolverOptions conduction;
    conduction.type = LinearSolverType::SimplicialLdlt;
    model.set_conduction_solver(conduction);
    return model;
  };
  // The heat entering through the bore, from the reactions of the prescribed
  // temperatures, against the exact flux k T'(a) over the quarter bore.
  const auto bore_heat_error = [&](const FemModel& model, const SolvedModel& solved) {
    const Scalar depth = model.dim() == 2 ? model.thickness() : model.mesh().bounding_box().extent().z();
    const Scalar entering = -solved.conduction.prescribed_heat;
    return std::abs(entering - heat_in_bore * depth) / std::abs(heat_in_bore * depth);
  };
  json::Value block = json::Value::make_object();
  const std::map<std::string, Ladder> ladders = run_sections(
      path_join(out_dir, "thermal_cylinder.csv"),
      {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet10}, a, b,
      build, [&](int) -> const AxisymmetricSolution& { return exact; }, &temperature, block,
      bore_heat_error);
  block.set("reference_self_check", json::Value::make_number(self_check));
  block.set("temperature_reference_check",
            json::Value::make_number(std::max({temperature_check, convection_mismatch,
                                               balance_mismatch})));
  block.set("exact_bore_heat_W_per_m", json::Value::make_number(heat_in_bore));
  block.set("exact_temperature_K",
            json::Value::make_string(app::format(temperature(a), 6) + " at the bore, " +
                                     app::format(temperature(b), 6) + " at the outer surface"));
  block.set("kind", json::Value::make_string("verification (exact conduction and elasticity)"));
  block.set("note",
            json::Value::make_string(
                "Quarter section a = 0.1 m, b = 0.2 m, k = 45 W/(m K), uniform generation "
                "1 MW/m^3, bore held at 400 K, convection h = 500 W/(m^2 K) to 300 K from "
                "the outer surface; steel with alpha = 1.2e-5 /K, stress free at 300 K, in "
                "plane strain with free surfaces. The temperature is the conduction "
                "solution on the same mesh; the extra error is that of the heat entering "
                "through the bore, from the reactions of the prescribed temperatures. "
                "Reference: T = -q r^2/(4k) + A ln(r/a) + B, and the thermoelastic solution "
                "of a long cylinder (Timoshenko and Goodier, section 150)."));
  summary.set("thermal_cylinder", block);

  std::ostringstream note;
  const Scalar shortfall_t = order_shortfall(
      ladders, [](const Ladder& l) -> const std::vector<Scalar>& { return l.temperature; },
      "T error", note);
  const Scalar shortfall_u = order_shortfall(
      ladders, [](const Ladder& l) -> const std::vector<Scalar>& { return l.displacement; },
      "u error", note);
  StudyOutcome outcome;
  outcome.name = "conduction and thermal stress in a thick cylinder vs exact";
  outcome.kind = "verification";
  outcome.metric = "largest shortfall of the temperature and displacement orders";
  outcome.value = std::max(shortfall_t, shortfall_u);
  outcome.tolerance = kOrderMargin;
  outcome.passed = outcome.value <= kOrderMargin && self_check < 1.0e-6 &&
                   std::max({temperature_check, convection_mismatch, balance_mismatch}) < 1.0e-12;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Bimetallic strip
// ---------------------------------------------------------------------------

StudyOutcome study_bimetal_strip(const std::string& out_dir, json::Value& summary) {
  // Layer 1 (bottom, y in [-t1, 0]) expands little, layer 2 (top, y in
  // [0, t2]) a lot.
  const Scalar length = 0.02;
  const Scalar t1 = 0.4e-3;
  const Scalar t2 = 0.6e-3;
  const Scalar width = 5.0e-3;
  const Scalar e1 = 140.0e9;
  const Scalar nu1 = 0.26;
  const Scalar alpha1 = 1.5e-6;
  const Scalar e2 = 100.0e9;
  const Scalar nu2 = 0.34;
  const Scalar alpha2 = 19.0e-6;
  const Scalar t_ref = 293.15;
  const Scalar delta_t = 50.0;

  // Exact curvature where the end effects have decayed: the axial strain
  // eps0 + kappa y makes sigma_xx = E_i (eps0 + kappa y - alpha_i dT) in each
  // layer with no resultant force or moment (sigma_yy = sigma_xy = 0 is then
  // an exact plane-stress state).
  Eigen::Matrix2d m;
  Eigen::Vector2d rhs;
  {
    const Scalar area[2] = {t1, t2};
    const Scalar first[2] = {-0.5 * t1 * t1, 0.5 * t2 * t2};
    const Scalar second[2] = {t1 * t1 * t1 / 3.0, t2 * t2 * t2 / 3.0};
    const Scalar young[2] = {e1, e2};
    const Scalar expansion[2] = {alpha1, alpha2};
    m.setZero();
    rhs.setZero();
    for (int i = 0; i < 2; ++i) {
      m(0, 0) += young[i] * area[i];
      m(0, 1) += young[i] * first[i];
      m(1, 0) += young[i] * first[i];
      m(1, 1) += young[i] * second[i];
      rhs(0) += young[i] * area[i] * expansion[i] * delta_t;
      rhs(1) += young[i] * first[i] * expansion[i] * delta_t;
    }
  }
  const Eigen::Vector2d solution = m.fullPivLu().solve(rhs);
  const Scalar kappa = solution(1);
  // Timoshenko's formula, 1/rho = 6 (a2 - a1) dT (1 + m)^2 /
  // (h [3 (1 + m)^2 + (1 + m n)(m^2 + 1/(m n))]), m = t1/t2, n = E1/E2.
  const Scalar mm = t1 / t2;
  const Scalar nn = e1 / e2;
  const Scalar h = t1 + t2;
  const Scalar kappa_timoshenko =
      6.0 * (alpha2 - alpha1) * delta_t * (1.0 + mm) * (1.0 + mm) /
      (h * (3.0 * (1.0 + mm) * (1.0 + mm) + (1.0 + mm * nn) * (mm * mm + 1.0 / (mm * nn))));
  const Scalar reference_check = std::abs(kappa_timoshenko - kappa) / std::abs(kappa);

  CsvWriter csv(path_join(out_dir, "bimetal_strip.csv"),
                {"nx", "ny", "h[m]", "num_dofs", "curvature[1/m]", "exact[1/m]", "error[-]",
                 "order[-]", "max_reaction[N]"});
  Ladder ladder;
  json::Value records = json::Value::make_array();
  for (const Index ny : {5, 10, 20, 40}) {
    const Index nx = 20 * ny;
    StructuredMeshSpec ms;
    ms.nx = nx;
    ms.ny = ny;
    ms.lx = length;
    ms.ly = h;
    ms.y0 = -t1;
    IsotropicMaterial low(e1, nu1, 8100.0, "low_expansion");
    low.set_thermal(alpha1, t_ref, 0.0);
    IsotropicMaterial high(e2, nu2, 8500.0, "high_expansion");
    high.set_thermal(alpha2, t_ref, 0.0);
    FemModel model(make_structured_quad_mesh(ms), low, width, StressState::PlaneStress,
                   IntegrationOptions());
    std::vector<Index> top;
    for (Index el = 0; el < model.mesh().num_elements(); ++el) {
      if (model.mesh().element_centroid(el).y() > 0.0) top.push_back(el);
    }
    model.assign_material(high, top);
    // The exact field has u_x = 0 on x = L/2: hold the interface node there in
    // x and y and the top node in x, which removes the rigid motions and
    // leaves the self-equilibrated thermal load no reaction to find.
    DisplacementConstraint centre;
    centre.region.name = "centre";
    Selector c;
    c.kind = SelectorKind::NearestNode;
    c.point = Vector3(0.5 * length, 0.0, 0.0);
    centre.region.members.push_back(c);
    centre.fix_x = centre.fix_y = true;
    model.constraints().push_back(centre);
    DisplacementConstraint upper;
    upper.region.name = "upper";
    Selector cu;
    cu.kind = SelectorKind::NearestNode;
    cu.point = Vector3(0.5 * length, t2, 0.0);
    upper.region.members.push_back(cu);
    upper.fix_x = true;
    model.constraints().push_back(upper);
    LoadCaseSpec load;
    load.name = "heated";
    load.temperature.source = TemperatureSpec::Source::Uniform;
    load.temperature.uniform = t_ref + delta_t;
    model.load_case_specs().push_back(load);
    const SolvedModel solved = solve_model(model);

    // Least-squares parabola u_y = c0 + c1 s + c2 s^2, s = x - L/2, through
    // the bottom-edge nodes of the middle half: kappa = -u_y'' = -2 c2.
    const Mesh& mesh = model.mesh();
    std::vector<Scalar> s_values;
    std::vector<Scalar> v_values;
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      const Vector3 x = mesh.node(n);
      if (std::abs(x.y() + t1) > 1.0e-3 * h / static_cast<Scalar>(ny)) continue;
      const Scalar s = x.x() - 0.5 * length;
      if (std::abs(s) > 0.25 * length * (1.0 + 1.0e-9)) continue;
      s_values.push_back(s);
      v_values.push_back(solved.solution.displacement(2 * n + 1));
    }
    Matrix design(static_cast<Eigen::Index>(s_values.size()), 3);
    Vector values(static_cast<Eigen::Index>(s_values.size()));
    for (std::size_t i = 0; i < s_values.size(); ++i) {
      const auto row = static_cast<Eigen::Index>(i);
      const Scalar s = s_values[i] / (0.25 * length);  // scaled for conditioning
      design(row, 0) = 1.0;
      design(row, 1) = s;
      design(row, 2) = s * s;
      values(row) = v_values[i];
    }
    const Vector coeff = design.colPivHouseholderQr().solve(values);
    const Scalar curvature = -2.0 * coeff(2) / (0.0625 * length * length);
    const Scalar error = std::abs(curvature - kappa) / std::abs(kappa);
    const Scalar max_reaction = solved.solution.reactions.cwiseAbs().maxCoeff();
    ladder.h.push_back(h / static_cast<Scalar>(ny));
    ladder.displacement.push_back(error);
    csv.raw_row({fmt(static_cast<Scalar>(nx)), fmt(static_cast<Scalar>(ny)),
                 fmt(ladder.h.back(), 8), fmt(static_cast<Scalar>(solved.num_dofs), 9),
                 fmt(curvature, 12), fmt(kappa, 12), fmt(error, 8),
                 fmt(ladder.running_order(ladder.displacement), 4), fmt(max_reaction, 4)});
    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(static_cast<Scalar>(nx)));
    rec.set("ny", json::Value::make_number(static_cast<Scalar>(ny)));
    rec.set("num_dofs", json::Value::make_number(static_cast<Scalar>(solved.num_dofs)));
    rec.set("curvature_per_m", json::Value::make_number(curvature));
    rec.set("relative_error", json::Value::make_number(error));
    rec.set("max_reaction_N", json::Value::make_number(max_reaction));
    records.push_back(rec);
  }
  csv.close();
  const Scalar order = ladder.order(ladder.displacement);
  json::Value block = json::Value::make_object();
  block.set("meshes", records);
  block.set("exact_curvature_per_m", json::Value::make_number(kappa));
  block.set("timoshenko_formula_per_m", json::Value::make_number(kappa_timoshenko));
  block.set("reference_vs_closed_form", json::Value::make_number(reference_check));
  block.set("order", json::Value::make_number(order));
  block.set("kind", json::Value::make_string("verification (exact plane-stress solution away "
                                              "from the ends)"));
  block.set("note",
            json::Value::make_string(
                "A 20 mm strip of two layers, 0.4 mm (E = 140 GPa, nu = 0.26, alpha = "
                "1.5e-6 /K) under 0.6 mm (E = 100 GPa, nu = 0.34, alpha = 1.9e-5 /K), "
                "heated uniformly by 50 K in plane stress with square Q4 cells and supports "
                "that react nothing. The curvature is fitted to the bottom-edge deflection "
                "of the middle half, five strip depths from the free ends, where the "
                "uniform-curvature state of Timoshenko's bimetal theory is an exact "
                "plane-stress solution."));
  summary.set("bimetal_strip", block);

  StudyOutcome outcome;
  outcome.name = "bimetallic strip under uniform heating vs Timoshenko";
  outcome.kind = "verification";
  outcome.metric = "finest-mesh relative curvature error";
  outcome.value = ladder.displacement.back();
  outcome.tolerance = 1.0e-3;
  outcome.passed = outcome.value <= outcome.tolerance && order >= 2.0 - kOrderMargin &&
                   reference_check < 1.0e-12;
  std::ostringstream note;
  note << "curvature " << app::format(kappa, 6) << " 1/m, order " << app::format(order, 3)
       << ", Timoshenko formula vs the force and moment balance "
       << app::format(reference_check, 2);
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// A bar hanging under its own weight
// ---------------------------------------------------------------------------

StudyOutcome study_self_weight(const std::string& out_dir, json::Value& summary) {
  // The bar occupies |x|, |y| <= w/2, 0 <= z <= L (a plate |x| <= w/2,
  // 0 <= y <= L in 2-D), its bottom free and its top pulled up by the traction
  // rho g L that carries its weight. The exact field, sigma_zz = rho g z and
  // every other stress zero, has
  //   u_x = -nu rho g x z / E,  u_y = -nu rho g y z / E,
  //   u_z = rho g (z^2 + nu (x^2 + y^2) - L^2) / (2 E),
  // a quadratic that the Tet10 space contains and the linear elements
  // approximate at O(h^2). Three-two-one supports at points where the exact
  // field is zero remove the rigid motions and react nothing.
  const Scalar length = 0.5;
  const Scalar width = 0.1;
  const Scalar rho = 7850.0;
  const Scalar g = 9.81;
  const Scalar e = 200.0e9;
  const Scalar nu = 0.3;
  const Scalar c = rho * g / e;

  CsvWriter csv(path_join(out_dir, "self_weight.csv"),
                {"element", "nx", "ny", "nz", "h[m]", "num_dofs", "u_rms_error[-]",
                 "u_rms_order[-]", "u_max_error[-]", "u_max_order[-]",
                 "u_max_interior_error[-]", "u_max_interior_order[-]",
                 "max_reaction_over_weight[-]", "solver"});
  std::map<std::string, Ladder> ladders;
  std::map<std::string, Scalar> worst_reaction;
  json::Value block = json::Value::make_object();
  const auto run = [&](ElementType type, Index nx) {
    const bool plane = type == ElementType::Quad4 || type == ElementType::Tri3;
    StructuredMeshSpec ms;
    ms.nx = nx;
    ms.ny = plane ? 5 * nx : nx;
    ms.nz = 5 * nx;
    ms.lx = width;
    ms.ly = plane ? length : width;
    ms.lz = length;
    ms.x0 = -0.5 * width;
    ms.y0 = plane ? 0.0 : -0.5 * width;
    Mesh mesh = type == ElementType::Quad4  ? make_structured_quad_mesh(ms)
                : type == ElementType::Tri3 ? make_structured_tri_mesh(ms)
                : type == ElementType::Hex8 ? make_structured_hex_mesh(ms)
                : type == ElementType::Tet4 ? make_structured_tet_mesh(ms)
                                            : make_structured_tet10_mesh(ms);
    const int dim = mesh.dim();
    const int up = dim - 1;  // the axial direction: y in 2-D, z in 3-D
    FemModel model(std::move(mesh), IsotropicMaterial(e, nu, rho, "steel"),
                   plane ? width : 1.0,
                   plane ? StressState::PlaneStress : StressState::ThreeDimensional,
                   IntegrationOptions());
    const auto point = [&](Scalar x, Scalar y, Scalar z) {
      Selector s;
      s.kind = SelectorKind::NearestNode;
      s.point = plane ? Vector3(x, z, 0.0) : Vector3(x, y, z);
      return s;
    };
    DisplacementConstraint top_centre;
    top_centre.region.name = "top_centre";
    top_centre.region.members.push_back(point(0.0, 0.0, length));
    top_centre.fix_x = top_centre.fix_y = true;
    top_centre.fix_z = !plane;
    model.constraints().push_back(top_centre);
    DisplacementConstraint bottom_centre;
    bottom_centre.region.name = "bottom_centre";
    bottom_centre.region.members.push_back(point(0.0, 0.0, 0.0));
    bottom_centre.fix_x = true;
    bottom_centre.fix_y = !plane;
    model.constraints().push_back(bottom_centre);
    if (!plane) {
      DisplacementConstraint top_edge;
      top_edge.region.name = "top_edge";
      top_edge.region.members.push_back(point(0.5 * width, 0.0, length));
      top_edge.fix_y = true;
      model.constraints().push_back(top_edge);
    }
    LoadCaseSpec load;
    load.name = "self_weight";
    load.gravity = Vector3::Zero();
    load.gravity(up) = -g;
    TractionLoadSpec pull;
    pull.region.name = "top";
    Selector top_face = box_selector();
    if (plane) {
      top_face.ymin = length;
    } else {
      top_face.zmin = length;
    }
    pull.region.members.push_back(top_face);
    pull.traction = Vector3::Zero();
    pull.traction(up) = rho * g * length;
    load.tractions.push_back(pull);
    model.load_case_specs().push_back(load);
    const SolvedModel solved = solve_model(model);

    const Mesh& m = model.mesh();
    Scalar err = 0.0;
    Scalar scale = 0.0;
    Scalar err_sq = 0.0;
    Scalar exact_sq = 0.0;
    Scalar interior = 0.0;  // largest error with 0.2 L < axial coordinate < 0.8 L
    for (Index n = 0; n < m.num_nodes(); ++n) {
      const Vector3 x = m.node(n);
      Vector3 exact = Vector3::Zero();
      if (plane) {
        exact.x() = -nu * c * x.x() * x.y();
        exact.y() = 0.5 * c * (x.y() * x.y() + nu * x.x() * x.x() - length * length);
      } else {
        exact.x() = -nu * c * x.x() * x.z();
        exact.y() = -nu * c * x.y() * x.z();
        exact.z() = 0.5 * c * (x.z() * x.z() + nu * (x.x() * x.x() + x.y() * x.y()) -
                               length * length);
      }
      Vector3 fe = Vector3::Zero();
      for (int k = 0; k < dim; ++k) fe(k) = solved.solution.displacement(dim * n + k);
      err = std::max(err, (fe - exact).norm());
      scale = std::max(scale, exact.norm());
      err_sq += (fe - exact).squaredNorm();
      exact_sq += exact.squaredNorm();
      const Scalar axial = x(up);
      if (axial > 0.2 * length && axial < 0.8 * length) {
        interior = std::max(interior, (fe - exact).norm());
      }
    }
    const Scalar weight = rho * g * length * width * width;
    const Scalar reaction = solved.solution.reactions.cwiseAbs().maxCoeff() / weight;
    const std::string name = to_string(type);
    Ladder& ladder = ladders[name];
    ladder.h.push_back(width / static_cast<Scalar>(nx));
    ladder.displacement.push_back(std::sqrt(err_sq / exact_sq));
    ladder.displacement_max.push_back(err / scale);
    ladder.extra.push_back(interior / scale);
    worst_reaction[name] = std::max(worst_reaction[name], reaction);
    csv.raw_row({name, fmt(static_cast<Scalar>(ms.nx)), fmt(static_cast<Scalar>(ms.ny)),
                 fmt(static_cast<Scalar>(plane ? 0 : ms.nz)), fmt(ladder.h.back(), 8),
                 fmt(static_cast<Scalar>(solved.num_dofs), 9),
                 fmt(ladder.displacement.back(), 8),
                 fmt(ladder.running_order(ladder.displacement), 4),
                 fmt(ladder.displacement_max.back(), 8),
                 fmt(ladder.running_order(ladder.displacement_max), 4),
                 fmt(ladder.extra.back(), 8), fmt(ladder.running_order(ladder.extra), 4),
                 fmt(reaction, 4), solved.solution.solver_name});
  };
  for (const Index nx : {2, 4, 8, 16}) run(ElementType::Quad4, nx);
  for (const Index nx : {2, 4, 8, 16}) run(ElementType::Tri3, nx);
  for (const Index nx : {2, 4, 8}) run(ElementType::Hex8, nx);
  for (const Index nx : {2, 4, 8, 16}) run(ElementType::Tet4, nx);
  for (const Index nx : {2, 4}) run(ElementType::Tet10, nx);
  csv.close();

  std::ostringstream note;
  Scalar shortfall = -1.0e300;
  for (const auto& kv : ladders) {
    json::Value entry = json::Value::make_object();
    json::Value rms = json::Value::make_array();
    json::Value largest = json::Value::make_array();
    for (Scalar v : kv.second.displacement) rms.push_back(json::Value::make_number(v));
    for (Scalar v : kv.second.displacement_max) largest.push_back(json::Value::make_number(v));
    entry.set("displacement_rms_errors", rms);
    entry.set("displacement_max_errors", largest);
    entry.set("max_reaction_over_weight", json::Value::make_number(worst_reaction[kv.first]));
    if (kv.first != "Tet10") {
      const Scalar p = kv.second.order(kv.second.displacement);
      entry.set("displacement_rms_order", json::Value::make_number(p));
      entry.set("displacement_max_order",
                json::Value::make_number(kv.second.order(kv.second.displacement_max)));
      entry.set("interior_max_order",
                json::Value::make_number(kv.second.order(kv.second.extra)));
      shortfall = std::max(shortfall, 2.0 - p);
      note << kv.first << " " << app::format(kv.second.displacement.back(), 3) << " (order "
           << app::format(p, 3) << "); ";
    }
    block.set(kv.first, entry);
  }
  const Scalar tet10_error = *std::max_element(ladders["Tet10"].displacement_max.begin(),
                                               ladders["Tet10"].displacement_max.end());
  Scalar reaction = 0.0;
  for (const auto& kv : worst_reaction) reaction = std::max(reaction, kv.second);
  note << "Tet10 " << app::format(tet10_error, 3) << " (exact space); largest reaction "
       << app::format(reaction, 2) << " of the weight";
  block.set("kind", json::Value::make_string("verification (exact elasticity solution)"));
  block.set("note",
            json::Value::make_string(
                "A 0.5 m steel bar of 0.1 m square section (a 0.1 m wide plate in plane stress "
                "for Q4 and Tri3) hanging under its own weight (g = 9.81 m/s^2), carried by "
                "the traction rho g L on its top face, with three-two-one supports at "
                "points where the exact field vanishes. Errors: the RMS nodal "
                "displacement error relative to the RMS exact displacement, whose order "
                "the study checks, and the largest nodal error over the largest exact "
                "displacement. The largest error sits at the four corners where the "
                "ends meet the free sides and approaches its O(h^2) rate slowly, as the "
                "maximum-norm estimate for bilinear elements (O(h^2 |log h|)) allows; "
                "the largest error away from the ends (0.2 L to 0.8 L, the interior "
                "columns) is reported beside it."));
  summary.set("self_weight", block);

  StudyOutcome outcome;
  outcome.name = "bar hanging under its own weight vs exact";
  outcome.kind = "verification";
  outcome.metric = "Tet10 displacement error (its space holds the exact field)";
  outcome.value = tet10_error;
  outcome.tolerance = 1.0e-10;
  outcome.passed = tet10_error <= outcome.tolerance && shortfall <= kOrderMargin &&
                   reaction < 1.0e-8;
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
