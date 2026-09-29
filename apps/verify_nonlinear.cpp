/// \file verify_nonlinear.cpp
/// \brief Verification of the geometrically non-linear analysis against exact
///        solutions of finite elasticity, and of its path-following through a
///        limit point.
///
/// Studies:
///   * `elastica`              a slender Tet10 cantilever under a dead force at
///                             its tip, up to k = P L^2 / (E I) = 10, against
///                             Euler's elastica (large rotation, small strain,
///                             Saint Venant-Kirchhoff);
///   * `hyperelastic-cylinder` a thick tube at finite strain against the exact
///                             axisymmetric solution, three ways: the
///                             compressible neo-Hookean material inflated by a
///                             pressure that follows its bore (hoop stretch
///                             1.48) and spinning, its centrifugal load at the
///                             deformed radius; a Saint Venant-Kirchhoff tube
///                             under a conducted temperature (the
///                             multiplicative thermal split). Quarter sections
///                             of Q4 and Tri3 in plane strain and one-layer Hex8
///                             and Tet10 sections held at u_z = 0;
///   * `arch-snap-through`     a clamped shallow arch under a crown load that
///                             snaps through: the arc-length path against
///                             displacement control at identical crown
///                             displacements, the tangent's inertia against the
///                             slope of the path, and load control stopping at the
///                             limit point with a bracket that encloses the limit
///                             load a fine displacement-controlled sweep finds.
///
/// The elastica is a beam theory: the continuum cantilever differs from it by
/// shear flexibility and by the curvature dependence of the Saint
/// Venant-Kirchhoff bending stiffness, both of order (h/L)^2 or (curvature
/// times h)^2. That study therefore separates the two errors: the
/// discretisation error, whose order of convergence three meshes measure by
/// Richardson's method, and the gap between the converged continuum and the
/// beam, which it bounds. The tube's solution is exact for the continuum the
/// model discretises, so there the error must vanish at the element's rate, as
/// in the linear studies of verify_loads.cpp.

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

Selector box() {
  Selector s;
  s.kind = SelectorKind::Box;
  return s;
}

SelectorGroup region_of(const std::string& name, const Selector& s) {
  SelectorGroup g;
  g.name = name;
  g.members.push_back(s);
  return g;
}

SelectorGroup node_region(const std::string& name, std::vector<Index> nodes) {
  Selector s;
  s.kind = SelectorKind::NodeIds;
  s.ids = std::move(nodes);
  return region_of(name, s);
}

NonlinearMonitor monitor(const std::string& name, const SelectorGroup& region, int component,
                         NonlinearMonitor::Quantity quantity =
                             NonlinearMonitor::Quantity::Displacement) {
  NonlinearMonitor m;
  m.name = name;
  m.region = region;
  m.component = component;
  m.quantity = quantity;
  return m;
}

/// The value of monitor `name` at a converged step.
Scalar monitor_value(const NonlinearResult& r, const NonlinearStep& s, const std::string& name) {
  for (std::size_t i = 0; i < r.monitor_names.size(); ++i) {
    if (r.monitor_names[i] == name) return s.monitors[i];
  }
  throw ConfigError("no monitor named '" + name + "'");
}

/// Observed order of three solutions on meshes refined by 2, from the
/// differences between them alone (Richardson): p = log2((q1 - q2) / (q2 - q3)).
/// NaN when the differences do not shrink monotonically.
Scalar three_mesh_order(Scalar q1, Scalar q2, Scalar q3) {
  const Scalar ratio = (q1 - q2) / (q2 - q3);
  return ratio > 0.0 ? std::log2(ratio) : std::nan("");
}

/// The limit the three solutions extrapolate to: q* = (q2^2 - q1 q3) /
/// (2 q2 - q1 - q3).
Scalar three_mesh_limit(Scalar q1, Scalar q2, Scalar q3) {
  return (q2 * q2 - q1 * q3) / (2.0 * q2 - q1 - q3);
}

// ---------------------------------------------------------------------------
// Euler's elastica
// ---------------------------------------------------------------------------

/// Tip of the elastica, in units of the length.
struct ElasticaTip {
  Scalar deflection = 0.0;  ///< v / L, in the direction of the force
  Scalar shortening = 0.0;  ///< (L - x_tip) / L
  Scalar rotation = 0.0;    ///< theta(1) [rad]
};

/// (theta, theta', x, v) at s = 1 of the elastica of a unit cantilever under
/// the load parameter k, from theta'(0) = w: theta'' = -k cos(theta),
/// x' = cos(theta), v' = sin(theta), by the classical Runge-Kutta method with
/// n steps.
std::array<Scalar, 4> elastica_end(Scalar k, Scalar w, int n) {
  const auto f = [k](const std::array<Scalar, 4>& y) {
    return std::array<Scalar, 4>{y[1], -k * std::cos(y[0]), std::cos(y[0]), std::sin(y[0])};
  };
  std::array<Scalar, 4> y{0.0, w, 0.0, 0.0};
  const Scalar h = 1.0 / static_cast<Scalar>(n);
  for (int i = 0; i < n; ++i) {
    const std::array<Scalar, 4> k1 = f(y);
    std::array<Scalar, 4> t{};
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k1[j];
    const std::array<Scalar, 4> k2 = f(t);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k2[j];
    const std::array<Scalar, 4> k3 = f(t);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + h * k3[j];
    const std::array<Scalar, 4> k4 = f(t);
    for (int j = 0; j < 4; ++j) y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
  }
  return y;
}

/// Euler's elastica of a cantilever under a dead force at its free end,
/// perpendicular to the undeformed axis, in k = P L^2 / (E I). The rotation
/// of the centre line satisfies theta'' = -k cos(theta), with theta(0) = 0
/// at the clamp and theta'(1) = 0 at the free end (no moment there). The
/// curvature at the clamp, theta'(0) = w, lies in (0, k] - it is k x_tip -
/// and is found by bisection on theta'(1), which is negative at w = 0 and
/// k (1 - x_tip) > 0 at w = k (Bisshopp and Drucker, 1945, give the solution
/// in elliptic integrals).
ElasticaTip elastica(Scalar k, int n) {
  Scalar lo = 0.0;
  Scalar hi = k;
  if (!(elastica_end(k, lo, n)[1] < 0.0) || !(elastica_end(k, hi, n)[1] > 0.0)) {
    throw SolverError("the elastica shooting interval does not bracket the solution");
  }
  for (int it = 0; it < 200 && hi - lo > 1.0e-16 * k; ++it) {
    const Scalar mid = 0.5 * (lo + hi);
    (elastica_end(k, mid, n)[1] < 0.0 ? lo : hi) = mid;
  }
  const std::array<Scalar, 4> y = elastica_end(k, 0.5 * (lo + hi), n);
  ElasticaTip tip;
  tip.deflection = y[3];
  tip.shortening = 1.0 - y[2];
  tip.rotation = y[0];
  return tip;
}

/// A slender cantilever of square section h x h, Tet10 cells (nx along it,
/// 2 x 2 across), clamped at x = 0 and loaded by a dead shear traction on its
/// tip face with resultant -force along y, solved by load control through the
/// load factors 0.1, 0.2, ..., 1.
struct CantileverRun {
  NonlinearResult result;
  Index num_dofs = 0;
  Scalar seconds = 0.0;
};

CantileverRun solve_cantilever(Index nx, Scalar length, Scalar h, Scalar youngs, Scalar force,
                               int levels) {
  StructuredMeshSpec ms;
  ms.nx = nx;
  ms.ny = 2;
  ms.nz = 2;
  ms.lx = length;
  ms.ly = h;
  ms.lz = h;
  FemModel model(make_structured_tet10_mesh(ms), IsotropicMaterial(youngs, 0.0, 7850.0, "steel"),
                 1.0, StressState::ThreeDimensional, IntegrationOptions());
  DisplacementConstraint clamp;
  Selector root = box();
  root.xmax = 0.0;
  clamp.region = region_of("clamped_root", root);
  clamp.fix_x = true;
  clamp.fix_y = true;
  clamp.fix_z = true;
  model.constraints().push_back(clamp);
  Selector tip = box();
  tip.xmin = length;
  LoadCaseSpec load;
  load.name = "tip_force";
  TractionLoadSpec traction;
  traction.region = region_of("tip_face", tip);
  traction.traction = Vector3(0.0, -force / (h * h), 0.0);
  load.tractions.push_back(traction);
  model.load_case_specs().push_back(load);
  model.finalize();
  Assembler assembler(model);

  NonlinearOptions options;
  options.law = HyperelasticModel::SaintVenantKirchhoff;
  options.steps = levels;
  for (int i = 1; i < levels; ++i) {
    options.load_factors.push_back(static_cast<Scalar>(i) / static_cast<Scalar>(levels));
  }
  Selector top = tip;
  top.ymin = h;
  Selector bottom = tip;
  bottom.ymax = 0.0;
  options.monitors = {monitor("tip_ux", region_of("tip_face", tip), 0),
                      monitor("tip_uy", region_of("tip_face", tip), 1),
                      monitor("top_ux", region_of("tip_top_edge", top), 0),
                      monitor("top_uy", region_of("tip_top_edge", top), 1),
                      monitor("bottom_ux", region_of("tip_bottom_edge", bottom), 0),
                      monitor("bottom_uy", region_of("tip_bottom_edge", bottom), 1)};
  Timer timer;
  NonlinearStaticAnalysis analysis(model, assembler, options);
  CantileverRun run;
  run.result = analysis.solve(0);
  run.seconds = timer.elapsed_seconds();
  run.num_dofs = model.dofs().num_dofs();
  return run;
}

// ---------------------------------------------------------------------------
// The thick tube at finite strain
// ---------------------------------------------------------------------------

/// A material law of a long tube in plane strain (l_z = 1), in its principal
/// strains e_r = l_r - 1 and e_theta = l_theta - 1 at the reference radius R:
/// the principal nominal (first Piola-Kirchhoff) stresses P_r, P_theta, the
/// derivatives of P_r in the strains and - for a law whose parameters vary
/// through the wall, a temperature field - in R at fixed strains. Written in
/// the strains so that it keeps full precision as they vanish.
class TubeMaterial {
 public:
  virtual ~TubeMaterial() = default;
  virtual Scalar radial(Scalar er, Scalar et, Scalar r) const = 0;
  virtual Scalar hoop(Scalar er, Scalar et, Scalar r) const = 0;
  virtual Scalar radial_by_er(Scalar er, Scalar et, Scalar r) const = 0;
  virtual Scalar radial_by_et(Scalar er, Scalar et, Scalar r) const = 0;
  virtual Scalar radial_by_radius(Scalar er, Scalar et, Scalar r) const {
    (void)er;
    (void)et;
    (void)r;
    return 0.0;
  }
};

/// The compressible neo-Hookean law of Hyperelastic.hpp:
/// P_i = [mu (l_i^2 - 1) + lambda ln J] / l_i with J = l_r l_theta.
class NeoHookeanTubeLaw final : public TubeMaterial {
 public:
  NeoHookeanTubeLaw(Scalar mu, Scalar lambda) : mu_(mu), lambda_(lambda) {}
  Scalar radial(Scalar er, Scalar et, Scalar) const override { return nominal(er, log_j(er, et)); }
  Scalar hoop(Scalar er, Scalar et, Scalar) const override { return nominal(et, log_j(er, et)); }
  Scalar radial_by_er(Scalar er, Scalar et, Scalar) const override {
    const Scalar l = 1.0 + er;
    return mu_ * (1.0 + 1.0 / (l * l)) + lambda_ * (1.0 - log_j(er, et)) / (l * l);
  }
  Scalar radial_by_et(Scalar er, Scalar et, Scalar) const override {
    return lambda_ / ((1.0 + er) * (1.0 + et));
  }

 private:
  static Scalar log_j(Scalar er, Scalar et) { return std::log1p(er + et + er * et); }
  Scalar nominal(Scalar ei, Scalar lnj) const {
    return (mu_ * ei * (2.0 + ei) + lambda_ * lnj) / (1.0 + ei);
  }
  Scalar mu_;
  Scalar lambda_;
};

/// The Saint Venant-Kirchhoff law with the multiplicative thermal split of
/// Hyperelastic.hpp under a temperature change dT(R): with the Green strains
/// E_i = e_i + e_i^2 / 2 (E_z = 0), theta = 1 + alpha dT and
/// E_th = (theta^2 - 1) / 2,
/// S_i = [lambda (E_r + E_theta - 3 E_th) + 2 mu (E_i - E_th)] / theta and
/// P_i = (1 + e_i) S_i.
class ThermalSvkTubeLaw final : public TubeMaterial {
 public:
  ThermalSvkTubeLaw(Scalar lambda, Scalar mu, Scalar alpha, std::function<Scalar(Scalar)> dt,
                    std::function<Scalar(Scalar)> dt_slope)
      : lambda_(lambda), mu_(mu), alpha_(alpha), dt_(std::move(dt)),
        dt_slope_(std::move(dt_slope)) {}
  Scalar radial(Scalar er, Scalar et, Scalar r) const override {
    return (1.0 + er) * pk2(er, er, et, r);
  }
  Scalar hoop(Scalar er, Scalar et, Scalar r) const override {
    return (1.0 + et) * pk2(et, er, et, r);
  }
  Scalar radial_by_er(Scalar er, Scalar et, Scalar r) const override {
    const Scalar l = 1.0 + er;
    return pk2(er, er, et, r) + l * l * (lambda_ + 2.0 * mu_) / theta(r);
  }
  Scalar radial_by_et(Scalar er, Scalar et, Scalar r) const override {
    return (1.0 + er) * lambda_ * (1.0 + et) / theta(r);
  }
  Scalar radial_by_radius(Scalar er, Scalar et, Scalar r) const override {
    // d/dR of N / theta at fixed strains, N = lambda (...) + 2 mu (...):
    // -(3 lambda + 2 mu) dE_th/dR / theta - S_r dtheta/dR / theta, with
    // dtheta/dR = alpha dT' and dE_th/dR = theta alpha dT'.
    const Scalar dtheta = alpha_ * dt_slope_(r);
    return (1.0 + er) * (-(3.0 * lambda_ + 2.0 * mu_) * dtheta -
                         pk2(er, er, et, r) * dtheta / theta(r));
  }

 private:
  Scalar theta(Scalar r) const { return 1.0 + alpha_ * dt_(r); }
  /// S_i of the principal strain ei.
  Scalar pk2(Scalar ei, Scalar er, Scalar et, Scalar r) const {
    const Scalar th = theta(r);
    const Scalar eth = 0.5 * (th * th - 1.0);
    const Scalar gr = er + 0.5 * er * er;
    const Scalar gt = et + 0.5 * et * et;
    const Scalar gi = ei + 0.5 * ei * ei;
    return (lambda_ * (gr + gt - 3.0 * eth) + 2.0 * mu_ * (gi - eth)) / th;
  }
  Scalar lambda_;
  Scalar mu_;
  Scalar alpha_;
  std::function<Scalar(Scalar)> dt_;
  std::function<Scalar(Scalar)> dt_slope_;
};

/// Exact radial deformation r = R + u(R) of a long thick tube (a <= R <= b,
/// plane strain) of a TubeMaterial under a pressure p on its bore that
/// follows the deforming surface and a steady rotation about its axis, whose
/// body force rho omega^2 r acts at the deformed radius. With the strains
/// e_r = u' and e_theta = u / R, the radial equilibrium in the reference
/// configuration, dP_r/dR + (P_r - P_theta)/R + rho omega^2 (R + u) = 0,
/// becomes
/// \f[
///   e_r' = \Big[\frac{P_\theta - P_r}{R} - \rho\omega^2 (R + u)
///          - \frac{\partial P_r}{\partial e_\theta}\,\frac{e_r - e_\theta}{R}
///          - \frac{\partial P_r}{\partial R}\Big]
///          \Big/ \frac{\partial P_r}{\partial e_r} .
/// \f]
/// The bore carries the pressure on its deformed area, so its nominal
/// traction is P_r(a) = -p (1 + e_theta(a)); the outer surface is free,
/// P_r(b) = 0. For a trial bore displacement u(a) the first condition gives
/// e_r(a) (by bisection: P_r grows with e_r), the classical Runge-Kutta method
/// integrates u and e_r from a to b in n steps, and bisection on u(a) makes
/// P_r(b) vanish.
class FiniteTube {
 public:
  FiniteTube(const TubeMaterial& law, Scalar a, Scalar b, Scalar p, Scalar rho_omega2, int n,
             Scalar u_scale, Scalar stress_scale)
      : law_(law), a_(a), b_(b), p_(p), rho_omega2_(rho_omega2), n_(n),
        stress_scale_(stress_scale) {
    // Bracket u(a) from zero, upwards and then downwards, in steps of a
    // tenth of the estimate u_scale, then bisect.
    Scalar lo = 0.0;
    Scalar hi = 0.0;
    const Scalar g0 = outer_traction_for(0.0);
    bool found = false;
    for (const Scalar direction : {1.0, -1.0}) {
      Scalar previous = 0.0;
      for (int i = 1; i <= 1000 && !found; ++i) {
        const Scalar trial = direction * 0.1 * u_scale * i;
        Scalar g = 0.0;
        try {
          g = outer_traction_for(trial);
        } catch (const SolverError&) {
          break;  // the law's range ends before a sign change
        }
        if ((g < 0.0) != (g0 < 0.0)) {
          lo = std::min(previous, trial);
          hi = std::max(previous, trial);
          found = true;
        }
        previous = trial;
      }
      if (found) break;
    }
    if (!found) {
      throw SolverError("the tube's shooting found no bore displacement that frees its outer "
                        "surface; the load may exceed the tube's limit");
    }
    Scalar g_lo = outer_traction_for(lo);
    for (int it = 0; it < 300 && hi - lo > 1.0e-16 * std::max(std::abs(lo), std::abs(hi)); ++it) {
      const Scalar mid = 0.5 * (lo + hi);
      const Scalar g = outer_traction_for(mid);
      if ((g < 0.0) == (g_lo < 0.0)) {
        lo = mid;
        g_lo = g;
      } else {
        hi = mid;
      }
    }
    integrate(0.5 * (lo + hi), &u_, &er_, outer_);
  }

  /// Radial displacement u(R), by cubic Hermite interpolation of the
  /// integration points with their slopes e_r.
  Scalar u(Scalar radius) const {
    const Scalar h = (b_ - a_) / static_cast<Scalar>(n_);
    const Scalar x = std::clamp((radius - a_) / h, 0.0, static_cast<Scalar>(n_));
    const int i = std::min(static_cast<int>(x), n_ - 1);
    const Scalar t = x - static_cast<Scalar>(i);
    const Scalar h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
    const Scalar h10 = t * (1.0 - t) * (1.0 - t);
    const Scalar h01 = t * t * (3.0 - 2.0 * t);
    const Scalar h11 = t * t * (t - 1.0);
    const auto k = static_cast<std::size_t>(i);
    return h00 * u_[k] + h10 * h * er_[k] + h01 * u_[k + 1] + h11 * h * er_[k + 1];
  }
  Scalar hoop_stretch(Scalar radius) const { return 1.0 + u(radius) / radius; }
  /// |P_r(b)| over the stress scale: the mismatch of the free outer surface.
  Scalar outer_mismatch() const { return std::abs(outer_) / stress_scale_; }
  /// Largest residual of the radial equilibrium over the wall by central
  /// differences of the integrated P_r, times (b - a) over the stress scale.
  Scalar equilibrium_residual() const {
    const Scalar h = (b_ - a_) / static_cast<Scalar>(n_);
    Scalar worst = 0.0;
    for (int i = 1; i < n_; ++i) {
      const auto k = static_cast<std::size_t>(i);
      const Scalar r = a_ + h * i;
      const Scalar et = u_[k] / r;
      const Scalar pr_plus = law_.radial(er_[k + 1], u_[k + 1] / (r + h), r + h);
      const Scalar pr_minus = law_.radial(er_[k - 1], u_[k - 1] / (r - h), r - h);
      const Scalar residual = (pr_plus - pr_minus) / (2.0 * h) +
                              (law_.radial(er_[k], et, r) - law_.hoop(er_[k], et, r)) / r +
                              rho_omega2_ * (r + u_[k]);
      worst = std::max(worst, std::abs(residual) * (b_ - a_) / stress_scale_);
    }
    return worst;
  }

 private:
  /// e_r at the bore for a bore displacement u_a: P_r(e_r, e_t) = -p (1 + e_t).
  Scalar bore_radial_strain(Scalar ua) const {
    const Scalar et = ua / a_;
    const Scalar target = -p_ * (1.0 + et);
    Scalar lo = -0.99;
    Scalar hi = 10.0;
    if (!(law_.radial(lo, et, a_) < target) || !(law_.radial(hi, et, a_) > target)) {
      throw SolverError("the tube's bore condition has no radial strain in (-0.99, 10)");
    }
    for (int it = 0; it < 200 && hi - lo > 1.0e-17; ++it) {
      const Scalar mid = 0.5 * (lo + hi);
      (law_.radial(mid, et, a_) < target ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
  }

  /// Integrate from the bore displacement ua; P_r at b in `outer`.
  void integrate(Scalar ua, std::vector<Scalar>* u, std::vector<Scalar>* er,
                 Scalar& outer) const {
    const auto f = [this](Scalar r, const std::array<Scalar, 2>& y) {
      const Scalar e_r = y[1];
      const Scalar e_t = y[0] / r;
      const Scalar slope = law_.radial_by_er(e_r, e_t, r);
      if (!(slope > 0.0) || !(e_r > -0.99)) {
        throw SolverError("the tube's integration left the range of the law");
      }
      const Scalar rhs = (law_.hoop(e_r, e_t, r) - law_.radial(e_r, e_t, r)) / r -
                         rho_omega2_ * (r + y[0]) -
                         law_.radial_by_et(e_r, e_t, r) * (e_r - e_t) / r -
                         law_.radial_by_radius(e_r, e_t, r);
      return std::array<Scalar, 2>{e_r, rhs / slope};
    };
    std::array<Scalar, 2> y{ua, bore_radial_strain(ua)};
    const Scalar h = (b_ - a_) / static_cast<Scalar>(n_);
    if (u != nullptr) {
      u->assign(static_cast<std::size_t>(n_) + 1, 0.0);
      er->assign(static_cast<std::size_t>(n_) + 1, 0.0);
      (*u)[0] = y[0];
      (*er)[0] = y[1];
    }
    for (int i = 0; i < n_; ++i) {
      const Scalar r = a_ + h * i;
      const std::array<Scalar, 2> k1 = f(r, y);
      const std::array<Scalar, 2> k2 = f(r + 0.5 * h, {y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]});
      const std::array<Scalar, 2> k3 = f(r + 0.5 * h, {y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]});
      const std::array<Scalar, 2> k4 = f(r + h, {y[0] + h * k3[0], y[1] + h * k3[1]});
      for (int j = 0; j < 2; ++j) y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
      if (u != nullptr) {
        (*u)[static_cast<std::size_t>(i) + 1] = y[0];
        (*er)[static_cast<std::size_t>(i) + 1] = y[1];
      }
    }
    outer = law_.radial(y[1], y[0] / b_, b_);
  }

  Scalar outer_traction_for(Scalar ua) const {
    Scalar outer = 0.0;
    integrate(ua, nullptr, nullptr, outer);
    return outer;
  }

  const TubeMaterial& law_;
  Scalar a_;
  Scalar b_;
  Scalar p_;
  Scalar rho_omega2_;
  int n_;
  Scalar stress_scale_;
  std::vector<Scalar> u_;
  std::vector<Scalar> er_;
  Scalar outer_ = 0.0;
};

/// Lame's small-strain bore displacement of the tube in plane strain.
Scalar lame_u(Scalar r, Scalar a, Scalar b, Scalar p, Scalar youngs, Scalar poisson) {
  const Scalar big_a = p * a * a / (b * b - a * a);
  const Scalar big_b = p * a * a * b * b / (b * b - a * a);
  return (1.0 + poisson) / youngs * ((1.0 - 2.0 * poisson) * big_a * r + big_b / r);
}

int element_order(ElementType type) { return type == ElementType::Tet10 ? 3 : 2; }

std::vector<Index> tube_ladder(ElementType type) {
  switch (type) {
    case ElementType::Quad4:
    case ElementType::Tri3: return {4, 8, 16, 32, 64};
    case ElementType::Hex8: return {4, 8, 16, 32};
    case ElementType::Tet4: return {4, 8, 16};
    case ElementType::Tet10: return {2, 4, 8, 16};
  }
  return {};
}

// ---------------------------------------------------------------------------
// The shallow arch
// ---------------------------------------------------------------------------

/// Half of a shallow circular arch of mid-surface radius `radius`, opening
/// angle 2 phi0 and thickness t: the unit box [0, 1]^2 of a structured Q4
/// mesh mapped by phi = phi0 x (x = 0 the crown, x = 1 the support) and
/// rho = radius - t/2 + t y onto (rho sin phi, rho cos phi - radius cos phi0).
Mesh half_arch(const Mesh& unit, Scalar radius, Scalar phi0, Scalar thickness) {
  Matrix x = unit.coordinates();
  for (Index n = 0; n < unit.num_nodes(); ++n) {
    const Scalar phi = phi0 * x(0, n);
    const Scalar rho = radius - 0.5 * thickness + thickness * x(1, n);
    x(0, n) = x(0, n) == 0.0 ? 0.0 : rho * std::sin(phi);
    x(1, n) = rho * std::cos(phi) - radius * std::cos(phi0);
  }
  Mesh mapped(std::move(x), unit.connectivity(), unit.element_type());
  mapped.validate();
  return mapped;
}

/// Nodes of a unit-box mesh on its side x = value.
std::vector<Index> unit_side_nodes(const Mesh& unit, Scalar value) {
  std::vector<Index> out;
  for (Index n = 0; n < unit.num_nodes(); ++n) {
    if (std::abs(unit.node(n).x() - value) < 1.0e-12) out.push_back(n);
  }
  return out;
}

Index unit_corner_node(const Mesh& unit, Scalar x, Scalar y) {
  for (Index n = 0; n < unit.num_nodes(); ++n) {
    const Vector3 p = unit.node(n);
    if (std::abs(p.x() - x) < 1.0e-12 && std::abs(p.y() - y) < 1.0e-12) return n;
  }
  throw ConfigError("the unit mesh has no node at the requested corner");
}

}  // namespace

// ---------------------------------------------------------------------------
// Euler's elastica
// ---------------------------------------------------------------------------

StudyOutcome study_elastica(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 1.0;
  const Scalar h = 0.01;
  const Scalar youngs = 210.0e9;
  const Scalar inertia = h * h * h * h / 12.0;
  const Scalar k_max = 10.0;
  const int levels = 10;  // k = 1, 2, ..., 10
  const Scalar force = k_max * youngs * inertia / (length * length);
  const std::vector<Index> ladder = {25, 50, 100};
  const Scalar tolerance = 1.0e-3;
  const Scalar order_margin = 0.3;
  const int expected_order = 3;
  const int rk_steps = 20000;

  // The reference, and its own error: the change when the Runge-Kutta steps
  // double.
  std::vector<ElasticaTip> reference(static_cast<std::size_t>(levels) + 1);
  Scalar reference_error = 0.0;
  for (int i = 1; i <= levels; ++i) {
    const Scalar k = k_max * i / levels;
    const ElasticaTip coarse = elastica(k, rk_steps);
    const ElasticaTip fine = elastica(k, 2 * rk_steps);
    reference[static_cast<std::size_t>(i)] = fine;
    reference_error = std::max({reference_error, std::abs(fine.deflection - coarse.deflection),
                                std::abs(fine.shortening - coarse.shortening),
                                std::abs(fine.rotation - coarse.rotation)});
  }

  CsvWriter csv(path_join(out_dir, "elastica.csv"),
                {"nx", "k[-]", "deflection[-]", "deflection_elastica[-]", "deflection_error[-]",
                 "shortening[-]", "shortening_elastica[-]", "shortening_error[-]",
                 "rotation[rad]", "rotation_elastica[rad]", "rotation_error[rad]",
                 "iterations", "negative_pivots"});
  // errors[m][i]: mesh m, level i (1-based)
  std::vector<std::vector<std::array<Scalar, 3>>> errors;
  json::Value meshes = json::Value::make_array();
  bool all_completed = true;
  Scalar worst_balance = 0.0;
  for (const Index nx : ladder) {
    const CantileverRun run = solve_cantilever(nx, length, h, youngs, force, levels);
    const NonlinearResult& r = run.result;
    all_completed = all_completed && r.completed;
    worst_balance = std::max(worst_balance, r.equilibrium.relative_force_error);
    std::vector<std::array<Scalar, 3>> e(static_cast<std::size_t>(levels) + 1,
                                         {std::nan(""), std::nan(""), std::nan("")});
    for (const NonlinearStep& s : r.steps) {
      const Scalar level = s.load_factor * levels;
      const int i = static_cast<int>(std::lround(level));
      if (std::abs(level - i) > 1.0e-9 || i < 1 || i > levels) continue;
      const ElasticaTip& ref = reference[static_cast<std::size_t>(i)];
      const Scalar v = -monitor_value(r, s, "tip_uy") / length;
      const Scalar u = -monitor_value(r, s, "tip_ux") / length;
      // The tip face's rotation: the chord from its bottom to its top edge.
      const Scalar dx = monitor_value(r, s, "top_ux") - monitor_value(r, s, "bottom_ux");
      const Scalar dy = h + monitor_value(r, s, "top_uy") - monitor_value(r, s, "bottom_uy");
      const Scalar theta = std::atan2(dx, dy);
      auto& slot = e[static_cast<std::size_t>(i)];
      slot = {v - ref.deflection, u - ref.shortening, theta - ref.rotation};
      csv.raw_row({fmt(static_cast<Scalar>(nx)), fmt(k_max * i / levels), fmt(v, 10),
                   fmt(ref.deflection, 10), fmt(slot[0], 6), fmt(u, 10), fmt(ref.shortening, 10),
                   fmt(slot[1], 6), fmt(theta, 10), fmt(ref.rotation, 10), fmt(slot[2], 6),
                   fmt(static_cast<Scalar>(s.iterations)),
                   fmt(static_cast<Scalar>(s.negative_pivots))});
    }
    errors.push_back(e);
    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(static_cast<Scalar>(nx)));
    rec.set("num_dofs", json::Value::make_number(static_cast<Scalar>(run.num_dofs)));
    rec.set("completed", json::Value::make_bool(r.completed));
    rec.set("steps", json::Value::make_number(static_cast<Scalar>(r.steps.size())));
    rec.set("iterations", json::Value::make_number(r.total_iterations));
    rec.set("cuts", json::Value::make_number(r.total_cuts));
    rec.set("max_green_strain", json::Value::make_number(r.max_green_strain));
    rec.set("force_balance_error", json::Value::make_number(r.equilibrium.relative_force_error));
    rec.set("moment_balance_error",
            json::Value::make_number(r.equilibrium.relative_moment_error));
    rec.set("seconds", json::Value::make_number(run.seconds));
    meshes.push_back(rec);
  }
  csv.close();

  // Finest-mesh errors, the observed orders and the extrapolated gaps.
  const std::vector<std::array<Scalar, 3>>& finest = errors.back();
  Scalar worst_finest = 0.0;
  Scalar worst_rotation = 0.0;
  Scalar lowest_order = 1.0e300;
  Scalar largest_gap = 0.0;
  json::Value levels_json = json::Value::make_array();
  CsvWriter order_csv(path_join(out_dir, "elastica_orders.csv"),
                      {"k[-]", "deflection_order[-]", "shortening_order[-]",
                       "deflection_gap[-]", "shortening_gap[-]"});
  for (int i = 1; i <= levels; ++i) {
    const auto k = static_cast<std::size_t>(i);
    for (int c = 0; c < 3; ++c) {
      if (std::isnan(finest[k][static_cast<std::size_t>(c)])) all_completed = false;
    }
    worst_finest = std::max({worst_finest, std::abs(finest[k][0]), std::abs(finest[k][1])});
    worst_rotation = std::max(worst_rotation, std::abs(finest[k][2]));
    Scalar orders[2];
    Scalar gaps[2];
    for (int c = 0; c < 2; ++c) {
      const auto cc = static_cast<std::size_t>(c);
      orders[c] = three_mesh_order(errors[0][k][cc], errors[1][k][cc], errors[2][k][cc]);
      gaps[c] = three_mesh_limit(errors[0][k][cc], errors[1][k][cc], errors[2][k][cc]);
      lowest_order = std::isnan(orders[c]) ? -1.0e300 : std::min(lowest_order, orders[c]);
      largest_gap = std::max(largest_gap, std::abs(gaps[c]));
    }
    order_csv.raw_row({fmt(k_max * i / levels), fmt(orders[0], 4), fmt(orders[1], 4),
                       fmt(gaps[0], 4), fmt(gaps[1], 4)});
    json::Value rec = json::Value::make_object();
    rec.set("k", json::Value::make_number(k_max * i / levels));
    rec.set("deflection_elastica", json::Value::make_number(reference[k].deflection));
    rec.set("shortening_elastica", json::Value::make_number(reference[k].shortening));
    rec.set("rotation_elastica_rad", json::Value::make_number(reference[k].rotation));
    rec.set("deflection_error_finest", json::Value::make_number(finest[k][0]));
    rec.set("shortening_error_finest", json::Value::make_number(finest[k][1]));
    rec.set("rotation_error_finest_rad", json::Value::make_number(finest[k][2]));
    rec.set("deflection_order", json::Value::make_number(orders[0]));
    rec.set("shortening_order", json::Value::make_number(orders[1]));
    rec.set("deflection_continuum_gap", json::Value::make_number(gaps[0]));
    rec.set("shortening_continuum_gap", json::Value::make_number(gaps[1]));
    levels_json.push_back(rec);
  }
  order_csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string(
                        "verification (the observed order of the discretisation error) + "
                        "validation (the gap between the continuum and beam theory)"));
  block.set("meshes", meshes);
  block.set("levels", levels_json);
  block.set("reference_rk4_change_on_doubling", json::Value::make_number(reference_error));
  block.set("tolerance_finest_error", json::Value::make_number(tolerance));
  block.set("lowest_observed_order", json::Value::make_number(lowest_order));
  block.set("expected_order", json::Value::make_number(expected_order));
  block.set("note",
            json::Value::make_string(
                "Cantilever L = 1 m of square section h = 0.01 m (E = 210 GPa, nu = 0), Tet10 "
                "cells nx x 2 x 2, clamped at x = 0, under a dead shear traction on its tip "
                "face whose resultant P gives k = P L^2 / EI = 1 ... 10 at the load factors "
                "0.1 ... 1 (Saint Venant-Kirchhoff, load control). Deflection, shortening "
                "(in units of L) and the rotation of the tip face against Euler's elastica, "
                "solved by shooting. The continuum differs from the elastica by shear "
                "(0.6 (h/L)^2 = 6e-5 in the linear range for nu = 0) and by the Saint "
                "Venant-Kirchhoff bending stiffness, M = EI c (1 - 0.3 (c h)^2) for the "
                "curvature c (c h <= 0.045 here): the observed order uses the three meshes "
                "alone (Richardson), and the gap is where they extrapolate to."));
  summary.set("elastica", block);

  std::ostringstream note;
  note << "finest mesh (nx = " << ladder.back() << ") largest |error| " << app::format(worst_finest, 3)
       << " L, rotation " << app::format(worst_rotation, 3) << " rad; lowest observed order "
       << app::format(lowest_order, 3) << "; continuum-elastica gap <= "
       << app::format(largest_gap, 3) << " L; reference change on doubling "
       << app::format(reference_error, 2);
  StudyOutcome outcome;
  outcome.name = "large-deflection cantilever vs Euler's elastica (Tet10, k <= 10)";
  outcome.kind = "verification + validation";
  outcome.metric = "largest tip-displacement error on the finest mesh [L]";
  outcome.value = worst_finest;
  outcome.tolerance = tolerance;
  outcome.passed = all_completed && worst_finest <= tolerance && worst_rotation <= tolerance &&
                   lowest_order >= expected_order - order_margin && reference_error < 1.0e-12 &&
                   worst_balance < 1.0e-8;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// The tube at finite strain: inflated, spinning, heated
// ---------------------------------------------------------------------------

StudyOutcome study_hyperelastic_cylinder(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 0.1;
  const Scalar b = 0.2;
  const int rk_steps = 4000;
  const Scalar order_margin = 0.3;

  // The three loadings of the quarter section.
  struct TubeCase {
    std::string name;
    HyperelasticModel law;
    Scalar youngs;
    Scalar poisson;
    Scalar density;
    Scalar alpha;
    Scalar pressure;   ///< follower pressure on the bore [Pa]
    Scalar omega;      ///< rotation about the axis [rad/s]
    Scalar bore_dt;    ///< bore temperature (outer 0, reference 0) [K]
  };
  const std::vector<TubeCase> cases = {
      {"inflation", HyperelasticModel::NeoHookean, 10.0e6, 0.3, 1100.0, 0.0, 1.5e6, 0.0, 0.0},
      {"spin", HyperelasticModel::NeoHookean, 10.0e6, 0.3, 1100.0, 0.0, 0.0, 200.0, 0.0},
      {"heating", HyperelasticModel::SaintVenantKirchhoff, 1.0e9, 0.3, 7800.0, 5.0e-4, 0.0, 0.0,
       100.0}};
  // The steady conducted temperature of the heated case, and its slope.
  const Scalar log_ratio = std::log(b / a);

  CsvWriter csv(path_join(out_dir, "hyperelastic_cylinder.csv"),
                {"case", "element", "n_r", "n_theta", "h[m]", "num_dofs", "u_rms_error[-]",
                 "u_rms_order[-]", "u_max_error[-]", "u_max_order[-]", "bore_hoop_stretch[-]",
                 "bore_hoop_stretch_exact[-]", "steps", "iterations", "linear_solver",
                 "seconds"});
  json::Value block = json::Value::make_object();
  Scalar worst_shortfall = -1.0e300;
  bool all_completed = true;
  Scalar worst_balance = 0.0;
  Scalar worst_self_check = 0.0;
  Scalar worst_reference_change = 0.0;
  std::ostringstream note;
  for (const TubeCase& tc : cases) {
    const Scalar lambda = tc.youngs * tc.poisson / ((1.0 + tc.poisson) * (1.0 - 2.0 * tc.poisson));
    const Scalar mu = tc.youngs / (2.0 * (1.0 + tc.poisson));
    const Scalar bore_dt = tc.bore_dt;
    const auto dt = [&, bore_dt](Scalar r) { return bore_dt * std::log(b / r) / log_ratio; };
    const auto dt_slope = [&, bore_dt](Scalar r) { return -bore_dt / (r * log_ratio); };
    std::unique_ptr<TubeMaterial> law;
    if (tc.law == HyperelasticModel::NeoHookean) {
      law = std::make_unique<NeoHookeanTubeLaw>(mu, lambda);
    } else {
      law = std::make_unique<ThermalSvkTubeLaw>(lambda, mu, tc.alpha, dt, dt_slope);
    }
    const Scalar rho_omega2 = tc.density * tc.omega * tc.omega;
    const Scalar stress_scale = std::max({tc.pressure, rho_omega2 * b * b,
                                          tc.youngs * tc.alpha * tc.bore_dt});
    const FiniteTube exact(*law, a, b, tc.pressure, rho_omega2, rk_steps, 0.01 * a,
                           stress_scale);
    const FiniteTube half_steps(*law, a, b, tc.pressure, rho_omega2, rk_steps / 2, 0.01 * a,
                                stress_scale);
    Scalar reference_change = 0.0;
    for (int i = 0; i <= 100; ++i) {
      const Scalar r = a + (b - a) * i / 100.0;
      reference_change = std::max(reference_change, std::abs(exact.u(r) - half_steps.u(r)));
    }
    reference_change /= std::abs(exact.u(a));
    const Scalar self_check = std::max(exact.outer_mismatch(), exact.equilibrium_residual());
    worst_self_check = std::max(worst_self_check, self_check);
    worst_reference_change = std::max(worst_reference_change, reference_change);

    IsotropicMaterial material(tc.youngs, tc.poisson, tc.density, "tube");
    if (tc.alpha != 0.0) material.set_thermal(tc.alpha, 0.0, 1.0);
    json::Value case_block = json::Value::make_object();
    for (const ElementType type :
         {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet10}) {
      const std::string name = to_string(type);
      std::vector<Scalar> hs;
      std::vector<Scalar> rms;
      std::vector<Scalar> maxes;
      json::Value records = json::Value::make_array();
      for (const Index nr : tube_ladder(type)) {
        const Index nt = 2 * nr;
        Mesh mesh = sector_mesh(type, nr, nt, a, b);
        const int dim = mesh.dim();
        FemModel model(std::move(mesh), material, 1.0,
                       dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional,
                       IntegrationOptions());
        add_quarter_supports(model);
        LoadCaseSpec load;
        load.name = tc.name;
        if (tc.pressure != 0.0) {
          PressureLoadSpec bore;
          bore.region = cylinder_surface("bore", a, true);
          bore.pressure = tc.pressure;
          load.pressures.push_back(bore);
        }
        if (tc.omega != 0.0) {
          load.centrifugal.enabled = true;
          load.centrifugal.angular_velocity = tc.omega;
          load.centrifugal.axis = Vector3::UnitZ();
        }
        if (tc.bore_dt != 0.0) {
          load.temperature.source = TemperatureSpec::Source::Conduction;
          RegionValue hot;
          hot.region = cylinder_surface("bore", a, true);
          hot.value = tc.bore_dt;
          RegionValue cold;
          cold.region = cylinder_surface("outer", b, false);
          cold.value = 0.0;
          load.temperature.conduction.prescribed = {hot, cold};
        }
        model.load_case_specs().push_back(load);
        model.finalize();
        Assembler assembler(model);
        NonlinearOptions options;
        options.law = tc.law;
        options.steps = 10;
        options.follower_pressure = true;
        Timer timer;
        NonlinearStaticAnalysis analysis(model, assembler, options);
        const NonlinearResult r = analysis.solve(0);
        const Scalar seconds = timer.elapsed_seconds();
        all_completed = all_completed && r.completed;
        worst_balance = std::max(worst_balance, r.equilibrium.relative_force_error);

        const Mesh& m = model.mesh();
        Scalar err_sq = 0.0;
        Scalar ref_sq = 0.0;
        Scalar err_max = 0.0;
        Scalar ref_max = 0.0;
        Scalar bore_stretch = 0.0;
        Index bore_nodes = 0;
        for (Index n = 0; n < m.num_nodes(); ++n) {
          const Vector3 x = m.node(n);
          const Scalar radius = std::hypot(x.x(), x.y());
          const Scalar ur = exact.u(radius);
          const Scalar ux = r.displacement(dim * n);
          const Scalar uy = r.displacement(dim * n + 1);
          const Scalar err = std::hypot(ux - ur * x.x() / radius, uy - ur * x.y() / radius);
          err_sq += err * err;
          ref_sq += ur * ur;
          err_max = std::max(err_max, err);
          ref_max = std::max(ref_max, std::abs(ur));
          if (radius <= a * (1.0 + 1.0e-9)) {
            bore_stretch += std::hypot(x.x() + ux, x.y() + uy) / radius;
            ++bore_nodes;
          }
        }
        hs.push_back((b - a) / static_cast<Scalar>(nr));
        rms.push_back(std::sqrt(err_sq / ref_sq));
        maxes.push_back(err_max / ref_max);
        bore_stretch /= static_cast<Scalar>(std::max<Index>(bore_nodes, 1));
        const std::size_t last = rms.size() - 1;
        const Scalar rms_order = last > 0 ? observed_order(hs[last - 1], rms[last - 1], hs[last],
                                                           rms[last])
                                          : std::nan("");
        const Scalar max_order = last > 0 ? observed_order(hs[last - 1], maxes[last - 1],
                                                           hs[last], maxes[last])
                                          : std::nan("");
        csv.raw_row({tc.name, name, fmt(static_cast<Scalar>(nr)), fmt(static_cast<Scalar>(nt)),
                     fmt(hs.back(), 8), fmt(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                     fmt(rms.back(), 8), fmt(rms_order, 4), fmt(maxes.back(), 8),
                     fmt(max_order, 4), fmt(bore_stretch, 10), fmt(exact.hoop_stretch(a), 10),
                     fmt(static_cast<Scalar>(r.steps.size())),
                     fmt(static_cast<Scalar>(r.total_iterations)), r.linear_solver,
                     fmt(seconds, 4)});
        json::Value rec = json::Value::make_object();
        rec.set("n_r", json::Value::make_number(static_cast<Scalar>(nr)));
        rec.set("num_dofs",
                json::Value::make_number(static_cast<Scalar>(model.dofs().num_dofs())));
        rec.set("displacement_rms_error", json::Value::make_number(rms.back()));
        rec.set("displacement_max_error", json::Value::make_number(maxes.back()));
        rec.set("bore_hoop_stretch", json::Value::make_number(bore_stretch));
        rec.set("completed", json::Value::make_bool(r.completed));
        rec.set("iterations", json::Value::make_number(r.total_iterations));
        rec.set("max_green_strain", json::Value::make_number(r.max_green_strain));
        rec.set("min_jacobian", json::Value::make_number(r.min_jacobian));
        rec.set("force_balance_error",
                json::Value::make_number(r.equilibrium.relative_force_error));
        rec.set("linear_solver", json::Value::make_string(r.linear_solver));
        records.push_back(rec);
      }
      const std::size_t last = rms.size() - 1;
      const Scalar order = observed_order(hs[last - 1], rms[last - 1], hs[last], rms[last]);
      worst_shortfall = std::max(worst_shortfall, element_order(type) - order);
      note << tc.name << " " << name << " " << app::format(rms.back(), 3) << " (order "
           << app::format(order, 3) << "); ";
      json::Value entry = json::Value::make_object();
      entry.set("meshes", records);
      entry.set("displacement_rms_order", json::Value::make_number(order));
      entry.set("expected_order", json::Value::make_number(element_order(type)));
      case_block.set(name, entry);
    }
    case_block.set("law", json::Value::make_string(to_string(tc.law)));
    case_block.set("bore_hoop_stretch_exact", json::Value::make_number(exact.hoop_stretch(a)));
    case_block.set("outer_hoop_stretch_exact", json::Value::make_number(exact.hoop_stretch(b)));
    case_block.set("reference_change_on_halving", json::Value::make_number(reference_change));
    case_block.set("reference_self_check", json::Value::make_number(self_check));
    block.set(tc.name, case_block);
  }
  csv.close();

  // The small-strain limit of the reference: at p1 and 2 p1 the bore
  // displacement of the inflated tube differs from Lame's by c p + O(p^2),
  // so 2 d(p1) - d(2 p1) must vanish to O(p1^2).
  const Scalar youngs = cases[0].youngs;
  const Scalar poisson = cases[0].poisson;
  const NeoHookeanTubeLaw small_law(youngs / (2.0 * (1.0 + poisson)),
                                    youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)));
  const auto lame_ratio = [&](Scalar p) {
    const Scalar u_lame = lame_u(a, a, b, p, youngs, poisson);
    const FiniteTube small(small_law, a, b, p, 0.0, rk_steps, u_lame, p);
    return small.u(a) / u_lame - 1.0;
  };
  const Scalar p1 = 1.0e-7 * youngs;
  const Scalar small_strain_limit = std::abs(2.0 * lame_ratio(p1) - lame_ratio(2.0 * p1));

  block.set("kind", json::Value::make_string("verification (exact finite-elasticity solutions)"));
  block.set("reference_small_strain_limit_vs_lame", json::Value::make_number(small_strain_limit));
  block.set("note",
            json::Value::make_string(
                "Quarter section of a tube a = 0.1 m, b = 0.2 m in plane strain, in 10 "
                "load-control steps, three ways: 'inflation' - compressible neo-Hookean "
                "(E = 10 MPa, nu = 0.3) under a bore pressure of 1.5 MPa that follows the "
                "deforming bore; 'spin' - the same material (rho = 1100 kg/m^3) spinning at "
                "200 rad/s about its axis, the centrifugal load at the deformed radius; "
                "'heating' - Saint Venant-Kirchhoff (E = 1 GPa, nu = 0.3, alpha = 5e-4 /K) "
                "under the conducted temperature of a bore at +100 K and an outer surface at "
                "0 K, the multiplicative thermal split. Symmetry supports; the 3-D sections "
                "one cell deep with u_z = 0. Error: the RMS nodal displacement error relative "
                "to the RMS exact displacement. Reference: the radial equilibrium of the "
                "finite deformation, a two-point boundary-value problem solved by shooting "
                "with fourth-order Runge-Kutta; the inflation reduces to Lame's solution as "
                "the pressure vanishes."));
  summary.set("hyperelastic_cylinder", block);

  note << "reference small-strain limit vs Lame " << app::format(small_strain_limit, 2);
  StudyOutcome outcome;
  outcome.name = "thick tube at finite strain vs exact: follower pressure, spin, heating";
  outcome.kind = "verification";
  outcome.metric = "largest shortfall of the displacement convergence order";
  outcome.value = worst_shortfall;
  outcome.tolerance = order_margin;
  outcome.passed = all_completed && worst_shortfall <= order_margin && worst_self_check < 1.0e-6 &&
                   worst_reference_change < 1.0e-9 && small_strain_limit < 1.0e-9 &&
                   worst_balance < 1.0e-8;
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Snap-through of a shallow arch
// ---------------------------------------------------------------------------

StudyOutcome study_arch_snap_through(const std::string& out_dir, json::Value& summary) {
  // Clamped shallow circular arch: half span S and rise H of the mid-surface,
  // thickness t, out-of-plane width w (plane stress), Q4 cells.
  const Scalar span = 1.0;
  const Scalar rise = 0.1;
  const Scalar thickness = 0.02;
  const Scalar width = 0.02;
  const Scalar youngs = 70.0e9;
  const Scalar poisson = 0.3;
  const Index n_along = 60;
  const Index n_through = 4;
  const Scalar radius = (span * span + rise * rise) / (2.0 * rise);
  const Scalar phi0 = std::asin(span / radius);
  const Scalar force = 20.0e3;  // on the half model at lambda = 1
  const Scalar tolerance = 1.0e-6;

  StructuredMeshSpec ms;
  ms.nx = n_along;
  ms.ny = n_through;
  ms.lx = 1.0;
  ms.ly = 1.0;
  const Mesh unit = make_structured_quad_mesh(ms);
  const std::vector<Index> support = unit_side_nodes(unit, 1.0);
  const Index crown = unit_corner_node(unit, 0.0, 1.0);

  // One model per run: the force-controlled one, or the crown driven down by
  // `drive` (> 0) with no force.
  const auto build = [&](Scalar drive) {
    FemModel model(half_arch(unit, radius, phi0, thickness),
                   IsotropicMaterial(youngs, poisson, 2700.0, "aluminium"), width,
                   StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint clamp;
    clamp.region = node_region("clamped_support", support);
    clamp.fix_x = true;
    clamp.fix_y = true;
    model.constraints().push_back(clamp);
    DisplacementConstraint symmetry;
    Selector axis = box();
    axis.xmax = 0.0;
    symmetry.region = region_of("crown_symmetry", axis);
    symmetry.fix_x = true;
    model.constraints().push_back(symmetry);
    LoadCaseSpec load;
    load.name = drive > 0.0 ? "crown_displacement" : "crown_force";
    if (drive > 0.0) {
      DisplacementConstraint push;
      push.region = node_region("crown", {crown});
      push.fix_y = true;
      push.value_y = -drive;
      model.constraints().push_back(push);
      load.prescribed_displacement_only = true;
    } else {
      PointLoadSpec p;
      p.region = node_region("crown", {crown});
      p.force = Vector3(0.0, -force, 0.0);
      load.point_loads.push_back(p);
    }
    model.load_case_specs().push_back(load);
    model.finalize();
    return model;
  };
  const std::vector<NonlinearMonitor> monitors = {
      monitor("crown_uy", node_region("crown", {crown}), 1),
      monitor("crown_ry", node_region("crown", {crown}), 1, NonlinearMonitor::Quantity::Reaction)};

  // 1. The arc-length path under the crown force.
  FemModel forced = build(0.0);
  const Assembler forced_assembler(forced);
  NonlinearOptions arc;
  arc.method = NonlinearOptions::Method::ArcLength;
  arc.steps = 40;
  arc.max_steps = 400;
  arc.target_load_factor = 1.0;
  arc.monitors = monitors;
  Timer arc_timer;
  const NonlinearResult path = NonlinearStaticAnalysis(forced, forced_assembler, arc).solve(0);
  const Scalar arc_seconds = arc_timer.elapsed_seconds();

  std::vector<Scalar> deflection;  // crown deflection (downwards) per arc-length step
  std::vector<Scalar> path_force;  // crown force lambda * P
  for (const NonlinearStep& s : path.steps) {
    deflection.push_back(-monitor_value(path, s, "crown_uy"));
    path_force.push_back(s.load_factor * force);
  }
  bool monotone = path.completed && !deflection.empty();
  for (std::size_t i = 1; i < deflection.size(); ++i) {
    monotone = monotone && deflection[i] > deflection[i - 1];
  }
  if (!monotone) {
    throw SolverError("the arch's arc-length path did not complete with a crown deflection "
                      "that grows at every step; the displacement-controlled comparison needs "
                      "one");
  }
  Scalar peak = 0.0;
  std::size_t peak_step = 0;
  for (std::size_t i = 0; i < path_force.size(); ++i) {
    if (path_force[i] > peak) {
      peak = path_force[i];
      peak_step = i;
    }
    if (i > 0 && path_force[i] < path_force[i - 1]) break;  // the first limit point
  }

  // 2. Displacement control through the same crown deflections.
  const Scalar drive = deflection.back();
  FemModel driven = build(drive);
  const Assembler driven_assembler(driven);
  NonlinearOptions control;
  control.steps = 40;
  control.max_steps = 1000;
  control.monitors = monitors;
  for (std::size_t i = 0; i + 1 < deflection.size(); ++i) {
    control.load_factors.push_back(deflection[i] / drive);
  }
  Timer control_timer;
  const NonlinearResult controlled =
      NonlinearStaticAnalysis(driven, driven_assembler, control).solve(0);
  const Scalar control_seconds = control_timer.elapsed_seconds();

  // 3. Load control under the crown force: it must stop at the limit point.
  NonlinearOptions load_control;
  load_control.steps = 40;
  load_control.monitors = monitors;
  const NonlinearResult stopped =
      NonlinearStaticAnalysis(forced, forced_assembler, load_control).solve(0);

  // The vertex of the parabola through three samples (x_i, y_i) of a path
  // around its maximum; the middle sample when they do not bend down.
  const auto parabola_vertex = [](Scalar x0, Scalar x1, Scalar x2, Scalar y0, Scalar y1,
                                  Scalar y2) {
    const Scalar f01 = (y1 - y0) / (x1 - x0);
    const Scalar f12 = (y2 - y1) / (x2 - x1);
    const Scalar f012 = (f12 - f01) / (x2 - x0);
    if (!(f012 < 0.0)) return y1;
    const Scalar xs = 0.5 * (x0 + x1) - f01 / (2.0 * f012);
    return y0 + f01 * (xs - x0) + f012 * (xs - x0) * (xs - x1);
  };
  if (peak_step == 0 || peak_step + 1 >= deflection.size()) {
    throw SolverError("the arch's arc-length path has no sample on each side of its first "
                      "limit point");
  }
  // A first estimate from the arc-length samples, which lie far apart.
  const Scalar coarse_limit =
      parabola_vertex(deflection[peak_step - 1], deflection[peak_step], deflection[peak_step + 1],
                      path_force[peak_step - 1], path_force[peak_step],
                      path_force[peak_step + 1]);

  // 4. The limit load itself: displacement control through `fine` stations
  // across the two sample intervals around the highest sample, and the
  // vertex of the parabola through the highest of them and its neighbours.
  // Its error estimate is the vertex's distance from the highest sample: the
  // parabola corrects the sample by far more than it errs itself.
  const int fine = 60;
  const Scalar fine_lo = deflection[peak_step - 1];
  const Scalar fine_hi = deflection[peak_step + 1];
  FemModel refined = build(fine_hi);
  const Assembler refined_assembler(refined);
  NonlinearOptions sweep;
  sweep.steps = 40;
  sweep.max_steps = 1000;
  sweep.monitors = monitors;
  for (int i = 0; i < fine; ++i) {
    sweep.load_factors.push_back((fine_lo + (fine_hi - fine_lo) * i / fine) / fine_hi);
  }
  const NonlinearResult swept =
      NonlinearStaticAnalysis(refined, refined_assembler, sweep).solve(0);
  std::vector<Scalar> sweep_deflection;
  std::vector<Scalar> sweep_force;
  for (const NonlinearStep& s : swept.steps) {
    const Scalar d = s.load_factor * fine_hi;
    if (d < fine_lo * (1.0 - 1.0e-9)) continue;
    sweep_deflection.push_back(d);
    sweep_force.push_back(-monitor_value(swept, s, "crown_ry"));
  }
  std::size_t top = 0;
  for (std::size_t i = 1; i < sweep_force.size(); ++i) {
    if (sweep_force[i] > sweep_force[top]) top = i;
  }
  if (!swept.completed || top == 0 || top + 1 >= sweep_force.size()) {
    throw SolverError("the displacement-controlled sweep around the arch's limit point did "
                      "not bracket its maximum");
  }
  const Scalar limit =
      parabola_vertex(sweep_deflection[top - 1], sweep_deflection[top], sweep_deflection[top + 1],
                      sweep_force[top - 1], sweep_force[top], sweep_force[top + 1]);
  const Scalar limit_error = std::abs(limit - sweep_force[top]);

  // Compare at the stations, which are the arc-length deflections.
  std::map<long long, const NonlinearStep*> at_station;
  for (const NonlinearStep& s : controlled.steps) {
    at_station[std::llround(s.load_factor * 1.0e12)] = &s;
  }
  std::vector<Scalar> control_force(deflection.size(), std::nan(""));
  std::vector<int> control_pivots(deflection.size(), -1);
  for (std::size_t i = 0; i < deflection.size(); ++i) {
    const Scalar station = i + 1 < deflection.size() ? deflection[i] / drive : 1.0;
    const auto it = at_station.find(std::llround(station * 1.0e12));
    if (it == at_station.end()) continue;
    control_force[i] = -monitor_value(controlled, *it->second, "crown_ry");
    control_pivots[i] = it->second->negative_pivots;
  }
  Scalar largest_difference = 0.0;
  bool all_matched = controlled.completed;
  for (std::size_t i = 0; i < deflection.size(); ++i) {
    if (std::isnan(control_force[i])) {
      all_matched = false;
      continue;
    }
    largest_difference = std::max(largest_difference, std::abs(control_force[i] - path_force[i]));
  }
  largest_difference /= limit;

  // Inertia: with the crown deflection prescribed the free set loses one DOF,
  // and by Haynsworth's inertia additivity the force-controlled tangent has
  // exactly one more negative eigenvalue than the displacement-controlled one
  // where the crown stiffness dF/d(deflection) - the Schur complement - is
  // negative. The slope is read off the displacement-controlled path by
  // central differences; points next to a sign change of the slope (the limit
  // points) are skipped.
  int checked = 0;
  int mismatches = 0;
  int unstable_steps = 0;
  CsvWriter csv(path_join(out_dir, "arch_snap_through.csv"),
                {"step", "crown_deflection[m]", "load_factor[-]", "force_arc_length[N]",
                 "force_displacement_control[N]", "difference[-]", "negative_pivots_arc_length",
                 "negative_pivots_displacement_control", "slope[N/m]", "inertia_consistent"});
  for (std::size_t i = 0; i < deflection.size(); ++i) {
    Scalar slope = std::nan("");
    std::string consistent;
    if (i > 0 && i + 1 < deflection.size() && !std::isnan(control_force[i - 1]) &&
        !std::isnan(control_force[i + 1])) {
      const Scalar left = (control_force[i] - control_force[i - 1]) / (deflection[i] - deflection[i - 1]);
      const Scalar right =
          (control_force[i + 1] - control_force[i]) / (deflection[i + 1] - deflection[i]);
      slope = (control_force[i + 1] - control_force[i - 1]) / (deflection[i + 1] - deflection[i - 1]);
      const int arc_pivots = path.steps[i].negative_pivots;
      if (arc_pivots > 0) ++unstable_steps;
      if ((left < 0.0) == (right < 0.0) && arc_pivots >= 0 && control_pivots[i] >= 0) {
        ++checked;
        const bool ok = arc_pivots == control_pivots[i] + (slope < 0.0 ? 1 : 0);
        if (!ok) ++mismatches;
        consistent = ok ? "yes" : "NO";
      }
    }
    csv.raw_row({fmt(static_cast<Scalar>(i + 1)), fmt(deflection[i], 10),
                 fmt(path.steps[i].load_factor, 10), fmt(path_force[i], 10),
                 fmt(control_force[i], 10),
                 fmt(std::abs(control_force[i] - path_force[i]) / limit, 4),
                 fmt(static_cast<Scalar>(path.steps[i].negative_pivots)),
                 fmt(static_cast<Scalar>(control_pivots[i])), fmt(slope, 8), consistent});
  }
  csv.close();

  // Load control stops at the limit point and brackets it: the last
  // converged load and the lowest rejected one enclose the limit load.
  const Scalar stop_force = stopped.load_factor * force;
  const Scalar bound_force = stopped.critical_bound * force;
  const bool stopped_at_limit = !stopped.completed && stop_force <= limit + limit_error &&
                                bound_force >= limit - limit_error &&
                                bound_force - stop_force <= 1.0e-3 * limit;

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (two solution methods, one path)"));
  block.set("arch", json::Value::make_string(
                        "clamped circular arch, half span 1 m, rise 0.1 m, thickness 0.02 m, "
                        "width 0.02 m (plane stress), E = 70 GPa, nu = 0.3; half model with a "
                        "symmetry support at the crown, 60 x 4 Q4 cells; crown force 20 kN on "
                        "the half model at lambda = 1 (Saint Venant-Kirchhoff)"));
  block.set("arc_length_steps", json::Value::make_number(static_cast<Scalar>(path.steps.size())));
  block.set("arc_length_iterations", json::Value::make_number(path.total_iterations));
  block.set("arc_length_cuts", json::Value::make_number(path.total_cuts));
  block.set("arc_length_seconds", json::Value::make_number(arc_seconds));
  block.set("displacement_control_steps",
            json::Value::make_number(static_cast<Scalar>(controlled.steps.size())));
  block.set("displacement_control_seconds", json::Value::make_number(control_seconds));
  block.set("highest_sampled_force_N", json::Value::make_number(peak));
  block.set("highest_sample_deflection_m", json::Value::make_number(deflection[peak_step]));
  block.set("limit_force_coarse_parabola_N", json::Value::make_number(coarse_limit));
  block.set("limit_sweep_stations", json::Value::make_number(static_cast<Scalar>(fine)));
  block.set("limit_force_N", json::Value::make_number(limit));
  block.set("limit_force_error_N", json::Value::make_number(limit_error));
  block.set("limit_step", json::Value::make_number(static_cast<Scalar>(peak_step + 1)));
  block.set("limit_deflection_m", json::Value::make_number(sweep_deflection[top]));
  block.set("final_crown_deflection_m", json::Value::make_number(drive));
  block.set("largest_force_difference_over_limit", json::Value::make_number(largest_difference));
  block.set("inertia_points_checked", json::Value::make_number(checked));
  block.set("inertia_mismatches", json::Value::make_number(mismatches));
  block.set("unstable_arc_length_steps", json::Value::make_number(unstable_steps));
  block.set("load_control_completed", json::Value::make_bool(stopped.completed));
  block.set("load_control_stop_force_N", json::Value::make_number(stop_force));
  block.set("load_control_stop_deflection_m",
            json::Value::make_number(stopped.steps.empty()
                                         ? 0.0
                                         : -monitor_value(stopped, stopped.steps.back(),
                                                          "crown_uy")));
  block.set("load_control_bound_force_N", json::Value::make_number(bound_force));
  block.set("load_control_brackets_limit", json::Value::make_bool(stopped_at_limit));
  block.set("load_control_termination", json::Value::make_string(stopped.termination));
  block.set("note",
            json::Value::make_string(
                "The crown force of the arc-length path at each converged step against the "
                "reaction of a displacement-controlled run whose load factors put the crown at "
                "exactly the same deflections (one discrete problem, two solution methods; "
                "the force difference is relative to the limit load). Inertia: the "
                "force-controlled tangent has one more negative pivot than the "
                "displacement-controlled one exactly where the path descends (Haynsworth). "
                "The limit load: displacement control through 60 stations across the two "
                "arc-length intervals around the highest sample, and the vertex of the "
                "parabola through the highest station and its neighbours (its error estimated "
                "by the vertex's distance from that station). Load control must stop at the "
                "limit point, say so, and bracket it: its last converged load and its lowest "
                "rejected one enclose the limit load, 0.1 % apart at most."));
  summary.set("arch_snap_through", block);

  std::ostringstream note;
  note << "limit force " << app::format(limit, 7) << " N at a crown deflection "
       << app::format(sweep_deflection[top], 4) << " m (coarse estimate "
       << app::format(coarse_limit, 6) << " N); " << path.steps.size()
       << " arc-length steps, " << unstable_steps << " of them unstable; inertia consistent at "
       << checked - mismatches << " of " << checked << " points; load control stopped at "
       << app::format(stop_force, 7) << " N, bracketing the limit below "
       << app::format(bound_force, 7) << " N";
  StudyOutcome outcome;
  outcome.name = "shallow-arch snap-through: arc length vs displacement control";
  outcome.kind = "verification";
  outcome.metric = "largest crown-force difference at equal deflection / limit force";
  outcome.value = largest_difference;
  outcome.tolerance = tolerance;
  outcome.passed = all_matched && largest_difference <= tolerance && mismatches == 0 &&
                   checked > 0 && unstable_steps > 0 && stopped_at_limit;
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
