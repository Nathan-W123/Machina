/// \file Buckling.hpp
/// \brief Linear (eigenvalue) buckling analysis of the constrained model.
///
/// A load case \f$f\f$ produces the displacement \f$u\f$ of \f$K u = f\f$ and
/// with it a stress field \f$\sigma = D B u\f$. Linearised buckling asks for
/// which multiple \f$\lambda\f$ of that stress state the stiffness loses
/// definiteness:
/// \f[
///   \left(K_{ff} + \lambda\, K_{G,ff}(u)\right)\phi = 0,
/// \f]
/// where \f$K_G\f$ is the geometric (initial-stress) stiffness
/// (`Element::geometric_stiffness`). The structure is predicted to buckle at
/// the load \f$\lambda f\f$: \f$\lambda > 1\f$ is a safety factor on the
/// reference load, \f$\lambda < 1\f$ means it buckles before reaching it. A
/// compressive stress makes \f$K_G\f$ negative in the direction of the
/// buckling mode; a tensile one stiffens, and a load case that only
/// stiffens the structure has no positive \f$\lambda\f$ at all (the
/// reversed load may still buckle it, at a negative \f$\lambda\f$).
///
/// **Assumptions.** Linear elasticity up to buckling, small pre-buckling
/// rotations (the pre-buckling state is the linear solution, scaled), loads
/// that keep their direction (no follower forces), and a perfect geometry:
/// \f$\lambda\f$ is the bifurcation load of the ideal structure and is an
/// upper bound on the collapse load of a real, imperfect one. It is a
/// stability *check*, not a post-buckling analysis.
///
/// **Algorithm.** With the positive definite \f$K\f$ on the right-hand side,
/// \f$\mu = 1/\lambda\f$ solves the symmetric-definite pencil
/// \f$(-K_G)\phi = \mu K \phi\f$, and the smallest positive load factors are
/// its largest eigenvalues. Subspace iteration on
/// \f$X \leftarrow K^{-1}(-K_G)X\f$ with a Rayleigh-Ritz projection onto the
/// \f$q = \min(n, \max(2m, m+8))\f$-dimensional subspace finds them with
/// the factorisation of \f$K_{ff}\f$ the static solve already made. But
/// \f$-K_G\f$ is indefinite: negative \f$\mu\f$ - load factors of the
/// reversed load - converge alongside, and when more of them than the spare
/// \f$q - m - 2\f$ slots exceed the wanted eigenvalues in magnitude (a
/// structure mostly in tension, or near-void SIMP material in tension) the
/// wanted ones cannot converge. The iteration then switches to the buckling
/// spectral transformation (Grimes, Lewis and Simon 1994),
/// \f$X \leftarrow (K + \sigma K_G)^{-1} K X\f$ with eigenvalues
/// \f$\nu = \lambda / (\lambda - \sigma)\f$: above 1 for the load factors
/// beyond \f$\sigma\f$ and in (0, 1) for every negative one, whatever its
/// size. For \f$0 < \sigma < \lambda_1\f$ the shifted matrix is positive
/// definite, and by Sylvester's law of inertia the number of negative pivots
/// of its \f$LDL^T\f$ factorisation equals the number of load factors in
/// \f$(0, \sigma)\f$; that count places \f$\sigma\f$ below \f$\lambda_1\f$
/// (lowered by halving from an estimate, or grown by factors of 8 from the
/// load case itself), and \f$\sigma\f$ is raised once to 0.8 of the settled
/// \f$\lambda_1\f$. A load case with no positive load factor below
/// \f$10^{15}\f$ is reported as only stiffening. The shifted factorisation
/// is a sparse \f$LDL^T\f$ whatever `linear` selects. Convergence requires
/// both the relative change of the requested load factors and every
/// eigenpair residual \f$\|K\phi + \lambda K_G\phi\| / \|K\phi\|\f$ below
/// their tolerances.
/// Systems with at most 400 free DOFs use a dense generalised eigensolve.
/// Modes are normalised to \f$\phi^T K \phi = 1\f$, the normalisation of the
/// sensitivity formula in BucklingConstraint.hpp, and returned full length
/// with zeros at prescribed DOFs.
///
/// **Density fields.** For a SIMP design `stiffness_scale` holds
/// \f$E_K(\rho)/E_0\f$ and `stress_scale` the factor of the stress that
/// enters \f$K_G\f$; with the same factor in both, void regions carry
/// artificial low-stiffness "pseudo" buckling modes, which is why the
/// topology optimiser passes \f$\rho^p\f$ without the \f$E_{min}\f$ floor
/// (see BucklingConstraint.hpp). The fraction of each mode's strain energy
/// in solid elements is reported so a void-localised mode is visible.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"

#include <functional>
#include <string>
#include <vector>

namespace sparlab {

/// When the subspace iteration uses the buckling spectral transformation.
enum class BucklingTransform {
  Auto,        ///< plain iteration, transformed once reversed-load modes crowd it
  ShiftInvert  ///< transformed from the start
};

struct BucklingOptions {
  int num_modes = 4;                 ///< number of smallest positive load factors
  int max_iterations = 400;          ///< subspace iteration cap
  Scalar tolerance = 1.0e-8;         ///< relative change of the requested load factors
  Scalar residual_tolerance = 1.0e-6;  ///< largest accepted eigenpair residual
  /// Fixed seed for the random starting subspace (bit-for-bit reproducible).
  unsigned int seed = 20240917u;
  /// Solver for K_ff when the caller supplies no factorisation.
  LinearSolverOptions linear;
  BucklingTransform transform = BucklingTransform::Auto;
  /// A positive estimate of the smallest load factor (the previous design's,
  /// in an optimisation) from which the transformation's shift starts; 0
  /// searches from the Ritz values or from 1.
  Scalar shift_hint = 0.0;
};

struct BucklingResult {
  std::string load_case;
  /// Smallest positive load factors, ascending [-]. Fewer than requested when
  /// the pencil has fewer positive eigenvalues; empty when the load case
  /// does not destabilise the structure at all.
  Vector load_factors;
  /// num_dofs x modes, \f$\phi^T K \phi = 1\f$, zero at prescribed DOFs [m].
  Matrix mode_shapes;
  Vector residuals;                  ///< per mode, \f$\|K\phi+\lambda K_G\phi\|/\|K\phi\|\f$
  /// Per mode, the fraction of \f$\phi^T K \phi\f$ in elements whose density
  /// is at or above 0.5 (1 for a plain analysis): a pseudo mode localised in
  /// void material scores near 0.
  Vector solid_energy_fraction;
  bool no_positive_load_factor = false;  ///< the load case only stiffens
  int iterations = 0;
  int subspace_size = 0;
  bool converged = false;
  bool transformed = false;          ///< the spectral transformation was used
  Scalar sigma = 0.0;                ///< its shift, a lower bound on lambda_1
  int factorizations = 0;            ///< LDL^T factorisations of K + sigma K_G
  Scalar final_change = 0.0;
  std::string linear_solver;
  Index linear_iterations = 0;
  std::vector<std::string> warnings;
  /// Reduced (free-DOF) modes of the whole converged subspace, the warm start
  /// of the next analysis of a slightly changed design.
  Matrix subspace;
};

/// Solves \f$K_{ff} x = b\f$ on the free DOFs (reduced vectors) with a
/// factorisation the caller already holds.
using FreeSolve = std::function<Vector(const Vector&)>;

/// Global geometric stiffness \f$K_G(u) = \sum_e s_e A_e^T K_{G,e}(u_e) A_e\f$
/// of the stress state of the full-length displacement `displacement`.
/// \param stress_scale optional per-element factors \f$s_e\f$ on the stress.
/// \param temperature optional nodal temperatures [K]: the prestress is then
///        \f$D(Bu - \varepsilon_0)\f$, so a heated, restrained structure can
///        buckle under its own thermal stress.
SparseMatrix assemble_geometric_stiffness(const FemModel& model, const Assembler& assembler,
                                          const Vector& displacement,
                                          const Vector* stress_scale = nullptr,
                                          const Vector* temperature = nullptr);

/// Solve the buckling eigenproblem of the stress state of `displacement`.
/// \param k_full the assembled global stiffness the displacement was solved
///        with (full size); its free-free block is the left matrix.
/// \param k_g_full the matching geometric stiffness (full size).
/// \param solve optional solver of \f$K_{ff}\f$ already factorised; when
///        empty a solver of `options.linear` is built and factorised here.
/// \param initial_subspace optional free-DOF vectors (n x k) placed first in
///        the starting subspace, typically `subspace` of a previous result.
/// \param stiffness_scale the per-element factors `k_full` was assembled
///        with, used to split each mode's strain energy over the elements.
/// \param solid_mask optional per-element flags (1 = solid) for
///        `solid_energy_fraction`; all elements count as solid when absent.
/// \throws ModelError for an ill-posed model, SolverError for a failed
///         eigensolve, ConvergenceError when the iteration budget is spent.
BucklingResult solve_buckling(const FemModel& model, const Assembler& assembler,
                              const SparseMatrix& k_full, const SparseMatrix& k_g_full,
                              const BucklingOptions& options, const FreeSolve& solve = {},
                              const Matrix* initial_subspace = nullptr,
                              const Vector* stiffness_scale = nullptr,
                              const std::vector<char>* solid_mask = nullptr);

/// Convenience driver for a plain (non-SIMP) analysis: assemble \f$K\f$,
/// solve the load case `load_case` of the model, assemble \f$K_G\f$ and
/// solve the buckling problem.
BucklingResult analyse_buckling(const FemModel& model, const Assembler& assembler,
                                std::size_t load_case, const BucklingOptions& options);

/// Euler's critical load of a column clamped at one end and free at the
/// other, \f$P = \pi^2 E I / (4 L^2)\f$ [N], and the Engesser shear
/// correction \f$1/P_s = 1/P + 1/(\kappa G A)\f$ with the rectangular
/// section's \f$\kappa = 5/6\f$, which a continuum model includes.
/// \{
Scalar euler_cantilever_load(Scalar e, Scalar second_moment, Scalar length);
Scalar engesser_cantilever_load(Scalar e, Scalar poisson, Scalar second_moment,
                                Scalar area, Scalar length);
/// \}

}  // namespace sparlab
