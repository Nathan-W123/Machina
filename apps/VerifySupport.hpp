/// \file VerifySupport.hpp
/// \brief What the studies of sparlab_verify share across its source files.
///
/// The structural studies live in sparlab_verify.cpp; the studies of the
/// volume, pressure and thermal loads live in verify_loads.cpp, those of the
/// geometrically non-linear analysis in verify_nonlinear.cpp, those of
/// plasticity in verify_plasticity.cpp and those of dynamics in
/// verify_dynamics.cpp, those of the incompatible-mode Hex8 in
/// verify_elements.cpp. All report a
/// `StudyOutcome` that the driver prints and writes to summary.json.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/Selector.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <cmath>
#include <string>

namespace sparlab {
namespace verify {

/// The verdict of one study.
struct StudyOutcome {
  std::string name;
  bool passed = true;
  std::string metric;
  Scalar value = 0.0;
  Scalar tolerance = 0.0;
  std::string kind;  ///< "verification" or "validation"
  std::string note;
};

/// Observed convergence order between two refinements of a quantity whose
/// error is e1 at mesh size h1 and e2 at h2: p = log(e1/e2) / log(h1/h2).
/// Zero when either error is not positive or the sizes coincide.
inline Scalar observed_order(Scalar h1, Scalar e1, Scalar h2, Scalar e2) {
  if (!(e1 > 0.0) || !(e2 > 0.0) || h1 == h2) return 0.0;
  return std::log(e1 / e2) / std::log(h1 / h2);
}

/// Studies of the volume, pressure and thermal loads (verify_loads.cpp).
/// \{
StudyOutcome study_lame_cylinder(const std::string& out_dir, json::Value& summary);
StudyOutcome study_rotating_disk(const std::string& out_dir, json::Value& summary);
StudyOutcome study_thermal_cylinder(const std::string& out_dir, json::Value& summary);
StudyOutcome study_bimetal_strip(const std::string& out_dir, json::Value& summary);
StudyOutcome study_self_weight(const std::string& out_dir, json::Value& summary);
/// \}

/// Studies of the geometrically non-linear analysis (verify_nonlinear.cpp).
/// \{
StudyOutcome study_elastica(const std::string& out_dir, json::Value& summary);
StudyOutcome study_hyperelastic_cylinder(const std::string& out_dir, json::Value& summary);
StudyOutcome study_arch_snap_through(const std::string& out_dir, json::Value& summary);
/// \}

/// Studies of the elastoplastic analysis (verify_plasticity.cpp).
/// \{
StudyOutcome study_plastic_cylinder(const std::string& out_dir, json::Value& summary);
StudyOutcome study_plastic_bending(const std::string& out_dir, json::Value& summary);
StudyOutcome study_plastic_cycle(const std::string& out_dir, json::Value& summary);
StudyOutcome study_chaboche_cycle(const std::string& out_dir, json::Value& summary);
StudyOutcome study_hill_directional(const std::string& out_dir, json::Value& summary);
/// \}

/// Studies of the logarithmic-strain finite plasticity (verify_logarithmic.cpp).
/// \{
StudyOutcome study_logarithmic_uniaxial(const std::string& out_dir, json::Value& summary);
StudyOutcome study_logarithmic_tube(const std::string& out_dir, json::Value& summary);
/// \}

/// Studies of the transient and harmonic analyses (verify_dynamics.cpp).
/// \{
StudyOutcome study_transient_modal(const std::string& out_dir, json::Value& summary);
StudyOutcome study_rod_harmonic(const std::string& out_dir, json::Value& summary);
StudyOutcome study_rod_transient(const std::string& out_dir, json::Value& summary);
StudyOutcome study_nonlinear_oscillator(const std::string& out_dir, json::Value& summary);
/// \}

/// Studies of the incompatible-mode Hex8 and the through-thickness rule in
/// bending (verify_elements.cpp).
/// \{
StudyOutcome study_im_cantilever(const std::string& out_dir, json::Value& summary);
StudyOutcome study_im_macneal_harder(const std::string& out_dir, json::Value& summary);
StudyOutcome study_im_plate(const std::string& out_dir, json::Value& summary);
StudyOutcome study_im_springback(const std::string& out_dir, json::Value& summary);
/// \}

/// Quarter sections of a cylinder (verify_loads.cpp).
/// \{
/// The unit box [0, 1]^2 of a structured mesh (times [0, depth] in 3-D) mapped
/// onto the annular sector a <= r <= b, 0 <= theta <= 90 deg by
/// r = a + (b - a) x, theta = 90 deg y. Every node moves, a Tet10's edge nodes
/// included, so a Tet10 cell follows the circles with curved faces while a
/// linear cell spans their chords. The mapped mesh carries no structured-grid
/// information: its cells differ.
Mesh quarter_annulus(const Mesh& unit, Scalar a, Scalar b);
/// Quarter section with n_r cells through the wall and n_theta round it. A
/// solid section is one layer of cells as deep as a cell is wide.
Mesh sector_mesh(ElementType type, Index nr, Index nt, Scalar a, Scalar b);
/// Symmetry supports of a quarter section, plus u_z = 0 at every node of a
/// solid one (plane strain).
void add_quarter_supports(FemModel& model);
/// The cylindrical surface r = radius: the bore (every node with r <= radius)
/// or the outer surface (every node with r >= radius).
SelectorGroup cylinder_surface(const std::string& name, Scalar radius, bool bore);
/// \}

}  // namespace verify
}  // namespace sparlab
