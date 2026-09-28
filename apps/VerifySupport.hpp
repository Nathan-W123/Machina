/// \file VerifySupport.hpp
/// \brief What the studies of sparlab_verify share across its source files.
///
/// The structural studies live in sparlab_verify.cpp; the studies of the
/// volume, pressure and thermal loads live in verify_loads.cpp. Both report a
/// `StudyOutcome` that the driver prints and writes to summary.json.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/io/Json.hpp"

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

}  // namespace verify
}  // namespace sparlab
