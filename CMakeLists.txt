cmake_minimum_required(VERSION 3.16)

project(SparLab
        VERSION 1.0.0
        DESCRIPTION "2-D/3-D finite-element structural solver with SIMP topology optimization"
        LANGUAGES CXX)

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
option(SPARLAB_BUILD_TESTS     "Build the Catch2 test suite"            ON)
option(SPARLAB_BUILD_APPS      "Build the command-line applications"    ON)
option(SPARLAB_ENABLE_WARNINGS "Enable a strict warning set"            ON)
option(SPARLAB_WARNINGS_AS_ERRORS "Treat compiler warnings as errors"   OFF)
option(SPARLAB_ENABLE_OPENMP   "Use OpenMP where the code supports it"  ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  message(STATUS "SparLab: CMAKE_BUILD_TYPE not set, defaulting to Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "SparLab: found Eigen ${Eigen3_VERSION}")

if(SPARLAB_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    message(STATUS "SparLab: OpenMP enabled (${OpenMP_CXX_VERSION})")
  else()
    message(STATUS "SparLab: OpenMP not found, building single-threaded")
  endif()
endif()

# ----------------------------------------------------------------------------
# Core library
# ----------------------------------------------------------------------------
add_library(sparlab STATIC
  src/core/Logging.cpp
  src/mesh/Mesh.cpp
  src/mesh/StructuredMesh.cpp
  src/mesh/SubMesh.cpp
  src/material/IsotropicMaterial.cpp
  src/material/Hyperelastic.cpp
  src/elements/Quadrature.cpp
  src/elements/Element.cpp
  src/elements/Quad4.cpp
  src/elements/Hex8.cpp
  src/elements/Tri3.cpp
  src/elements/Tet4.cpp
  src/elements/Tet10.cpp
  src/elements/FaceGeometry.cpp
  src/fem/DofManager.cpp
  src/fem/Selector.cpp
  src/fem/BoundaryConditions.cpp
  src/fem/FemModel.cpp
  src/fem/Assembler.cpp
  src/fem/LinearSolver.cpp
  src/fem/Multigrid.cpp
  src/fem/ModelDiagnostics.cpp
  src/fem/StaticAnalysis.cpp
  src/fem/StressRecovery.cpp
  src/fem/ModalAnalysis.cpp
  src/fem/Buckling.cpp
  src/fem/Loads.cpp
  src/fem/HeatConduction.cpp
  src/fem/TotalLagrangian.cpp
  src/fem/NonlinearStatic.cpp
  src/topopt/SimpInterpolation.cpp
  src/topopt/DensityFilter.cpp
  src/topopt/DesignDomain.cpp
  src/topopt/Sensitivity.cpp
  src/topopt/OptimalityCriteria.cpp
  src/topopt/Mma.cpp
  src/topopt/Projection.cpp
  src/topopt/StressConstraint.cpp
  src/topopt/BucklingConstraint.cpp
  src/topopt/OverhangFilter.cpp
  src/topopt/LengthScale.cpp
  src/topopt/TopologyOptimizer.cpp
  src/io/Json.cpp
  src/io/Config.cpp
  src/io/CsvWriter.cpp
  src/io/VtkWriter.cpp
  src/io/StlWriter.cpp
  src/io/CalculixWriter.cpp
  src/io/MeshReader.cpp
  src/io/ResultWriter.cpp
)

target_include_directories(sparlab PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(sparlab PUBLIC Eigen3::Eigen)

if(SPARLAB_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(sparlab PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(sparlab PUBLIC SPARLAB_HAVE_OPENMP=1)
endif()

target_compile_definitions(sparlab PUBLIC
  SPARLAB_VERSION_STRING="${PROJECT_VERSION}"
  SPARLAB_BUILD_TYPE="$<CONFIG>"
)

# Eigen's own asserts are valuable during development but very slow; keep them
# in Debug only.
target_compile_definitions(sparlab PUBLIC
  $<$<NOT:$<CONFIG:Debug>>:EIGEN_NO_DEBUG>
)

if(SPARLAB_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(sparlab PRIVATE /W4 /permissive-)
  else()
    target_compile_options(sparlab PRIVATE
      -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Wcast-align
      -Wunused -Woverloaded-virtual -Wdouble-promotion -Wformat=2)
  endif()
  if(SPARLAB_WARNINGS_AS_ERRORS)
    if(MSVC)
      target_compile_options(sparlab PRIVATE /WX)
    else()
      target_compile_options(sparlab PRIVATE -Werror)
    endif()
  endif()
endif()

# ----------------------------------------------------------------------------
# Applications
# ----------------------------------------------------------------------------
if(SPARLAB_BUILD_APPS)
  add_library(sparlab_app_support STATIC apps/AppSupport.cpp)
  target_link_libraries(sparlab_app_support PUBLIC sparlab)
  target_include_directories(sparlab_app_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/apps)

  foreach(app sparlab_solve sparlab_topopt sparlab_verify sparlab_bench)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE sparlab_app_support)
  endforeach()
  target_sources(sparlab_verify PRIVATE apps/verify_loads.cpp apps/verify_nonlinear.cpp)
endif()

# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------
if(SPARLAB_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
message(STATUS "SparLab configuration summary")
message(STATUS "  build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ std    : ${CMAKE_CXX_STANDARD}")
message(STATUS "  tests      : ${SPARLAB_BUILD_TESTS}")
message(STATUS "  apps       : ${SPARLAB_BUILD_APPS}")
