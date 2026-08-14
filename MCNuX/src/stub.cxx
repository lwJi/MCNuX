#include "mcnux_units.hxx"

#include <cctk.h>

// Including mcnux_units.hxx here keeps the compile-time verification of
// specs/conventions-and-units.md [MCNX-CNV-02..06] (its static_asserts) wired
// into every build of the thorn without a dedicated translation unit.

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The full gate additionally asserts
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 once the thorn
// consumes AMReX headers.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
