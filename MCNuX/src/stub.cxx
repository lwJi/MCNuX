#include "mcnux_packets.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

// This translation unit is the aggregation point for the thorn's compile-time
// self-tests: every shared header is included here, so its static_asserts are
// re-verified by every build even if no other TU happens to need it.
// mcnux_units.hxx carries specs/conventions-and-units.md [MCNX-CNV-02..06];
// mcnux_rng.hxx carries specs/rng-and-statistical-acceptance.md
// [MCNX-RNG-01..06] (Philox4x32-10 KAT vectors, the [0, 1) construction, and
// purity); mcnux_packets.hxx carries the layout pin of
// specs/particle-container-and-gpu.md [MCNX-GPU-01..03].

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. These are its language-level legs; the AMReX legs
// (sizeof(amrex::Real) and sizeof(amrex::ParticleReal)) live -- permanently --
// in mcnux_packets.hxx, the header family that consumes AMReX, and are
// re-verified here through the include above.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
