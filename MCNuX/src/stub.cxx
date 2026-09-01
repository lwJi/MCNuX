#include "mcnux_coefficients.hxx"
#include "mcnux_deposit.hxx"
#include "mcnux_emission.hxx"
#include "mcnux_escape.hxx"
#include "mcnux_fluid.hxx"
#include "mcnux_gather.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_interactions.hxx"
#include "mcnux_opacity.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_selftest.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_stats.hxx"
#include "mcnux_stats_beam.hxx"
#include "mcnux_stats_emission.hxx"
#include "mcnux_stats_scatterbox.hxx"
#include "mcnux_stats_writer.hxx"
#include "mcnux_table_coeffs.hxx"
#include "mcnux_table_range.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_trp.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

// Compile-time selftest aggregation point: every shared mcnux_*.hxx header is
// #included here so its static_asserts run on every build, independent of
// which production TUs happen to include it. Each header documents the spec
// requirements its own asserts pin.

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The AMReX legs of the gate —
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 — live in
// mcnux_particles.hxx (headers carry their own selftests); the plain
// double/CCTK_REAL legs stay here.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
