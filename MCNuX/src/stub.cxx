#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

// Including mcnux_units.hxx here keeps the compile-time verification of
// specs/conventions-and-units.md [MCNX-CNV-02..06] (its static_asserts) wired
// into every build of the thorn without a dedicated translation unit; likewise
// mcnux_rng.hxx for specs/rng-and-statistical-acceptance.md [MCNX-RNG-01..06]
// (Philox4x32-10 KAT vectors, the [0, 1) construction, and purity), and
// mcnux_particles.hxx for the pure-SoA layout pin of
// specs/particle-container-and-gpu.md [MCNX-GPU-01..03].

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The AMReX legs of the gate —
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 — live in
// mcnux_particles.hxx (headers carry their own selftests); the plain
// double/CCTK_REAL legs stay here.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
