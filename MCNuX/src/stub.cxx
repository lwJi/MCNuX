#include "mcnux_opacity.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

// Including mcnux_units.hxx here keeps the compile-time verification of
// specs/conventions-and-units.md [MCNX-CNV-02..06] (its static_asserts) wired
// into every build of the thorn without a dedicated translation unit; likewise
// mcnux_rng.hxx for specs/rng-and-statistical-acceptance.md [MCNX-RNG-01..06]
// (Philox4x32-10 KAT vectors, the [0, 1) construction, and purity),
// mcnux_particles.hxx for the pure-SoA layout pin of
// specs/particle-container-and-gpu.md [MCNX-GPU-01..03], and
// mcnux_opacity.hxx for the WeakLibInterp call-boundary contract of
// specs/opacity-eos-evaluation.md [MCNX-OPA-01/02/03/08] (type boundary,
// species-axis grounding, Brem Alpha, AxisBounds copyability), and
// mcnux_tetrad.hxx for the tetrad/null-closure helpers of
// specs/geodesic-propagation.md [MCNX-GEO-02] and
// specs/packet-representation-and-sampling.md [MCNX-PKT-04] (analytic 3x3
// inverse, flat-static tetrad identity, flat-limit transform/closure), and
// mcnux_srcterms.hxx for the source-term ledger arithmetic of
// specs/hydro-coupling-source-terms.md [MCNX-HYD-02] (sign fixtures, exact
// fixture values, nu_x zero-lepton, no-degeneracy-in-ledger).

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The AMReX legs of the gate —
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 — live in
// mcnux_particles.hxx (headers carry their own selftests); the plain
// double/CCTK_REAL legs stay here.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
