#include "mcnux_coefficients.hxx"
#include "mcnux_opacity.hxx"
#include "mcnux_opacity_analytic.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_trp.hxx"
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
// fixture values, nu_x zero-lepton, no-degeneracy-in-ledger), and
// mcnux_trp.hxx for the trapped-regime opacity-relabeling pure map of
// specs/trapped-regime-treatment.md [MCNX-TRP-02/03] (the two exact
// invariants, the 2021-convention monotonicity tripwire, the alpha = 1
// identity endpoint), and
// mcnux_opacity_analytic.hxx for the analytic (gray) opacity mode of
// specs/verification-suite-design.md [MCNX-VER-05] feeding the coefficient
// interface of specs/opacity-eos-evaluation.md:107-126 (per-species
// pass-through of kappa_a0/kappa_s0, eta_scale = 0 -> eta = 0, the
// Kirchhoff J_eq identity, E^3 scaling of 4 pi E^3/(hc)^3, nu_x emission),
// and
// mcnux_coefficients.hxx for the source-agnostic runtime dispatch of the same
// interface (specs/opacity-eos-evaluation.md:118-126, [MCNX-VER-05]:
// pinned OpacitySource values, the param.ccl keyword map, CoefficientSource
// copyability).

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The AMReX legs of the gate —
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 — live in
// mcnux_particles.hxx (headers carry their own selftests); the plain
// double/CCTK_REAL legs stay here.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
