#include "mcnux_coefficients.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_opacity.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_table_coeffs.hxx"
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
// identity endpoint), and mcnux_coefficients.hxx for the source-agnostic
// coefficient interface and its analytic (gray) source of
// specs/verification-suite-design.md [MCNX-VER-05] (exact kappa_a/kappa_s
// pass-through, source-enumerator distinctness, device-capturable parameter
// layout; the hc pin itself is asserted in mcnux_units.hxx), and
// mcnux_table_coeffs.hxx for the baseline table-source coefficient assembly
// and nu_x mapping of specs/opacity-eos-evaluation.md [MCNX-OPA-04/05]
// (TableEval-contract conformance, trivially-copyable view bundle/functor,
// nu_x has-no-slot special-case pin, lepton-number mu_nu sign map), and
// mcnux_geodesic.hxx for the flat-spacetime limit of the geodesic push of
// specs/geodesic-propagation.md [MCNX-GEO-01/04] (dp_i/dt == 0 and |dx/dt|
// == 1 exactly on the flat snapshot, one RK4 step an exact straight line
// with bitwise-constant momentum, trivially-copyable gather functor).

// Precision gate of specs/build-and-integration.md [MCNX-BLD-03]: all reals
// are IEEE-754 binary64. The AMReX legs of the gate —
// sizeof(amrex::Real) == sizeof(amrex::ParticleReal) == 8 — live in
// mcnux_particles.hxx (headers carry their own selftests); the plain
// double/CCTK_REAL legs stay here.
static_assert(sizeof(double) == 8, "MCNuX requires IEEE-754 binary64 doubles");
static_assert(sizeof(CCTK_REAL) == 8, "MCNuX requires CCTK_REAL == binary64");
