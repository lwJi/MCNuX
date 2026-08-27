#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace MCNuX {

// Single-rate pin of specs/carpetx-thorn-integration.md [MCNX-CTX-06].
//
// MCNuX transport is operator-split and advances every packet on every
// refinement level exactly once per coarsest-level Delta t ([MCNX-CTX-05]).
// With CarpetX::use_subcycling = yes the refined levels advance in per-level
// batches with different Delta t, and neither the transport cadence nor the
// packet redistribution point is specified for that regime — it is an open
// question of the spec, and running unspecified is not permitted. The guard
// below keeps that region unreachable: it aborts at CCTK_PARAMCHECK, before
// any evolution can start.
//
// use_subcycling enters scope unqualified via DECLARE_CCTK_PARAMETERS from
// the SHARES: Driver / USES BOOLEAN block of MCNuX/param.ccl.

extern "C" void MCNuX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (use_subcycling)
    CCTK_VERROR(
        "MCNuX requires CarpetX::use_subcycling = no, but this run sets it to "
        "yes. MCNuX transport is single-rate: it advances all packets on all "
        "refinement levels exactly once per coarsest-level dt, which is "
        "well-defined only when every level shares that timestep. Transport "
        "under subcycling is unspecified (see the open question of "
        "specs/carpetx-thorn-integration.md [MCNX-CTX-06]). Set "
        "CarpetX::use_subcycling = no, or deactivate the MCNuX thorn.");

  // Emission-loop guards ([MCNX-CTX-06] pattern: abort at PARAMCHECK, before
  // any evolution, in exactly the configurations the loop cannot yet serve).
  if (enable_emission && CCTK_EQUALS(opacity_source, "table"))
    CCTK_VERROR(
        "MCNuX::enable_emission = yes requires MCNuX::opacity_source = "
        "\"analytic\": no table-residency layer exists yet, so the emission "
        "loop's RangedTableCoefficients slot has no live table views to "
        "evaluate ([MCNX-OPA-04]). Use the analytic coefficient source, or "
        "disable emission.");
  // Episode-driver guard (same no-table-residency reason as the emission
  // guard above; enable_interactions WITH test_synthetic_packets is allowed
  // — the fixture is the test population of `interactions-fixedseed`).
  if (enable_interactions && CCTK_EQUALS(opacity_source, "table"))
    CCTK_VERROR(
        "MCNuX::enable_interactions = yes requires MCNuX::opacity_source = "
        "\"analytic\": no table-residency layer exists yet, so the episode "
        "driver's RangedTableCoefficients slot has no live table views to "
        "evaluate ([MCNX-OPA-04]). Use the analytic coefficient source, or "
        "disable interactions.");
  if (enable_emission && test_synthetic_packets)
    CCTK_VERROR(
        "MCNuX::enable_emission = yes is incompatible with "
        "MCNuX::test_synthetic_packets = yes: the synthetic fixture seeds "
        "hard-coded packet ids 1..8 without touching the AMReX id counter, so "
        "production packets created by the emission loop would collide with "
        "them ([MCNX-GPU-04] id uniqueness, "
        "specs/particle-container-and-gpu.md).");
}

} // namespace MCNuX
