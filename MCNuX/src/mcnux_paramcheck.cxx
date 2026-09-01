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
  // Trapped-regime guard (specs/trapped-regime-treatment.md [MCNX-TRP-02]):
  // the relabeled scheme's alpha selection reads the unprimed kappa_a off
  // the coefficient interface, and no table-residency layer exists — the
  // same no-live-table-views reason as the emission/interactions guards
  // above.
  if (CCTK_EQUALS(trapped_scheme, "relabeled") &&
      CCTK_EQUALS(opacity_source, "table"))
    CCTK_VERROR(
        "MCNuX::trapped_scheme = \"relabeled\" requires "
        "MCNuX::opacity_source = \"analytic\": no table-residency layer "
        "exists yet, so the relabeling's alpha selection has no live table "
        "views from which to read the unprimed kappa_a ([MCNX-OPA-04]). Use "
        "the analytic coefficient source, or keep the explicit scheme.");
  // Ledger-closure auditability guard ([MCNX-HYD-05]): the synthetic-deposit
  // fixture normalizes with synthfix::dV * synthfix::dt = 0.5 while the real
  // emission contributor normalizes with the grid's cell volume times the
  // run's dt. Both land in the same grid sum, so with both active no single
  // scalar multiplier can reconstruct Sum_cells (source)*dV*dt and the
  // closure is un-auditable. Forbid the mix rather than silently mis-close.
  if (test_ledger_closure && test_synthetic_deposit && enable_emission)
    CCTK_VERROR(
        "MCNuX::test_ledger_closure = yes cannot audit a run with BOTH "
        "MCNuX::test_synthetic_deposit = yes AND MCNuX::enable_emission = "
        "yes: the two contributors deposit with different dV*dt "
        "normalization pairs (the fixture's synthetic constants vs the "
        "grid's cell volume times the run's dt) into the same source-term "
        "sum, so no single grid-side multiplier closes the ledger "
        "([MCNX-HYD-05], specs/hydro-coupling-source-terms.md:177-188). "
        "Disable one of the two contributors, or turn the closure audit "
        "off.");
  // The same auditability guard for the episode driver: its scattering and
  // absorption deposits also normalize with the grid's cell volume times the
  // run's dt, so mixing it with the synthetic fixture's 0.5 pair under the
  // closure is equally un-auditable.
  if (test_ledger_closure && test_synthetic_deposit && enable_interactions)
    CCTK_VERROR(
        "MCNuX::test_ledger_closure = yes cannot audit a run with BOTH "
        "MCNuX::test_synthetic_deposit = yes AND MCNuX::enable_interactions = "
        "yes: the two contributors deposit with different dV*dt "
        "normalization pairs (the fixture's synthetic constants vs the "
        "grid's cell volume times the run's dt) into the same source-term "
        "sum, so no single grid-side multiplier closes the ledger "
        "([MCNX-HYD-05], specs/hydro-coupling-source-terms.md:177-188). "
        "Disable one of the two contributors, or turn the closure audit "
        "off.");

  // Statistical-reduction writer guard ([MCNX-VER-07]): the stats-emission
  // rows reduce the emission step, so without the emission loop there is
  // nothing to reduce and every packet row would abort on its zero-count
  // guard at the first transport step.
  if (test_stats_emission && !enable_emission)
    CCTK_VERROR(
        "MCNuX::test_stats_emission = yes requires MCNuX::enable_emission = "
        "yes: the statistical-reduction writer MCNuX_StatsEmission reduces "
        "the production emission step ([MCNX-PKT-06]/[MCNX-VER-07], the "
        "`stats-emission` benchmark of "
        "specs/verification-suite-design.md:255-279); without the emission "
        "loop there is no step to reduce. Enable emission, or turn the "
        "statistical writer off.");

  // stats-beam guards ([MCNX-VER-07]; the test_stats_emission precedent
  // above and the id-collision precedent below): the beam benchmark reduces
  // the episode driver's absorption/escape resolution, so the driver must
  // run; and its seeder writes hard-coded packet ids 1..beam_num_packets
  // with cpu 0 without touching the AMReX id counter, so any other packet
  // producer/seeder in the same run is an id-collision or duplicate-seeder
  // hazard.
  if (test_stats_beam && !enable_interactions)
    CCTK_VERROR(
        "MCNuX::test_stats_beam = yes requires MCNuX::enable_interactions = "
        "yes: the beam-attenuation writer MCNuX_StatsBeam reduces the episode "
        "driver's pure-absorber resolution of the beam "
        "([MCNX-INT-01]/[MCNX-INT-03], the `stats-beam` benchmark of "
        "specs/verification-suite-design.md:265); without the driver no "
        "packet is ever absorbed or escapes. Enable interactions, or turn "
        "the beam benchmark off.");
  if (test_stats_beam && enable_emission)
    CCTK_VERROR(
        "MCNuX::test_stats_beam = yes is incompatible with "
        "MCNuX::enable_emission = yes: the beam seeder writes hard-coded "
        "packet ids 1..beam_num_packets without touching the AMReX id "
        "counter, so production packets created by the emission loop would "
        "collide with them ([MCNX-GPU-04] id uniqueness), and emitted "
        "packets would contaminate the beam's transmitted/absorbed count "
        "closure.");
  if (test_stats_beam && test_synthetic_packets)
    CCTK_VERROR(
        "MCNuX::test_stats_beam = yes is incompatible with "
        "MCNuX::test_synthetic_packets = yes: both are initial-time packet "
        "seeders writing hard-coded id ranges starting at 1, so their ids "
        "would collide ([MCNX-GPU-04] id uniqueness), and the fixture "
        "packets would contaminate the beam's transmitted/absorbed count "
        "closure.");

  // stats-scatterbox guards (the stats-beam block above, cloned: same
  // driver-must-run reasoning and the same hard-coded-ids-1..N seeder
  // hazard; additionally the beam seeder itself collides — both seeders
  // write ids 1..N with cpu 0 and would double-seed the population).
  if (test_stats_scatterbox && !enable_interactions)
    CCTK_VERROR(
        "MCNuX::test_stats_scatterbox = yes requires "
        "MCNuX::enable_interactions = yes: the collision-statistics writer "
        "MCNuX_StatsScatterbox reduces the episode driver's scattering "
        "resolution of the box ([MCNX-INT-02]/[MCNX-INT-04], the "
        "`stats-scatterbox` benchmark of "
        "specs/verification-suite-design.md:266); without the driver no "
        "packet ever scatters. Enable interactions, or turn the scatterbox "
        "benchmark off.");
  if (test_stats_scatterbox && enable_emission)
    CCTK_VERROR(
        "MCNuX::test_stats_scatterbox = yes is incompatible with "
        "MCNuX::enable_emission = yes: the scatterbox seeder writes "
        "hard-coded packet ids 1..scatterbox_num_packets without touching "
        "the AMReX id counter, so production packets created by the emission "
        "loop would collide with them ([MCNX-GPU-04] id uniqueness), and "
        "emitted packets would contaminate the per-packet Poisson count "
        "statistics.");
  if (test_stats_scatterbox && test_synthetic_packets)
    CCTK_VERROR(
        "MCNuX::test_stats_scatterbox = yes is incompatible with "
        "MCNuX::test_synthetic_packets = yes: both are initial-time packet "
        "seeders writing hard-coded id ranges starting at 1, so their ids "
        "would collide ([MCNX-GPU-04] id uniqueness), and the fixture "
        "packets would contaminate the per-packet Poisson count statistics.");
  if (test_stats_scatterbox && test_stats_beam)
    CCTK_VERROR(
        "MCNuX::test_stats_scatterbox = yes is incompatible with "
        "MCNuX::test_stats_beam = yes: both are initial-time packet seeders "
        "writing hard-coded id ranges starting at 1, so their ids would "
        "collide ([MCNX-GPU-04] id uniqueness), and each fixture would "
        "contaminate the other benchmark's statistics.");

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
