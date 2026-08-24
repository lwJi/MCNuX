// The MCNuX unit self-test battery — the runtime half of the `unit-selftest`
// benchmark of specs/verification-suite-design.md [MCNX-VER-01],
// [MCNX-VER-06].
//
// The battery re-runs, at runtime, the checks that the shared headers already
// assert at compile time, and writes them out as a table of rows so that the
// Cactus regression harness can diff them against committed golden data
// (MCNuX/test/unit-selftest/, produced by CarpetX's TSV writer at
// ABSTOL=RELTOL=1e-12). The spec's binding shape is "one row per named check
// (measured error, pass flag)"; here a row is one index of the 1d grid array
// MCNuX::mcnux_selftest, carrying `selftest_error` and `selftest_pass`.
//
// Why both compile-time and runtime? The static_asserts prove the predicates
// hold for the *compiler's* constant evaluator; these rows additionally pin
// the *measured* errors as golden numbers, so that a change which keeps a
// check passing but moves its error (a re-derived factor, a different
// evaluation order) is still caught at 1e-12.
//
// Layout: this file owns the authoritative index -> name row table below and
// the ordered appender call sequence in run_battery() (that sequence IS the
// row ordering); the check code itself lives in one per-domain
// mcnux_selftest_<domain>.cxx per appender, with the shared Battery machinery
// and the appender declarations in mcnux_selftest.hxx.
//
// ---------------------------------------------------------------------------
// Row table (index -> name).  APPEND-ONLY.
// ---------------------------------------------------------------------------
// Golden data is keyed by row index, so rows are never reordered, renamed in
// meaning, or removed: later tasks append their checks at the end and
// regenerate the golden file (and bump SIZE in interface.ccl). To append,
// add a new append_*_rows function — in the matching existing
// mcnux_selftest_<domain>.cxx, or a new one (added to make.code.defn) if the
// domain is new — declare it in mcnux_selftest.hxx, and call it at the END
// of run_battery() below; never insert rows into an existing appender. The
// full battery the spec eventually requires (tetrad/null identities, Tmunu
// closed forms, nu_x mapping, relabeling identities, sign/range/inversion
// fixtures, sizeof(amrex::Real)) is built up this way.
//
//   0..2    rng.philox.kat{1,2,3}                 [MCNX-RNG-01], [MCNX-RNG-02]
//   3..6    rng.uniform.{u64_zero, u64_2p11_minus_1,
//                        u64_2p11, u64_max}       [MCNX-RNG-04]
//   7       rng.purity                            [MCNX-RNG-05], [MCNX-RNG-06]
//   8       rng.pairing                           [MCNX-RNG-03]
//   9..20   units.factor.<name>                   [MCNX-CNV-03], [MCNX-CNV-04]
//  21..32   units.roundtrip.<name>                [MCNX-CNV-02]
//  33       units.species_table                   [MCNX-CNV-06]
//  34..35   precision.sizeof_{double, cctk_real}  [MCNX-BLD-03]
//  36..37   precision.sizeof_{amrex_real,
//                             amrex_particle_real} [MCNX-BLD-03]
//  38..46   opacity.parity.{eos_evaluate,
//                eos_invert_dey_noguess, eos_invert_dey_guess,
//                emab, iso_imom1, nes_pair, nes_detailed_balance,
//                pair_crossing, brem_summed}      [MCNX-OPA-01..03]
//  47..49   opacity.units.{log10_kelvin_once,
//                nes_pair_kelvin_axis,
//                nes_db_dual_temperature}         [MCNX-OPA-02]
//  50..51   opacity.species.{emab,iso}_slots      [MCNX-OPA-08]
//  52       opacity.bounds.endpoints              (T11 range-policy hook)
//  53       geo.inverse_metric.identity           [MCNX-GEO-02]
//  54       geo.null_closure                      [MCNX-GEO-02]
//  55..57   tetrad.{orthonormality, null_norm,
//                   nu_recovery}                  [MCNX-PKT-04]
//  58..59   tetrad.flat_static.{momentum,
//                               identity}         [MCNX-PKT-04]
//  60..61   srcterm.emission.{lepton,energy}_sign [MCNX-HYD-02]
//  62..63   srcterm.absorption.{lepton,
//                               energy}_sign      [MCNX-HYD-02]
//  64..65   srcterm.scattering.{lepton_zero,
//                               momentum_only}    [MCNX-HYD-02]
//  66..67   trp.relabel.{eta_over_kappa_a,
//                        total_rate}              [MCNX-TRP-02]
//  68..69   trp.relabel.monotonicity_{ka, ks}     [MCNX-TRP-03]
//  70       trp.relabel.alpha_one_identity        [MCNX-TRP-02] (alpha = 1)
//  71       units.hc                              hc pin (opacity-eos-
//                                                 evaluation.md:356-361)
//  72..75   coef.analytic.{kappa_passthrough,
//                          kirchhoff_identity,
//                          eta_zero, nux_emission} [MCNX-VER-05]
//  76..78   coef.dispatch.{analytic, table,
//                          parameter}             [MCNX-VER-05] (source
//                                                 switch invisible below the
//                                                 coefficient interface)
//  79..81   geo.gather.{linear_values,
//                      linear_derivs,
//                      constant_derivs_zero}      [MCNX-GEO-03]
//  82       geo.rhs.null_closure_identity         [MCNX-GEO-02]
//  83..84   geo.push.flat_{straight_line,
//                         unit_speed}             [MCNX-GEO-04] (flat limit)
//  85       geo.push.order_ge2                    [MCNX-GEO-04] (convergence)
//  86..87   tblcoef.kirchhoff.{nue, nuebar}       [MCNX-OPA-04] (Kirchhoff
//                                                 closure, synthetic tables)
//  88..89   tblcoef.nux.{kappa_a_zero, eta_zero}  [MCNX-OPA-05] (exact)
//  90       tblcoef.nux.kappa_s_halfmean          [MCNX-OPA-05]
//  91       tblcoef.dispatch.table_plumb          [MCNX-OPA-04] (plugs into
//                                                 the [MCNX-VER-05] dispatch)
//  92..93   tblrange.floor.{zero_no_touch, nux}   [MCNX-OPA-06] (transparency
//                                                 floor; NaN-poisoned table
//                                                 data prove no table touch)
//  94..96   tblrange.clamp.{equality, counters,
//                           inrange_noop}         [MCNX-OPA-06] (component-
//                                                 wise clamping + per-axis
//                                                 counters, exact)
//  97       tblrange.nonan                        [MCNX-OPA-06] (no-NaN
//                                                 guarantee, hostile probe)
//  98       tblrange.rho_min                      [MCNX-OPA-06] (floor
//                                                 density = max of the
//                                                 tabulated lower bounds)
//  99       tblrange.inversion_protocol           [MCNX-OPA-07] (pinned
//                                                 failing inversion input)
// 100       stats.z.exact                         [MCNX-VER-07] (z closed
//                                                 form, exact fixture)
// 101       stats.z.sign                          [MCNX-VER-07] (negative z,
//                                                 exact, in band)
// 102       stats.z.zero                          [MCNX-VER-07] (estimate ==
//                                                 expected -> z == 0, pass)
// 103       stats.z.boundary_pass                 [MCNX-RNG-07] (z exactly 4
//                                                 passes: non-strict <=)
// 104       stats.z.boundary_fail                 [MCNX-RNG-07] (one 2^-20
//                                                 step beyond 4 sigma fails)
// 105       stats.constants                       [MCNX-RNG-07] (pinned N_p,
//                                                 seeds, 4 sigma bound)
//
// Tiers, per specs/README.md and the tolerance discussion in mcnux_units.hxx:
//   * exact   — pass requires a measured error of identically 0 (KATs, the
//               [0, 1) closed forms, purity/pairing, the species table, the
//               precision rows).
//   * machine — the round-trip rows; MCNuX::detail::round_trips judges at
//               detail::rtol_machine = 1e-14 internally, so the row is
//               boolean and its measured error is 0 or 1.
//   * pinned  — the factor-recomputation rows, judged at
//               detail::rtol_pinned = 1e-9. Their measured errors are ~5e-10
//               by construction: the spec's table values carry 10 significant
//               digits, so agreement cannot be asserted at the machine tier
//               (see the rtol_pinned comment in mcnux_units.hxx).

#include "mcnux_selftest.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>

namespace MCNuX {

namespace {

// The battery itself. The ordered calls below ARE the row ordering (each
// appender's own row range is documented in its mcnux_selftest_<domain>.cxx);
// new appenders go at the END, per the append-only contract above.
void run_battery(Battery &b) {
  // Rows 0..37 — the RNG contract of
  // specs/rng-and-statistical-acceptance.md [MCNX-RNG-01..06], the unit
  // factors/round-trips/species table of specs/conventions-and-units.md
  // [MCNX-CNV-02..04, -06], and the binary64 precision gate of
  // specs/build-and-integration.md [MCNX-BLD-03].
  append_rng_units_rows(b);

  // Rows 38..52 — the WeakLibInterp call-boundary wrapper checks of
  // specs/opacity-eos-evaluation.md (kernel parity on synthetic tables,
  // unit-boundary conversions, species-axis grounding, T11 bounds hooks).
  append_opacity_rows(b);

  // Rows 53..59 — the tetrad/null-closure helper checks of
  // specs/geodesic-propagation.md [MCNX-GEO-02] and
  // specs/packet-representation-and-sampling.md [MCNX-PKT-04] on synthetic
  // flat-static and curved inputs.
  append_tetrad_rows(b);

  // Rows 60..65 — the source-term ledger sign fixtures of
  // specs/hydro-coupling-source-terms.md [MCNX-HYD-02] through the
  // mcnux_srcterms.hxx helpers.
  append_srcterm_rows(b);

  // Rows 66..70 — the trapped-regime relabeling invariants, convention
  // tripwire, and alpha = 1 identity of specs/trapped-regime-treatment.md
  // [MCNX-TRP-02/03] through the mcnux_trp.hxx pure map.
  append_trp_rows(b);

  // Rows 71..78 — the hc pin, the analytic (gray) opacity mode, and the
  // source-agnostic coefficient dispatch of
  // specs/verification-suite-design.md [MCNX-VER-05] through
  // mcnux_coefficients.hxx and its parameter glue.
  append_coefficient_rows(b);

  // Rows 79..85 — the geodesic gather/RHS/integrator legs of
  // specs/geodesic-propagation.md [MCNX-GEO-02/03/04] through
  // mcnux_geodesic.hxx on synthetic vertex data and analytic metrics.
  append_geodesic_rows(b);

  // Rows 86..91 — the baseline table-source coefficient assembly (Kirchhoff
  // closure) and nu_x mapping of specs/opacity-eos-evaluation.md
  // [MCNX-OPA-04/05] through the mcnux_table_coeffs.hxx callable on
  // synthetic in-memory tables, plus its plumb-through of the [MCNX-VER-05]
  // dispatch.
  append_baseline_table_rows(b);

  // Rows 92..99 — the table-range enforcement policy (transparency floor,
  // intersection clamping, per-axis clamp counters, no-NaN guarantee) and
  // the EOS-inversion error protocol of specs/opacity-eos-evaluation.md
  // [MCNX-OPA-06/07] through mcnux_table_range.hxx.
  append_table_range_rows(b);

  // Rows 100..105 — the statistical-acceptance z-score reduction and pinned
  // constants of specs/rng-and-statistical-acceptance.md [MCNX-RNG-07] and
  // specs/verification-suite-design.md [MCNX-VER-07] through mcnux_stats.hxx.
  append_stats_rows(b);
}

// Size of the MCNuX::mcnux_selftest array as declared in interface.ccl.
int selftest_array_size() {
  const int gi = CCTK_GroupIndex("MCNuX::mcnux_selftest");
  if (gi < 0)
    CCTK_VERROR("MCNuX::mcnux_selftest is not a known group");
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("MCNuX::mcnux_selftest must be a 1-dimensional grid array");
  return int(*sizes[0]);
}

} // namespace

extern "C" void MCNuX_SelfTest(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SelfTest;
  DECLARE_CCTK_PARAMETERS;

  Battery battery;
  run_battery(battery);

  const int nrows = battery.size();
  const int nslots = selftest_array_size();
  if (nrows != nslots)
    CCTK_VERROR("MCNuX self-test battery has %d rows, but "
                "MCNuX::mcnux_selftest is declared with SIZE=%d in "
                "MCNuX/interface.ccl; the two must agree (and the golden data "
                "in MCNuX/test/unit-selftest/ must be regenerated whenever "
                "rows are appended)",
                nrows, nslots);

  int nfailed = 0;
  for (int i = 0; i < nrows; ++i) {
    const CheckRow &row = battery[i];
    selftest_error[i] = row.error;
    selftest_pass[i] = row.pass ? 1.0 : 0.0;

    // The TSV table is numeric only (it is diffed column-wise by the harness),
    // so the index -> name mapping is echoed to the log instead.
    CCTK_VINFO("MCNuX selftest row %3d  %-40s  error = %.17g  %s", i, row.name,
               row.error, row.pass ? "PASS" : "FAIL");
    if (!row.pass)
      ++nfailed;
  }

  if (nfailed > 0)
    CCTK_VWARN(CCTK_WARN_ALERT,
               "MCNuX unit self-test battery: %d of %d checks FAILED (see the "
               "per-row lines above; the golden diff of the `unit-selftest` "
               "benchmark is the scoring mechanism)",
               nfailed, nrows);
  else
    CCTK_VINFO("MCNuX unit self-test battery: all %d checks passed", nrows);
}

} // namespace MCNuX
