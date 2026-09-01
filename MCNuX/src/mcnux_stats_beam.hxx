#ifndef MCNUX_STATS_BEAM_HXX
#define MCNUX_STATS_BEAM_HXX

// Benchmark-owned standard-error (sigma) formulas of the `stats-beam`
// statistical benchmark ([MCNX-INT-01]/[MCNX-INT-03], the beam-attenuation
// check of specs/neutrino-matter-interactions.md:197-200; [MCNX-VER-07],
// specs/verification-suite-design.md:265: transmitted fraction exp(-kappa_a L)
// through an eta_scale = 0 slab).
//
// Sigma ownership is per-benchmark, NEVER mcnux_stats.hxx (that header's
// binding non-goal: it carries the zscore() reduction and the pinned
// constants only) — the mcnux_stats_emission.hxx precedent, one small header
// per stats-* benchmark family member.
//
// Estimator convention (documented methodology of this benchmark,
// verification-suite-design.md:255-263):
//   * estimate  = the transmitted packet COUNT (the [MCNX-GPU-05] escape
//     tally of the beam species — every non-absorbed packet of a
//     pure-absorber run reaches the far face);
//   * expected  = N_p * p with p = p_transmission(kappa_code, L_code)
//                 = exp(-kappa_a L) (Beer-Lambert, kappa and L both in
//     geometrized code units — the [MCNX-INT-01] draw law's own unit
//     system);
//   * sigma     = sqrt(N_p p (1 - p)) — the binomial standard error of the
//     transmitted count: each of the N_p independent monoenergetic packets
//     survives the slab with probability p (the per-cell episode restart is
//     statistically exact by memorylessness, the [MCNX-INT-02] framing), so
//     the count is Binomial(N_p, p). The count form matches the
//     sigma_bernoulli_sum style of the stats-emission rows.
//
// The rational variance form is a separate constexpr function so the formula
// shape is pinned by exact compile-time fixtures (the mcnux_stats_emission.hxx
// precedent: production sqrt is std::sqrt, compile-time predicates stay
// rational); the sqrt leg and the exp transmission law get runtime selftest
// rows instead (mcnux_selftest_stats_beam.cxx, rows 166..169).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes.

#include <cassert>
#include <cmath>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Transmission probability (Beer-Lambert)  [MCNX-INT-01]/[MCNX-INT-03]
// ---------------------------------------------------------------------------

// p = exp(-kappa L), kappa the code-unit absorption opacity (the caller
// converts cm^-1 -> code via opacity_cgs_to_code, never a retyped literal)
// and L the code-unit slab thickness along the beam. Non-constexpr by design
// (std::exp; the interaction_time/std::log precedent) — covered by the
// runtime selftest rows, not the static_asserts below.
inline double p_transmission(double kappa_code, double L_code) noexcept {
  assert(kappa_code >= 0.0);
  assert(L_code >= 0.0);
  return std::exp(-kappa_code * L_code);
}

// ---------------------------------------------------------------------------
// Rational variance form (constexpr, exact-fixture-testable)
// ---------------------------------------------------------------------------

// Variance of the transmitted count of N_p independent packets each
// surviving with probability p: Binomial(N_p, p) -> N_p p (1 - p).
constexpr double var_beam_transmission(double N_p, double p) noexcept {
  assert(N_p > 0.0);
  assert(p >= 0.0 && p <= 1.0);
  return N_p * p * (1.0 - p);
}

// ---------------------------------------------------------------------------
// Standard error (std::sqrt on the production path)
// ---------------------------------------------------------------------------

// sigma of the transmitted count: sqrt(N_p p (1 - p)). Exact on
// perfect-square variances; 0 at the deterministic endpoints p in {0, 1}
// (an all-absorbed or all-transmitted beam has no binomial spread — the
// benchmark never runs there: its writer aborts on a degenerate p).
inline double sigma_beam_transmission(double N_p, double p) noexcept {
  return std::sqrt(var_beam_transmission(N_p, p));
}

// ---------------------------------------------------------------------------
// Compile-time verification (rational form; the sqrt/exp legs are runtime
// rows 166..169 of the unit-selftest battery, mcnux_selftest_stats_beam.cxx)
// ---------------------------------------------------------------------------

// The formula shape on exactly-representable binary fractions, including the
// pinned N_p = 2^20 of [MCNX-RNG-07] at p = 1/2 (0.5 * 0.5 = 0.25 and
// 2^20 * 0.25 = 2^18 are exact).
static_assert(var_beam_transmission(1048576.0, 0.5) == 262144.0 &&
                  var_beam_transmission(4.0, 0.25) == 0.75,
              "[MCNX-VER-07] transmitted-count variance is N_p p (1 - p) "
              "(binomial; the stats-beam estimator)");

// Deterministic endpoints: p = 0 (all absorbed) and p = 1 (all transmitted)
// have zero variance, exactly.
static_assert(var_beam_transmission(1048576.0, 0.0) == 0.0 &&
                  var_beam_transmission(1048576.0, 1.0) == 0.0,
              "a deterministic beam outcome has no binomial spread");

// p <-> 1 - p symmetry on exact quarters, and monotone growth with N_p.
static_assert(var_beam_transmission(8.0, 0.25) ==
                      var_beam_transmission(8.0, 0.75) &&
                  var_beam_transmission(2.0, 0.5) <
                      var_beam_transmission(4.0, 0.5),
              "binomial variance is symmetric in p <-> 1 - p and grows "
              "with N_p");

} // namespace MCNuX

#endif // MCNUX_STATS_BEAM_HXX
