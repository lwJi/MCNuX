#ifndef MCNUX_STATS_EMISSION_HXX
#define MCNUX_STATS_EMISSION_HXX

// Benchmark-owned standard-error (sigma) formulas of the `stats-emission`
// statistical benchmark ([MCNX-PKT-06], [MCNX-VER-07],
// specs/verification-suite-design.md:255-279) and the RNG
// uniformity-at-the-pinned-bar rows of
// specs/rng-and-statistical-acceptance.md:118.
//
// Sigma ownership is per-benchmark, NEVER mcnux_stats.hxx (that header's
// binding non-goal: it carries the zscore() reduction and the pinned
// constants only; verification-suite-design.md:457-458 grants each benchmark
// its own documented methodology). These formulas are shared analytic
// content for the stats-* benchmark family, so they live in their own small
// header rather than file-local in mcnux_emission.cxx — second consumer
// MCNuX_StatsScatterbox (mcnux_interactions.cxx), whose post-scatter
// isotropy rows reuse the sigma_isotropy_* pair below verbatim: the
// [MCNX-INT-04] scatter draw is the identical cos theta = 2 u1 - 1,
// phi = 2 pi u2 map as emission's [MCNX-PKT-04] direction draw:
//
//   * sigma_bernoulli_sum — the count-law estimator's standard error. The
//     [MCNX-PKT-02] floor + Bernoulli-remainder construction makes each
//     (cell, species, bin) realized count floor(N_p) + Bernoulli(p) with
//     p = frac(N_p), so the total count's variance is Sum p(1 - p) over the
//     independent per-cell draws and sigma = sqrt(Sum p(1 - p)).
//   * sigma_uniform_mean — the spec-stated sigma of the RNG uniformity mean
//     row: Var(U(0,1)) = 1/12, sigma = 1/sqrt(12 N)
//     (rng-and-statistical-acceptance.md:118, spec-pinned).
//   * sigma_uniform_second_moment — the benchmark-documented sigma of the
//     second-moment row: Var(u^2) = E[u^4] - E[u^2]^2 = 1/5 - 1/9 = 4/45,
//     sigma = sqrt(4/(45 N)).
//   * sigma_isotropy_costheta — the cos(theta) first-moment row: cos(theta)
//     = 2 u1 - 1 is U(-1, 1) with variance 1/3 ([MCNX-PKT-04] isotropic
//     sampling), so the sample mean over n packets has
//     sigma = sqrt(1/(3 n)).
//   * sigma_isotropy_phi — the phi first-moment row: phi = 2 pi u2 is
//     U(0, 2 pi) with variance (2 pi)^2/12 = pi^2/3, so the sample mean has
//     sigma = pi sqrt(1/(3 n)) — implemented literally as pi times the
//     cos(theta) sigma so the proportionality is bitwise.
//
// The rational variance forms are separate constexpr functions so the
// formula shape is pinned by exact compile-time fixtures (the tetrad-header
// precedent: production sqrt is std::sqrt, compile-time predicates stay
// rational). Runtime measured-error rows live in
// mcnux_selftest_stats_emission.cxx (rows 150..155).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes.

#include "mcnux_tetrad.hxx" // detail::pi (the single in-repo pi literal)

#include <cassert>
#include <cmath>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Rational variance forms (constexpr, exact-fixture-testable)
// ---------------------------------------------------------------------------

// Variance of the sample mean of N iid U(0, 1) draws: 1/(12 N).
constexpr double var_uniform_mean(double N) noexcept {
  assert(N > 0.0);
  return 1.0 / (12.0 * N);
}

// Variance of the sample second moment of N iid U(0, 1) draws:
// Var(u^2)/N = (1/5 - 1/9)/N = 4/(45 N).
constexpr double var_uniform_second_moment(double N) noexcept {
  assert(N > 0.0);
  return 4.0 / (45.0 * N);
}

// Variance of the sample mean of n iid cos(theta) = 2 u1 - 1 ~ U(-1, 1)
// draws: (1/3)/n.
constexpr double var_isotropy_costheta(double n) noexcept {
  assert(n > 0.0);
  return 1.0 / (3.0 * n);
}

// ---------------------------------------------------------------------------
// Standard errors (std::sqrt on the production path)
// ---------------------------------------------------------------------------

// sigma of the total realized packet count: sqrt(Sum p(1 - p)) over the
// independent per-(cell, species, bin) Bernoulli remainders, p = frac(N_p).
// The caller accumulates Sum p(1 - p) from the PRE-floor real-valued N_p of
// the count pass. Exact on perfect-square sums; 0 at 0 (every p in {0} —
// a deterministic count has no Bernoulli spread).
inline double sigma_bernoulli_sum(double sum_p_one_minus_p) noexcept {
  assert(sum_p_one_minus_p >= 0.0);
  return std::sqrt(sum_p_one_minus_p);
}

// sigma of the RNG uniformity sample mean: 1/sqrt(12 N) (spec-stated,
// rng-and-statistical-acceptance.md:118).
inline double sigma_uniform_mean(double N) noexcept {
  return std::sqrt(var_uniform_mean(N));
}

// sigma of the RNG uniformity sample second moment: sqrt(4/(45 N))
// (benchmark-documented; Var(u^2) = 4/45).
inline double sigma_uniform_second_moment(double N) noexcept {
  return std::sqrt(var_uniform_second_moment(N));
}

// sigma of the per-species cos(theta) sample mean over n realized packets:
// sqrt(1/(3 n)) (documented sample-size choice: n is the realized count).
inline double sigma_isotropy_costheta(double n) noexcept {
  return std::sqrt(var_isotropy_costheta(n));
}

// sigma of the per-species phi sample mean over n realized packets:
// pi sqrt(1/(3 n)) (Var(U(0, 2 pi)) = pi^2/3), literally pi times the
// cos(theta) sigma so the two are bitwise proportional.
inline double sigma_isotropy_phi(double n) noexcept {
  return detail::pi * sigma_isotropy_costheta(n);
}

// ---------------------------------------------------------------------------
// Compile-time verification (rational forms; the sqrt legs are runtime rows
// 150..155 of the unit-selftest battery, mcnux_selftest_stats_emission.cxx)
// ---------------------------------------------------------------------------

// The formula shapes at exactly-representable multiplications (12 * N,
// 45 * N, 3 * n exact), against the identically-rounded literal quotients.
static_assert(var_uniform_mean(1.0) == 1.0 / 12.0 &&
                  var_uniform_mean(0.5) == 1.0 / 6.0,
              "[MCNX-RNG-07] uniformity-mean variance is 1/(12 N) "
              "(rng-and-statistical-acceptance.md:118)");
static_assert(var_uniform_second_moment(1.0) == 4.0 / 45.0 &&
                  var_uniform_second_moment(2.0) == 4.0 / 90.0,
              "[MCNX-VER-07] uniformity second-moment variance is 4/(45 N) "
              "(Var(u^2) = 1/5 - 1/9 = 4/45)");
static_assert(var_isotropy_costheta(1.0) == 1.0 / 3.0 &&
                  var_isotropy_costheta(4.0) == 1.0 / 12.0,
              "[MCNX-PKT-04]/[MCNX-VER-07] cos(theta) sample-mean variance "
              "is 1/(3 n) (Var(U(-1, 1)) = 1/3)");

// Monotone shrinkage with sample size (exact comparisons on exact halves).
static_assert(var_uniform_mean(2.0) < var_uniform_mean(1.0) &&
                  var_uniform_second_moment(2.0) <
                      var_uniform_second_moment(1.0) &&
                  var_isotropy_costheta(2.0) < var_isotropy_costheta(1.0),
              "standard errors must shrink with sample size");

} // namespace MCNuX

#endif // MCNUX_STATS_EMISSION_HXX
