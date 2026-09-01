#ifndef MCNUX_STATS_SCATTERBOX_HXX
#define MCNUX_STATS_SCATTERBOX_HXX

// Benchmark-owned standard-error (sigma) formulas of the `stats-scatterbox`
// statistical benchmark ([MCNX-INT-02]/[MCNX-INT-04], the collision-
// statistics check of specs/neutrino-matter-interactions.md:201-204;
// [MCNX-VER-07], specs/verification-suite-design.md:266: per-packet
// scattering-event counts Poisson(lambda = kappa_s l_path) in a uniform
// pure-scattering box, sample mean AND sample variance within 4 sigma).
//
// Sigma ownership is per-benchmark, NEVER mcnux_stats.hxx (that header's
// binding non-goal: it carries the zscore() reduction and the pinned
// constants only) — the mcnux_stats_emission.hxx / mcnux_stats_beam.hxx
// precedent, one small header per stats-* benchmark family member. The
// post-scatter isotropy rows of the same benchmark (spec :205-207) REUSE
// sigma_isotropy_costheta / sigma_isotropy_phi of mcnux_stats_emission.hxx
// (identical distribution: the [MCNX-INT-04] scatter draw is
// cos theta = 2 u1 - 1, phi = 2 pi u2, the same map as emission's
// [MCNX-PKT-04] direction draw) — only the Poisson formulas below are new
// content.
//
// Estimator convention (documented methodology of this benchmark,
// verification-suite-design.md:255-263):
//   * l_path := c * Delta t of the SINGLE transport step = cctk_delta_time
//     in code units (c = 1). A null packet on the at-rest Minkowski
//     background moves at coordinate speed 1 (the stats-beam precedent), so
//     one transport step gives every packet path length exactly dt —
//     including across scattering direction changes (speed stays 1) and
//     cell crossings.
//   * The per-cell episode restart ([MCNX-INT-02]: a cell crossing discards
//     the undrawn remainder and redraws the exponential interaction time) is
//     statistically exact by memorylessness, so in a UNIFORM medium the
//     per-packet event count over path length l is exactly
//     Poisson(kappa_s_code * l) — kappa_s in code units (cm^-1 -> code via
//     opacity_cgs_to_code, never a retyped literal).
//   * Row-0 estimate = the sample mean kbar = Sum k / N over the N
//     independent per-packet counts; Var(kbar) = lambda/N (Poisson variance
//     equals the mean), sigma = sqrt(lambda/N).
//   * Row-1 estimate = the UNBIASED sample variance
//     S^2 = (Sum k^2 - N kbar^2)/(N - 1) (the N - 1 divisor is the pinned
//     estimator choice; the divisor difference is O(1/N), below the
//     large-N approximation already made in Var(S^2) below).
//     Var(S^2) ~ (mu_4 - sigma^4)/N to leading order in 1/N; for Poisson
//     the central moments are mu_2 = lambda and mu_4 = lambda(1 + 3 lambda),
//     so mu_4 - mu_2^2 = lambda + 3 lambda^2 - lambda^2 = lambda +
//     2 lambda^2 and sigma = sqrt((lambda + 2 lambda^2)/N). This is a
//     documented large-N approximation (exactly the leading order at the
//     pinned N = 2^20; NOT the Gaussian 2 lambda^2/N — Poisson kurtosis
//     contributes the extra lambda/N term at the same order).
//
// The rational variance forms are separate constexpr functions so the
// formula shapes are pinned by exact compile-time fixtures (the
// mcnux_stats_beam.hxx precedent: production sqrt is std::sqrt, compile-time
// predicates stay rational); the sqrt legs get runtime selftest rows instead
// (mcnux_selftest_stats_scatterbox.cxx, rows 170..173).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes.

#include <cassert>
#include <cmath>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Rational variance forms (constexpr, exact-fixture-testable)
// ---------------------------------------------------------------------------

// Variance of the sample mean of n iid Poisson(lambda) counts: lambda/n
// (Poisson variance equals the mean).
constexpr double var_poisson_mean(double lambda, double n) noexcept {
  assert(lambda >= 0.0);
  assert(n > 0.0);
  return lambda / n;
}

// Leading-order variance of the (unbiased) sample variance of n iid
// Poisson(lambda) counts: (mu_4 - mu_2^2)/n = (lambda + 2 lambda^2)/n
// (mu_2 = lambda, mu_4 = lambda(1 + 3 lambda); large-n approximation, see
// the header doc — valid at the pinned n = 2^20).
constexpr double var_poisson_variance(double lambda, double n) noexcept {
  assert(lambda >= 0.0);
  assert(n > 0.0);
  return (lambda + 2.0 * lambda * lambda) / n;
}

// ---------------------------------------------------------------------------
// Standard errors (std::sqrt on the production path)
// ---------------------------------------------------------------------------

// sigma of the sample mean: sqrt(lambda/n). Exact on perfect-square
// quotients; 0 at lambda = 0 (a zero-rate medium scatters deterministically
// never — the benchmark never runs there: its writer aborts on lambda = 0).
inline double sigma_poisson_mean(double lambda, double n) noexcept {
  return std::sqrt(var_poisson_mean(lambda, n));
}

// sigma of the unbiased sample variance: sqrt((lambda + 2 lambda^2)/n).
inline double sigma_poisson_variance(double lambda, double n) noexcept {
  return std::sqrt(var_poisson_variance(lambda, n));
}

// ---------------------------------------------------------------------------
// Compile-time verification (rational forms; the sqrt legs are runtime rows
// 170..173 of the unit-selftest battery, mcnux_selftest_stats_scatterbox.cxx)
// ---------------------------------------------------------------------------

// The formula shapes on exactly-representable binary fractions: at
// lambda = 1/2, n = 4 the mean variance is (1/2)/4 = 1/8 exactly and the
// variance-of-variance is (1/2 + 2/4)/4 = 1/4 exactly.
static_assert(var_poisson_mean(0.5, 4.0) == 0.125 &&
                  var_poisson_mean(4.0, 1048576.0) == 3.814697265625e-06,
              "[MCNX-VER-07] Poisson sample-mean variance is lambda/n "
              "(the stats-scatterbox row-0 estimator)");

// The lambda + 2 lambda^2 numerator (NEVER the Gaussian 2 lambda^2: the
// Poisson kurtosis term lambda matters at the same order), on exact binary
// fractions incl. the lambda = 1 point where the two differ by 50%.
static_assert(var_poisson_variance(0.5, 4.0) == 0.25 &&
                  var_poisson_variance(1.0, 1.0) == 3.0 &&
                  var_poisson_variance(2.0, 1.0) == 10.0,
              "[MCNX-VER-07] Poisson sample-variance variance is "
              "(lambda + 2 lambda^2)/n (mu_4 - mu_2^2 = lambda + 2 lambda^2)");

// Degenerate endpoint and monotone shrinkage with sample size.
static_assert(var_poisson_mean(0.0, 1048576.0) == 0.0 &&
                  var_poisson_variance(0.0, 1048576.0) == 0.0,
              "a zero-rate Poisson count is deterministic (no spread)");
static_assert(var_poisson_mean(0.5, 8.0) < var_poisson_mean(0.5, 4.0) &&
                  var_poisson_variance(0.5, 8.0) <
                      var_poisson_variance(0.5, 4.0),
              "standard errors must shrink with sample size");

} // namespace MCNuX

#endif // MCNUX_STATS_SCATTERBOX_HXX
