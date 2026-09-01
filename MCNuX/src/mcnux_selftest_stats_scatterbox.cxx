// ---------------------------------------------------------------------------
// Rows 170..173 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the benchmark-owned Poisson standard-error formulas
// of the stats-scatterbox statistical benchmark
// ([MCNX-INT-02]/[MCNX-INT-04] collision statistics of
// specs/neutrino-matter-interactions.md:201-207; [MCNX-VER-07],
// specs/verification-suite-design.md:266) through
// mcnux_stats_scatterbox.hxx. Sigma ownership stays per-benchmark
// (mcnux_stats.hxx binding non-goal); these rows pin the sqrt legs the
// header's compile-time asserts cannot (the rational variance forms are
// asserted there): exact perfect-square sigma fixtures for both the
// sample-mean and sample-variance estimators, incl. one at the pinned
// N = 2^20. The isotropy sigmas the benchmark also consumes are the reused
// sigma_isotropy_* pair of mcnux_stats_emission.hxx, already pinned by
// rows 154..155 — not duplicated here.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_stats.hxx"
#include "mcnux_stats_scatterbox.hxx"

namespace MCNuX {

void append_stats_scatterbox_rows(Battery &b) {
  // Row 170 — the rational sample-mean variance lambda/n on the exact
  // binary-fraction fixture lambda = 1/2, n = 4 -> 1/8 exactly.
  b.add_exact("statsscat.var.mean_exact", var_poisson_mean(0.5, 4.0), 0.125);

  // Row 171 — the rational sample-variance variance (lambda + 2 lambda^2)/n
  // on the same fixture: (1/2 + 1/2)/4 = 1/4 exactly. The lambda term is
  // load-bearing (the Gaussian-only 2 lambda^2/n would give 1/8 here).
  b.add_exact("statsscat.var.variance_exact", var_poisson_variance(0.5, 4.0),
              0.25);

  // Row 172 — sample-mean sigma on an exact perfect square at the pinned
  // N_p = 2^20 ([MCNX-RNG-07]): var(4, 2^20) = 2^-18, sigma = 2^-9 bitwise
  // (sqrt is correctly rounded).
  b.add_exact("statsscat.sigma.mean_exact",
              sigma_poisson_mean(4.0, double(stats_default_num_packets)),
              0.001953125);

  // Row 173 — sample-variance sigma on an exact perfect square:
  // var(1/2, 4) = 1/4, sigma = 1/2 bitwise.
  b.add_exact("statsscat.sigma.variance_exact",
              sigma_poisson_variance(0.5, 4.0), 0.5);
}

} // namespace MCNuX
