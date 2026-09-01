// ---------------------------------------------------------------------------
// Rows 166..169 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the benchmark-owned transmission-probability and
// standard-error formulas of the stats-beam statistical benchmark
// ([MCNX-INT-01]/[MCNX-INT-03] beam attenuation of
// specs/neutrino-matter-interactions.md:197-200; [MCNX-VER-07],
// specs/verification-suite-design.md:265) through mcnux_stats_beam.hxx.
// Sigma ownership stays per-benchmark (mcnux_stats.hxx binding non-goal);
// these rows pin the exp/sqrt legs the header's compile-time asserts cannot
// (the rational binomial variance form is asserted there): the exact
// zero-depth transmission endpoint, the machine-tier Beer-Lambert
// multiplicativity, an exact perfect-square sigma fixture at the pinned
// N_p = 2^20, and the sigma = 0 degenerate endpoints p in {0, 1}.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_stats.hxx"
#include "mcnux_stats_beam.hxx"

namespace MCNuX {

void append_stats_beam_rows(Battery &b) {
  // Row 166 — the zero-depth endpoint: exp(-kappa * 0) == 1 exactly
  // (IEEE 754 exp(±0) is exact), any kappa.
  b.add_exact("statsbeam.p.zero_depth_one", p_transmission(0.53125, 0.0), 1.0);

  // Row 167 — Beer-Lambert multiplicativity at machine tier on exact binary
  // depths: p(kappa, L1 + L2) == p(kappa, L1) * p(kappa, L2) (the
  // memorylessness that makes the per-cell episode restart statistically
  // exact, [MCNX-INT-02]).
  {
    const double kappa = 0.5;
    const double L1 = 0.75, L2 = 1.125; // L1 + L2 = 1.875 exactly
    b.add_boolean("statsbeam.p.multiplicative",
                  detail::approx_eq(p_transmission(kappa, L1 + L2),
                                    p_transmission(kappa, L1) *
                                        p_transmission(kappa, L2),
                                    detail::rtol_machine));
  }

  // Row 168 — binomial count sigma on an exact perfect square at the pinned
  // N_p = 2^20 ([MCNX-RNG-07]): var(2^20, 1/2) = 2^18, sigma = 512 bitwise
  // (sqrt is correctly rounded).
  b.add_exact("statsbeam.sigma.exact",
              sigma_beam_transmission(double(stats_default_num_packets), 0.5),
              512.0);

  // Row 169 — the degenerate deterministic endpoints: p = 0 (all absorbed)
  // and p = 1 (all transmitted) give sigma exactly 0.
  b.add_boolean("statsbeam.sigma.degenerate_zero",
                sigma_beam_transmission(1048576.0, 0.0) == 0.0 &&
                    sigma_beam_transmission(1048576.0, 1.0) == 0.0);
}

} // namespace MCNuX
