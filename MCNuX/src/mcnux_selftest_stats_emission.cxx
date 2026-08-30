// ---------------------------------------------------------------------------
// Rows 150..155 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the benchmark-owned standard-error formulas of the
// stats-emission statistical benchmark ([MCNX-PKT-06], [MCNX-VER-07],
// specs/verification-suite-design.md:255-279; the RNG uniformity sigmas of
// specs/rng-and-statistical-acceptance.md:118) through
// mcnux_stats_emission.hxx. Sigma ownership stays per-benchmark
// (mcnux_stats.hxx binding non-goal); these rows pin the sqrt legs the
// header's compile-time asserts cannot (the rational variance forms are
// asserted there): exact perfect-square fixtures, the sigma = 0 degenerate
// point of a deterministic count, machine-tier squared identities
// sigma^2 == variance at the pinned N = 2^20, and the bitwise
// pi-proportionality of the phi sigma.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_stats.hxx"
#include "mcnux_stats_emission.hxx"
#include "mcnux_tetrad.hxx"

namespace MCNuX {

void append_stats_emission_rows(Battery &b) {
  // Row 150 — Bernoulli-sum sigma on exact perfect squares: sqrt is
  // correctly rounded (IEEE 754), so sqrt(0.25) == 0.5 and sqrt(2.25) == 1.5
  // bitwise.
  b.add_boolean("statsem.sigma.bernoulli_exact",
                sigma_bernoulli_sum(0.25) == 0.5 &&
                    sigma_bernoulli_sum(2.25) == 1.5);

  // Row 151 — the degenerate deterministic-count point: Sum p(1 - p) = 0
  // (every remainder p in {0}) gives sigma exactly 0.
  b.add_exact("statsem.sigma.bernoulli_zero", sigma_bernoulli_sum(0.0), 0.0);

  // Row 152 — uniform-mean sigma at the pinned N = 2^20: squared identity
  // 12 N sigma^2 == 1 at machine tier ([MCNX-RNG-07] spec-stated
  // sigma = 1/sqrt(12 N)).
  {
    const double N = double(stats_default_num_packets);
    const double s = sigma_uniform_mean(N);
    b.add_boolean("statsem.sigma.uniform_mean",
                  detail::approx_eq(12.0 * N * s * s, 1.0,
                                    detail::rtol_machine));
  }

  // Row 153 — uniform-second-moment sigma at the pinned N = 2^20: squared
  // identity 45 N sigma^2 / 4 == 1 at machine tier (Var(u^2) = 4/45).
  {
    const double N = double(stats_default_num_packets);
    const double s = sigma_uniform_second_moment(N);
    b.add_boolean("statsem.sigma.uniform_m2",
                  detail::approx_eq(45.0 * N * s * s / 4.0, 1.0,
                                    detail::rtol_machine));
  }

  // Row 154 — cos(theta) sample-mean sigma: squared identity
  // 3 n sigma^2 == 1 at machine tier on a non-square fixture
  // (Var(U(-1, 1)) = 1/3).
  {
    const double n = 12345.0;
    const double s = sigma_isotropy_costheta(n);
    b.add_boolean("statsem.sigma.costheta",
                  detail::approx_eq(3.0 * n * s * s, 1.0,
                                    detail::rtol_machine));
  }

  // Row 155 — the phi sigma is bitwise pi times the cos(theta) sigma
  // (Var(U(0, 2 pi)) = pi^2/3; implemented literally as the product).
  b.add_exact("statsem.sigma.phi_pi_proportional", sigma_isotropy_phi(12345.0),
              detail::pi * sigma_isotropy_costheta(12345.0));
}

} // namespace MCNuX
