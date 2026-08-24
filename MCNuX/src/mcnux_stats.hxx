#ifndef MCNUX_STATS_HXX
#define MCNUX_STATS_HXX

#include "mcnux_units.hxx" // detail::cabs

#include <cstdint>
#include <type_traits>

// The MCNuX statistical-acceptance machinery: the pure 4-sigma z-score
// reduction plus the pinned packet count and seeds.
//
// Governing specs:
//   * specs/rng-and-statistical-acceptance.md [MCNX-RNG-07] — the acceptance
//     bar |estimate - expected| <= 4 sigma (NON-STRICT), the pinned default
//     packet count N_p = 1048576 (2^20), the primary seed 1296518744, and
//     the secondary seed 20260730 for seed-robustness cross-checks.
//     NOTE (spec defect, recorded): the spec glosses the primary seed as
//     "(0x4D434E58)", but 0x4D434E58 == 1296256600 != 1296518744; the
//     decimal literal governs — it is what the spec's own validator pins
//     (rng-and-statistical-acceptance.md:123) and what its Verification
//     section restates (:118). 1296518744 == 0x4D474E58. The hex gloss
//     needs a spec correction, not a code change.
//   * specs/verification-suite-design.md [MCNX-VER-07] — the per-check output
//     shape: estimate, expected, sigma, z = (estimate - expected)/sigma, and
//     acceptance |z| <= 4; archived z values later become golden numbers.
//
// Deliberate non-goals (binding scope of this header):
//   * NO sigma formula lives here. Sigma ownership is per-benchmark
//     (verification-suite-design.md:453-454): each statistical benchmark
//     documents its own analytically propagated (or documented sample)
//     standard error and passes it in. This machinery never computes one,
//     and carries no estimator or N_p-dependence.
//   * The pinned values below are DEFAULTS: a benchmark pinning different
//     values documents them in its parfile
//     (verification-suite-design.md:277-279). Golden statistical data is
//     always generated with the primary seed. A statistical check may not be
//     re-run with fresh seeds to pass ("seed shopping", [MCNX-RNG-07] last
//     sentence) — a prohibition on process, documented here, not enforceable
//     in code.
//
// Everything is constexpr noexcept (the mcnux_rng.hxx/mcnux_units.hxx
// convention), so this header carries its own compile-time verification block
// and is includable from host code, device code (AMReX builds with
// --expt-relaxed-constexpr), and standalone unit tests alike.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes.

namespace MCNuX {

// ---------------------------------------------------------------------------
// Pinned statistical-acceptance constants  [MCNX-RNG-07]
// ---------------------------------------------------------------------------
// rng-and-statistical-acceptance.md:104. These literals appear in MCNuX
// exactly once — here. The seed type matches the S argument of u(S, q, e, k)
// (mcnux_rng.hxx).

// Default packet count for statistical checks: 2^20.
inline constexpr std::int64_t stats_default_num_packets = 1048576;

// Primary seed: all golden statistical data uses this seed. (Decimal is the
// normative spelling; the spec's hex gloss is wrong — see the header note.)
inline constexpr std::uint64_t stats_seed_primary = 1296518744;

// Secondary seed: seed-robustness cross-checks only, never golden data.
inline constexpr std::uint64_t stats_seed_secondary = 20260730;

// The acceptance bar in units of sigma  ([MCNX-RNG-07]/[MCNX-VER-07]: 4).
inline constexpr double stats_sigma_bound = 4.0;

// ---------------------------------------------------------------------------
// The z-score reduction  [MCNX-VER-07]
// ---------------------------------------------------------------------------
// The four-field per-check output shape of verification-suite-design.md
// :255-279, so the stats-* benchmarks can emit the fields directly, plus the
// judged acceptance flag.

struct ZScore {
  double estimate; // the Monte Carlo estimate
  double expected; // the analytic/expected value
  double sigma;    // the caller-supplied standard error (per-benchmark)
  double z;        // (estimate - expected) / sigma
  bool pass;       // |estimate - expected| <= stats_sigma_bound * sigma
};

// Build the ZScore for one check. `pass` is judged in the spec's primitive
// multiplication form |estimate - expected| <= 4 * sigma — non-strict <=
// (rng-and-statistical-acceptance.md:101), no division — so it stays
// well-defined at the degenerate point sigma = 0 (estimate == expected
// passes; any difference fails). The reported z is only meaningful for
// sigma > 0; no benchmark may rely on sigma = 0 behavior of z.
constexpr ZScore zscore(double estimate, double expected,
                        double sigma) noexcept {
  return ZScore{estimate, expected, sigma, (estimate - expected) / sigma,
                detail::cabs(estimate - expected) <= stats_sigma_bound * sigma};
}

// ---------------------------------------------------------------------------
// Compile-time verification  ([MCNX-RNG-07], [MCNX-VER-07])
// ---------------------------------------------------------------------------
// Exact binary-fraction fixtures only, so every equality below is bitwise.
// The runtime mirror rows live in mcnux_selftest.cxx (rows 100..105).

// The pinned constants, against their spec literals — including the 2^20
// identity the spec states in both forms. The primary seed is asserted
// against the normative decimal literal and its true hex spelling
// 0x4D474E58 (the spec's "(0x4D434E58)" gloss is arithmetically wrong; see
// the header note).
static_assert(stats_default_num_packets == 1048576 &&
                  stats_default_num_packets == (std::int64_t(1) << 20),
              "[MCNX-RNG-07] default statistical packet count is 2^20");
static_assert(stats_seed_primary == 1296518744 &&
                  stats_seed_primary == 0x4D474E58,
              "[MCNX-RNG-07] primary statistical seed pin");
static_assert(stats_seed_secondary == 20260730,
              "[MCNX-RNG-07] secondary statistical seed pin");
static_assert(stats_sigma_bound == 4.0,
              "[MCNX-RNG-07]/[MCNX-VER-07] acceptance bar is 4 sigma");

// z on an exact binary fraction: (2.5 - 1.0)/0.5 == 3.0 bitwise, in band.
static_assert(zscore(2.5, 1.0, 0.5).z == 3.0 && zscore(2.5, 1.0, 0.5).pass,
              "[MCNX-VER-07] z = (estimate - expected)/sigma, exact fixture");

// Negative side: (0.25 - 1.0)/0.25 == -3.0 bitwise, in band (|z| judged).
static_assert(zscore(0.25, 1.0, 0.25).z == -3.0 &&
                  zscore(0.25, 1.0, 0.25).pass,
              "[MCNX-VER-07] z carries the sign of estimate - expected");

// The non-strict boundary: |diff| = 0.5 = 4 * 0.125 exactly -> pass.
static_assert(zscore(1.5, 1.0, 0.125).z == 4.0 && zscore(1.5, 1.0, 0.125).pass,
              "[MCNX-RNG-07] the 4 sigma bound is non-strict (<=)");

// One ulp-scale step beyond the boundary fails: diff = 0.5 + 2^-20.
static_assert(!zscore(1.5 + 0x1p-20, 1.0, 0.125).pass,
              "[MCNX-RNG-07] beyond 4 sigma fails");

// The four-field shape echoes its inputs unchanged.
static_assert(zscore(2.5, 1.0, 0.5).estimate == 2.5 &&
                  zscore(2.5, 1.0, 0.5).expected == 1.0 &&
                  zscore(2.5, 1.0, 0.5).sigma == 0.5,
              "[MCNX-VER-07] estimate/expected/sigma pass through unchanged");

// Dumb, device-capturable aggregate (the BaselineTableViews precedent).
static_assert(std::is_trivially_copyable<ZScore>::value,
              "ZScore must be trivially copyable");

} // namespace MCNuX

#endif // MCNUX_STATS_HXX
