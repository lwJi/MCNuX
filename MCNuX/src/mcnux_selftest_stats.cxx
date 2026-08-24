// ---------------------------------------------------------------------------
// Rows 100..105 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the statistical-acceptance machinery of
// specs/rng-and-statistical-acceptance.md [MCNX-RNG-07] and
// specs/verification-suite-design.md [MCNX-VER-07] (T26) through
// mcnux_stats.hxx. All rows are exact tier on exact binary-fraction fixtures
// (halves/quarters/2^-20 steps), so every equality is bitwise: the z closed
// form and its sign, the z = 0 degenerate agreement point, the NON-STRICT
// 4 sigma boundary (z exactly 4 passes; one 2^-20 step beyond fails), and
// the pinned constants (N_p = 2^20, primary seed 0x4D434E58, secondary seed,
// the 4 sigma bound). Sigma is a caller input throughout — no sigma formula
// is exercised or owned here (sigma ownership is per-benchmark,
// verification-suite-design.md:453-454).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_stats.hxx"

#include <cstdint>

namespace MCNuX {

void append_stats_rows(Battery &b) {
  // Row 100 — z closed form, exact: (2.5 - 1.0)/0.5 == 3.0 bitwise.
  b.add_exact("stats.z.exact", zscore(2.5, 1.0, 0.5).z, 3.0);

  // Row 101 — z carries the sign of estimate - expected and the four-field
  // shape echoes its inputs: (0.25 - 1.0)/0.25 == -3.0 exactly, in band.
  {
    const ZScore r = zscore(0.25, 1.0, 0.25);
    b.add_boolean("stats.z.sign", r.z == -3.0 && r.pass && r.estimate == 0.25 &&
                                      r.expected == 1.0 && r.sigma == 0.25);
  }

  // Row 102 — estimate == expected: z is exactly 0 and the check passes.
  {
    const ZScore r = zscore(1.5, 1.5, 0.25);
    b.add_boolean("stats.z.zero", r.z == 0.0 && r.pass);
  }

  // Row 103 — the boundary is non-strict: |diff| = 0.5 = 4 * 0.125 exactly,
  // z exactly 4, PASS ([MCNX-RNG-07]: |estimate - expected| <= 4 sigma).
  {
    const ZScore r = zscore(1.5, 1.0, 0.125);
    b.add_boolean("stats.z.boundary_pass", r.z == 4.0 && r.pass);
  }

  // Row 104 — one 2^-20 step beyond the boundary fails: diff = 0.5 + 2^-20,
  // z > 4 strictly.
  {
    const ZScore r = zscore(1.5 + 0x1p-20, 1.0, 0.125);
    b.add_boolean("stats.z.boundary_fail", !r.pass && r.z > 4.0);
  }

  // Row 105 — the pinned constants against their spec literals, including
  // the N_p == 2^20 identity. The primary seed is judged against the
  // normative decimal 1296518744 == 0x4D474E58; the spec's "(0x4D434E58)"
  // gloss is arithmetically wrong (see the mcnux_stats.hxx header note).
  b.add_boolean("stats.constants",
                stats_default_num_packets == 1048576 &&
                    stats_default_num_packets == (std::int64_t(1) << 20) &&
                    stats_seed_primary == 1296518744 &&
                    stats_seed_primary == 0x4D474E58 &&
                    stats_seed_secondary == 20260730 &&
                    stats_sigma_bound == 4.0);
}

} // namespace MCNuX
