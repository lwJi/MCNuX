// ---------------------------------------------------------------------------
// Rows 66..70 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the [MCNX-TRP-02/03] trapped-regime relabeling pure
// map of specs/trapped-regime-treatment.md (T21a), through mcnux_trp.hxx.
// Synthetic pinned (eta, kappa_a, kappa_s) tuples (kappa_a > 0 so the ratio
// invariant is defined) swept over alpha in (0, 1] — no live opacity tables.
// The two invariant rows are machine-tier booleans (detail::approx_eq at
// detail::rtol_machine, the spec's ~1e-14 tier); the monotonicity rows are
// the exact strict-inequality convention tripwire of [MCNX-TRP-03] (the
// flipped 2020-letter alpha convention reverses both); the alpha = 1 row is
// the exact bitwise identity endpoint the [MCNX-TRP-05] limit benchmark
// later builds on. Alpha SELECTION (the xi control rule) is T21b, not here.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_trp.hxx"
#include "mcnux_units.hxx"

namespace MCNuX {

void append_trp_rows(Battery &b) {
  struct Unprimed {
    double eta, kappa_a, kappa_s;
  };
  // Pinned sweep states: exact binary fractions plus generic magnitudes
  // spanning ~13 decades in the coefficients.
  const Unprimed states[] = {
      {1.5, 0.5, 0.25},
      {3.7e-3, 1.2e2, 4.9e-1},
      {6.02e5, 9.1e3, 2.2e4},
      {1.0e-8, 3.3e-6, 7.7e-7},
  };
  // Pinned alpha sweep, strictly increasing within (0, 1].
  const double alphas[] = {0.0625, 0.25, 0.5, 0.75, 0.9, 1.0};
  constexpr int nstates = int(sizeof(states) / sizeof(states[0]));
  constexpr int nalphas = int(sizeof(alphas) / sizeof(alphas[0]));

  // Rows 66..67 — the two exact-consequence invariants of [MCNX-TRP-02],
  // machine tier at every (state, alpha) of the sweep.
  bool ratio_ok = true, total_ok = true;
  for (int s = 0; s < nstates; ++s)
    for (int a = 0; a < nalphas; ++a) {
      const RelabeledCoefficients r =
          relabel(states[s].eta, states[s].kappa_a, states[s].kappa_s,
                  alphas[a]);
      ratio_ok = ratio_ok &&
                 detail::approx_eq(r.eta_p / r.kappa_a_p,
                                   states[s].eta / states[s].kappa_a,
                                   detail::rtol_machine);
      total_ok = total_ok &&
                 detail::approx_eq(r.kappa_a_p + r.kappa_s_p,
                                   states[s].kappa_a + states[s].kappa_s,
                                   detail::rtol_machine);
    }
  b.add_boolean("trp.relabel.eta_over_kappa_a", ratio_ok);
  b.add_boolean("trp.relabel.total_rate", total_ok);

  // Rows 68..69 — the [MCNX-TRP-03] convention tripwire, exact and strict:
  // over every adjacent alpha_1 < alpha_2 pair of the sweep, kappa_a' is
  // strictly increasing and kappa_s' strictly decreasing in alpha. A
  // flipped (2020-letter) implementation reverses both.
  bool ka_ok = true, ks_ok = true;
  for (int s = 0; s < nstates; ++s)
    for (int a = 0; a + 1 < nalphas; ++a) {
      const RelabeledCoefficients lo =
          relabel(states[s].eta, states[s].kappa_a, states[s].kappa_s,
                  alphas[a]);
      const RelabeledCoefficients hi =
          relabel(states[s].eta, states[s].kappa_a, states[s].kappa_s,
                  alphas[a + 1]);
      ka_ok = ka_ok && lo.kappa_a_p < hi.kappa_a_p;
      ks_ok = ks_ok && lo.kappa_s_p > hi.kappa_s_p;
    }
  b.add_boolean("trp.relabel.monotonicity_ka", ka_ok);
  b.add_boolean("trp.relabel.monotonicity_ks", ks_ok);

  // Row 70 — alpha = 1 identity endpoint, bitwise on every state: the
  // relabeling is exactly the identity (the (1 - alpha) kappa_a term
  // vanishes exactly at alpha = 1).
  bool id_ok = true;
  for (int s = 0; s < nstates; ++s) {
    const RelabeledCoefficients r =
        relabel(states[s].eta, states[s].kappa_a, states[s].kappa_s, 1.0);
    id_ok = id_ok && r.eta_p == states[s].eta &&
            r.kappa_a_p == states[s].kappa_a &&
            r.kappa_s_p == states[s].kappa_s;
  }
  b.add_boolean("trp.relabel.alpha_one_identity", id_ok);
}

} // namespace MCNuX
