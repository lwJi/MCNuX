// ---------------------------------------------------------------------------
// Rows 162..165 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the [MCNX-TRP-04] alpha-selection rule of
// specs/trapped-regime-treatment.md:145-171 (T21b), through mcnux_trp.hxx:
// alpha = min(1, xi/(kappa_a Dt_c)) on the UNPRIMED kappa_a, in mutually
// consistent (code) units. Synthetic pinned (kappa_a, Dt_c, xi) tuples — no
// live opacity tables, no grid. The clamp row is exact (alpha == 1.0
// bitwise, including the kappa_a = 0 transparent cell and the exact
// kappa_a Dt_c == xi boundary); the scaling row reproduces xi/(kappa_a Dt_c)
// bitwise (the "exact" tier of verification-suite-design.md:374); the bound
// row composes alpha_select -> relabel and checks the constructive
// invariant kappa_a' Dt_c <= xi over the whole sweep (the ONLY binding
// stability control — the beta-dependent Eq. 50 bound is a deferred open
// question, spec :450-455); the override row pins the [MCNX-TRP-05]
// fixed-alpha path (alpha_fixed > 0 bypasses the rule; alpha_fixed = 1 is
// the bitwise identity endpoint; the -1 sentinel falls through to the
// rule). Relabeling map invariants themselves are rows 66..70
// (mcnux_selftest_trp.cxx), not re-tested here.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_trp.hxx"
#include "mcnux_units.hxx"

namespace MCNuX {

namespace {

// The runtime override composition (mcnux_interactions.cxx /
// mcnux_emission.cxx wiring, restated as a pure expression): a positive
// alpha_fixed replaces the selection rule; the -1 sentinel selects it.
constexpr double effective_alpha(double alpha_fixed, double kappa_a,
                                 double dt_c, double xi) noexcept {
  return (alpha_fixed > 0.0) ? alpha_fixed
                             : alpha_select(kappa_a, dt_c, xi);
}

} // namespace

void append_trp_select_rows(Battery &b) {
  struct SelectFix {
    double kappa_a, dt_c, xi;
  };

  // Clamp-branch sweep: kappa_a Dt_c <= xi, spanning exact binary fractions,
  // generic magnitudes, the transparent cell kappa_a = 0, and the exact
  // boundary kappa_a Dt_c == xi (<= is the clamp side, alpha = 1).
  const SelectFix clamp_states[] = {
      {0.5, 0.25, 1.0},   // well inside, exact fractions
      {0.0, 0.25, 1.0},   // transparent cell: no division, alpha = 1
      {0.0, 3.7e-2, 0.5}, // transparent cell, generic Dt_c
      {2.0, 0.5, 1.0},    // exact boundary kappa_a Dt_c == xi
      {1.3e2, 4.0e-3, 0.75},
  };
  // Scaling-branch sweep: kappa_a Dt_c > xi, exact fractions and generic
  // magnitudes spanning ~8 decades in kappa_a Dt_c / xi.
  const SelectFix scale_states[] = {
      {16.0, 0.5, 1.0},     // x = 8, alpha = 0.125 exactly
      {16.0, 0.5, 2.0},     // xi-scaling: alpha doubles to 0.25
      {3.3e4, 1.7e-2, 1.0},
      {9.1e2, 6.25e-2, 0.5},
      {1.0e6, 2.5e-1, 1.0},
  };
  constexpr int nclamp = int(sizeof(clamp_states) / sizeof(clamp_states[0]));
  constexpr int nscale = int(sizeof(scale_states) / sizeof(scale_states[0]));

  // Row 162 — clamp-at-1 branch, exact: alpha_select returns EXACTLY 1.0
  // (the [MCNX-TRP-05] identity endpoint) whenever kappa_a Dt_c <= xi,
  // including kappa_a = 0 and the equality boundary.
  bool clamp_ok = true;
  for (int i = 0; i < nclamp; ++i)
    clamp_ok = clamp_ok &&
               alpha_select(clamp_states[i].kappa_a, clamp_states[i].dt_c,
                            clamp_states[i].xi) == 1.0;
  b.add_boolean("trp.select.clamp_at_one", clamp_ok);

  // Row 163 — scaling branch, bitwise: alpha_select reproduces
  // xi/(kappa_a Dt_c) exactly (==, the pinned evaluation order xi / (kappa_a
  // * Dt_c)), including the exact-fraction anchors 0.125 and 0.25 of the
  // xi-scaling pair.
  bool scale_ok = true;
  for (int i = 0; i < nscale; ++i) {
    const SelectFix &f = scale_states[i];
    scale_ok = scale_ok &&
               alpha_select(f.kappa_a, f.dt_c, f.xi) ==
                   f.xi / (f.kappa_a * f.dt_c);
  }
  scale_ok = scale_ok && alpha_select(16.0, 0.5, 1.0) == 0.125 &&
             alpha_select(16.0, 0.5, 2.0) == 0.25;
  b.add_boolean("trp.select.xi_scaling_exact", scale_ok);

  // Row 164 — the constructive invariant of the rule ([MCNX-TRP-04], the
  // only binding bound): composing alpha_select -> relabel gives
  // kappa_a' Dt_c <= xi at EVERY state of BOTH sweeps (kappa_s and eta are
  // spectators; map by name, never positionally).
  bool bound_ok = true;
  for (int i = 0; i < nclamp + nscale; ++i) {
    const SelectFix &f = i < nclamp ? clamp_states[i]
                                    : scale_states[i - nclamp];
    const double alpha = alpha_select(f.kappa_a, f.dt_c, f.xi);
    const RelabeledCoefficients rc = relabel(1.5, f.kappa_a, 0.25, alpha);
    bound_ok = bound_ok && rc.kappa_a_p * f.dt_c <= f.xi;
  }
  b.add_boolean("trp.select.kappa_ap_dtc_bound", bound_ok);

  // Row 165 — the fixed-alpha override path ([MCNX-TRP-04] "may replace the
  // selection rule", [MCNX-TRP-05] prep): alpha_fixed = 0.5 bypasses the
  // rule on a state where the rule would pick 0.125 (relabel then yields the
  // exact by-name primed values); alpha_fixed = 1 through relabel is the
  // bitwise identity endpoint; the -1 sentinel falls through to the rule.
  bool override_ok = true;
  {
    // Rule would give 0.125 here; the override must win.
    override_ok = override_ok &&
                  effective_alpha(0.5, 16.0, 0.5, 1.0) == 0.5;
    const RelabeledCoefficients rc = relabel(1.5, 16.0, 0.25, 0.5);
    override_ok = override_ok && rc.eta_p == 0.75 && rc.kappa_a_p == 8.0 &&
                  rc.kappa_s_p == 8.25;
    // alpha_fixed = 1: exact identity through the map (the [MCNX-TRP-05]
    // bitwise-limit endpoint).
    override_ok = override_ok &&
                  effective_alpha(1.0, 16.0, 0.5, 1.0) == 1.0;
    const RelabeledCoefficients ri = relabel(1.5, 16.0, 0.25, 1.0);
    override_ok = override_ok && ri.eta_p == 1.5 && ri.kappa_a_p == 16.0 &&
                  ri.kappa_s_p == 0.25;
    // -1 sentinel: the rule is used (both branches).
    override_ok = override_ok &&
                  effective_alpha(-1.0, 16.0, 0.5, 1.0) == 0.125 &&
                  effective_alpha(-1.0, 0.5, 0.25, 1.0) == 1.0;
  }
  b.add_boolean("trp.select.fixed_alpha_override", override_ok);
}

} // namespace MCNuX
