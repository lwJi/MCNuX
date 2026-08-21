#ifndef MCNUX_TRP_HXX
#define MCNUX_TRP_HXX

// Trapped-regime opacity relabeling — the pure post-map of
// specs/trapped-regime-treatment.md [MCNX-TRP-02] (provenance
// arXiv:2103.16588v2 Eq. 43, the BINDING 2021 alpha convention):
//
//   eta' = alpha * eta
//   kappa_a' = alpha * kappa_a
//   kappa_s' = kappa_s + (1 - alpha) * kappa_a ,       alpha in (0, 1]
//
// applied per (cell, species, energy bin) to the unprimed coefficients of
// specs/opacity-eos-evaluation.md, and substituted for (eta, kappa_a,
// kappa_s) in the emission count law and the event-sampling law — and
// nowhere else (tallies remain the actual packet events).
//
// [MCNX-TRP-03] ALPHA-CONVENTION GUARD — the corpus's highest-consequence
// trap: the 2020 letter arXiv:2008.08089 (its Eq. 4) uses the SAME symbol
// alpha with the OPPOSITE meaning — the letter's alpha is this corpus's
// 1 - alpha. This file (and every alpha-dependent quantity in MCNuX) pins
// the 2021 convention: alpha is the RETAINED fraction of absorption, so
// alpha -> 1 is the explicit scheme (no relabeling) and increasing alpha
// means LESS relabeling. Checkable monotonicity signature (asserted below):
// for fixed unprimed coefficients and alpha_1 < alpha_2 in (0, 1], kappa_a'
// is strictly increasing and kappa_s' strictly decreasing in alpha — a
// flipped (2020-letter) implementation reverses both and fails the asserts.
//
// Two identities are exact consequences of the map and must hold to the
// machine tier (~1e-14 relative) at every evaluated state ([MCNX-TRP-02]):
//
//   eta'/kappa_a' = eta/kappa_a          (equilibrium energy density
//                                         unmodified; kappa_a > 0)
//   kappa_a' + kappa_s' = kappa_a + kappa_s   (total interaction rate /
//                                              diffusion timescale unmodified)
//
// Out of scope here (later tasks): the per-cell alpha SELECTION rule
// alpha = min(1, xi/(kappa_a Dt_c)) and the kappa_a' Dt_c <= xi diagnostic
// ([MCNX-TRP-04] control plumbing), and the alpha -> 1 bitwise-limit
// benchmark ([MCNX-TRP-05]). The alpha in (0, 1] range is the corpus's
// recorded assumption (trapped-regime-treatment.md, "Open questions");
// enforcement is paramcheck's job, not this pure map's.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// constexpr C++ so that host code, device code, and the self-test battery
// can all include it (the mcnux_units.hxx pattern).

#include "mcnux_units.hxx"

namespace MCNuX {

// The relabeled (primed) coefficient triple, in the same units as the
// unprimed inputs (the map is dimensionless in alpha).
struct RelabeledCoefficients {
  double eta_p, kappa_a_p, kappa_s_p;
};

// The one normative map, [MCNX-TRP-02] / arXiv:2103.16588 Eq. 43 (2021
// convention). Pure function of its arguments; alpha in (0, 1] is assumed,
// never checked here.
constexpr RelabeledCoefficients relabel(double eta, double kappa_a,
                                        double kappa_s,
                                        double alpha) noexcept {
  return {alpha * eta, alpha * kappa_a, kappa_s + (1.0 - alpha) * kappa_a};
}

// ---------------------------------------------------------------------------
// Compile-time verification — the [MCNX-TRP-02] invariants, the
// [MCNX-TRP-03] monotonicity tripwire, and the alpha = 1 identity endpoint,
// on exact binary-fraction fixtures so every non-tolerance comparison below
// is exact.
// ---------------------------------------------------------------------------

namespace detail {

// Fixture: eta = 1.5, kappa_a = 0.5, kappa_s = 0.25 (all exact in binary;
// kappa_a > 0 so the ratio invariant is well defined).
inline constexpr double trpfix_eta = 1.5;
inline constexpr double trpfix_ka = 0.5;
inline constexpr double trpfix_ks = 0.25;

inline constexpr RelabeledCoefficients trpfix_a25 =
    relabel(trpfix_eta, trpfix_ka, trpfix_ks, 0.25);
inline constexpr RelabeledCoefficients trpfix_a75 =
    relabel(trpfix_eta, trpfix_ka, trpfix_ks, 0.75);
inline constexpr RelabeledCoefficients trpfix_a100 =
    relabel(trpfix_eta, trpfix_ka, trpfix_ks, 1.0);

} // namespace detail

// Exact map values at alpha = 0.25: eta' = 0.375, kappa_a' = 0.125,
// kappa_s' = 0.25 + 0.75*0.5 = 0.625 (all exact binary fractions).
static_assert(detail::trpfix_a25.eta_p == 0.375 &&
                  detail::trpfix_a25.kappa_a_p == 0.125 &&
                  detail::trpfix_a25.kappa_s_p == 0.625,
              "[MCNX-TRP-02] relabeling map values at alpha = 0.25");

// Invariant 1: eta'/kappa_a' = eta/kappa_a (equilibrium energy density
// unmodified), machine tier.
static_assert(detail::approx_eq(detail::trpfix_a25.eta_p /
                                    detail::trpfix_a25.kappa_a_p,
                                detail::trpfix_eta / detail::trpfix_ka,
                                detail::rtol_machine) &&
                  detail::approx_eq(detail::trpfix_a75.eta_p /
                                        detail::trpfix_a75.kappa_a_p,
                                    detail::trpfix_eta / detail::trpfix_ka,
                                    detail::rtol_machine),
              "[MCNX-TRP-02] eta'/kappa_a' must equal eta/kappa_a");

// Invariant 2: kappa_a' + kappa_s' = kappa_a + kappa_s (total interaction
// rate unmodified), machine tier.
static_assert(detail::approx_eq(detail::trpfix_a25.kappa_a_p +
                                    detail::trpfix_a25.kappa_s_p,
                                detail::trpfix_ka + detail::trpfix_ks,
                                detail::rtol_machine) &&
                  detail::approx_eq(detail::trpfix_a75.kappa_a_p +
                                        detail::trpfix_a75.kappa_s_p,
                                    detail::trpfix_ka + detail::trpfix_ks,
                                    detail::rtol_machine),
              "[MCNX-TRP-02] kappa_a' + kappa_s' must equal kappa_a + kappa_s");

// [MCNX-TRP-03] monotonicity tripwire (exact, strict): alpha_1 = 0.25 <
// alpha_2 = 0.75 gives kappa_a' strictly increasing and kappa_s' strictly
// decreasing. The flipped 2020-letter convention reverses BOTH and fails
// here.
static_assert(detail::trpfix_a25.kappa_a_p < detail::trpfix_a75.kappa_a_p,
              "[MCNX-TRP-03] kappa_a' must be strictly increasing in alpha "
              "(2021 convention; a flipped 2020-letter map fails this)");
static_assert(detail::trpfix_a25.kappa_s_p > detail::trpfix_a75.kappa_s_p,
              "[MCNX-TRP-03] kappa_s' must be strictly decreasing in alpha "
              "(2021 convention; a flipped 2020-letter map fails this)");

// alpha = 1 identity endpoint, exact: the relabeling is the identity
// (eta' = eta, kappa_a' = kappa_a, kappa_s' = kappa_s, bitwise — the
// (1 - alpha) kappa_a term vanishes exactly). The [MCNX-TRP-05] bitwise-run
// limit benchmark builds on this.
static_assert(detail::trpfix_a100.eta_p == detail::trpfix_eta &&
                  detail::trpfix_a100.kappa_a_p == detail::trpfix_ka &&
                  detail::trpfix_a100.kappa_s_p == detail::trpfix_ks,
              "alpha = 1 must reproduce the unprimed coefficients exactly");

} // namespace MCNuX

#endif // MCNUX_TRP_HXX
