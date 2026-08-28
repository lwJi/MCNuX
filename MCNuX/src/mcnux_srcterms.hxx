#ifndef MCNUX_SRCTERMS_HXX
#define MCNUX_SRCTERMS_HXX

// Source-term ledger arithmetic of specs/hydro-coupling-source-terms.md
// [MCNX-HYD-02].
//
// Per cell c and transport step, with dP^t, dP_i, dL the net change of total
// packet content in c due to interaction events, the MCNuX-owned source-term
// grid variables hold
//
//   G^t = -dP^t/(dV dt) ,  G_i = -dP_i/(dV dt) ,  S_l = -dL/(dV dt)
//
// with dL = l_s * N per event and the event-contribution table of the spec
// (hydro-coupling-source-terms.md:104-120,133-149). The minus sign makes the
// variables the FLUID-SIDE sources: what packets gain, the fluid loses
// (positive = gained by the fluid). All quantities are geometrized code
// units, per unit coordinate volume per unit coordinate time, undensitized —
// no unit conversion happens here ([MCNX-HYD-02]'s normalization; the
// partner's densitization is its own business, [MCNX-HYD-04]).
//
// The species lepton numbers come from mcnux_units.hxx (lepton_number(),
// [MCNX-CNV-06]); the {+1, -1, 0} table is never re-authored here. The
// degeneracy factor g = 4 enters exactly once, at packet creation — it
// appears NOWHERE in this ledger (N already counts all lumped states,
// hydro-coupling-source-terms.md:119-120).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// constexpr C++ so that host code, device code, and the self-test battery can
// all include it (the mcnux_units.hxx pattern).

#include "mcnux_units.hxx"

namespace MCNuX {

// Net change of total packet content in one cell over one transport step (or
// the contribution of a single event to it): dP^t, dP_i (time component index
// up, spatial components index down, matching the packet state (p^t, p_i)),
// and the lepton-number change dL.
struct LedgerDelta {
  double dPt, dPx, dPy, dPz, dL;
};

// The fluid-side source contribution a LedgerDelta induces: {G^t, G_x, G_y,
// G_z, S_l}, laid out to match the grid variables {rf_t, rf_x, rf_y, rf_z,
// lep_src} of MCNuX/interface.ccl.
struct SourceContribution {
  double Gt, Gx, Gy, Gz, Sl;
};

// The one normative formula, applied channel-wise: source = -delta/(dV dt).
constexpr double source_rate(double delta, double dV, double dt) noexcept {
  return -delta / (dV * dt);
}

constexpr SourceContribution source_from_delta(const LedgerDelta &d, double dV,
                                               double dt) noexcept {
  return {source_rate(d.dPt, dV, dt), source_rate(d.dPx, dV, dt),
          source_rate(d.dPy, dV, dt), source_rate(d.dPz, dV, dt),
          source_rate(d.dL, dV, dt)};
}

// Event-contribution table of hydro-coupling-source-terms.md:110-116.

// Packet creation (emission): +N p^t, +N p_i, +l_s N.
constexpr LedgerDelta emission_delta(Species s, double N, double pt, double px,
                                     double py, double pz) noexcept {
  return {N * pt, N * px, N * py, N * pz,
          static_cast<double>(lepton_number(s)) * N};
}

// Absorption (discrete or diffusion-test): -N p^t, -N p_i, -l_s N.
constexpr LedgerDelta absorption_delta(Species s, double N, double pt,
                                       double px, double py,
                                       double pz) noexcept {
  return {-N * pt, -N * px, -N * py, -N * pz,
          -static_cast<double>(lepton_number(s)) * N};
}

// Elastic scattering / diffusion advection: +N (p_out - p_in) per component,
// zero lepton-number change (species-independent by the table).
constexpr LedgerDelta scattering_delta(double N, double pt_in, double px_in,
                                       double py_in, double pz_in,
                                       double pt_out, double px_out,
                                       double py_out, double pz_out) noexcept {
  return {N * (pt_out - pt_in), N * (px_out - px_in), N * (py_out - py_in),
          N * (pz_out - pz_in), 0.0};
}

// ---------------------------------------------------------------------------
// Conservation-ledger closure audit — [MCNX-HYD-05]
// (hydro-coupling-source-terms.md:177-188)
// ---------------------------------------------------------------------------

// Event-side audit accumulator: per channel X in {P^t, P_x, P_y, P_z, L},
// the signed net packet change Sum_events dX and the gross event activity
// Sum_events |dX| (the closure denominator — quiet steps with near-total
// cancellation are judged against activity, not the net).
//
// The audit MUST be built from LedgerDelta values (the pre-negation
// packet-side deltas of emission_delta/absorption_delta/scattering_delta),
// never from the deposited SourceContribution: deriving the event side from
// the already-negated grid values would trivially self-cancel instead of
// independently verifying the deposit.
struct LedgerAudit {
  LedgerDelta net{0.0, 0.0, 0.0, 0.0, 0.0};   // Sum_events dX (signed)
  LedgerDelta gross{0.0, 0.0, 0.0, 0.0, 0.0}; // Sum_events |dX|

  constexpr void accumulate(const LedgerDelta &d) noexcept {
    net.dPt += d.dPt;
    net.dPx += d.dPx;
    net.dPy += d.dPy;
    net.dPz += d.dPz;
    net.dL += d.dL;
    gross.dPt += detail::cabs(d.dPt);
    gross.dPx += detail::cabs(d.dPx);
    gross.dPy += detail::cabs(d.dPy);
    gross.dPz += detail::cabs(d.dPz);
    gross.dL += detail::cabs(d.dL);
  }

  // Fold another audit in whole (net and gross already event-summed by its
  // owner): channel-wise += on both accumulators. Used by
  // MCNuX_LedgerClosure to merge the per-contributor step audits
  // (emission_step_audit(), later interaction_step_audit()) into its local
  // event-side total. No absolute values here — the folded audit's gross is
  // already Sum_events |dX|.
  constexpr void accumulate(const LedgerAudit &a) noexcept {
    net.dPt += a.net.dPt;
    net.dPx += a.net.dPx;
    net.dPy += a.net.dPy;
    net.dPz += a.net.dPz;
    net.dL += a.net.dL;
    gross.dPt += a.gross.dPt;
    gross.dPx += a.gross.dPx;
    gross.dPy += a.gross.dPy;
    gross.dPz += a.gross.dPz;
    gross.dL += a.gross.dL;
  }
};

// The per-step event-side audit of the production emission loop
// ([MCNX-HYD-05] event side): reset and refilled once per transport step by
// MCNuX_Emission (mcnux_emission.cxx, which defines it) from the SAME
// pre-negation LedgerDelta values it deposits, with the identical code-unit
// dV/dt normalization pair. Consumed by MCNuX_LedgerClosure
// (mcnux_srcterms.cxx), which folds it into its event-side total each step
// (a provable no-op when MCNuX_Emission has not run: the accumulator is
// zero-initialized and only MCNuX_Emission writes it).
LedgerAudit &emission_step_audit();

// The per-step event-side audit of the episode driver ([MCNX-HYD-05] event
// side, interaction events): reset and refilled once per transport step by
// MCNuX_EpisodeDriver (mcnux_interactions.cxx, which defines it) from the
// SAME pre-negation LedgerDelta values it deposits, with the identical
// code-unit dV/dt normalization pair. Separate from emission_step_audit()
// above — the two contributors own their accumulators independently (never
// mutate the other's). Consumed by MCNuX_LedgerClosure (mcnux_srcterms.cxx),
// which folds it into its event-side total each step alongside the emission
// audit (a provable no-op when MCNuX_EpisodeDriver has not run: the
// accumulator is zero-initialized and only MCNuX_EpisodeDriver writes it).
LedgerAudit &interaction_step_audit();

// The [MCNX-HYD-05] closure tolerance: relative 1e-13 per step (the
// specs/README.md conservation tolerance).
inline constexpr double ledger_rtol = 1e-13;

// The documented tiny floor eps0 guarding the all-zero case
// (hydro-coupling-source-terms.md:188 requires the floor be documented):
// pinned to 1e-30, matching the absolute floor of the approx_eq closeness
// predicate (mcnux_units.hxx, conventions-and-units.md:129). With zero gross
// activity the bound becomes ledger_rtol * 1e-30, so an exactly-balanced
// all-zero step passes (residual 0) while any nonzero residual on a zero-
// activity step fails loudly.
inline constexpr double ledger_eps0 = 1e-30;

// One channel's closure verdict. grid_total is Sum_cells (source value) *
// dV * dt for that channel; net and gross come from the LedgerAudit. Since
// the grid holds G = -dX/(dV dt) ([MCNX-HYD-02]), grid_total and net are
// opposite-signed and must cancel: residual = |grid_total + net|.
struct LedgerClosure {
  double residual; // |grid_total + net|
  double bound;    // ledger_rtol * max(gross, ledger_eps0)
  bool pass;       // residual <= bound
};

constexpr LedgerClosure ledger_closure(double grid_total, double net,
                                       double gross) noexcept {
  const double residual = detail::cabs(grid_total + net);
  const double bound =
      ledger_rtol * (gross > ledger_eps0 ? gross : ledger_eps0);
  return {residual, bound, residual <= bound};
}

// ---------------------------------------------------------------------------
// Shared synthetic-event fixture (the `source-zero-add` benchmark)
// ---------------------------------------------------------------------------

// One source of truth for the deterministic, iteration-scaled event list:
// MCNuX_SyntheticDeposit (grid side), MCNuX_LedgerClosure's event-side audit,
// and the ledger selftest rows all evaluate exactly these values — they must
// stay bit-identical or the committed golden data of `source-zero-add` and
// `unit-selftest` moves. dV and dt are fixed synthetic denominators
// (dV dt = 0.5, exact in binary) — NOT any grid's cell size or the run's
// time step; the closure multiply must use the same dV dt as the deposit.
namespace synthfix {

inline constexpr double dV = 2.0;
inline constexpr double dt = 0.25;

// nu_e emission: N = 3 scale, p^t = 2.5, p_i = (0.5, -1.25, 0.75).
constexpr LedgerDelta emission(double scale) noexcept {
  return emission_delta(Species::NuE, 3.0 * scale, 2.5, 0.5, -1.25, 0.75);
}

// nu_e absorption: N = 1.5 scale, p^t = 3.0, p_i = (-0.5, 1.0, 0.25).
constexpr LedgerDelta absorption(double scale) noexcept {
  return absorption_delta(Species::NuE, 1.5 * scale, 3.0, -0.5, 1.0, 0.25);
}

// Elastic scattering: N = 2 scale, p^t unchanged (1.5), momentum rotated
// (1.0, 0.25, -0.5) -> (-0.25, 0.75, 0.5).
constexpr LedgerDelta scattering(double scale) noexcept {
  return scattering_delta(2.0 * scale, 1.5, 1.0, 0.25, -0.5, 1.5, -0.25, 0.75,
                          0.5);
}

// The event-side audit of the whole list at one scale.
constexpr LedgerAudit audit(double scale) noexcept {
  LedgerAudit a{};
  a.accumulate(emission(scale));
  a.accumulate(absorption(scale));
  a.accumulate(scattering(scale));
  return a;
}

} // namespace synthfix

// ---------------------------------------------------------------------------
// Compile-time verification — the [MCNX-HYD-02] sign fixtures
// (hydro-coupling-source-terms.md:143-145,260-262) on exact binary-fraction
// inputs, so every comparison below is exact (no tolerance).
// ---------------------------------------------------------------------------

namespace detail {

// Fixture: a nu_e packet with N = 3, p^t = 2.5, p_i = (0.5, -1.25, 0.75),
// deposited over dV = 2, dt = 0.25 (dV dt = 0.5; all values exact in binary).
inline constexpr double srcfix_N = 3.0;
inline constexpr double srcfix_pt = 2.5;
inline constexpr double srcfix_px = 0.5;
inline constexpr double srcfix_py = -1.25;
inline constexpr double srcfix_pz = 0.75;
inline constexpr double srcfix_dV = 2.0;
inline constexpr double srcfix_dt = 0.25;

inline constexpr SourceContribution srcfix_emission = source_from_delta(
    emission_delta(Species::NuE, srcfix_N, srcfix_pt, srcfix_px, srcfix_py,
                   srcfix_pz),
    srcfix_dV, srcfix_dt);

inline constexpr SourceContribution srcfix_absorption = source_from_delta(
    absorption_delta(Species::NuE, srcfix_N, srcfix_pt, srcfix_px, srcfix_py,
                     srcfix_pz),
    srcfix_dV, srcfix_dt);

// Elastic scattering fixture: N = 2, p^t unchanged (1.5), momentum rotated
// (1.0, 0.25, -0.5) -> (-0.25, 0.75, 0.5).
inline constexpr SourceContribution srcfix_scattering =
    source_from_delta(scattering_delta(2.0, 1.5, 1.0, 0.25, -0.5, 1.5, -0.25,
                                       0.75, 0.5),
                      srcfix_dV, srcfix_dt);

} // namespace detail

// Worked signs of [MCNX-HYD-02]: emission of a nu_e packet gives S_l < 0 and
// G^t < 0 (the fluid loses a lepton and energy); absorption reverses both.
static_assert(detail::srcfix_emission.Sl < 0.0,
              "nu_e emission must yield S_l < 0 (the fluid loses a lepton)");
static_assert(detail::srcfix_emission.Gt < 0.0,
              "nu_e emission must yield G^t < 0 (the fluid loses energy)");
static_assert(detail::srcfix_absorption.Sl > 0.0,
              "nu_e absorption must yield S_l > 0");
static_assert(detail::srcfix_absorption.Gt > 0.0,
              "nu_e absorption must yield G^t > 0");

// Exact values: G = -dP/(dV dt) with dV dt = 0.5, so G^t = -(3*2.5)/0.5 = -15
// and S_l = -(+1*3)/0.5 = -6; absorption is the exact negation.
static_assert(detail::srcfix_emission.Gt == -15.0, "emission G^t value");
static_assert(detail::srcfix_emission.Gx == -3.0, "emission G_x value");
static_assert(detail::srcfix_emission.Gy == 7.5, "emission G_y value");
static_assert(detail::srcfix_emission.Gz == -4.5, "emission G_z value");
static_assert(detail::srcfix_emission.Sl == -6.0, "emission S_l value");
static_assert(detail::srcfix_absorption.Gt == -detail::srcfix_emission.Gt &&
                  detail::srcfix_absorption.Gx == -detail::srcfix_emission.Gx &&
                  detail::srcfix_absorption.Gy == -detail::srcfix_emission.Gy &&
                  detail::srcfix_absorption.Gz == -detail::srcfix_emission.Gz &&
                  detail::srcfix_absorption.Sl == -detail::srcfix_emission.Sl,
              "absorption must be the exact negation of emission");

// Elastic scattering changes momentum components only: G^t and S_l are
// exactly zero, the spatial components are not.
static_assert(detail::srcfix_scattering.Gt == 0.0,
              "elastic scattering must not touch G^t");
static_assert(detail::srcfix_scattering.Sl == 0.0,
              "elastic scattering must not touch S_l");
static_assert(detail::srcfix_scattering.Gx == 5.0 &&
                  detail::srcfix_scattering.Gy == -2.0 &&
                  detail::srcfix_scattering.Gz == -4.0,
              "scattering momentum transfer values");

// nu_x carries zero net lepton number: its emission/absorption never
// contributes to S_l ([MCNX-CNV-06] via lepton_number()).
static_assert(source_from_delta(emission_delta(Species::NuX, 3.0, 2.5, 0.5,
                                               -1.25, 0.75),
                                2.0, 0.25)
                      .Sl == 0.0,
              "nu_x packets never contribute lepton number");

// nu_e_bar emission gains the fluid a lepton: S_l > 0.
static_assert(source_from_delta(emission_delta(Species::NuEBar, 3.0, 2.5, 0.5,
                                               -1.25, 0.75),
                                2.0, 0.25)
                      .Sl > 0.0,
              "nu_e_bar emission must yield S_l > 0");

// The degeneracy factor g appears nowhere in the ledger: for equal N the
// energy/momentum deltas are species-independent (g = 4 would make NuX differ
// from NuE by a factor 4 if it leaked in).
static_assert(emission_delta(Species::NuX, 3.0, 2.5, 0.5, -1.25, 0.75).dPt ==
                  emission_delta(Species::NuE, 3.0, 2.5, 0.5, -1.25, 0.75).dPt,
              "the degeneracy factor g must not enter the ledger");

// ---------------------------------------------------------------------------
// Compile-time verification — the [MCNX-HYD-05] closure audit on the shared
// synthetic-event fixture list at scale 1 (exact binary-fraction inputs, so
// every comparison below is exact).
// ---------------------------------------------------------------------------

namespace detail {

// Event side: net and gross of the whole synthetic list at scale 1.
inline constexpr LedgerAudit audfix = synthfix::audit(1.0);

// Grid side of the same list: each event deposits into exactly one cell, so
// Sum_cells (source value) * dV * dt is the sum of the three deposited
// SourceContribution values times dV dt.
inline constexpr SourceContribution audfix_emis =
    source_from_delta(synthfix::emission(1.0), synthfix::dV, synthfix::dt);
inline constexpr SourceContribution audfix_absn =
    source_from_delta(synthfix::absorption(1.0), synthfix::dV, synthfix::dt);
inline constexpr SourceContribution audfix_scat =
    source_from_delta(synthfix::scattering(1.0), synthfix::dV, synthfix::dt);
inline constexpr double audfix_dVdt = synthfix::dV * synthfix::dt;

inline constexpr double audfix_grid_Pt =
    (audfix_emis.Gt + audfix_absn.Gt + audfix_scat.Gt) * audfix_dVdt;
inline constexpr double audfix_grid_Px =
    (audfix_emis.Gx + audfix_absn.Gx + audfix_scat.Gx) * audfix_dVdt;
inline constexpr double audfix_grid_Py =
    (audfix_emis.Gy + audfix_absn.Gy + audfix_scat.Gy) * audfix_dVdt;
inline constexpr double audfix_grid_Pz =
    (audfix_emis.Gz + audfix_absn.Gz + audfix_scat.Gz) * audfix_dVdt;
inline constexpr double audfix_grid_L =
    (audfix_emis.Sl + audfix_absn.Sl + audfix_scat.Sl) * audfix_dVdt;

} // namespace detail

// Exact event-side values at scale 1: emission d = (7.5, 1.5, -3.75, 2.25,
// +3), absorption d = (-4.5, 0.75, -1.5, -0.375, -1.5), scattering
// d = (0, -2.5, 1, 2, 0).
static_assert(detail::audfix.net.dPt == 3.0 && detail::audfix.net.dPx == -0.25 &&
                  detail::audfix.net.dPy == -4.25 &&
                  detail::audfix.net.dPz == 3.875 &&
                  detail::audfix.net.dL == 1.5,
              "synthetic-list audit net values");
static_assert(detail::audfix.gross.dPt == 12.0 &&
                  detail::audfix.gross.dPx == 4.75 &&
                  detail::audfix.gross.dPy == 6.25 &&
                  detail::audfix.gross.dPz == 4.625 &&
                  detail::audfix.gross.dL == 4.5,
              "synthetic-list audit gross activity values");

// The worked identity: Sum_cells G * dV * dt = -dX, so grid total + net = 0
// exactly on every channel (dV dt = 0.5 is a power of two, so the divide/
// multiply round-trips bitwise).
static_assert(detail::audfix_grid_Pt + detail::audfix.net.dPt == 0.0 &&
                  detail::audfix_grid_Px + detail::audfix.net.dPx == 0.0 &&
                  detail::audfix_grid_Py + detail::audfix.net.dPy == 0.0 &&
                  detail::audfix_grid_Pz + detail::audfix.net.dPz == 0.0 &&
                  detail::audfix_grid_L + detail::audfix.net.dL == 0.0,
              "grid total and event net must cancel exactly on the fixture");

// The closure verdicts: the consistent fixture passes with residual exactly
// 0; a doubled grid total (the signature of a missing zero point) fails; the
// all-zero case is judged against the documented eps0 floor and passes.
static_assert(ledger_closure(detail::audfix_grid_Pt, detail::audfix.net.dPt,
                             detail::audfix.gross.dPt)
                      .pass &&
                  ledger_closure(detail::audfix_grid_Pt,
                                 detail::audfix.net.dPt,
                                 detail::audfix.gross.dPt)
                          .residual == 0.0,
              "consistent fixture must close with residual exactly 0");
static_assert(!ledger_closure(2.0 * detail::audfix_grid_Pt,
                              detail::audfix.net.dPt, detail::audfix.gross.dPt)
                   .pass,
              "a doubled grid total must fail closure");
static_assert(ledger_closure(0.0, 0.0, 0.0).pass &&
                  ledger_closure(0.0, 0.0, 0.0).bound ==
                      ledger_rtol * ledger_eps0,
              "the all-zero case must pass via the documented eps0 floor");

} // namespace MCNuX

#endif // MCNUX_SRCTERMS_HXX
