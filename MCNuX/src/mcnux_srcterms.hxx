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

} // namespace MCNuX

#endif // MCNUX_SRCTERMS_HXX
