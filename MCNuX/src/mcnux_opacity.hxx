#ifndef MCNUX_OPACITY_HXX
#define MCNUX_OPACITY_HXX

// The WeakLibInterp call-boundary wrapper layer: one thin wrapper per WLI
// `_Point` entry point of the five families (EOS evaluate, EOS inversion,
// EmAb, Iso, NES/Pair, Brem), pinning MCNuX's argument/units/log10/index
// contract at the boundary.
//
// Governing spec: specs/opacity-eos-evaluation.md, requirements
// [MCNX-OPA-01] (complete, exact consumption surface — wrap the `_Point`
// kernels, never reimplement the interpolation; family table at
// specs/opacity-eos-evaluation.md:136-157), [MCNX-OPA-02] (Kelvin at EVERY
// temperature argument, including NES/Pair whose table axis is Kelvin per
// WeakLibInterp/specs/opacity-nes-pair.md), [MCNX-OPA-03] (log10 ownership:
// the EOS families receive raw coordinates, the four opacity families receive
// coordinates already log10'd by MCNuX; Ye always raw/linear), and the
// species-axis grounding of [MCNX-OPA-08].
//
// Unit contract at this boundary (specs/opacity-eos-evaluation.md:159-170):
//   * Temperature is MeV everywhere inside MCNuX; every wrapper below takes
//     T in MeV and applies MCNuX::mev_to_kelvin (mcnux_units.hxx) exactly
//     once per temperature argument. The factor's literal lives only in
//     mcnux_units.hxx (single-literal sweep, [MCNX-CNV-03]).
//   * Energies in MeV, rho in g/cm^3, Ye dimensionless. log10 is taken of
//     the value in exactly these units (Kelvin for temperature axes).
//   * All wrapper internals and WLI calls are wli::Real (= double
//     unconditionally, wli_real.H); the boundary static_asserts below pin
//     that no narrowing can hide here ([MCNX-BLD-03] context).
//
// Index contract ([MCNX-OPA-01]): all integer indices (iSpecies, iMom, iEp,
// iE, kernel/moment) are zero-based and pass through unreinterpreted; the
// energy and moment axes are discrete and never interpolated; every
// [nOpacities, nMoments] offset lookup goes through wli::IsoOffset
// (species-major, iSpecies + nOpacities*iMom) — offset indexing is never
// hand-rolled here.
//
// Scope boundary (deliberate):
//   * NO range policy lives here. WeakLibInterp clamps only the bracket
//     index, never the delta — out-of-range queries silently extrapolate
//     from the edge cell, and a non-positive value reaching a log10 makes a
//     silent NaN. The [MCNX-OPA-06] transparency floor / clamping / clamp
//     counters that prevent this are the T11 range-policy layer; this header
//     only provides the minimal AxisBounds surface T11 builds on.
//   * NO kappa/eta assembly, no Kirchhoff closure, no nu_x mapping
//     ([MCNX-OPA-04/05], T10), and no degeneracy factor anywhere (g = 4
//     enters exactly once at packet creation, never in per-species transport
//     coefficients).
//   * NO eta_e assembly helper (the baseline never calls NES/Pair; the EOS
//     evaluate wrapper suffices to fetch mu_e when a consumer needs eta_e).
//
// Species-axis grounding [MCNX-OPA-08]: the production EmAb/Iso tables carry
// exactly TWO species datasets — `Electron Neutrino` (slot 0) and `Electron
// Antineutrino` (slot 1), Fortran shape (40,185,81,30) each — and no
// heavy-lepton dataset; see
// WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls:22-27 and
// specs/conventions-and-units.md:131,148. The structural constants below pin
// that mapping; live-fixture drift detection stays with
// specs/tools/validate_specs.sh:280-287.
//
// WLI kernels are not constexpr, so the kernel-parity checks are runtime
// rows of the unit-selftest battery (mcnux_selftest.cxx); only the
// compile-time-checkable structure is static_asserted here.

#include "mcnux_units.hxx"

#include "wli_eos.H"
#include "wli_eos_inversion.H"
#include "wli_opacity_brem.H"
#include "wli_opacity_emab_iso.H"
#include "wli_opacity_nes_pair.H"
#include "wli_real.H"

#include <cmath>
#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Type boundary  ([MCNX-BLD-03] context)
// ---------------------------------------------------------------------------
// wli::Real is double unconditionally (wli_real.H) and MCNuX's constants are
// double; the asserts make any future divergence a build failure rather than
// an implicit conversion.
static_assert(std::is_same_v<wli::Real, double>,
              "the WLI boundary requires wli::Real == double");
static_assert(sizeof(wli::Real) == 8, "wli::Real must be binary64");

// ---------------------------------------------------------------------------
// Temperature unit boundary  [MCNX-OPA-02]
// ---------------------------------------------------------------------------
// The ONLY place a wrapper converts a temperature. Applied exactly once per
// temperature argument of every family — including NES/Pair, whose
// temperature axis is Kelvin (assuming MeV there is the documented
// easy-to-make error, specs/opacity-eos-evaluation.md:161-165).

// T[K] = T[MeV] * mev_to_kelvin (raw Kelvin, for the EOS families that log
// internally).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
temperature_kelvin_from_mev(wli::Real T_MeV) noexcept {
  return T_MeV * wli::Real(mev_to_kelvin);
}

// T[MeV] = T[K] / mev_to_kelvin — the exact inverse operation, used to map an
// inversion result (Kelvin) back into MCNuX's internal MeV. 0/f == 0, so the
// inversion-failure signal T == 0 survives the conversion unchanged.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
temperature_mev_from_kelvin(wli::Real T_K) noexcept {
  return T_K / wli::Real(mev_to_kelvin);
}

// log10(T[K]) for the pre-logged opacity families ([MCNX-OPA-03]: MCNuX owns
// the log10 of every opacity-family coordinate).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
log10_temperature_kelvin_from_mev(wli::Real T_MeV) noexcept {
  return std::log10(temperature_kelvin_from_mev(T_MeV));
}

// ---------------------------------------------------------------------------
// EOS evaluate  (library owns the log10s; raw rho, T[K], Ye)
// ---------------------------------------------------------------------------
// Ds/Ts/Ys are the table's raw physical axis nodes (Ts in KELVIN — table
// native); table is the flat column-major log-stored sub-table
// (wli_eos.H:43-57). The wrapper's sole job is the MeV->K conversion.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real eos_evaluate(
    wli::Real rho, wli::Real T_MeV, wli::Real Ye, wli::Real const *Ds, int nD,
    wli::Real const *Ts, int nT, wli::Real const *Ys, int nY, wli::Real OS,
    wli::Real const *table) noexcept {
  return wli::EosInterpolateSingleVariable3DPoint(
      rho, temperature_kelvin_from_mev(T_MeV), Ye, Ds, nD, Ts, nT, Ys, nY, OS,
      table);
}

// ---------------------------------------------------------------------------
// EOS inversion  (raw rho, X, Ye; result mapped K -> MeV; [MCNX-OPA-07])
// ---------------------------------------------------------------------------
// Thin pass-throughs over the six wli::ComputeTemperatureWith_* wrappers
// (wli_eos_inversion.H:265-355). The full {T, Error} result is surfaced —
// error != 0 means T == 0 is NOT a temperature ([MCNX-OPA-07]); no error
// policy is applied here (that protocol layer is T11's concern). The
// EosInversionBounds struct passes through untouched. The returned table
// temperature (Kelvin) is mapped to MCNuX-internal MeV; the _Guess variants'
// T_guess_MeV argument is the one temperature argument and receives the
// single MeV->K conversion.

struct EosInversionResultMeV {
  wli::Real T_MeV; // recovered temperature in MeV; 0 whenever error != 0
  int error;       // WLI error code, forwarded unreinterpreted
};

namespace detail {
// Shared K->MeV result mapping (keeps the six wrappers textually identical).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
inversion_result_mev(wli::EosInversionResult const &r) noexcept {
  return EosInversionResultMeV{temperature_mev_from_kelvin(r.T), r.Error};
}
} // namespace detail

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dey_noguess(wli::Real rho, wli::Real E, wli::Real Ye,
                       wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                       wli::Real const *Ys, int nY, wli::Real OS,
                       wli::Real const *Es,
                       wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DEY_NoGuess(
      rho, E, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Es, b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dey_guess(wli::Real rho, wli::Real E, wli::Real Ye,
                     wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                     wli::Real const *Ys, int nY, wli::Real OS,
                     wli::Real const *Es, wli::Real T_guess_MeV,
                     wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DEY_Guess(
      rho, E, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Es,
      temperature_kelvin_from_mev(T_guess_MeV), b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dpy_noguess(wli::Real rho, wli::Real P, wli::Real Ye,
                       wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                       wli::Real const *Ys, int nY, wli::Real OS,
                       wli::Real const *Ps,
                       wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DPY_NoGuess(
      rho, P, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Ps, b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dpy_guess(wli::Real rho, wli::Real P, wli::Real Ye,
                     wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                     wli::Real const *Ys, int nY, wli::Real OS,
                     wli::Real const *Ps, wli::Real T_guess_MeV,
                     wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DPY_Guess(
      rho, P, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Ps,
      temperature_kelvin_from_mev(T_guess_MeV), b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dsy_noguess(wli::Real rho, wli::Real S, wli::Real Ye,
                       wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                       wli::Real const *Ys, int nY, wli::Real OS,
                       wli::Real const *Ss,
                       wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DSY_NoGuess(
      rho, S, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Ss, b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE EosInversionResultMeV
eos_invert_dsy_guess(wli::Real rho, wli::Real S, wli::Real Ye,
                     wli::Real const *Ds, int nD, wli::Real const *Ts, int nT,
                     wli::Real const *Ys, int nY, wli::Real OS,
                     wli::Real const *Ss, wli::Real T_guess_MeV,
                     wli::EosInversionBounds const &b) noexcept {
  return detail::inversion_result_mev(wli::ComputeTemperatureWith_DSY_Guess(
      rho, S, Ye, Ds, nD, Ts, nT, Ys, nY, OS, Ss,
      temperature_kelvin_from_mev(T_guess_MeV), b));
}

// ---------------------------------------------------------------------------
// Species-axis grounding  [MCNX-OPA-08]
// ---------------------------------------------------------------------------
// The EmAb/Iso species axis of the production tables holds exactly TWO
// electron-type dataset slots, in the binding order nu_e = 0, nu_e_bar = 1
// (WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls:22-27:
// exactly `Electron Antineutrino` and `Electron Neutrino`, Fortran shape
// (40,185,81,30) each; specs/conventions-and-units.md:131,148). nu_x has NO
// dataset in either table (its mapping is the T10 [MCNX-OPA-05] layer, not a
// table slot). These constants are the hardcoded structural mirror; drift of
// the live fixtures is detected by specs/tools/validate_specs.sh:280-287.

inline constexpr int num_tabulated_species = 2; // {nu_e, nu_e_bar}, no nu_x
inline constexpr int no_table_slot = -1;        // "this species has no dataset"

constexpr int emab_species_slot(Species s) noexcept {
  return s == Species::NuE      ? 0
         : s == Species::NuEBar ? 1
                                : no_table_slot;
}

constexpr int iso_species_slot(Species s) noexcept {
  return s == Species::NuE      ? 0
         : s == Species::NuEBar ? 1
                                : no_table_slot;
}

namespace detail {

constexpr bool species_axis_holds_emab() noexcept {
  return num_tabulated_species == 2 &&
         emab_species_slot(Species::NuE) == 0 &&
         emab_species_slot(Species::NuEBar) == 1 &&
         emab_species_slot(Species::NuX) == no_table_slot && no_table_slot < 0;
}

constexpr bool species_axis_holds_iso() noexcept {
  return num_tabulated_species == 2 && iso_species_slot(Species::NuE) == 0 &&
         iso_species_slot(Species::NuEBar) == 1 &&
         iso_species_slot(Species::NuX) == no_table_slot && no_table_slot < 0;
}

} // namespace detail

static_assert(detail::species_axis_holds_emab() &&
                  detail::species_axis_holds_iso(),
              "[MCNX-OPA-08] exactly two electron-type EmAb/Iso table slots, "
              "ordered nu_e = 0, nu_e_bar = 1, and no nu_x dataset");

// ---------------------------------------------------------------------------
// EmAb / Iso  (caller-side log10 of E, rho, T[K]; Ye raw; [MCNX-OPA-03])
// ---------------------------------------------------------------------------
// LogEs/LogDs/LogTs are the table's already-log10'd axis-node grids (log10 of
// MeV, g/cm^3, Kelvin); Ys is raw. The wrapper computes the three log10
// query coordinates in exactly those units; Ye is never logged. OS is the
// scalar per-dataset offset the caller pre-selected (EmAb: 1D Offsets by
// species slot; Iso: via iso_offset below).

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real emab_evaluate(
    wli::Real E_MeV, wli::Real rho, wli::Real T_MeV, wli::Real Ye,
    wli::Real const *LogEs, int nE, wli::Real const *LogDs, int nD,
    wli::Real const *LogTs, int nT, wli::Real const *Ys, int nY, wli::Real OS,
    wli::Real const *table) noexcept {
  return wli::EmAbInterpolateSingleVariable4DPoint(
      std::log10(E_MeV), std::log10(rho),
      log10_temperature_kelvin_from_mev(T_MeV), Ye, LogEs, nE, LogDs, nD,
      LogTs, nT, Ys, nY, OS, table);
}

// The one sanctioned [nOpacities, nMoments] offset lookup (species-major,
// iSpecies + nOpacities*iMom): forwards to wli::IsoOffset, never hand-rolled.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
iso_offset(wli::Real const *offsets, int nOpacities, int nMoments,
           int iSpecies, int iMom) noexcept {
  return wli::IsoOffset(offsets, nOpacities, nMoments, iSpecies, iMom);
}

// Iso at a FIXED moment index iMom (discrete, never interpolated); OS is the
// scalar the caller selected via iso_offset for the same (species, moment).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real iso_evaluate(
    wli::Real E_MeV, wli::Real rho, wli::Real T_MeV, wli::Real Ye,
    wli::Real const *LogEs, int nE, wli::Real const *LogDs, int nD,
    wli::Real const *LogTs, int nT, wli::Real const *Ys, int nY, int iMom,
    int nMom, wli::Real OS, wli::Real const *table) noexcept {
  return wli::IsoInterpolateSingleVariable5DPoint(
      std::log10(E_MeV), std::log10(rho),
      log10_temperature_kelvin_from_mev(T_MeV), Ye, LogEs, nE, LogDs, nD,
      LogTs, nT, Ys, nY, iMom, nMom, OS, table);
}

// ---------------------------------------------------------------------------
// NES / Pair  (caller-side log10 of T[K] and eta_e; Kelvin axis!)
// ---------------------------------------------------------------------------
// LogTs is the table's log10-KELVIN temperature grid
// (WeakLibInterp/specs/opacity-nes-pair.md pins the Kelvin axis; the MeV
// assumption is the documented easy-to-make error), LogEtas the log10 grid of
// the dimensionless eta_e = mu_e/(k_B T). (iEp, iE, kernel) are discrete
// pass-through indices into the (nEp, nE, nMom, nT, nEta) table.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real nes_pair_evaluate(
    wli::Real T_MeV, wli::Real eta, wli::Real const *LogTs, int nT,
    wli::Real const *LogEtas, int nEta, int iEp, int iE, int nEp, int nE,
    int kernel, int nMom, wli::Real OS, wli::Real const *table) noexcept {
  return wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
      log10_temperature_kelvin_from_mev(T_MeV), std::log10(eta), LogTs, nT,
      LogEtas, nEta, iEp, iE, nEp, nE, kernel, nMom, OS, table);
}

// NES detailed-balance fill (upper triangle iEp >= iE). CRITICAL
// dual-temperature split (specs/opacity-eos-evaluation.md:142,161-165;
// wli_opacity_nes_pair.H:175-179,201-203): the underlying kernel takes TWO
// temperature slots in DIFFERENT units — the log10-KELVIN table coordinate
// AND the bare energy-unit (MeV) temperature used unlogged in the Boltzmann
// factor exp((E[iE]-E[iEp])/T). The wrapper takes the single MCNuX-internal
// T_MeV and derives both slots as two DISTINCT variables so the converted
// Kelvin value can never be reused for the Boltzmann slot (or vice versa).
// E_MeV is the physical (unlogged) tabulated energy grid in MeV, strictly
// increasing.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real nes_detailed_balance_fill(
    wli::Real T_MeV, wli::Real eta, wli::Real const *LogTs, int nT,
    wli::Real const *LogEtas, int nEta, int iEp, int iE, int nEp, int nE,
    int kernel, int nMom, wli::Real OS, wli::Real const *table,
    wli::Real const *E_MeV) noexcept {
  // Slot 1: the Kelvin table axis, log10'd — MeV->K applied exactly once.
  const wli::Real log10_T_kelvin = log10_temperature_kelvin_from_mev(T_MeV);
  // Slot 2: the Boltzmann temperature, unlogged, in the SAME energy unit as
  // E_MeV (MeV) — deliberately NOT converted.
  const wli::Real T_boltzmann_MeV = T_MeV;
  return wli::NESDetailedBalanceFillPoint(log10_T_kelvin, std::log10(eta),
                                          LogTs, nT, LogEtas, nEta, iEp, iE,
                                          nEp, nE, kernel, nMom, OS, table,
                                          E_MeV, T_boltzmann_MeV);
}

// Pair crossing-symmetry fill (upper triangle iEp > iE only; the diagonal is
// NOT a fixed point). Exact relabeling — no Boltzmann factor, no E/T beyond
// the table coordinate. The caller supplies the already-swapped kernel index
// and the OS matching THAT swapped kernel.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real pair_crossing_fill(
    wli::Real T_MeV, wli::Real eta, wli::Real const *LogTs, int nT,
    wli::Real const *LogEtas, int nEta, int iEp, int iE, int nEp, int nE,
    int kernelSwapped, int nMom, wli::Real OS,
    wli::Real const *table) noexcept {
  return wli::PairCrossingSymmetryFillPoint(
      log10_temperature_kelvin_from_mev(T_MeV), std::log10(eta), LogTs, nT,
      LogEtas, nEta, iEp, iE, nEp, nE, kernelSwapped, nMom, OS, table);
}

// ---------------------------------------------------------------------------
// Brem  (caller-side log10 of the 3 effective densities and T[K])
// ---------------------------------------------------------------------------
// CRITICAL axis order (the copy-paste trap, wli_opacity_brem.H:45-53): the
// Brem thermodynamic plane is (rho, T) — density BEFORE temperature, the
// TRANSPOSE of NES/Pair's (T, eta) — with the 5D extents (nEp,nE,nMom,nD,nT).
// The three effective densities rho*xp, rho*xn, rho*sqrt(xp*xn) carry the
// pinned recovered-value weights Alpha = [1, 1, 28/3]
// (specs/opacity-eos-evaluation.md:148-151).

inline constexpr int brem_num_effective_densities = 3;

// Alpha(l) for l in {0, 1, 2} = [1, 1, 28/3]; constexpr so device wrappers
// fold it into literals (no namespace-scope array crosses the GPU boundary).
constexpr wli::Real brem_alpha(int l) noexcept {
  return l == 2 ? wli::Real(28) / wli::Real(3) : wli::Real(1);
}

// Summed Brem evaluate: raw effective densities [g/cm^3] in, log10'd
// MCNuX-side ([MCNX-OPA-03]); one MeV->K conversion for the single shared
// temperature; Alpha applied by the WLI kernel to RECOVERED values.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real brem_evaluate_summed(
    wli::Real const (&rho_eff)[brem_num_effective_densities], wli::Real T_MeV,
    wli::Real const *LogDs, int nD, wli::Real const *LogTs, int nT, int iEp,
    int iE, int nEp, int nE, int moment, int nMom, wli::Real OS,
    wli::Real const *table) noexcept {
  wli::Real LogD[brem_num_effective_densities];
  for (int l = 0; l < brem_num_effective_densities; ++l)
    LogD[l] = std::log10(rho_eff[l]);
  const wli::Real Alpha[brem_num_effective_densities] = {
      brem_alpha(0), brem_alpha(1), brem_alpha(2)};
  return wli::BremInterpolateSingleVariable2D2DAlignedSummedPoint(
      LogD, Alpha, brem_num_effective_densities,
      log10_temperature_kelvin_from_mev(T_MeV), LogDs, nD, LogTs, nT, iEp, iE,
      nEp, nE, moment, nMom, OS, table);
}

static_assert(brem_alpha(0) == wli::Real(1) && brem_alpha(1) == wli::Real(1) &&
                  brem_alpha(2) == wli::Real(28) / wli::Real(3),
              "Brem recovered-value weights must be Alpha = [1, 1, 28/3]");

// ---------------------------------------------------------------------------
// Table-range awareness hooks  (bounds surface only — NOT [MCNX-OPA-06])
// ---------------------------------------------------------------------------
// Only the bounds-reading surface the range-policy layer builds on: no
// clamping, no transparency floor, no clamp counters, no NaN guards live
// here — those are mcnux_table_range.hxx ([MCNX-OPA-06], T11). Axis nodes
// are strictly ascending, so the endpoints ARE the bounds.

struct AxisBounds {
  wli::Real lo, hi;
};

static_assert(std::is_trivially_copyable_v<AxisBounds>,
              "AxisBounds must be trivially copyable (device-capturable)");

// Endpoints of a strictly-ascending axis-node array, in the axis's own
// stored units (raw or log10 — whatever the table stores).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
axis_bounds(wli::Real const *xs, int n) noexcept {
  return AxisBounds{xs[0], xs[n - 1]};
}

// Temperature-axis variant for RAW Kelvin node arrays (the EOS families),
// reported in MCNuX-internal MeV so range comparisons happen in MCNuX units.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
axis_bounds_temperature_mev(wli::Real const *Ts_kelvin, int n) noexcept {
  return AxisBounds{temperature_mev_from_kelvin(Ts_kelvin[0]),
                    temperature_mev_from_kelvin(Ts_kelvin[n - 1])};
}

} // namespace MCNuX

#endif // MCNUX_OPACITY_HXX
