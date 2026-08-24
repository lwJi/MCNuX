#ifndef MCNUX_TABLE_RANGE_HXX
#define MCNUX_TABLE_RANGE_HXX

// The table-range enforcement policy of specs/opacity-eos-evaluation.md
// [MCNX-OPA-06] (T11) — WeakLibInterp clamps only the bracket index, never
// the interpolation delta, so out-of-range queries extrapolate linearly from
// the edge cell, silently; range enforcement is MCNuX's responsibility.
// This header is the range-policy layer that the scope comments of
// mcnux_opacity.hxx (AxisBounds hooks) and mcnux_table_coeffs.hxx (no range
// policy in the assembly) defer to. Binding policy, in order:
//
//   1. TRANSPARENCY FLOOR: if rho < rho_min(table) — the table's lowest
//      tabulated density — the cell is transparent for the table source:
//      kappa_a = kappa_s = eta = 0 for all species and energies, and NO
//      table evaluation occurs (the floor precedes even the nu_x branch,
//      which otherwise evaluates the Iso table twice).
//   2. CLAMPING ELSEWHERE: for rho >= rho_min the query coordinates
//      (E, rho, T, Ye) are clamped component-wise into the tabulated axis
//      ranges before evaluation, in the same log/linear space the axis is
//      stored in. (eta_e is not clamped here because the baseline assembly
//      never calls NES/Pair — spec's "as applicable".) No evaluated
//      coordinate ever lies outside its table axis range; combined with the
//      strictly-positive table axes this also closes the [MCNX-OPA-03]
//      no-negative-into-log10 guarantee, so no WeakLibInterp NaN can
//      propagate into transport.
//   3. DIAGNOSTICS: per-axis clamp counters (a diagnostic, not a physics
//      output) make silent saturation observable.
//
// Intersection clamping: the baseline assembly evaluates up to three
// co-resident tables per query (EmAb, Iso, EOS chemical potentials), each
// with its own axis grids. Each coordinate is clamped into the INTERSECTION
// of every range that consumes it ([max of lows, min of highs]), so one
// clamped triple is valid for every downstream call — mu_nu's three
// eos_evaluate calls, the EmAb call, and the Iso calls all see identical,
// in-range coordinates. Clamping the raw value against pow(10, log-bound)
// is equivalent to log-space clamping (pow10 is monotone and the wrappers
// own the log10s, [MCNX-OPA-03]); clamping T in MeV against MeV-mapped
// bounds is likewise equivalent (the MeV->K conversion is monotone and
// applied once inside the wrappers, [MCNX-OPA-02]).
//
// rho_min decision (recorded per the T11 brief): the spec's "the table's
// lowest tabulated density" is singular, but the view bundle carries three
// density axes (EmAb/Iso log10-stored, EOS raw). rho_min is the MAX of the
// three lower bounds — conservative: below the floor, NO table can be
// touched. In the production table set the three coincide (the .h5ls
// snapshots, opacity-eos-evaluation.md:78-89). rho_min is table-derived at
// runtime: no parameter, no new physical-constant literal.
//
// Also here: the [MCNX-OPA-07] EOS-inversion error-protocol consumption
// helper (inversion_usable). The wrappers of mcnux_opacity.hxx already
// surface the full EosInversionResultMeV{T_MeV, error}; the protocol is
// that error != 0 implies the returned T_MeV == 0 is NOT a temperature and
// must not be used, and that in test configurations an inversion error on
// the baseline path is a hard failure. The hard-failure plumbing belongs to
// the future consumer (the baseline transport path receives T in MeV
// externally and performs no inversion); this header pins only the
// pass/fail predicate so any future consumer inherits the protocol.
//
// Counters and device use: ClampCounters is caller-owned host state behind
// a nullable pointer (nullptr => skip counting), so RangedTableCoefficients
// stays trivially copyable with const noexcept methods. The plain
// increments are single-thread correct (the selftest battery is
// single-threaded); a device-side consumer must point `counters` at
// device-visible memory and replace the increments with atomics (or pass
// nullptr and count separately) — that adaptation is deliberately left to
// the first device consumer.
//
// Deliberately free of cctk.h and CarpetX includes, like the layers below
// it; the static_asserts run on every build via inclusion from stub.cxx.

#include "mcnux_opacity.hxx"
#include "mcnux_table_coeffs.hxx"

#include <cmath>
#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Clamp diagnostics  ([MCNX-OPA-06] "Diagnostics")
// ---------------------------------------------------------------------------
// One counter per clamped axis of the baseline assembly. eta_e has no
// counter because the baseline never evaluates NES/Pair (spec's "as
// applicable"); an inelastic extension adds it alongside its assembly.
struct ClampCounters {
  long long E = 0;   // fluid-frame energy clamps
  long long rho = 0; // density clamps
  long long T = 0;   // temperature clamps
  long long Ye = 0;  // electron-fraction clamps
};

// ---------------------------------------------------------------------------
// Bounds assembly from the view bundle (reuses the AxisBounds hook surface
// of mcnux_opacity.hxx — never re-derives xs[0]/xs[n-1] inline).
// ---------------------------------------------------------------------------

// Intersection of two ranges: [max of lows, min of highs].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
intersect_bounds(AxisBounds a, AxisBounds b) noexcept {
  return AxisBounds{std::fmax(a.lo, b.lo), std::fmin(a.hi, b.hi)};
}

// Un-log a log10-stored axis range into raw physical units (monotone, so
// the endpoints map to endpoints).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
pow10_bounds(AxisBounds b) noexcept {
  return AxisBounds{std::pow(wli::Real(10), b.lo),
                    std::pow(wli::Real(10), b.hi)};
}

// Temperature range of a log10-KELVIN axis, reported in MCNuX-internal MeV
// (the mcnux_opacity.hxx axis_bounds_temperature_mev counterpart for the
// pre-logged opacity-family axes).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
log_kelvin_bounds_temperature_mev(wli::Real const *LogTs_K, int n) noexcept {
  const AxisBounds k = pow10_bounds(axis_bounds(LogTs_K, n));
  return AxisBounds{temperature_mev_from_kelvin(k.lo),
                    temperature_mev_from_kelvin(k.hi)};
}

// E [MeV]: intersection of the EmAb and Iso energy ranges.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
energy_bounds_mev(const BaselineTableViews &v) noexcept {
  return intersect_bounds(pow10_bounds(axis_bounds(v.emab_LogEs, v.emab_nE)),
                          pow10_bounds(axis_bounds(v.iso_LogEs, v.iso_nE)));
}

// rho [g/cm^3]: intersection of the EmAb, Iso (log10-stored), and EOS (raw)
// density ranges.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
density_bounds_cgs(const BaselineTableViews &v) noexcept {
  return intersect_bounds(
      intersect_bounds(pow10_bounds(axis_bounds(v.emab_LogDs, v.emab_nD)),
                       pow10_bounds(axis_bounds(v.iso_LogDs, v.iso_nD))),
      axis_bounds(v.eos_Ds, v.eos_nD));
}

// T [MeV]: intersection of the EmAb, Iso (log10-Kelvin), and EOS (raw
// Kelvin) temperature ranges, all mapped into MCNuX-internal MeV.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
temperature_bounds_mev(const BaselineTableViews &v) noexcept {
  return intersect_bounds(
      intersect_bounds(
          log_kelvin_bounds_temperature_mev(v.emab_LogTs, v.emab_nT),
          log_kelvin_bounds_temperature_mev(v.iso_LogTs, v.iso_nT)),
      axis_bounds_temperature_mev(v.eos_Ts_K, v.eos_nT));
}

// Ye (dimensionless, always raw/linear — never logged, [MCNX-OPA-03]).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AxisBounds
ye_bounds(const BaselineTableViews &v) noexcept {
  return intersect_bounds(
      intersect_bounds(axis_bounds(v.emab_Ys, v.emab_nY),
                       axis_bounds(v.iso_Ys, v.iso_nY)),
      axis_bounds(v.eos_Ys, v.eos_nY));
}

// The transparency-floor density: the max of the three tabulated density
// lower bounds (see the header comment for the singular-vs-three decision).
// Identically density_bounds_cgs(v).lo — one derivation, no duplicate.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
rho_min_cgs(const BaselineTableViews &v) noexcept {
  return density_bounds_cgs(v).lo;
}

// ---------------------------------------------------------------------------
// The counting clamp
// ---------------------------------------------------------------------------
// min(max(x, lo), hi) via fmin/fmax, so a NaN input degrades to the lower
// bound instead of propagating (fmax(NaN, lo) == lo) — part of the no-NaN
// guarantee. The counter increments only when the value actually changed
// (a NaN input "changed": !(y == x) holds).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE wli::Real
clamp_count(wli::Real x, AxisBounds b, long long &ctr) noexcept {
  const wli::Real y = std::fmin(std::fmax(x, b.lo), b.hi);
  if (!(y == x))
    ++ctr;
  return y;
}

// ---------------------------------------------------------------------------
// The range-enforced table-source callable  ([MCNX-OPA-06])
// ---------------------------------------------------------------------------
// Satisfies the frozen TableEval contract of evaluate_coefficients()
// exactly like the BaselineTableCoefficients it wraps, so it drop-in
// replaces the bare assembly in the `table` slot of the [MCNX-VER-05]
// dispatch. Policy order: floor first (before any species dispatch, so the
// nu_x Iso calls are also guarded), then component-wise intersection
// clamping with per-axis counting, then ONE delegated call on the clamped
// triple (all downstream table calls see identical coordinates).
struct RangedTableCoefficients {
  BaselineTableCoefficients base;
  ClampCounters *counters; // caller-owned; nullptr => skip counting

  Coefficients operator()(Species s, double E_MeV,
                          const FluidState &st) const noexcept {
    // 1. Transparency floor: below the lowest tabulated density the cell is
    //    transparent — zero triple, no table touch, no counting (the floor
    //    is a policy outcome, not a clamp).
    if (st.rho_cgs < rho_min_cgs(base.v))
      return Coefficients{0.0, 0.0, 0.0};

    // 2./3. Component-wise intersection clamping with per-axis counters.
    ClampCounters scratch; // discarded when counters == nullptr
    ClampCounters &c = counters ? *counters : scratch;
    const double E = clamp_count(E_MeV, energy_bounds_mev(base.v), c.E);
    const double rho =
        clamp_count(st.rho_cgs, density_bounds_cgs(base.v), c.rho);
    const double T = clamp_count(st.T_MeV, temperature_bounds_mev(base.v), c.T);
    const double Ye = clamp_count(st.Ye, ye_bounds(base.v), c.Ye);

    // 4. One delegated evaluation on the clamped triple.
    return base(s, E, FluidState{rho, T, Ye});
  }
};

// ---------------------------------------------------------------------------
// EOS-inversion error protocol  ([MCNX-OPA-07])
// ---------------------------------------------------------------------------
// The one sanctioned pass/fail read of an inversion result: error != 0
// means the returned T_MeV == 0 is NOT a temperature and must not be used.
// In test configurations an inversion error on the baseline path is a hard
// failure — that hard-fail is the (future) consumer's obligation at its
// call site; no consumer exists on the baseline path (T arrives in MeV from
// the GRMHD partner, hydro-coupling-source-terms.md).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
inversion_usable(const EosInversionResultMeV &r) noexcept {
  return r.error == 0;
}

// ---------------------------------------------------------------------------
// Compile-time verification — structure only (the value-level policy rows
// are runtime rows of the unit-selftest battery, mcnux_selftest.cxx).
// ---------------------------------------------------------------------------

// The ranged functor satisfies the frozen TableEval contract, so it plugs
// into evaluate_coefficients() wherever the bare assembly does.
static_assert(
    std::is_invocable_r<Coefficients, const RangedTableCoefficients &, Species,
                        double, const FluidState &>::value,
    "[MCNX-OPA-06] RangedTableCoefficients must satisfy the TableEval "
    "contract (Species, double E_MeV, const FluidState&) -> Coefficients");

// Dumb, device-capturable aggregates (the BaselineTableViews precedent).
static_assert(std::is_trivially_copyable<ClampCounters>::value,
              "ClampCounters must be trivially copyable");
static_assert(std::is_trivially_copyable<RangedTableCoefficients>::value,
              "RangedTableCoefficients must be trivially copyable (captured "
              "by value into device kernels; counters is a caller-owned "
              "pointer, never owned state)");

} // namespace MCNuX

#endif // MCNUX_TABLE_RANGE_HXX
