// The MCNuX unit self-test battery — the runtime half of the `unit-selftest`
// benchmark of specs/verification-suite-design.md [MCNX-VER-01],
// [MCNX-VER-06].
//
// The battery re-runs, at runtime, the checks that the shared headers already
// assert at compile time, and writes them out as a table of rows so that the
// Cactus regression harness can diff them against committed golden data
// (MCNuX/test/unit-selftest/, produced by CarpetX's TSV writer at
// ABSTOL=RELTOL=1e-12). The spec's binding shape is "one row per named check
// (measured error, pass flag)"; here a row is one index of the 1d grid array
// MCNuX::mcnux_selftest, carrying `selftest_error` and `selftest_pass`.
//
// Why both compile-time and runtime? The static_asserts prove the predicates
// hold for the *compiler's* constant evaluator; these rows additionally pin
// the *measured* errors as golden numbers, so that a change which keeps a
// check passing but moves its error (a re-derived factor, a different
// evaluation order) is still caught at 1e-12.
//
// ---------------------------------------------------------------------------
// Row table (index -> name).  APPEND-ONLY.
// ---------------------------------------------------------------------------
// Golden data is keyed by row index, so rows are never reordered, renamed in
// meaning, or removed: later tasks append their checks at the end and
// regenerate the golden file (and bump SIZE in interface.ccl). The full
// battery the spec eventually requires (tetrad/null identities, Tmunu closed
// forms, nu_x mapping, relabeling identities, sign/range/inversion fixtures,
// sizeof(amrex::Real)) is built up this way.
//
//   0..2    rng.philox.kat{1,2,3}                 [MCNX-RNG-01], [MCNX-RNG-02]
//   3..6    rng.uniform.{u64_zero, u64_2p11_minus_1,
//                        u64_2p11, u64_max}       [MCNX-RNG-04]
//   7       rng.purity                            [MCNX-RNG-05], [MCNX-RNG-06]
//   8       rng.pairing                           [MCNX-RNG-03]
//   9..20   units.factor.<name>                   [MCNX-CNV-03], [MCNX-CNV-04]
//  21..32   units.roundtrip.<name>                [MCNX-CNV-02]
//  33       units.species_table                   [MCNX-CNV-06]
//  34..35   precision.sizeof_{double, cctk_real}  [MCNX-BLD-03]
//  36..37   precision.sizeof_{amrex_real,
//                             amrex_particle_real} [MCNX-BLD-03]
//  38..46   opacity.parity.{eos_evaluate,
//                eos_invert_dey_noguess, eos_invert_dey_guess,
//                emab, iso_imom1, nes_pair, nes_detailed_balance,
//                pair_crossing, brem_summed}      [MCNX-OPA-01..03]
//  47..49   opacity.units.{log10_kelvin_once,
//                nes_pair_kelvin_axis,
//                nes_db_dual_temperature}         [MCNX-OPA-02]
//  50..51   opacity.species.{emab,iso}_slots      [MCNX-OPA-08]
//  52       opacity.bounds.endpoints              (T11 range-policy hook)
//  53       geo.inverse_metric.identity           [MCNX-GEO-02]
//  54       geo.null_closure                      [MCNX-GEO-02]
//  55..57   tetrad.{orthonormality, null_norm,
//                   nu_recovery}                  [MCNX-PKT-04]
//  58..59   tetrad.flat_static.{momentum,
//                               identity}         [MCNX-PKT-04]
//  60..61   srcterm.emission.{lepton,energy}_sign [MCNX-HYD-02]
//  62..63   srcterm.absorption.{lepton,
//                               energy}_sign      [MCNX-HYD-02]
//  64..65   srcterm.scattering.{lepton_zero,
//                               momentum_only}    [MCNX-HYD-02]
//  66..67   trp.relabel.{eta_over_kappa_a,
//                        total_rate}              [MCNX-TRP-02]
//  68..69   trp.relabel.monotonicity_{ka, ks}     [MCNX-TRP-03]
//  70       trp.relabel.alpha_one_identity        [MCNX-TRP-02] (alpha = 1)
//  71       units.hc                              hc pin (opacity-eos-
//                                                 evaluation.md:356-361)
//  72..75   coef.analytic.{kappa_passthrough,
//                          kirchhoff_identity,
//                          eta_zero, nux_emission} [MCNX-VER-05]
//  76..78   coef.dispatch.{analytic, table,
//                          parameter}             [MCNX-VER-05] (source
//                                                 switch invisible below the
//                                                 coefficient interface)
//  79..81   geo.gather.{linear_values,
//                      linear_derivs,
//                      constant_derivs_zero}      [MCNX-GEO-03]
//  82       geo.rhs.null_closure_identity         [MCNX-GEO-02]
//  83..84   geo.push.flat_{straight_line,
//                         unit_speed}             [MCNX-GEO-04] (flat limit)
//  85       geo.push.order_ge2                    [MCNX-GEO-04] (convergence)
//
// Tiers, per specs/README.md and the tolerance discussion in mcnux_units.hxx:
//   * exact   — pass requires a measured error of identically 0 (KATs, the
//               [0, 1) closed forms, purity/pairing, the species table, the
//               precision rows).
//   * machine — the round-trip rows; MCNuX::detail::round_trips judges at
//               detail::rtol_machine = 1e-14 internally, so the row is
//               boolean and its measured error is 0 or 1.
//   * pinned  — the factor-recomputation rows, judged at
//               detail::rtol_pinned = 1e-9. Their measured errors are ~5e-10
//               by construction: the spec's table values carry 10 significant
//               digits, so agreement cannot be asserted at the machine tier
//               (see the rtol_pinned comment in mcnux_units.hxx).

#include "mcnux_coefficients.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_opacity.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_trp.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace MCNuX {

namespace {

// One row of the battery.
struct CheckRow {
  const char *name;
  double error; // measured error, in the row's own tier
  bool pass;
};

// Relative error |got/want - 1|, with an absolute fallback where the anchor is
// zero. This is the reported quantity for the `pinned` tier; the boolean tiers
// report 0 (pass) or 1 (fail), and the exact double comparisons report the
// absolute difference (which is 0 exactly when the check passes).
constexpr double rel_error(double got, double want) noexcept {
  const double d = detail::cabs(got - want);
  return want == 0.0 ? d : d / detail::cabs(want);
}

// Ordered, fixed-capacity accumulator; `capacity` only bounds the builder, the
// authoritative length is the number of rows actually appended, which is
// cross-checked against the interface.ccl array size below.
class Battery {
public:
  static constexpr int capacity = 256;

  void add_boolean(const char *name, bool ok) { add(name, ok ? 0.0 : 1.0, ok); }

  // Exact double comparison: no tolerance at all, and the measured error is
  // the absolute difference so a near miss is still visible in the golden data.
  void add_exact(const char *name, double got, double want) {
    add(name, detail::cabs(got - want), got == want);
  }

  // Factor recomputation against a pinned 10-significant-digit table anchor.
  void add_pinned(const char *name, double got, double want) {
    add(name, rel_error(got, want),
        detail::approx_eq(got, want, detail::rtol_pinned));
  }

  int size() const { return n_; }
  const CheckRow &operator[](int i) const { return rows_[i]; }

private:
  void add(const char *name, double error, bool pass) {
    if (n_ >= capacity)
      CCTK_VERROR("MCNuX self-test battery exceeded its builder capacity of "
                  "%d rows while adding \"%s\"",
                  capacity, name);
    rows_[n_].name = name;
    rows_[n_].error = error;
    rows_[n_].pass = pass;
    ++n_;
  }

  CheckRow rows_[capacity] = {};
  int n_ = 0;
};

// ---------------------------------------------------------------------------
// Rows 38..52 — the WeakLibInterp call-boundary wrapper layer of
// specs/opacity-eos-evaluation.md ([MCNX-OPA-01/02/03/08]; T9).
//
// Kernel parity is checked on small synthetic in-memory tables built inline
// below (monotone axes, smooth log-stored values — the pattern of
// WeakLibInterp/test/test_parallelfor_wrappers.cpp; NO HDF5, no file I/O:
// the battery stays table-file-free). For each family the MCNuX wrapper
// (mcnux_opacity.hxx) and the raw WLI `_Point` kernel are called on identical
// inputs and compared bitwise (exact tier, [MCNX-OPA-01]). The unit-boundary
// rows pin the MeV->Kelvin/log10 conversion applied exactly once per
// temperature argument — including the NES/Pair Kelvin axis and the
// detailed-balance dual-temperature split ([MCNX-OPA-02],
// specs/opacity-eos-evaluation.md:161-165,182-184). The species rows are the
// hardcoded structural mirror of [MCNX-OPA-08] (exactly two electron-type
// EmAb/Iso slots, nu_e = 0 / nu_e_bar = 1, no nu_x dataset; grounded on
// WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls:22-27 and
// specs/conventions-and-units.md:131,148 — live-fixture drift detection
// remains with specs/tools/validate_specs.sh:280-287).
// ---------------------------------------------------------------------------

// A bounded, finite, smooth log-stored synthetic table value keyed on the
// flat column-major offset (the test_parallelfor_wrappers.cpp generator).
double synth_tval(std::size_t k) {
  return 0.35 + 0.4 * std::sin(0.6 * double(k) + 0.2);
}

void append_opacity_rows(Battery &b) {
  using wli::Real;

  // --- Shared EOS 3D synthetic table: raw axes (rho [g/cm^3], T [KELVIN —
  // table native], Ye), flat column-major, rho fastest. ---
  constexpr int nD = 4, nT = 5, nY = 3;
  const Real Ds[nD] = {1.0e3, 5.0e5, 2.0e8, 6.0e11};
  const Real Ts_K[nT] = {1.0e9, 3.0e9, 1.0e10, 5.0e10, 1.0e11};
  const Real Ys[nY] = {0.05, 0.30, 0.55};
  const Real OS_eos = 3.5;
  Real eos_table[nD * nT * nY];
  for (std::size_t k = 0; k < std::size_t(nD) * nT * nY; ++k)
    eos_table[k] = synth_tval(k);

  // In-range query state, MCNuX units (T in MeV; 1.7 MeV ~ 1.97e10 K).
  const Real rho_q = 7.3e6, T_q_MeV = 1.7, Ye_q = 0.22;

  // Row 38 — EOS evaluate parity: wrapper (MeV in, one MeV->K conversion)
  // vs the raw kernel fed the independently converted Kelvin temperature.
  {
    const Real got = eos_evaluate(rho_q, T_q_MeV, Ye_q, Ds, nD, Ts_K, nT, Ys,
                                  nY, OS_eos, eos_table);
    const Real want = wli::EosInterpolateSingleVariable3DPoint(
        rho_q, T_q_MeV * mev_to_kelvin, Ye_q, Ds, nD, Ts_K, nT, Ys, nY, OS_eos,
        eos_table);
    b.add_exact("opacity.parity.eos_evaluate", got, want);
  }

  // --- EOS inversion synthetic sub-table: stored value strictly increasing
  // in the T index so a unique in-range root exists. ---
  Real inv_table[nD * nT * nY];
  for (int iY = 0; iY < nY; ++iY)
    for (int iT = 0; iT < nT; ++iT)
      for (int iD = 0; iD < nD; ++iD)
        inv_table[iD + nD * (iT + nT * iY)] =
            0.08 * iT + 0.01 * iD + 0.005 * iY;
  const Real OS_inv = 2.0;
  const Real X_q = -0.3; // between the recovered face values at Ts_K[0]/[nT-1]
  const wli::EosInversionBounds bounds{1.0e3, 6.0e11, -1.0, 1.0,
                                       0.05,  0.55,   true};

  // Row 39 — inversion parity, DEY no-guess: identical {T, Error} up to the
  // wrapper's single K->MeV result mapping, and a successful (Error == 0)
  // recovery so the row exercises the real bisection path.
  {
    const EosInversionResultMeV got = eos_invert_dey_noguess(
        rho_q, X_q, Ye_q, Ds, nD, Ts_K, nT, Ys, nY, OS_inv, inv_table, bounds);
    const wli::EosInversionResult raw = wli::ComputeTemperatureWith_DEY_NoGuess(
        rho_q, X_q, Ye_q, Ds, nD, Ts_K, nT, Ys, nY, OS_inv, inv_table, bounds);
    b.add_boolean("opacity.parity.eos_invert_dey_noguess",
                  got.error == 0 && raw.Error == 0 &&
                      got.T_MeV == raw.T / mev_to_kelvin);
  }

  // Row 40 — inversion parity, DEY guess-driven: T_Guess is a TEMPERATURE
  // ARGUMENT, so the wrapper takes it in MeV and converts once ([MCNX-OPA-02]);
  // the raw kernel receives the independently converted Kelvin guess.
  {
    const Real T_guess_MeV = 2.0;
    const EosInversionResultMeV got =
        eos_invert_dey_guess(rho_q, X_q, Ye_q, Ds, nD, Ts_K, nT, Ys, nY,
                             OS_inv, inv_table, T_guess_MeV, bounds);
    const wli::EosInversionResult raw = wli::ComputeTemperatureWith_DEY_Guess(
        rho_q, X_q, Ye_q, Ds, nD, Ts_K, nT, Ys, nY, OS_inv, inv_table,
        T_guess_MeV * mev_to_kelvin, bounds);
    b.add_boolean("opacity.parity.eos_invert_dey_guess",
                  got.error == 0 && raw.Error == 0 &&
                      got.T_MeV == raw.T / mev_to_kelvin);
  }

  // --- Shared EmAb/Iso synthetic grids: already-log10'd E [MeV], rho
  // [g/cm^3], T [KELVIN]; Ye raw ([MCNX-OPA-03]: MCNuX owns these log10s). ---
  constexpr int oE = 4, oD = 5, oT = 4, oY = 3, oMom = 2;
  const Real LogEs[oE] = {-1.0, 0.5, 1.5, 2.0};
  const Real LogDs[oD] = {6.0, 8.0, 10.0, 12.0, 14.0};
  const Real LogTs[oT] = {9.5, 10.0, 10.7, 11.3};
  const Real Ys2[oY] = {0.05, 0.25, 0.50};
  const Real E_q = 10.5, rho_q2 = 3.3e11, T_q2_MeV = 2.3, Ye_q2 = 0.31;

  // Row 41 — EmAb parity: wrapper logs (E, rho, T[K]) itself; the raw kernel
  // is fed independently computed log10 coordinates.
  {
    Real emab_table[oE * oD * oT * oY];
    for (std::size_t k = 0; k < std::size_t(oE) * oD * oT * oY; ++k)
      emab_table[k] = synth_tval(k);
    const Real OS = 1.25;
    const Real got = emab_evaluate(E_q, rho_q2, T_q2_MeV, Ye_q2, LogEs, oE,
                                   LogDs, oD, LogTs, oT, Ys2, oY, OS,
                                   emab_table);
    const Real want = wli::EmAbInterpolateSingleVariable4DPoint(
        std::log10(E_q), std::log10(rho_q2),
        std::log10(T_q2_MeV * mev_to_kelvin), Ye_q2, LogEs, oE, LogDs, oD,
        LogTs, oT, Ys2, oY, OS, emab_table);
    b.add_exact("opacity.parity.emab", got, want);
  }

  // Row 42 — Iso parity at a NONZERO moment index (iMom = 1), with the OS
  // selected through the species-major [nOpacities, nMoments] offset lookup
  // (iso_offset -> wli::IsoOffset; never hand-rolled).
  {
    Real iso_table[oE * oMom * oD * oT * oY];
    for (std::size_t k = 0; k < std::size_t(oE) * oMom * oD * oT * oY; ++k)
      iso_table[k] = synth_tval(k + 7);
    const Real offsets[num_tabulated_species * oMom] = {0.5, 1.0, 1.5, 2.0};
    const int iSpecies = iso_species_slot(Species::NuEBar), iMom = 1;
    const Real os_got =
        iso_offset(offsets, num_tabulated_species, oMom, iSpecies, iMom);
    const Real os_want =
        wli::IsoOffset(offsets, num_tabulated_species, oMom, iSpecies, iMom);
    const Real got =
        iso_evaluate(E_q, rho_q2, T_q2_MeV, Ye_q2, LogEs, oE, LogDs, oD, LogTs,
                     oT, Ys2, oY, iMom, oMom, os_got, iso_table);
    const Real want = wli::IsoInterpolateSingleVariable5DPoint(
        std::log10(E_q), std::log10(rho_q2),
        std::log10(T_q2_MeV * mev_to_kelvin), Ye_q2, LogEs, oE, LogDs, oD,
        LogTs, oT, Ys2, oY, iMom, oMom, os_want, iso_table);
    b.add_boolean("opacity.parity.iso_imom1",
                  os_got == os_want && got == want);
  }

  // --- Shared NES/Pair synthetic table: (nEp, nE, nMom, nT, nEta), axes
  // log10 KELVIN and log10 eta_e. T = 1.9 MeV ~ 2.2e10 K, in range. ---
  constexpr int pEp = 3, pE = 3, pMom = 4, pT = 4, pEta = 3;
  const Real LogTs_np[pT] = {10.0, 10.5, 11.0, 11.5};
  const Real LogEtas[pEta] = {-1.0, 0.0, 1.0};
  const Real OS_np = 0.75;
  Real np_table[pEp * pE * pMom * pT * pEta];
  for (std::size_t k = 0; k < std::size_t(pEp) * pE * pMom * pT * pEta; ++k)
    np_table[k] = synth_tval(k + 3);
  const Real T_np_MeV = 1.9, eta_q = 2.4;
  const Real E_np_MeV[pEp] = {4.0, 10.0, 24.0}; // physical grid, increasing

  // Row 43 — NES/Pair aligned-point parity (Kelvin axis).
  {
    const Real got =
        nes_pair_evaluate(T_np_MeV, eta_q, LogTs_np, pT, LogEtas, pEta, 1, 2,
                          pEp, pE, 2, pMom, OS_np, np_table);
    const Real want = wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
        std::log10(T_np_MeV * mev_to_kelvin), std::log10(eta_q), LogTs_np, pT,
        LogEtas, pEta, 1, 2, pEp, pE, 2, pMom, OS_np, np_table);
    b.add_exact("opacity.parity.nes_pair", got, want);
  }

  // Row 44 — NES detailed-balance fill parity: the raw kernel receives the
  // log10-KELVIN axis coordinate AND the bare MeV Boltzmann temperature.
  {
    const int iEp = 2, iE = 0, kernel = 1;
    const Real got = nes_detailed_balance_fill(T_np_MeV, eta_q, LogTs_np, pT,
                                               LogEtas, pEta, iEp, iE, pEp, pE,
                                               kernel, pMom, OS_np, np_table,
                                               E_np_MeV);
    const Real want = wli::NESDetailedBalanceFillPoint(
        std::log10(T_np_MeV * mev_to_kelvin), std::log10(eta_q), LogTs_np, pT,
        LogEtas, pEta, iEp, iE, pEp, pE, kernel, pMom, OS_np, np_table,
        E_np_MeV, T_np_MeV);
    b.add_exact("opacity.parity.nes_detailed_balance", got, want);
  }

  // Row 45 — Pair crossing-symmetry fill parity (upper triangle iEp > iE,
  // caller-swapped kernel component; exact relabeling, no Boltzmann factor).
  {
    const int iEp = 2, iE = 1, kernelSwapped = 3;
    const Real got =
        pair_crossing_fill(T_np_MeV, eta_q, LogTs_np, pT, LogEtas, pEta, iEp,
                           iE, pEp, pE, kernelSwapped, pMom, OS_np, np_table);
    const Real want = wli::PairCrossingSymmetryFillPoint(
        std::log10(T_np_MeV * mev_to_kelvin), std::log10(eta_q), LogTs_np, pT,
        LogEtas, pEta, iEp, iE, pEp, pE, kernelSwapped, pMom, OS_np, np_table);
    b.add_exact("opacity.parity.pair_crossing", got, want);
  }

  // Row 46 — Brem summed parity: (rho, T) plane with rho BEFORE T (nD != nT
  // here, so a transposed wrapper cannot pass), Alpha = [1, 1, 28/3] applied
  // to recovered values, three caller-logged effective densities.
  {
    constexpr int bEp = 3, bE = 3, bMom = 1, bD = 4, bT = 5;
    const Real LogDs_b[bD] = {8.0, 10.0, 12.0, 14.0};
    const Real LogTs_b[bT] = {9.5, 10.2, 10.8, 11.2, 11.6};
    const Real OS_b = 0.4;
    Real brem_table[bEp * bE * bMom * bD * bT];
    for (std::size_t k = 0; k < std::size_t(bEp) * bE * bMom * bD * bT; ++k)
      brem_table[k] = synth_tval(k + 11);
    const Real rho_eff[brem_num_effective_densities] = {2.1e11, 3.4e11,
                                                        2.7e11};
    const Real T_b_MeV = 2.1;
    const int iEp = 0, iE = 2, moment = 0;
    const Real got =
        brem_evaluate_summed(rho_eff, T_b_MeV, LogDs_b, bD, LogTs_b, bT, iEp,
                             iE, bEp, bE, moment, bMom, OS_b, brem_table);
    Real LogD_eff[brem_num_effective_densities];
    for (int l = 0; l < brem_num_effective_densities; ++l)
      LogD_eff[l] = std::log10(rho_eff[l]);
    const Real Alpha[brem_num_effective_densities] = {1.0, 1.0, 28.0 / 3.0};
    const Real want = wli::BremInterpolateSingleVariable2D2DAlignedSummedPoint(
        LogD_eff, Alpha, brem_num_effective_densities,
        std::log10(T_b_MeV * mev_to_kelvin), LogDs_b, bD, LogTs_b, bT, iEp, iE,
        bEp, bE, moment, bMom, OS_b, brem_table);
    b.add_exact("opacity.parity.brem_summed", got, want);
  }

  // Row 47 — unit boundary, exact: the wrapper layer's internally-used
  // Kelvin/log10 coordinate equals log10(T[MeV] * mev_to_kelvin) computed
  // independently — the conversion is applied exactly once, in this order.
  b.add_exact("opacity.units.log10_kelvin_once",
              log10_temperature_kelvin_from_mev(T_q_MeV),
              std::log10(T_q_MeV * mev_to_kelvin));

  // Row 48 — the NES/Pair temperature axis is KELVIN ([MCNX-OPA-02]'s
  // "including NES/Pair", the documented easy-to-make error): the wrapper
  // result equals the raw kernel at the Kelvin coordinate and DIFFERS from
  // the raw kernel fed the mistaken log10(T[MeV]) coordinate.
  {
    const Real got =
        nes_pair_evaluate(T_np_MeV, eta_q, LogTs_np, pT, LogEtas, pEta, 1, 2,
                          pEp, pE, 2, pMom, OS_np, np_table);
    const Real want_kelvin =
        wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
            std::log10(T_np_MeV * mev_to_kelvin), std::log10(eta_q), LogTs_np,
            pT, LogEtas, pEta, 1, 2, pEp, pE, 2, pMom, OS_np, np_table);
    const Real wrong_mev = wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
        std::log10(T_np_MeV), std::log10(eta_q), LogTs_np, pT, LogEtas, pEta,
        1, 2, pEp, pE, 2, pMom, OS_np, np_table);
    b.add_boolean("opacity.units.nes_pair_kelvin_axis",
                  got == want_kelvin && got != wrong_mev);
  }

  // Row 49 — the detailed-balance dual-temperature split: LogT is the
  // KELVIN table coordinate while the unlogged Boltzmann temperature stays
  // in MeV (the energy unit of E). The wrapper must match the split form and
  // differ from the misuse that feeds the Kelvin value into the Boltzmann
  // factor.
  {
    const int iEp = 2, iE = 0, kernel = 1;
    const Real got = nes_detailed_balance_fill(T_np_MeV, eta_q, LogTs_np, pT,
                                               LogEtas, pEta, iEp, iE, pEp, pE,
                                               kernel, pMom, OS_np, np_table,
                                               E_np_MeV);
    // The symmetric lower-triangle read Phi(iE, iEp): energy indices swapped.
    const Real lower = wli::NESPairInterpolateSingleVariable2D2DAlignedPoint(
        std::log10(T_np_MeV * mev_to_kelvin), std::log10(eta_q), LogTs_np, pT,
        LogEtas, pEta, iE, iEp, pEp, pE, kernel, pMom, OS_np, np_table);
    const Real want =
        lower * std::exp((E_np_MeV[iE] - E_np_MeV[iEp]) / T_np_MeV);
    const Real wrong_kelvin_boltzmann =
        lower * std::exp((E_np_MeV[iE] - E_np_MeV[iEp]) /
                         (T_np_MeV * mev_to_kelvin));
    b.add_boolean("opacity.units.nes_db_dual_temperature",
                  got == want && got != wrong_kelvin_boltzmann);
  }

  // Rows 50..51 — [MCNX-OPA-08] species-axis grounding, structural constants
  // (no file I/O; the species_table_holds pattern): exactly two electron-type
  // slots ordered nu_e = 0, nu_e_bar = 1, and nu_x has NO EmAb/Iso dataset.
  b.add_boolean("opacity.species.emab_slots", detail::species_axis_holds_emab());
  b.add_boolean("opacity.species.iso_slots", detail::species_axis_holds_iso());

  // Row 52 — the T11 range-policy hook surface: AxisBounds reads the axis
  // endpoints, and the temperature variant reports them in MCNuX-internal MeV.
  {
    const AxisBounds bd = axis_bounds(Ds, nD);
    const AxisBounds bt = axis_bounds_temperature_mev(Ts_K, nT);
    b.add_boolean("opacity.bounds.endpoints",
                  bd.lo == Ds[0] && bd.hi == Ds[nD - 1] &&
                      bt.lo == Ts_K[0] / mev_to_kelvin &&
                      bt.hi == Ts_K[nT - 1] / mev_to_kelvin);
  }
}

// ---------------------------------------------------------------------------
// Rows 53..59 — tetrad construction + null-closure helpers of
// specs/geodesic-propagation.md [MCNX-GEO-02] and
// specs/packet-representation-and-sampling.md [MCNX-PKT-04] (T13).
//
// Inputs are synthetic and self-contained: a flat-static pair (Minkowski,
// u^mu = (1,0,0,0)) probing the closed-form limit of
// packet-representation-and-sampling.md:187-190, and a curved pair (SPD
// gamma_ij with off-diagonal terms, alpha != 1, beta^i != 0, and a boosted
// u^mu whose u^t is solved at runtime from g_munu u^mu u^nu = -1) probing the
// generic identities. All rows are machine-tier booleans (error 0 or 1), so
// the golden data is bitwise-stable across libm ulp differences. The
// direction draw takes fixed literal uniforms — the helper deliberately
// performs no u(S, q, e, k) call ([MCNX-PKT-05] draw-map bookkeeping belongs
// to the consuming code, not to these helpers).
// ---------------------------------------------------------------------------

void append_tetrad_rows(Battery &b) {
  // Absolute machine-tier closeness for O(1) quantities whose target may be
  // exactly 0 (orthonormality off-diagonals, the null norm), where a relative
  // comparison is undefined.
  const auto near_abs = [](double got, double want) {
    return detail::cabs(got - want) <= detail::rtol_machine;
  };

  // --- Synthetic curved input. ---
  const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
  const double alp = 0.83;
  const double betau[3] = {0.06, -0.11, 0.04};
  const SpacetimeMetric m = spacetime_metric_from_adm(alp, betau, gamma);

  // Boosted u^mu: spatial components chosen, u^t fixed by the normalization
  // g_munu u^mu u^nu = -1 (future-pointing root of the quadratic
  // g_tt (u^t)^2 + 2 g_ti u^i u^t + (gamma_ij u^i u^j + 1) = 0).
  double uvec[4] = {0.0, 0.15, -0.05, 0.10};
  {
    const double a = m.g[0][0];
    double bb = 0.0, cc = 1.0;
    for (int i = 1; i < 4; ++i) {
      bb += m.g[0][i] * uvec[i];
      for (int j = 1; j < 4; ++j)
        cc += m.g[i][j] * uvec[i] * uvec[j];
    }
    uvec[0] = (bb + std::sqrt(bb * bb - a * cc)) / (-a);
  }

  // Row 53 — [MCNX-GEO-02] analytic inverse: gamma^{ij} gamma_jk = delta^i_k
  // (machine) and det(gamma_ij) > 0 on the SPD test metric.
  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
  {
    const double lo[3][3] = {{gamma.xx, gamma.xy, gamma.xz},
                             {gamma.xy, gamma.yy, gamma.yz},
                             {gamma.xz, gamma.yz, gamma.zz}};
    const double up[3][3] = {
        {gu.xx, gu.xy, gu.xz}, {gu.xy, gu.yy, gu.yz}, {gu.xz, gu.yz, gu.zz}};
    bool ok = gu.det > 0.0;
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j)
          s += up[i][j] * lo[j][k];
        ok = ok && near_abs(s, i == k ? 1.0 : 0.0);
      }
    b.add_boolean("geo.inverse_metric.identity", ok);
  }

  // One drawn packet on the curved background: fixed literal uniforms into
  // the isotropic draw, then the pinned x->y->z Gram-Schmidt tetrad and the
  // [MCNX-PKT-04] transform.
  const double nu = 3.25;
  const UnitVector3 nhat = fluid_frame_direction(0.35, 0.62);
  const Tetrad tet = build_tetrad(alp, betau, gamma, uvec);
  const FourVectorUp pup = fluid_momentum_to_coordinate(tet, nu, nhat);
  const CoordinateMomentum pm = transform_to_coordinate_frame(tet, m, nu, nhat);

  // Row 54 — [MCNX-GEO-02] null closure on the transformed momentum:
  // p^t = sqrt(gamma^{ij} p_i p_j)/alpha recovers the tetrad's p^t (the
  // physical nullness content) and alpha^2 (p^t)^2 = gamma^{ij} p_i p_j
  // (geodesic-propagation.md:96-97), both machine tier.
  {
    const double pt = p_t_closure(pm.px, pm.py, pm.pz, gu, alp);
    const double n2 = momentum_norm2(gu, pm.px, pm.py, pm.pz);
    b.add_boolean("geo.null_closure",
                  detail::approx_eq(pt, pm.pt, detail::rtol_machine) &&
                      detail::approx_eq(alp * alp * pt * pt, n2,
                                        detail::rtol_machine));
  }

  // Row 55 — tetrad orthonormality: g_munu e^mu_(a) e^nu_(b) = eta_ab for
  // all 10 independent pairs, machine tier.
  {
    bool ok = true;
    for (int a = 0; a < 4; ++a)
      for (int c = a; c < 4; ++c) {
        const double eta = a != c ? 0.0 : (a == 0 ? -1.0 : 1.0);
        ok = ok && near_abs(metric_dot(m, tet.e[a], tet.e[c]), eta);
      }
    b.add_boolean("tetrad.orthonormality", ok);
  }

  // Row 56 — transformed momentum is null: g_munu p^mu p^nu = 0 normalized
  // by (alpha p^t)^2 (packet-representation-and-sampling.md:138-140).
  {
    const double norm = metric_dot(m, pup.v, pup.v);
    const double scale = (alp * pm.pt) * (alp * pm.pt);
    b.add_boolean("tetrad.null_norm", near_abs(norm / scale, 0.0));
  }

  // Row 57 — nu recovery: -p_mu u^mu equals the drawn nu, machine tier.
  {
    double pu = 0.0;
    for (int mu = 0; mu < 4; ++mu)
      for (int nu4 = 0; nu4 < 4; ++nu4)
        pu += m.g[mu][nu4] * pup.v[nu4] * uvec[mu];
    b.add_boolean("tetrad.nu_recovery",
                  detail::approx_eq(-pu, nu, detail::rtol_machine));
  }

  // --- Flat-static limit: Minkowski, u^mu = (1, 0, 0, 0). ---
  const SpatialMetric flat{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
  const SpacetimeMetric mf = spacetime_metric_from_adm(1.0, beta0, flat);
  const Tetrad tf = build_tetrad(1.0, beta0, flat, ustat);
  const CoordinateMomentum pf = transform_to_coordinate_frame(tf, mf, nu, nhat);

  // Row 58 — flat-static momentum: p^t = nu and p_i = nu n_hat_i, machine
  // tier (packet-representation-and-sampling.md:187-190).
  b.add_boolean("tetrad.flat_static.momentum",
                detail::approx_eq(pf.pt, nu, detail::rtol_machine) &&
                    detail::approx_eq(pf.px, nu * nhat.x,
                                      detail::rtol_machine) &&
                    detail::approx_eq(pf.py, nu * nhat.y,
                                      detail::rtol_machine) &&
                    detail::approx_eq(pf.pz, nu * nhat.z,
                                      detail::rtol_machine));

  // Row 59 — flat-static tetrad legs equal the coordinate basis, exactly
  // (the x->y->z Gram-Schmidt touches only exact values on this input).
  {
    bool ok = true;
    for (int a = 0; a < 4; ++a)
      for (int mu = 0; mu < 4; ++mu)
        ok = ok && tf.e[a][mu] == (a == mu ? 1.0 : 0.0);
    b.add_boolean("tetrad.flat_static.identity", ok);
  }
}

// ---------------------------------------------------------------------------
// Rows 60..65 — the [MCNX-HYD-02] source-term ledger arithmetic of
// specs/hydro-coupling-source-terms.md (T17). Pure-function sign fixtures
// through the mcnux_srcterms.hxx helpers on the header's own compile-time
// fixture set (exact binary-fraction inputs, dV dt = 0.5), so every row is
// exact: the worked signs of the spec (nu_e emission -> S_l < 0, G^t < 0;
// nu_e absorption reverses both; elastic scattering changes momentum
// components only) plus the exact fixture values, judged with no tolerance.
// ---------------------------------------------------------------------------

void append_srcterm_rows(Battery &b) {
  constexpr SourceContribution emis = detail::srcfix_emission;
  constexpr SourceContribution absn = detail::srcfix_absorption;
  constexpr SourceContribution scat = detail::srcfix_scattering;

  // Rows 60..61 — nu_e emission: the fluid loses a lepton and energy
  // (S_l < 0, G^t < 0), at the exact fixture values -6 and -15.
  b.add_boolean("srcterm.emission.lepton_sign",
                emis.Sl < 0.0 && emis.Sl == -6.0);
  b.add_boolean("srcterm.emission.energy_sign",
                emis.Gt < 0.0 && emis.Gt == -15.0);

  // Rows 62..63 — nu_e absorption reverses both signs and is the exact
  // negation of emission, component-wise.
  b.add_boolean("srcterm.absorption.lepton_sign",
                absn.Sl > 0.0 && absn.Sl == -emis.Sl);
  b.add_boolean("srcterm.absorption.energy_sign",
                absn.Gt > 0.0 && absn.Gt == -emis.Gt && absn.Gx == -emis.Gx &&
                    absn.Gy == -emis.Gy && absn.Gz == -emis.Gz);

  // Rows 64..65 — elastic scattering changes momentum components only:
  // S_l and G^t are exactly zero, the momentum transfer is the exact,
  // nonzero fixture value.
  b.add_exact("srcterm.scattering.lepton_zero", scat.Sl, 0.0);
  b.add_boolean("srcterm.scattering.momentum_only",
                scat.Gt == 0.0 && scat.Gx == 5.0 && scat.Gy == -2.0 &&
                    scat.Gz == -4.0);
}

// ---------------------------------------------------------------------------
// Rows 66..70 — the [MCNX-TRP-02/03] trapped-regime relabeling pure map of
// specs/trapped-regime-treatment.md (T21a), through mcnux_trp.hxx. Synthetic
// pinned (eta, kappa_a, kappa_s) tuples (kappa_a > 0 so the ratio invariant
// is defined) swept over alpha in (0, 1] — no live opacity tables. The two
// invariant rows are machine-tier booleans (detail::approx_eq at
// detail::rtol_machine, the spec's ~1e-14 tier); the monotonicity rows are
// the exact strict-inequality convention tripwire of [MCNX-TRP-03] (the
// flipped 2020-letter alpha convention reverses both); the alpha = 1 row is
// the exact bitwise identity endpoint the [MCNX-TRP-05] limit benchmark
// later builds on. Alpha SELECTION (the xi control rule) is T21b, not here.
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Rows 71..78 — the analytic (gray) opacity mode and the source-agnostic
// coefficient dispatch of specs/verification-suite-design.md [MCNX-VER-05]
// (T12) through mcnux_coefficients.hxx, plus the hc pin those formulas rest
// on (specs/opacity-eos-evaluation.md:356-361).
//
// The analytic rows sweep a synthetic AnalyticOpacityParams (exact binary
// fractions, distinct per species, kappa_a0(NuX) > 0) over a literal (E, T)
// grid spanning E/T from ~1e-2 to ~1e2: kappa pass-through is bitwise
// (exact tier); the Kirchhoff identity eta/kappa_a = eta_scale * c * 4 pi
// E^3/(hc)^3/(exp(E/T) + 1) is written out independently here and judged at
// detail::rtol_machine (machine tier, verification-suite-design.md:320-321);
// eta vanishes exactly when eta_scale = 0 or kappa_a0 = 0; and nu_x emission
// is strictly positive (the nu_x coverage mechanism, :225-227). The dispatch
// rows pin that evaluate_coefficients() returns the analytic triple bitwise
// under Analytic and a synthetic table callable's triple bitwise under
// Table — and that the two differ, so the switch is observable — and that
// the parameter glue reads the parfile literals of unit-selftest.par back
// bitwise (the parfile sets opacity_source = "analytic" and the three arrays
// to the exact binary fractions restated below).
// ---------------------------------------------------------------------------

void append_coefficient_rows(Battery &b) {
  // Row 71 — the hc pin, recomputed from h, c, e against the 10-digit
  // anchor of specs/opacity-eos-evaluation.md:359-360.
  b.add_pinned("units.hc", hc_MeV_cm, pinned::hc_MeV_cm);

  const AnalyticOpacityParams ap = {
      {0.5, 0.25, 0.125}, {0.0625, 0.03125, 0.015625}, {1.0, 0.5, 1.0}};
  const Species species[NUM_SPECIES] = {Species::NuE, Species::NuEBar,
                                        Species::NuX};
  const double energies[] = {0.125, 1.0, 3.7, 12.5, 80.0};
  const double temperatures[] = {0.5, 2.0, 11.0};
  constexpr int nE = int(sizeof(energies) / sizeof(energies[0]));
  constexpr int nT = int(sizeof(temperatures) / sizeof(temperatures[0]));

  // Row 72 — kappa_a = kappa_a0(s), kappa_s = kappa_s0(s), bitwise over the
  // whole sweep (energy- and temperature-independent).
  bool pass_ok = true;
  for (int is = 0; is < NUM_SPECIES; ++is)
    for (int ie = 0; ie < nE; ++ie)
      for (int it = 0; it < nT; ++it) {
        const Coefficients c =
            analytic_coefficients(ap, species[is], energies[ie],
                                  temperatures[it]);
        pass_ok = pass_ok && c.kappa_a == ap.kappa_a0[is] &&
                  c.kappa_s == ap.kappa_s0[is];
      }
  b.add_boolean("coef.analytic.kappa_passthrough", pass_ok);

  // Row 73 — the Kirchhoff identity, independently restated from
  // verification-suite-design.md:216, machine tier wherever kappa_a0 > 0
  // (every species of this fixture).
  bool kirch_ok = true;
  for (int is = 0; is < NUM_SPECIES; ++is)
    for (int ie = 0; ie < nE; ++ie)
      for (int it = 0; it < nT; ++it) {
        const double E = energies[ie], T = temperatures[it];
        const Coefficients c =
            analytic_coefficients(ap, species[is], E, T);
        const double want = ap.eta_scale[is] * c_cgs * 4.0 * detail::pi *
                            E * E * E / (hc_MeV_cm * hc_MeV_cm * hc_MeV_cm) /
                            (std::exp(E / T) + 1.0);
        kirch_ok = kirch_ok && detail::approx_eq(c.eta / c.kappa_a, want,
                                                 detail::rtol_machine);
      }
  b.add_boolean("coef.analytic.kirchhoff_identity", kirch_ok);

  // Row 74 — eta = 0 exactly when eta_scale = 0 (absorption-only medium) and
  // when kappa_a0 = 0 (no absorber, hence no Kirchhoff emitter).
  AnalyticOpacityParams ap_noeta = ap;
  AnalyticOpacityParams ap_noabs = ap;
  for (int is = 0; is < NUM_SPECIES; ++is) {
    ap_noeta.eta_scale[is] = 0.0;
    ap_noabs.kappa_a0[is] = 0.0;
  }
  double eta_zero_err = 0.0;
  for (int is = 0; is < NUM_SPECIES; ++is)
    for (int ie = 0; ie < nE; ++ie)
      for (int it = 0; it < nT; ++it) {
        const double e1 = analytic_coefficients(ap_noeta, species[is],
                                                energies[ie], temperatures[it])
                              .eta;
        const double e2 = analytic_coefficients(ap_noabs, species[is],
                                                energies[ie], temperatures[it])
                              .eta;
        eta_zero_err += detail::cabs(e1) + detail::cabs(e2);
      }
  b.add_exact("coef.analytic.eta_zero", eta_zero_err, 0.0);

  // Row 75 — nu_x emission is reachable: eta(NuX) > 0 with kappa_a0(NuX) > 0
  // and eta_scale(NuX) = 1, at every (E, T) of the sweep (E > 0 throughout).
  bool nux_ok = ap.kappa_a0[species_index(Species::NuX)] > 0.0 &&
                ap.eta_scale[species_index(Species::NuX)] == 1.0;
  for (int ie = 0; ie < nE; ++ie)
    for (int it = 0; it < nT; ++it)
      nux_ok = nux_ok && analytic_coefficients(ap, Species::NuX, energies[ie],
                                               temperatures[it])
                                 .eta > 0.0;
  b.add_boolean("coef.analytic.nux_emission", nux_ok);

  // Synthetic table-source callable with the dispatch contract
  // (Species, double E_MeV, const FluidState&) -> Coefficients: a smooth
  // function of all three inputs that cannot coincide with the analytic
  // triple (its kappa_a carries the species index and rho).
  const auto table = [](Species s, double E, const FluidState &st) {
    const double k = 1.0 + species_index(s);
    return Coefficients{k * 1e-3 * st.rho_cgs, k * 2e-4 * E * st.Ye,
                        k * 7.0 * st.T_MeV * E};
  };
  const FluidState states[] = {
      {1.0e10, 2.0, 0.25}, {3.0e12, 11.0, 0.375}, {5.0e14, 0.5, 0.1}};
  constexpr int nS = int(sizeof(states) / sizeof(states[0]));

  // Rows 76..77 — dispatch exactness: Analytic returns the analytic triple
  // bitwise and differs from the table triple; Table returns the table
  // triple bitwise. Both through the identical Coefficients type.
  bool disp_an_ok = true, disp_tb_ok = true;
  for (int is = 0; is < NUM_SPECIES; ++is)
    for (int ie = 0; ie < nE; ++ie)
      for (int ist = 0; ist < nS; ++ist) {
        const Species s = species[is];
        const double E = energies[ie];
        const FluidState &st = states[ist];
        const Coefficients an = analytic_coefficients(ap, s, E, st.T_MeV);
        const Coefficients tb = table(s, E, st);
        const Coefficients got_an = evaluate_coefficients(
            CoefficientSource::Analytic, ap, table, s, E, st);
        const Coefficients got_tb =
            evaluate_coefficients(CoefficientSource::Table, ap, table, s, E, st);
        disp_an_ok = disp_an_ok && got_an.kappa_a == an.kappa_a &&
                     got_an.kappa_s == an.kappa_s && got_an.eta == an.eta &&
                     (got_an.kappa_a != tb.kappa_a ||
                      got_an.kappa_s != tb.kappa_s || got_an.eta != tb.eta);
        disp_tb_ok = disp_tb_ok && got_tb.kappa_a == tb.kappa_a &&
                     got_tb.kappa_s == tb.kappa_s && got_tb.eta == tb.eta;
      }
  b.add_boolean("coef.dispatch.analytic", disp_an_ok);
  b.add_boolean("coef.dispatch.table", disp_tb_ok);

  // Row 78 — the parameter glue: MCNuX/test/unit-selftest.par sets
  // opacity_source = "analytic" and kappa_a0 = {0.5, 0.25, 0.125},
  // kappa_s0 = {0.0625, 0.03125, 0.015625}, eta_scale = {1.0, 0.5, 1.0};
  // the same literals as the fixture `ap` above. Bitwise read-back.
  const AnalyticOpacityParams from_par = analytic_params_from_parameters();
  bool par_ok = selected_coefficient_source() == CoefficientSource::Analytic;
  for (int is = 0; is < NUM_SPECIES; ++is)
    par_ok = par_ok && from_par.kappa_a0[is] == ap.kappa_a0[is] &&
             from_par.kappa_s0[is] == ap.kappa_s0[is] &&
             from_par.eta_scale[is] == ap.eta_scale[is];
  b.add_boolean("coef.dispatch.parameter", par_ok);
}

// ---------------------------------------------------------------------------
// Rows 79..85 — the geodesic push of specs/geodesic-propagation.md (T15)
// through mcnux_geodesic.hxx: the vertex-centered trilinear gather
// [MCNX-GEO-03] on synthetic in-memory vertex data (host buffers wrapped as
// amrex::Array4 views over a small vertex box — the in-memory synthetic-table
// discipline of the opacity rows; no grid, no driver), the null-closure
// identity [MCNX-GEO-02] at every RK stage of a curved analytic metric, the
// flat-spacetime exactness of [MCNX-GEO-04] (straight line to rtol_machine,
// momentum bitwise constant, |dx/dt| = 1), and the convergence order of the
// integrator on a smooth analytic metric (RK4 is used; the row asserts the
// spec's order >= 2 bound, ratio >= 3.9 on step halving). Booleans only.
// ---------------------------------------------------------------------------

// Synthetic vertex box: 5 x 5 x 5 vertices (indices 0..4), anisotropic
// spacing and an off-origin prob_lo so that index/coordinate mix-ups cannot
// pass. Ten fields in ADMBaseX component order (metric 6, lapse 1, shift 3).
struct SyntheticVertexData {
  static constexpr int n = 5;
  static constexpr double prob_lo[3] = {-1.0, -0.5, 0.25};
  static constexpr double dx[3] = {0.25, 0.5, 0.125};
  std::vector<double> metric, lapse, shift; // component-major buffers

  SyntheticVertexData()
      : metric(std::size_t(n) * n * n * 6), lapse(std::size_t(n) * n * n),
        shift(std::size_t(n) * n * n * 3) {}

  static std::size_t at(int i, int j, int k, int c) {
    return std::size_t(i) + std::size_t(n) * (std::size_t(j) + std::size_t(n) * (std::size_t(k) + std::size_t(n) * c));
  }

  // Fill every field as the affine function a + b . (x, y, z) of the vertex
  // coordinates (a, b per field from `coef`, ten rows of four).
  void fill_affine(const double coef[10][4]) {
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          const double x = prob_lo[0] + i * dx[0];
          const double y = prob_lo[1] + j * dx[1];
          const double z = prob_lo[2] + k * dx[2];
          for (int c = 0; c < 10; ++c) {
            const double v = coef[c][0] + coef[c][1] * x + coef[c][2] * y +
                             coef[c][3] * z;
            if (c < 6)
              metric[at(i, j, k, c)] = v;
            else if (c == 6)
              lapse[at(i, j, k, 0)] = v;
            else
              shift[at(i, j, k, c - 7)] = v;
          }
        }
  }

  VertexMetricGather gather() const {
    const amrex::Dim3 lo{0, 0, 0}, hi{n, n, n}; // hi exclusive
    VertexMetricGather g{};
    g.metric = amrex::Array4<const amrex::Real>(metric.data(), lo, hi, 6);
    g.lapse = amrex::Array4<const amrex::Real>(lapse.data(), lo, hi, 1);
    g.shift = amrex::Array4<const amrex::Real>(shift.data(), lo, hi, 3);
    for (int d = 0; d < 3; ++d) {
      g.prob_lo[d] = prob_lo[d];
      g.dx[d] = dx[d];
    }
    return g;
  }
};

// A smooth curved analytic metric with exact derivatives (polynomial, SPD
// and diagonally dominant on |x^i| <~ 1.5): the fixture for the closure
// identity and convergence rows. Host-only test scaffolding.
struct CurvedAnalyticMetric {
  MetricSnapshot operator()(double x, double y, double z) const noexcept {
    MetricSnapshot m{};
    m.alpha = 1.0 + 0.2 * x + 0.1 * y * z;
    m.dalpha[0] = 0.2;
    m.dalpha[1] = 0.1 * z;
    m.dalpha[2] = 0.1 * y;
    m.beta[0] = 0.1 * y;
    m.beta[1] = 0.05 * z;
    m.beta[2] = 0.1 * x * y;
    m.dbeta[1][0] = 0.1;      // d_y beta^x
    m.dbeta[2][1] = 0.05;     // d_z beta^y
    m.dbeta[0][2] = 0.1 * y;  // d_x beta^z
    m.dbeta[1][2] = 0.1 * x;  // d_y beta^z
    m.g = SpatialMetric{1.0 + 0.1 * x * x, 0.1 * z,  0.05 * y,
                        1.0 + 0.2 * y * y, 0.1 * x,  1.0 + 0.1 * z * z};
    m.dg[0] = SpatialMetric{0.2 * x, 0.0, 0.0, 0.0, 0.1, 0.0};
    m.dg[1] = SpatialMetric{0.0, 0.0, 0.05, 0.4 * y, 0.0, 0.0};
    m.dg[2] = SpatialMetric{0.0, 0.1, 0.0, 0.0, 0.0, 0.2 * z};
    return m;
  }
};

void append_geodesic_rows(Battery &b) {
  // --- Rows 79..80: affine fields, values and derivatives at off-vertex
  // points (machine tier). Distinct, exact-binary coefficients per field. ---
  const double coef[10][4] = {
      {1.25, 0.125, -0.0625, 0.03125},  // gxx
      {0.0625, -0.03125, 0.125, 0.25},  // gxy
      {-0.125, 0.25, 0.0625, -0.125},   // gxz
      {1.5, -0.25, 0.125, 0.0625},      // gyy
      {0.03125, 0.0625, -0.125, 0.25},  // gyz
      {1.125, 0.0625, 0.25, -0.03125},  // gzz
      {0.875, 0.25, -0.125, 0.0625},    // alp
      {0.125, -0.0625, 0.25, 0.125},    // betax
      {-0.25, 0.125, 0.0625, -0.25},    // betay
      {0.0625, 0.03125, -0.25, 0.125},  // betaz
  };
  SyntheticVertexData data;
  data.fill_affine(coef);
  const VertexMetricGather gather = data.gather();

  // Query points strictly inside the 4 x 4 x 4 cell box, at generic
  // (non-binary) fractions of a cell, plus one exactly on a vertex.
  const double queries[][3] = {{-0.63, 0.37, 0.41},
                               {-0.07, 1.23, 0.66},
                               {-0.9, -0.45, 0.27},
                               {-0.5, 0.5, 0.5},
                               {-0.26, 1.49, 0.74}};
  bool values_ok = true, derivs_ok = true;
  for (const auto &q : queries) {
    const MetricSnapshot m = gather(q[0], q[1], q[2]);
    double want[10], got[10], dgot[10][3];
    for (int c = 0; c < 10; ++c)
      want[c] = coef[c][0] + coef[c][1] * q[0] + coef[c][2] * q[1] +
                coef[c][3] * q[2];
    const double *const gv[6] = {&m.g.xx, &m.g.xy, &m.g.xz,
                                 &m.g.yy, &m.g.yz, &m.g.zz};
    for (int c = 0; c < 6; ++c) {
      got[c] = *gv[c];
      for (int d = 0; d < 3; ++d) {
        const double *const dgv[6] = {&m.dg[d].xx, &m.dg[d].xy, &m.dg[d].xz,
                                      &m.dg[d].yy, &m.dg[d].yz, &m.dg[d].zz};
        dgot[c][d] = *dgv[c];
      }
    }
    got[6] = m.alpha;
    for (int d = 0; d < 3; ++d)
      dgot[6][d] = m.dalpha[d];
    for (int k = 0; k < 3; ++k) {
      got[7 + k] = m.beta[k];
      for (int d = 0; d < 3; ++d)
        dgot[7 + k][d] = m.dbeta[d][k];
    }
    for (int c = 0; c < 10; ++c) {
      values_ok =
          values_ok && detail::approx_eq(got[c], want[c], detail::rtol_machine);
      for (int d = 0; d < 3; ++d)
        derivs_ok = derivs_ok && detail::approx_eq(dgot[c][d], coef[c][1 + d],
                                                   detail::rtol_machine);
    }
  }
  b.add_boolean("geo.gather.linear_values", values_ok);
  b.add_boolean("geo.gather.linear_derivs", derivs_ok);

  // --- Row 81: constant fields — derivatives identically zero and values
  // reproduced bitwise (exact tier; the mechanism behind flat exactness). ---
  double const_coef[10][4] = {};
  for (int c = 0; c < 10; ++c)
    const_coef[c][0] = coef[c][0];
  SyntheticVertexData cdata;
  cdata.fill_affine(const_coef);
  const VertexMetricGather cgather = cdata.gather();
  bool const_ok = true;
  for (const auto &q : queries) {
    const MetricSnapshot m = cgather(q[0], q[1], q[2]);
    const double *const gv[6] = {&m.g.xx, &m.g.xy, &m.g.xz,
                                 &m.g.yy, &m.g.yz, &m.g.zz};
    for (int c = 0; c < 6; ++c)
      const_ok = const_ok && *gv[c] == coef[c][0];
    const_ok = const_ok && m.alpha == coef[6][0];
    for (int k = 0; k < 3; ++k)
      const_ok = const_ok && m.beta[k] == coef[7 + k][0];
    for (int d = 0; d < 3; ++d) {
      const_ok = const_ok && m.dalpha[d] == 0.0 && m.dg[d].xx == 0.0 &&
                 m.dg[d].xy == 0.0 && m.dg[d].xz == 0.0 && m.dg[d].yy == 0.0 &&
                 m.dg[d].yz == 0.0 && m.dg[d].zz == 0.0;
      for (int k = 0; k < 3; ++k)
        const_ok = const_ok && m.dbeta[d][k] == 0.0;
    }
  }
  b.add_boolean("geo.gather.constant_derivs_zero", const_ok);

  // --- Row 82: null-closure identity alpha^2 (p^t)^2 = gamma^{ij} p_i p_j
  // at EVERY RHS evaluation (all four RK stages of every step) along a
  // curved trajectory, machine tier. ---
  {
    const CurvedAnalyticMetric curved;
    bool closure_ok = true;
    int nstages = 0;
    const auto observe = [&](const MetricSnapshot &m, const double *p,
                             const GeodesicRhs &r) {
      const InverseSpatialMetric gu = spatial_metric_inverse(m.g);
      const double lhs = (m.alpha * r.pt) * (m.alpha * r.pt);
      const double rhs = momentum_norm2(gu, p[0], p[1], p[2]);
      closure_ok = closure_ok && detail::approx_eq(lhs, rhs, detail::rtol_machine);
      ++nstages;
    };
    double x[3] = {0.1, -0.2, 0.15};
    double p[3] = {0.6, 0.3, -0.4};
    for (int step = 0; step < 8; ++step)
      geodesic_step_with(curved, x, p, 0.125, detail::SqrtStd{}, observe);
    b.add_boolean("geo.rhs.null_closure_identity", closure_ok && nstages == 32);
  }

  // --- Rows 83..84: flat spacetime. Straight line x(t) = x0 + t p/|p| to
  // rtol_machine over the benchmark's four steps with p bitwise unchanged;
  // |dx/dt| = 1 to rtol_machine for several momenta. ---
  {
    const detail::FlatGather flat;
    const double x0[3] = {0.125, -0.25, 0.0625};
    const double p0[3] = {0.5, -0.25, 0.125};
    double x[3] = {x0[0], x0[1], x0[2]};
    double p[3] = {p0[0], p0[1], p0[2]};
    const double dt = 0.125;
    const int nsteps = 4;
    bool line_ok = true;
    const double pnorm = std::sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]);
    for (int step = 1; step <= nsteps; ++step) {
      geodesic_step(flat, x, p, dt);
      const double t = step * dt;
      for (int i = 0; i < 3; ++i)
        line_ok = line_ok &&
                  detail::approx_eq(x[i], x0[i] + t * p0[i] / pnorm,
                                    detail::rtol_machine) &&
                  p[i] == p0[i];
    }
    b.add_boolean("geo.push.flat_straight_line", line_ok);

    const double momenta[][3] = {{1.0, 0.0, 0.0},     {0.0, -2.0, 0.0},
                                 {0.0, 0.0, 0.25},    {0.5, -0.25, 0.125},
                                 {-0.375, 0.75, 0.5}, {3.0, 4.0, 12.0}};
    bool speed_ok = true;
    for (const auto &pm : momenta) {
      const GeodesicRhs r = geodesic_rhs(flat_metric_snapshot(), pm);
      const double v2 = r.dxdt[0] * r.dxdt[0] + r.dxdt[1] * r.dxdt[1] +
                        r.dxdt[2] * r.dxdt[2];
      speed_ok = speed_ok && detail::approx_eq(std::sqrt(v2), 1.0, detail::rtol_machine) &&
                 r.dpdt[0] == 0.0 && r.dpdt[1] == 0.0 && r.dpdt[2] == 0.0;
    }
    b.add_boolean("geo.push.flat_unit_speed", speed_ok);
  }

  // --- Row 85: convergence order >= 2 on the curved analytic metric:
  // errors against a 32x finer reference at dt and dt/2 must shrink by at
  // least 3.9 (RK4 gives ~16). ---
  {
    const CurvedAnalyticMetric curved;
    const double x0[3] = {0.1, -0.2, 0.15};
    const double p0[3] = {0.6, 0.3, -0.4};
    const double T = 1.0;
    const auto integrate = [&](int nsteps, double out[6]) {
      double x[3] = {x0[0], x0[1], x0[2]};
      double p[3] = {p0[0], p0[1], p0[2]};
      for (int step = 0; step < nsteps; ++step)
        geodesic_step(curved, x, p, T / nsteps);
      for (int i = 0; i < 3; ++i) {
        out[i] = x[i];
        out[3 + i] = p[i];
      }
    };
    double ref[6], coarse[6], fine[6];
    integrate(128, ref);
    integrate(4, coarse);
    integrate(8, fine);
    double e_coarse = 0.0, e_fine = 0.0;
    for (int i = 0; i < 6; ++i) {
      e_coarse = std::max(e_coarse, detail::cabs(coarse[i] - ref[i]));
      e_fine = std::max(e_fine, detail::cabs(fine[i] - ref[i]));
    }
    const bool order_ok = e_fine > 0.0 && e_coarse / e_fine >= 3.9;
    CCTK_VINFO("MCNuX selftest geo.push.order_ge2: error(dt) = %.3e, "
               "error(dt/2) = %.3e, ratio = %.3f",
               e_coarse, e_fine, e_fine > 0.0 ? e_coarse / e_fine : 0.0);
    b.add_boolean("geo.push.order_ge2", order_ok);
  }
}

// ---------------------------------------------------------------------------
// The battery itself. Every check is a call into the shared headers; no
// constant, fixture, or tolerance is restated here.
// ---------------------------------------------------------------------------
void run_battery(Battery &b) {
  // Rows 0..2 — Philox4x32-10 known-answer vectors, exact word equality.
  // rng-and-statistical-acceptance.md [MCNX-RNG-01], [MCNX-RNG-02].
  const char *const kat_names[detail::num_kat_vectors] = {
      "rng.philox.kat1", "rng.philox.kat2", "rng.philox.kat3"};
  for (int i = 0; i < detail::num_kat_vectors; ++i)
    b.add_boolean(kat_names[i], detail::kat_holds(i));

  // Rows 3..6 — the [0, 1) construction's closed forms, exact.
  // [MCNX-RNG-04].
  const char *const cf_names[detail::num_uniform_closed_forms] = {
      "rng.uniform.u64_zero", "rng.uniform.u64_2p11_minus_1",
      "rng.uniform.u64_2p11", "rng.uniform.u64_max"};
  for (int i = 0; i < detail::num_uniform_closed_forms; ++i) {
    const detail::UniformClosedForm f = detail::uniform_closed_form(i);
    b.add_exact(cf_names[i], uniform_from_u64(f.u64), f.r);
  }

  // Row 7 — u(S, q, e, k) is a pure function of its arguments: repeated,
  // permuted, reversed, and batched evaluation all agree bitwise.
  // [MCNX-RNG-05], [MCNX-RNG-06].
  b.add_boolean("rng.purity", detail::purity_holds());

  // Row 8 — draws k = 2b and k = 2b+1 share Philox block b as word pairs 0
  // and 1, over the counter/key packing of [MCNX-RNG-03].
  b.add_boolean("rng.pairing", detail::pairing_holds());

  // Rows 9..20 — each derived conversion factor, recomputed from the five
  // defining constants, against the spec's pinned table anchor.
  // conventions-and-units.md [MCNX-CNV-03], [MCNX-CNV-04].
  b.add_pinned("units.factor.length_code_to_cgs", length_code_to_cgs,
               pinned::length_code_to_cgs);
  b.add_pinned("units.factor.time_code_to_cgs", time_code_to_cgs,
               pinned::time_code_to_cgs);
  b.add_pinned("units.factor.mass_code_to_cgs", mass_code_to_cgs,
               pinned::mass_code_to_cgs);
  b.add_pinned("units.factor.energy_code_to_cgs", energy_code_to_cgs,
               pinned::energy_code_to_cgs);
  b.add_pinned("units.factor.energy_code_to_MeV", energy_code_to_MeV,
               pinned::energy_code_to_MeV);
  b.add_pinned("units.factor.mass_density_code_to_cgs",
               mass_density_code_to_cgs, pinned::mass_density_code_to_cgs);
  b.add_pinned("units.factor.energy_density_code_to_cgs",
               energy_density_code_to_cgs, pinned::energy_density_code_to_cgs);
  b.add_pinned("units.factor.opacity_code_to_cgs", opacity_code_to_cgs,
               pinned::opacity_code_to_cgs);
  b.add_pinned("units.factor.rate_code_to_cgs", rate_code_to_cgs,
               pinned::rate_code_to_cgs);
  b.add_pinned("units.factor.emissivity_code_to_cgs", emissivity_code_to_cgs,
               pinned::emissivity_code_to_cgs);
  b.add_pinned("units.factor.mev_to_kelvin", mev_to_kelvin,
               pinned::mev_to_kelvin);
  b.add_pinned("units.factor.k_B_MeV_per_K", k_B_MeV_per_K,
               pinned::k_B_MeV_per_K);

  // Rows 21..32 — cgs -> code -> cgs (and back) is the identity to the machine
  // tier for the same factors, probed over several decades. [MCNX-CNV-02].
  b.add_boolean("units.roundtrip.length_code_to_cgs",
                detail::round_trips(length_code_to_cgs));
  b.add_boolean("units.roundtrip.time_code_to_cgs",
                detail::round_trips(time_code_to_cgs));
  b.add_boolean("units.roundtrip.mass_code_to_cgs",
                detail::round_trips(mass_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_code_to_cgs",
                detail::round_trips(energy_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_code_to_MeV",
                detail::round_trips(energy_code_to_MeV));
  b.add_boolean("units.roundtrip.mass_density_code_to_cgs",
                detail::round_trips(mass_density_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_density_code_to_cgs",
                detail::round_trips(energy_density_code_to_cgs));
  b.add_boolean("units.roundtrip.opacity_code_to_cgs",
                detail::round_trips(opacity_code_to_cgs));
  b.add_boolean("units.roundtrip.rate_code_to_cgs",
                detail::round_trips(rate_code_to_cgs));
  b.add_boolean("units.roundtrip.emissivity_code_to_cgs",
                detail::round_trips(emissivity_code_to_cgs));
  b.add_boolean("units.roundtrip.mev_to_kelvin",
                detail::round_trips(mev_to_kelvin));
  b.add_boolean("units.roundtrip.k_B_MeV_per_K",
                detail::round_trips(k_B_MeV_per_K));

  // Row 33 — the species enumeration, degeneracies (including nu_x's g = 4),
  // and lepton numbers, exact. [MCNX-CNV-06].
  b.add_boolean("units.species_table", detail::species_table_holds());

  // Rows 34..35 — the binary64 precision gate of
  // specs/build-and-integration.md [MCNX-BLD-03], measured at runtime rather
  // than only asserted in stub.cxx.
  b.add_exact("precision.sizeof_double", double(sizeof(double)), 8.0);
  b.add_exact("precision.sizeof_cctk_real", double(sizeof(CCTK_REAL)), 8.0);

  // Rows 36..37 — the AMReX legs of the same gate, per the coverage matrix of
  // specs/verification-suite-design.md (both sizeof(amrex::Real) and
  // sizeof(amrex::ParticleReal)). Compile-time asserted in
  // mcnux_particles.hxx; measured here as golden rows.
  b.add_exact("precision.sizeof_amrex_real", double(sizeof(amrex::Real)), 8.0);
  b.add_exact("precision.sizeof_amrex_particle_real",
              double(sizeof(amrex::ParticleReal)), 8.0);

  // Rows 38..52 — the WeakLibInterp call-boundary wrapper checks of
  // specs/opacity-eos-evaluation.md (kernel parity on synthetic tables,
  // unit-boundary conversions, species-axis grounding, T11 bounds hooks).
  append_opacity_rows(b);

  // Rows 53..59 — the tetrad/null-closure helper checks of
  // specs/geodesic-propagation.md [MCNX-GEO-02] and
  // specs/packet-representation-and-sampling.md [MCNX-PKT-04] on synthetic
  // flat-static and curved inputs.
  append_tetrad_rows(b);

  // Rows 60..65 — the source-term ledger sign fixtures of
  // specs/hydro-coupling-source-terms.md [MCNX-HYD-02] through the
  // mcnux_srcterms.hxx helpers.
  append_srcterm_rows(b);

  // Rows 66..70 — the trapped-regime relabeling invariants, convention
  // tripwire, and alpha = 1 identity of specs/trapped-regime-treatment.md
  // [MCNX-TRP-02/03] through the mcnux_trp.hxx pure map.
  append_trp_rows(b);

  // Rows 71..78 — the hc pin, the analytic (gray) opacity mode, and the
  // source-agnostic coefficient dispatch of
  // specs/verification-suite-design.md [MCNX-VER-05] through
  // mcnux_coefficients.hxx and its parameter glue.
  append_coefficient_rows(b);

  // Rows 79..85 — the geodesic gather/RHS/integrator legs of
  // specs/geodesic-propagation.md [MCNX-GEO-02/03/04] through
  // mcnux_geodesic.hxx on synthetic vertex data and analytic metrics.
  append_geodesic_rows(b);
}

// Size of the MCNuX::mcnux_selftest array as declared in interface.ccl.
int selftest_array_size() {
  const int gi = CCTK_GroupIndex("MCNuX::mcnux_selftest");
  if (gi < 0)
    CCTK_VERROR("MCNuX::mcnux_selftest is not a known group");
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("MCNuX::mcnux_selftest must be a 1-dimensional grid array");
  return int(*sizes[0]);
}

} // namespace

extern "C" void MCNuX_SelfTest(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SelfTest;
  DECLARE_CCTK_PARAMETERS;

  Battery battery;
  run_battery(battery);

  const int nrows = battery.size();
  const int nslots = selftest_array_size();
  if (nrows != nslots)
    CCTK_VERROR("MCNuX self-test battery has %d rows, but "
                "MCNuX::mcnux_selftest is declared with SIZE=%d in "
                "MCNuX/interface.ccl; the two must agree (and the golden data "
                "in MCNuX/test/unit-selftest/ must be regenerated whenever "
                "rows are appended)",
                nrows, nslots);

  int nfailed = 0;
  for (int i = 0; i < nrows; ++i) {
    const CheckRow &row = battery[i];
    selftest_error[i] = row.error;
    selftest_pass[i] = row.pass ? 1.0 : 0.0;

    // The TSV table is numeric only (it is diffed column-wise by the harness),
    // so the index -> name mapping is echoed to the log instead.
    CCTK_VINFO("MCNuX selftest row %3d  %-40s  error = %.17g  %s", i, row.name,
               row.error, row.pass ? "PASS" : "FAIL");
    if (!row.pass)
      ++nfailed;
  }

  if (nfailed > 0)
    CCTK_VWARN(CCTK_WARN_ALERT,
               "MCNuX unit self-test battery: %d of %d checks FAILED (see the "
               "per-row lines above; the golden diff of the `unit-selftest` "
               "benchmark is the scoring mechanism)",
               nfailed, nrows);
  else
    CCTK_VINFO("MCNuX unit self-test battery: all %d checks passed", nrows);
}

} // namespace MCNuX
