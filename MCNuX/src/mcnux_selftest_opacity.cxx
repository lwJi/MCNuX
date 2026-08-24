// ---------------------------------------------------------------------------
// Rows 38..52 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the WeakLibInterp call-boundary wrapper layer of
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

#include "mcnux_selftest.hxx"

#include "mcnux_opacity.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <cstddef>

namespace MCNuX {

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

} // namespace MCNuX
