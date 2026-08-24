// ---------------------------------------------------------------------------
// Rows 92..99 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the table-range enforcement policy of
// specs/opacity-eos-evaluation.md [MCNX-OPA-06] and the EOS-inversion error
// protocol [MCNX-OPA-07] (T11) through the RangedTableCoefficients wrapper
// of mcnux_table_range.hxx.
//
// Fixture: the T10 synthetic view bundle (same axes and generator phases as
// rows 86..91) built TWICE — once with real data for the clamp rows, once
// with NaN-POISONED table data (axis arrays real) for the transparency-floor
// rows, so any accidental table evaluation below the floor surfaces as NaN
// and fails the row (this is how "no table evaluation occurs" is made
// falsifiable, opacity-eos-evaluation.md:305-306). Fixture bounds: E in
// [10^-1, 10^2] MeV, rho intersection [10^6, 6e11] g/cm^3 (EmAb/Iso floor
// 10^6 dominates the EOS 1e3), T intersection [10^9.5 K, 1e11 K] in MeV,
// Ye intersection [0.05, 0.50]. All rows exact/boolean tier. The inversion
// row reuses the row-39 synthetic inversion fixture with a pinned FAILING
// input (X above MaxX => error 2, T = 0).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_coefficients.hxx"
#include "mcnux_table_coeffs.hxx"
#include "mcnux_table_range.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <cstddef>
#include <limits>

namespace MCNuX {

void append_table_range_rows(Battery &b) {
  using wli::Real;

  // --- Axes (the T10 grids verbatim). ---
  constexpr int oE = 4, oD = 5, oT = 4, oY = 3, oMom = 2;
  const Real LogEs[oE] = {-1.0, 0.5, 1.5, 2.0};
  const Real LogDs[oD] = {6.0, 8.0, 10.0, 12.0, 14.0};
  const Real LogTs[oT] = {9.5, 10.0, 10.7, 11.3};
  const Real Ys2[oY] = {0.05, 0.25, 0.50};
  constexpr int nD = 4, nT = 5, nY = 3;
  const Real Ds[nD] = {1.0e3, 5.0e5, 2.0e8, 6.0e11};
  const Real Ts_K[nT] = {1.0e9, 3.0e9, 1.0e10, 5.0e10, 1.0e11};
  const Real Ys[nY] = {0.05, 0.30, 0.55};

  // --- Real and NaN-poisoned datasets (same shapes; the poisoned bundle
  // shares the REAL axis arrays, so rho_min_cgs still works on it). ---
  constexpr std::size_t n_emab = std::size_t(oE) * oD * oT * oY;
  constexpr std::size_t n_iso = std::size_t(oE) * oMom * oD * oT * oY;
  constexpr std::size_t n_eos = std::size_t(nD) * nT * nY;
  static Real emab_nue[n_emab], emab_nueb[n_emab];
  static Real iso_nue[n_iso], iso_nueb[n_iso];
  static Real mu_e_tab[n_eos], mu_p_tab[n_eos], mu_n_tab[n_eos];
  static Real emab_nan[n_emab], iso_nan[n_iso], eos_nan[n_eos];
  const Real qnan = std::numeric_limits<Real>::quiet_NaN();
  for (std::size_t k = 0; k < n_emab; ++k) {
    emab_nue[k] = synth_tval(k);
    emab_nueb[k] = synth_tval(k + 5);
    emab_nan[k] = qnan;
  }
  for (std::size_t k = 0; k < n_iso; ++k) {
    iso_nue[k] = synth_tval(k + 7);
    iso_nueb[k] = synth_tval(k + 13);
    iso_nan[k] = qnan;
  }
  for (std::size_t k = 0; k < n_eos; ++k) {
    mu_e_tab[k] = synth_tval(k + 2);
    mu_p_tab[k] = synth_tval(k + 17);
    mu_n_tab[k] = synth_tval(k + 23);
    eos_nan[k] = qnan;
  }
  const Real iso_offsets[num_tabulated_species * oMom] = {0.5, 1.0, 1.5, 2.0};
  const Real eos_OS[num_chemical_potentials] = {1.5, 2.5, 3.0};

  // View-bundle builder: axes fixed, datasets swappable.
  const auto make_views = [&](Real const *ea0, Real const *ea1, Real const *is0,
                              Real const *is1, Real const *me, Real const *mp,
                              Real const *mn) {
    BaselineTableViews v{};
    v.emab_LogEs = LogEs;
    v.emab_nE = oE;
    v.emab_LogDs = LogDs;
    v.emab_nD = oD;
    v.emab_LogTs = LogTs;
    v.emab_nT = oT;
    v.emab_Ys = Ys2;
    v.emab_nY = oY;
    v.emab_OS[emab_species_slot(Species::NuE)] = 1.25;
    v.emab_OS[emab_species_slot(Species::NuEBar)] = 0.75;
    v.emab_table[emab_species_slot(Species::NuE)] = ea0;
    v.emab_table[emab_species_slot(Species::NuEBar)] = ea1;
    v.iso_LogEs = LogEs;
    v.iso_nE = oE;
    v.iso_LogDs = LogDs;
    v.iso_nD = oD;
    v.iso_LogTs = LogTs;
    v.iso_nT = oT;
    v.iso_Ys = Ys2;
    v.iso_nY = oY;
    v.iso_nMom = oMom;
    v.iso_offsets = iso_offsets;
    v.iso_table[iso_species_slot(Species::NuE)] = is0;
    v.iso_table[iso_species_slot(Species::NuEBar)] = is1;
    v.eos_Ds = Ds;
    v.eos_nD = nD;
    v.eos_Ts_K = Ts_K;
    v.eos_nT = nT;
    v.eos_Ys = Ys;
    v.eos_nY = nY;
    for (int m = 0; m < num_chemical_potentials; ++m)
      v.eos_OS[m] = eos_OS[m];
    v.eos_table[chemical_potential_index(ChemicalPotential::Mu_e)] = me;
    v.eos_table[chemical_potential_index(ChemicalPotential::Mu_p)] = mp;
    v.eos_table[chemical_potential_index(ChemicalPotential::Mu_n)] = mn;
    return v;
  };
  const BaselineTableViews views = make_views(
      emab_nue, emab_nueb, iso_nue, iso_nueb, mu_e_tab, mu_p_tab, mu_n_tab);
  const BaselineTableViews views_poisoned = make_views(
      emab_nan, emab_nan, iso_nan, iso_nan, eos_nan, eos_nan, eos_nan);
  const BaselineTableCoefficients base{views};

  // Pinned in-range probe (the T10 values) and the floor probe (below the
  // 10^6 g/cm^3 EmAb/Iso density floor).
  const double E_q = 10.5;
  const FluidState st_in{3.3e11, 2.3, 0.31};
  const FluidState st_floor{1.0e5, 2.3, 0.31};

  // Rows 92..93 — transparency floor, exact zero with NO table touch: the
  // poisoned bundle turns any accidental evaluation into NaN (which would
  // fail the exact comparison). NuX proves the floor precedes the nu_x
  // branch and its two unconditional Iso evaluations.
  {
    const RangedTableCoefficients ranged{
        BaselineTableCoefficients{views_poisoned}, nullptr};
    const Coefficients ce = ranged(Species::NuE, E_q, st_floor);
    const Coefficients cx = ranged(Species::NuX, E_q, st_floor);
    b.add_exact("tblrange.floor.zero_no_touch",
                detail::cabs(ce.kappa_a) + detail::cabs(ce.kappa_s) +
                    detail::cabs(ce.eta),
                0.0);
    b.add_exact("tblrange.floor.nux",
                detail::cabs(cx.kappa_a) + detail::cabs(cx.kappa_s) +
                    detail::cabs(cx.eta),
                0.0);
  }

  // All-axes-out probe: E above max (250 > 100), rho above max
  // (1e13 > 6e11), T below min (0.1 MeV < 10^9.5 K in MeV), Ye above max
  // (0.9 > 0.50).
  const double E_out = 250.0;
  const FluidState st_out{1.0e13, 0.1, 0.9};

  // Row 94 — clamp equality, bitwise (same code path): the ranged functor at
  // the out-of-bounds probe equals the ranged functor at the hand-clamped
  // coordinates, for an electron species and for nu_x. counters == nullptr
  // exercises the nullable skip-counting leg.
  {
    const RangedTableCoefficients ranged{base, nullptr};
    const FluidState st_cl{density_bounds_cgs(views).hi,
                           temperature_bounds_mev(views).lo,
                           ye_bounds(views).hi};
    const double E_cl = energy_bounds_mev(views).hi;
    const Coefficients g_e = ranged(Species::NuE, E_out, st_out);
    const Coefficients w_e = ranged(Species::NuE, E_cl, st_cl);
    const Coefficients g_x = ranged(Species::NuX, E_out, st_out);
    const Coefficients w_x = ranged(Species::NuX, E_cl, st_cl);
    b.add_boolean("tblrange.clamp.equality",
                  g_e.kappa_a == w_e.kappa_a && g_e.kappa_s == w_e.kappa_s &&
                      g_e.eta == w_e.eta && g_x.kappa_a == w_x.kappa_a &&
                      g_x.kappa_s == w_x.kappa_s && g_x.eta == w_x.eta);
  }

  // Row 95 — per-axis counters: ONE evaluation at the all-axes-out probe
  // increments each of the four counters exactly once.
  {
    ClampCounters cc;
    const RangedTableCoefficients ranged{base, &cc};
    (void)ranged(Species::NuE, E_out, st_out);
    b.add_boolean("tblrange.clamp.counters",
                  cc.E == 1 && cc.rho == 1 && cc.T == 1 && cc.Ye == 1);
  }

  // Row 96 — in-range no-op: at the pinned in-range probe the ranged result
  // is bitwise the plain BaselineTableCoefficients result for all three
  // species, and every counter stays 0.
  {
    ClampCounters cc;
    const RangedTableCoefficients ranged{base, &cc};
    const Species species[NUM_SPECIES] = {Species::NuE, Species::NuEBar,
                                          Species::NuX};
    bool ok = true;
    for (int is = 0; is < NUM_SPECIES; ++is) {
      const Coefficients got = ranged(species[is], E_q, st_in);
      const Coefficients want = base(species[is], E_q, st_in);
      ok = ok && got.kappa_a == want.kappa_a && got.kappa_s == want.kappa_s &&
           got.eta == want.eta;
    }
    ok = ok && cc.E == 0 && cc.rho == 0 && cc.T == 0 && cc.Ye == 0;
    b.add_boolean("tblrange.clamp.inrange_noop", ok);
  }

  // Row 97 — no-NaN guarantee: a hostile probe (E <= 0 and T <= 0, which
  // unguarded would reach a log10; huge rho; Ye > 1) yields finite, non-NaN
  // coefficients for an electron species and for nu_x ([MCNX-OPA-03]'s
  // no-negative-into-log10 closed by the clamp).
  {
    const RangedTableCoefficients ranged{base, nullptr};
    const FluidState st_bad{1.0e30, -2.0, 1.5};
    const Coefficients ce = ranged(Species::NuE, -3.0, st_bad);
    const Coefficients cx = ranged(Species::NuX, -3.0, st_bad);
    b.add_boolean("tblrange.nonan",
                  !std::isnan(ce.kappa_a) && !std::isnan(ce.kappa_s) &&
                      !std::isnan(ce.eta) && !std::isnan(cx.kappa_a) &&
                      !std::isnan(cx.kappa_s) && !std::isnan(cx.eta));
  }

  // Row 98 — the transparency-floor density is the max of the tabulated
  // density lower bounds (EmAb/Iso 10^6 dominates the EOS 1e3), computed
  // independently here; bitwise (same pow/fmax calls).
  b.add_exact("tblrange.rho_min", rho_min_cgs(views),
              std::fmax(std::pow(10.0, LogDs[0]), Ds[0]));

  // Row 99 — [MCNX-OPA-07] inversion protocol on a pinned FAILING input:
  // the row-39 synthetic inversion fixture with X = 5.0 above MaxX = 1.0
  // => error != 0 (code 2), returned T_MeV == 0 exactly (not a
  // temperature), and inversion_usable rejects it.
  {
    Real inv_table[nD * nT * nY];
    for (int iY = 0; iY < nY; ++iY)
      for (int iT = 0; iT < nT; ++iT)
        for (int iD = 0; iD < nD; ++iD)
          inv_table[iD + nD * (iT + nT * iY)] =
              0.08 * iT + 0.01 * iD + 0.005 * iY;
    const Real OS_inv = 2.0;
    const wli::EosInversionBounds bounds{1.0e3, 6.0e11, -1.0, 1.0,
                                         0.05,  0.55,   true};
    const EosInversionResultMeV r =
        eos_invert_dey_noguess(7.3e6, 5.0, 0.22, Ds, nD, Ts_K, nT, Ys, nY,
                               OS_inv, inv_table, bounds);
    b.add_boolean("tblrange.inversion_protocol",
                  r.error != 0 && r.T_MeV == 0.0 && !inversion_usable(r));
  }
}

} // namespace MCNuX
