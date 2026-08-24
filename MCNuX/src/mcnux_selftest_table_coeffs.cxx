// ---------------------------------------------------------------------------
// Rows 86..91 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the baseline table-source coefficient assembly and
// nu_x mapping of specs/opacity-eos-evaluation.md [MCNX-OPA-04/05] (T10)
// through the BaselineTableCoefficients callable of mcnux_table_coeffs.hxx.
//
// Fixture: synthetic in-memory tables only, the T9 pattern verbatim
// (synth_tval generator, flat column-major arrays; already-log10'd axes for
// EmAb/Iso, raw axes — Ts in KELVIN — for the EOS chemical-potential
// sub-tables; NO HDF5, no file I/O). One EmAb dataset + offset per
// electron-type slot, one Iso 5D dataset per slot plus the species-major
// [nOpacities, nMoments] offsets array, and three EOS sub-tables sharing one
// raw axis set (one per mu_e/mu_p/mu_n). The pinned query point reuses the
// T9 in-range values (E = 10.5 MeV; rho = 3.3e11 g/cm^3, T = 2.3 MeV,
// Ye = 0.31 — in range for both the EmAb/Iso grids and the EOS axes).
//
// Rows: the Kirchhoff closure eta = c kappa_a f_eq 4 pi E^3/(hc)^3 with
// mu_nu = +-(mu_e + mu_p - mu_n) is independently restated here (mu's via
// direct eos_evaluate calls, the analytic-leg idiom of row 73) and judged at
// detail::rtol_machine for nu_e and nu_e_bar; the nu_x zeros are exact
// (kappa_a = eta = 0 from tabulated physics, [MCNX-OPA-05]); the nu_x
// half-mean equals 1/2 (iso_evaluate(nu_e) + iso_evaluate(nu_e_bar))
// computed directly, machine tier; and evaluate_coefficients() under
// CoefficientSource::Table returns the callable's triple bitwise for all
// three species (the table source plugs into the T12 dispatch). No range
// policy is exercised (T11) and no NES/Pair/Brem is assembled
// (opacity-eos-evaluation.md:340-349).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_coefficients.hxx"
#include "mcnux_table_coeffs.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <cstddef>

namespace MCNuX {

void append_baseline_table_rows(Battery &b) {
  using wli::Real;

  // --- EmAb/Iso synthetic grids: already-log10'd E [MeV], rho [g/cm^3],
  // T [KELVIN]; Ye raw (the T9 grids of rows 41..42). ---
  constexpr int oE = 4, oD = 5, oT = 4, oY = 3, oMom = 2;
  const Real LogEs[oE] = {-1.0, 0.5, 1.5, 2.0};
  const Real LogDs[oD] = {6.0, 8.0, 10.0, 12.0, 14.0};
  const Real LogTs[oT] = {9.5, 10.0, 10.7, 11.3};
  const Real Ys2[oY] = {0.05, 0.25, 0.50};

  // One EmAb dataset + scalar offset per electron-type slot (distinct
  // generator phases so the two species cannot alias).
  Real emab_nue[oE * oD * oT * oY], emab_nueb[oE * oD * oT * oY];
  for (std::size_t k = 0; k < std::size_t(oE) * oD * oT * oY; ++k) {
    emab_nue[k] = synth_tval(k);
    emab_nueb[k] = synth_tval(k + 5);
  }

  // One Iso 5D dataset per slot; species-major [nOpacities, nMoments]
  // offsets (selected only through iso_offset inside the callable).
  Real iso_nue[oE * oMom * oD * oT * oY], iso_nueb[oE * oMom * oD * oT * oY];
  for (std::size_t k = 0; k < std::size_t(oE) * oMom * oD * oT * oY; ++k) {
    iso_nue[k] = synth_tval(k + 7);
    iso_nueb[k] = synth_tval(k + 13);
  }
  const Real iso_offsets[num_tabulated_species * oMom] = {0.5, 1.0, 1.5, 2.0};

  // --- EOS chemical-potential sub-tables: raw axes (rho, T [KELVIN — table
  // native], Ye), the row-38 fixture axes; one (offset, dataset) per mu. ---
  constexpr int nD = 4, nT = 5, nY = 3;
  const Real Ds[nD] = {1.0e3, 5.0e5, 2.0e8, 6.0e11};
  const Real Ts_K[nT] = {1.0e9, 3.0e9, 1.0e10, 5.0e10, 1.0e11};
  const Real Ys[nY] = {0.05, 0.30, 0.55};
  Real mu_e_tab[nD * nT * nY], mu_p_tab[nD * nT * nY], mu_n_tab[nD * nT * nY];
  for (std::size_t k = 0; k < std::size_t(nD) * nT * nY; ++k) {
    mu_e_tab[k] = synth_tval(k + 2);
    mu_p_tab[k] = synth_tval(k + 17);
    mu_n_tab[k] = synth_tval(k + 23);
  }
  const Real eos_OS[num_chemical_potentials] = {1.5, 2.5, 3.0};

  // --- The view bundle + callable. ---
  BaselineTableViews views{};
  views.emab_LogEs = LogEs;
  views.emab_nE = oE;
  views.emab_LogDs = LogDs;
  views.emab_nD = oD;
  views.emab_LogTs = LogTs;
  views.emab_nT = oT;
  views.emab_Ys = Ys2;
  views.emab_nY = oY;
  views.emab_OS[emab_species_slot(Species::NuE)] = 1.25;
  views.emab_OS[emab_species_slot(Species::NuEBar)] = 0.75;
  views.emab_table[emab_species_slot(Species::NuE)] = emab_nue;
  views.emab_table[emab_species_slot(Species::NuEBar)] = emab_nueb;
  views.iso_LogEs = LogEs;
  views.iso_nE = oE;
  views.iso_LogDs = LogDs;
  views.iso_nD = oD;
  views.iso_LogTs = LogTs;
  views.iso_nT = oT;
  views.iso_Ys = Ys2;
  views.iso_nY = oY;
  views.iso_nMom = oMom;
  views.iso_offsets = iso_offsets;
  views.iso_table[iso_species_slot(Species::NuE)] = iso_nue;
  views.iso_table[iso_species_slot(Species::NuEBar)] = iso_nueb;
  views.eos_Ds = Ds;
  views.eos_nD = nD;
  views.eos_Ts_K = Ts_K;
  views.eos_nT = nT;
  views.eos_Ys = Ys;
  views.eos_nY = nY;
  for (int m = 0; m < num_chemical_potentials; ++m)
    views.eos_OS[m] = eos_OS[m];
  views.eos_table[chemical_potential_index(ChemicalPotential::Mu_e)] = mu_e_tab;
  views.eos_table[chemical_potential_index(ChemicalPotential::Mu_p)] = mu_p_tab;
  views.eos_table[chemical_potential_index(ChemicalPotential::Mu_n)] = mu_n_tab;
  const BaselineTableCoefficients table{views};

  // Pinned in-range query point (the T9 values of rows 38/41).
  const double E_q = 10.5;
  const FluidState st{3.3e11, 2.3, 0.31};

  // Rows 86..87 — Kirchhoff closure, independently restated: mu's via three
  // DIRECT eos_evaluate calls, f_eq and the 4 pi E^3/(hc)^3 normalization
  // written out here, judged at detail::rtol_machine.
  {
    const double mu_e =
        eos_evaluate(st.rho_cgs, st.T_MeV, st.Ye, Ds, nD, Ts_K, nT, Ys, nY,
                     eos_OS[0], mu_e_tab);
    const double mu_p =
        eos_evaluate(st.rho_cgs, st.T_MeV, st.Ye, Ds, nD, Ts_K, nT, Ys, nY,
                     eos_OS[1], mu_p_tab);
    const double mu_n =
        eos_evaluate(st.rho_cgs, st.T_MeV, st.Ye, Ds, nD, Ts_K, nT, Ys, nY,
                     eos_OS[2], mu_n_tab);
    const double mu_nu_nue = mu_e + mu_p - mu_n; // mu_nu(nu_e_bar) = -this
    const struct {
      const char *name;
      Species s;
      double mu;
    } kirch[2] = {{"tblcoef.kirchhoff.nue", Species::NuE, mu_nu_nue},
                  {"tblcoef.kirchhoff.nuebar", Species::NuEBar, -mu_nu_nue}};
    for (const auto &kc : kirch) {
      const Coefficients c = table(kc.s, E_q, st);
      const double f_eq =
          1.0 / (std::exp((E_q - kc.mu) / st.T_MeV) + 1.0);
      const double want = c_cgs * c.kappa_a * f_eq * 4.0 * detail::pi * E_q *
                          E_q * E_q /
                          (hc_MeV_cm * hc_MeV_cm * hc_MeV_cm);
      b.add_boolean(kc.name,
                    c.kappa_a > 0.0 &&
                        detail::approx_eq(c.eta, want, detail::rtol_machine));
    }
  }

  // Rows 88..90 — the nu_x mapping: kappa_a and eta exactly zero (no table
  // is touched for them), kappa_s the half-mean of the two directly-evaluated
  // electron-type Iso values (machine tier, 1e-14 relative).
  {
    const Coefficients cx = table(Species::NuX, E_q, st);
    b.add_exact("tblcoef.nux.kappa_a_zero", cx.kappa_a, 0.0);
    b.add_exact("tblcoef.nux.eta_zero", cx.eta, 0.0);

    const int s0 = iso_species_slot(Species::NuE);
    const int s1 = iso_species_slot(Species::NuEBar);
    const Real os0 = iso_offset(iso_offsets, num_tabulated_species, oMom, s0,
                                iso_baseline_moment);
    const Real os1 = iso_offset(iso_offsets, num_tabulated_species, oMom, s1,
                                iso_baseline_moment);
    const double ks0 =
        iso_evaluate(E_q, st.rho_cgs, st.T_MeV, st.Ye, LogEs, oE, LogDs, oD,
                     LogTs, oT, Ys2, oY, iso_baseline_moment, oMom, os0,
                     iso_nue);
    const double ks1 =
        iso_evaluate(E_q, st.rho_cgs, st.T_MeV, st.Ye, LogEs, oE, LogDs, oD,
                     LogTs, oT, Ys2, oY, iso_baseline_moment, oMom, os1,
                     iso_nueb);
    b.add_boolean("tblcoef.nux.kappa_s_halfmean",
                  cx.kappa_s > 0.0 &&
                      detail::approx_eq(cx.kappa_s, 0.5 * (ks0 + ks1),
                                        detail::rtol_machine));
  }

  // Row 91 — dispatch plumb-through: evaluate_coefficients() under
  // CoefficientSource::Table returns the callable's triple bitwise for all
  // three species (the [MCNX-OPA-04] source plugs into the [MCNX-VER-05]
  // dispatch). The analytic parameter set is irrelevant under Table.
  {
    const AnalyticOpacityParams ap_unused = {};
    const Species species[NUM_SPECIES] = {Species::NuE, Species::NuEBar,
                                          Species::NuX};
    bool plumb_ok = true;
    for (int is = 0; is < NUM_SPECIES; ++is) {
      const Coefficients direct = table(species[is], E_q, st);
      const Coefficients via = evaluate_coefficients(
          CoefficientSource::Table, ap_unused, table, species[is], E_q, st);
      plumb_ok = plumb_ok && via.kappa_a == direct.kappa_a &&
                 via.kappa_s == direct.kappa_s && via.eta == direct.eta;
    }
    b.add_boolean("tblcoef.dispatch.table_plumb", plumb_ok);
  }
}

} // namespace MCNuX
