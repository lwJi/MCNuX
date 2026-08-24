// ---------------------------------------------------------------------------
// Rows 71..78 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the analytic (gray) opacity mode and the
// source-agnostic coefficient dispatch of specs/verification-suite-design.md
// [MCNX-VER-05] (T12) through mcnux_coefficients.hxx, plus the hc pin those
// formulas rest on (specs/opacity-eos-evaluation.md:356-361).
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

#include "mcnux_selftest.hxx"

#include "mcnux_coefficients.hxx"
#include "mcnux_units.hxx"

#include <cmath>

namespace MCNuX {

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

} // namespace MCNuX
