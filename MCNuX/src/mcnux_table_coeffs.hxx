#ifndef MCNUX_TABLE_COEFFS_HXX
#define MCNUX_TABLE_COEFFS_HXX

// The baseline table-source transport-coefficient assembly of
// specs/opacity-eos-evaluation.md [MCNX-OPA-04] (EmAb + Iso only, Kirchhoff
// closure for eta) and the nu_x dataset mapping of [MCNX-OPA-05] (exact
// tier), packaged as the table-source callable of the source-agnostic
// dispatch in mcnux_coefficients.hxx:
//
//   BaselineTableCoefficients::operator()(Species, double E_MeV,
//                                         const FluidState&) -> Coefficients
//
// is exactly the frozen TableEval contract of evaluate_coefficients(), so
// the functor plugs into the `table` slot under CoefficientSource::Table.
//
// Assembly ([MCNX-OPA-04], opacity-eos-evaluation.md:193-220), electron-type
// species s in {nu_e, nu_e_bar}:
//
//   kappa_a(s, E) = EmAb(s; E, rho, T, Ye)                        [cm^-1]
//   kappa_s(s, E) = Iso(s; E, rho, T, Ye; iMom = 0)               [cm^-1]
//   eta(s, E)     = c * kappa_a * f_eq(E; s) * 4 pi E^3/(hc)^3
//                                              [MeV cm^-3 s^-1 MeV^-1]
//   f_eq(E; s)    = 1/(exp((E - mu_nu(s))/T) + 1)   (E, mu_nu, T in MeV)
//   mu_nu(nu_e)   = mu_e + mu_p - mu_n,  mu_nu(nu_e_bar) = -mu_nu(nu_e)
//
// with mu_e, mu_p, mu_n the EOS table's chemical-potential datasets (MeV)
// evaluated at the same (rho, T, Ye). The Iso moment index is 0 (the l = 0
// Legendre coefficient, the total elastic rate); the l = 1 slice is not
// consumed by the baseline.
//
// nu_x mapping ([MCNX-OPA-05], opacity-eos-evaluation.md:222-239): nu_x has
// NO dataset in the production EmAb/Iso tables (emab_species_slot /
// iso_species_slot return no_table_slot for it — never a valid index), so it
// is special-cased BEFORE any table call:
//
//   kappa_a(nu_x, E) = 0        (from tabulated physics)
//   eta(nu_x, E)     = 0        (from tabulated physics)
//   kappa_s(nu_x, E) = 1/2 [kappa_s(nu_e, E) + kappa_s(nu_e_bar, E)]
//
// The half-mean is the exact lump-average of the production chirality map
// over the 2 nu + 2 nu_bar heavy species that nu_x aggregates. The
// degeneracy factor g = 4 enters exactly once, at packet creation — NEVER
// here (per-species coefficients are g = 1 by construction; mcnux_units.hxx
// species_degeneracy, mcnux_coefficients.hxx conventions block).
//
// Everything evaluates through the mcnux_opacity.hxx call-boundary wrappers
// (emab_evaluate / iso_evaluate / eos_evaluate / iso_offset and the species
// slot maps) — this header re-authors no WLI call, no unit conversion, and
// no log10 ([MCNX-OPA-01/02/03] stay owned by the wrapper layer). Constants
// come exclusively from mcnux_units.hxx (c_cgs, hc_MeV_cm) and
// mcnux_tetrad.hxx (detail::pi): zero new physical-constant literals
// (single-literal sweep, [MCNX-CNV-03]).
//
// Scope boundary (deliberate):
//   * NO range policy: inputs are assumed in-range; the [MCNX-OPA-06]
//     transparency floor / clamping / clamp counters are the T11 layer
//     (mcnux_opacity.hxx AxisBounds is its hook surface, not this header).
//   * NO NES/Pair/Brem assembly (structured open question,
//     opacity-eos-evaluation.md:340-349): the baseline is EmAb + Iso only,
//     and tabulated nu_x emission is identically zero.
//   * NO table residency/loading layer: the view bundle below is a dumb,
//     trivially-copyable aggregate of pointers and extents, filled by the
//     caller (selftest fixtures today; a residency layer later).
//   * NO caching: the nu_x kappa_s calls iso_evaluate twice per query
//     (opacity-eos-evaluation.md:325-329 allows "none").
//
// Deliberately free of cctk.h and CarpetX includes: plain portable C++ (the
// AMReX GPU qualifiers arrive via mcnux_opacity.hxx exactly as they do for
// the wrapper layer itself), so host code, device kernels, and the self-test
// battery can all include it. Its static_asserts run on every build via
// inclusion from stub.cxx.

#include "mcnux_coefficients.hxx"
#include "mcnux_opacity.hxx"
#include "mcnux_tetrad.hxx" // detail::pi
#include "mcnux_units.hxx"

#include <cmath>
#include <type_traits>

namespace MCNuX {

// The three EOS chemical-potential sub-tables the Kirchhoff closure consumes
// (the `Electron/Proton/Neutron Chemical Potential` datasets, MeV, in the
// table's own rest-mass convention — opacity-eos-evaluation.md:370-373).
enum class ChemicalPotential : int { Mu_e = 0, Mu_p = 1, Mu_n = 2 };
inline constexpr int num_chemical_potentials = 3;

constexpr int chemical_potential_index(ChemicalPotential m) noexcept {
  return static_cast<int>(m);
}

// The Iso moment index consumed by the baseline: the l = 0 Legendre
// coefficient (total elastic rate). An index pin, not a physical constant.
inline constexpr int iso_baseline_moment = 0;

// ---------------------------------------------------------------------------
// The table-view bundle: every pointer/extent the assembly needs, and
// nothing else. Pointers are caller-owned inputs (the T9 synthetic-fixture
// pattern); no residency, no I/O, no ownership here. Trivially copyable so
// the whole bundle can be captured by value into device kernels.
// ---------------------------------------------------------------------------
struct BaselineTableViews {
  // EmAb (4D per species): shared already-log10'd axes (log10 of MeV,
  // g/cm^3, KELVIN; Ye raw), one dataset pointer + scalar offset per
  // electron-type slot (indexed by emab_species_slot).
  wli::Real const *emab_LogEs;
  int emab_nE;
  wli::Real const *emab_LogDs;
  int emab_nD;
  wli::Real const *emab_LogTs;
  int emab_nT;
  wli::Real const *emab_Ys;
  int emab_nY;
  wli::Real emab_OS[num_tabulated_species];
  wli::Real const *emab_table[num_tabulated_species];

  // Iso (5D per species): shared axes as EmAb, one dataset pointer per
  // electron-type slot (indexed by iso_species_slot), the species-major
  // [nOpacities, nMoments] offsets array (selected only through
  // iso_offset -> wli::IsoOffset, never hand-indexed), and the moment count.
  wli::Real const *iso_LogEs;
  int iso_nE;
  wli::Real const *iso_LogDs;
  int iso_nD;
  wli::Real const *iso_LogTs;
  int iso_nT;
  wli::Real const *iso_Ys;
  int iso_nY;
  int iso_nMom;
  wli::Real const *iso_offsets; // [num_tabulated_species * iso_nMom]
  wli::Real const *iso_table[num_tabulated_species];

  // EOS chemical potentials (3D): shared RAW axes (Ds [g/cm^3], Ts [KELVIN —
  // table native], Ys), one (offset, dataset) pair per mu (indexed by
  // chemical_potential_index).
  wli::Real const *eos_Ds;
  int eos_nD;
  wli::Real const *eos_Ts_K;
  int eos_nT;
  wli::Real const *eos_Ys;
  int eos_nY;
  wli::Real eos_OS[num_chemical_potentials];
  wli::Real const *eos_table[num_chemical_potentials];
};

// ---------------------------------------------------------------------------
// The table-source callable ([MCNX-OPA-04/05]) — satisfies the frozen
// TableEval contract of evaluate_coefficients() in mcnux_coefficients.hxx.
// ---------------------------------------------------------------------------
struct BaselineTableCoefficients {
  BaselineTableViews v;

  // One EOS chemical potential [MeV] at the fluid state (raw coordinates in,
  // the wrapper owns the single MeV->K conversion).
  double chemical_potential(ChemicalPotential m,
                            const FluidState &st) const noexcept {
    const int i = chemical_potential_index(m);
    return eos_evaluate(st.rho_cgs, st.T_MeV, st.Ye, v.eos_Ds, v.eos_nD,
                        v.eos_Ts_K, v.eos_nT, v.eos_Ys, v.eos_nY, v.eos_OS[i],
                        v.eos_table[i]);
  }

  // mu_nu(s) = l(s) * (mu_e + mu_p - mu_n): +1 for nu_e, -1 for nu_e_bar
  // (opacity-eos-evaluation.md:207), via the pinned lepton numbers of
  // mcnux_units.hxx (no re-enumeration).
  double mu_nu(Species s, const FluidState &st) const noexcept {
    return lepton_number(s) * (chemical_potential(ChemicalPotential::Mu_e, st) +
                               chemical_potential(ChemicalPotential::Mu_p, st) -
                               chemical_potential(ChemicalPotential::Mu_n, st));
  }

  // kappa_s for one tabulated electron-type slot, iMom = 0, offset through
  // the sanctioned iso_offset lookup.
  double kappa_s_slot(int slot, double E_MeV,
                      const FluidState &st) const noexcept {
    const wli::Real OS = iso_offset(v.iso_offsets, num_tabulated_species,
                                    v.iso_nMom, slot, iso_baseline_moment);
    return iso_evaluate(E_MeV, st.rho_cgs, st.T_MeV, st.Ye, v.iso_LogEs,
                        v.iso_nE, v.iso_LogDs, v.iso_nD, v.iso_LogTs, v.iso_nT,
                        v.iso_Ys, v.iso_nY, iso_baseline_moment, v.iso_nMom, OS,
                        v.iso_table[slot]);
  }

  // The frozen TableEval contract: (Species, E [MeV], FluidState) ->
  // Coefficients {kappa_a [cm^-1], kappa_s [cm^-1],
  // eta [MeV cm^-3 s^-1 MeV^-1]}, g = 1 per species.
  Coefficients operator()(Species s, double E_MeV,
                          const FluidState &st) const noexcept {
    // nu_x FIRST — it has no EmAb/Iso dataset (emab_species_slot /
    // iso_species_slot return no_table_slot = -1, not a valid index):
    // kappa_a = eta = 0 exactly, kappa_s the exact half-mean of the two
    // electron-type slots ([MCNX-OPA-05]).
    if (s == Species::NuX) {
      const double ks_nue =
          kappa_s_slot(iso_species_slot(Species::NuE), E_MeV, st);
      const double ks_nueb =
          kappa_s_slot(iso_species_slot(Species::NuEBar), E_MeV, st);
      return {0.0, 0.5 * (ks_nue + ks_nueb), 0.0};
    }

    // Electron-type assembly ([MCNX-OPA-04]).
    const int ea = emab_species_slot(s);
    const double kappa_a = emab_evaluate(
        E_MeV, st.rho_cgs, st.T_MeV, st.Ye, v.emab_LogEs, v.emab_nE,
        v.emab_LogDs, v.emab_nD, v.emab_LogTs, v.emab_nT, v.emab_Ys, v.emab_nY,
        v.emab_OS[ea], v.emab_table[ea]);
    const double kappa_s = kappa_s_slot(iso_species_slot(s), E_MeV, st);

    // Kirchhoff closure: eta = c kappa_a f_eq 4 pi E^3/(hc)^3, with the
    // equilibrium occupancy at mu_nu(s); E, mu_nu, T all in MeV (the Kelvin
    // conversion happens only inside the wrappers, never here).
    const double f_eq =
        1.0 / (std::exp((E_MeV - mu_nu(s, st)) / st.T_MeV) + 1.0);
    const double eta = c_cgs * kappa_a * f_eq * 4.0 * detail::pi *
                       (E_MeV * E_MeV * E_MeV) /
                       (hc_MeV_cm * hc_MeV_cm * hc_MeV_cm);
    return {kappa_a, kappa_s, eta};
  }
};

// ---------------------------------------------------------------------------
// Compile-time verification — structure only (the WLI kernels are not
// constexpr; the value-level identities are runtime rows of the
// unit-selftest battery, mcnux_selftest.cxx).
// ---------------------------------------------------------------------------

// The functor satisfies the frozen TableEval contract of
// evaluate_coefficients() exactly.
static_assert(
    std::is_invocable_r<Coefficients, const BaselineTableCoefficients &,
                        Species, double, const FluidState &>::value,
    "[MCNX-OPA-04] BaselineTableCoefficients must satisfy the TableEval "
    "contract (Species, double E_MeV, const FluidState&) -> Coefficients");

// Dumb, device-capturable view bundle and functor (the AxisBounds /
// AnalyticOpacityParams precedent).
static_assert(std::is_trivially_copyable<BaselineTableViews>::value,
              "BaselineTableViews must be trivially copyable (captured by "
              "value into device kernels)");
static_assert(std::is_trivially_copyable<BaselineTableCoefficients>::value,
              "BaselineTableCoefficients must be trivially copyable");

// The nu_x special case is mandatory, not an optimization: the slot maps
// yield no valid table index for it.
static_assert(emab_species_slot(Species::NuX) == no_table_slot &&
                  iso_species_slot(Species::NuX) == no_table_slot,
              "[MCNX-OPA-05] nu_x has no EmAb/Iso dataset slot and must be "
              "special-cased before any table call");

// The mu_nu sign map rides on the pinned lepton numbers: +1 (nu_e),
// -1 (nu_e_bar), 0 (nu_x) — so mu_nu(nu_e_bar) = -mu_nu(nu_e) exactly.
static_assert(lepton_number(Species::NuE) == 1 &&
                  lepton_number(Species::NuEBar) == -1 &&
                  lepton_number(Species::NuX) == 0,
              "[MCNX-OPA-04] mu_nu(s) = l(s) (mu_e + mu_p - mu_n) requires "
              "the pinned lepton numbers (+1, -1, 0)");

// Three chemical-potential slots, distinct and zero-based.
static_assert(chemical_potential_index(ChemicalPotential::Mu_e) == 0 &&
                  chemical_potential_index(ChemicalPotential::Mu_p) == 1 &&
                  chemical_potential_index(ChemicalPotential::Mu_n) == 2 &&
                  num_chemical_potentials == 3,
              "exactly three EOS chemical-potential sub-tables, indexed "
              "mu_e = 0, mu_p = 1, mu_n = 2");

// The baseline consumes only the l = 0 Iso moment.
static_assert(iso_baseline_moment == 0,
              "[MCNX-OPA-04] the baseline elastic rate is the l = 0 Legendre "
              "coefficient (iMom = 0)");

} // namespace MCNuX

#endif // MCNUX_TABLE_COEFFS_HXX
