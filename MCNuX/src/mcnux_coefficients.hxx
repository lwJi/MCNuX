#ifndef MCNUX_COEFFICIENTS_HXX
#define MCNUX_COEFFICIENTS_HXX

// The source-agnostic per-species transport-coefficient interface of
// specs/opacity-eos-evaluation.md:105-126 and its analytic (gray) source,
// specs/verification-suite-design.md [MCNX-VER-05] (formulas normative at
// :213-217, units and constants at :219-220).
//
// Everything downstream of this interface (sampling, events, tallies, RNG)
// consumes one triple per (species, energy, fluid state),
//
//   Coefficients { kappa_a, kappa_s, eta }
//
// in microphysics units — kappa in cm^-1, eta in MeV cm^-3 s^-1 MeV^-1 —
// and never learns which source produced it ([MCNX-VER-05]: "the mode switch
// is invisible below the interface", verification-suite-design.md:322-324).
// Conversion to geometrized units is the consumer's job
// (specs/neutrino-matter-interactions.md:96-104), not this header's.
//
// Two sources feed the interface, selected by the runtime parameter
// MCNuX::opacity_source (param.ccl):
//
//   Table    — the production weaklib tables, assembled by the
//              [MCNX-OPA-04] layer (a later task); it plugs into the
//              `table` callable slot of evaluate_coefficients() below with
//              the contract (Species, double E_MeV, const FluidState&)
//              -> Coefficients.
//   Analytic — the gray Kirchhoff-form coefficients of this header, per
//              species s, E and T in MeV (verification-suite-design.md:213-217,
//              verbatim):
//
//                kappa_a(s, E) = kappa_a0(s)
//                kappa_s(s, E) = kappa_s0(s)
//                eta(s, E)     = eta_scale(s) * c * kappa_a0(s)
//                                * 4 pi E^3/(hc)^3 * 1/(exp(E/T) + 1)
//
//              with c = c_cgs and hc = hc_MeV_cm from mcnux_units.hxx (the
//              same pins as specs/opacity-eos-evaluation.md:356-361) and pi
//              from mcnux_tetrad.hxx (its single literal). eta_scale = 1 is
//              exact Kirchhoff equilibrium at mu_nu = 0 (J_eq = eta/kappa_a/c
//              = 4 pi E^3/(hc)^3 f_eq); eta_scale = 0 is an absorption-only
//              medium. The analytic mode is the nu_x coverage mechanism:
//              kappa_a0(NuX) > 0 with eta_scale(NuX) = 1 exercises nu_x
//              emission and absorption, which the tables cannot reach.
//
// Conventions pinned here:
//   * g = 1 per species in the coefficient ([MCNX-VER-05]); the nu_x
//     degeneracy g = 4 enters exactly once, downstream at packet creation,
//     never in a coefficient (specs/conventions-and-units.md [MCNX-CNV-06],
//     `species_degeneracy` in mcnux_units.hxx).
//   * T_MeV > 0 is assumed and NOT checked: E = 0 with T = 0 is 0/0 in
//     exp(E/T), and a silent guard would hide a caller defect (the alpha
//     policy of mcnux_trp.hxx). Range enforcement belongs to the fluid-state
//     range policy, not this pure map. A large E/T overflows std::exp to
//     +inf and gives f_eq = 0 exactly — the intended limit, not an error.
//   * Trapped-regime relabeling (mcnux_trp.hxx) is a pure post-map on the
//     Coefficients triple and applies to the analytic coefficients
//     identically (verification-suite-design.md:228-230).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// C++ so that host code, device kernels, and the self-test battery can all
// include it. The parameter glue (selected_coefficient_source(),
// analytic_params_from_parameters()) is declared here but defined in
// mcnux_coefficients.cxx, the only translation unit that touches
// cctk_Parameters.h for this interface.

#include "mcnux_tetrad.hxx" // detail::pi
#include "mcnux_units.hxx"

#include <cmath>
#include <type_traits>

namespace MCNuX {

// The source-agnostic triple, microphysics units (cm^-1, cm^-1,
// MeV cm^-3 s^-1 MeV^-1), one species, g = 1.
struct Coefficients {
  double kappa_a, kappa_s, eta;
};

// The fluid-state input of the interface (opacity-eos-evaluation.md:118-119):
// rest-mass density [g cm^-3], temperature [MeV], electron fraction.
struct FluidState {
  double rho_cgs, T_MeV, Ye;
};

// The analytic-mode parameter set, indexed by species_index(Species)
// (NuE = 0, NuEBar = 1, NuX = 2; mcnux_units.hxx, never re-enumerated here).
// kappa_a0, kappa_s0 in cm^-1; eta_scale dimensionless. Plain aggregate of
// doubles so it is trivially copyable into device kernels by value.
struct AnalyticOpacityParams {
  double kappa_a0[NUM_SPECIES];
  double kappa_s0[NUM_SPECIES];
  double eta_scale[NUM_SPECIES];
};

// Which source feeds the interface.
enum class CoefficientSource : int { Table = 0, Analytic = 1 };

// ---------------------------------------------------------------------------
// Analytic (gray) source  [MCNX-VER-05]
// ---------------------------------------------------------------------------

// The Kirchhoff-form ratio eta/kappa_a at eta_scale = 1 and mu_nu = 0,
//   c * 4 pi E^3/(hc)^3 * 1/(exp(E/T) + 1)     [MeV cm^-2 s^-1 MeV^-1],
// i.e. c times the equilibrium spectral energy density J_eq(E) at occupancy
// f_eq = 1/(exp(E/T) + 1). Not constexpr (std::exp).
inline double kirchhoff_eta_over_kappa_a(double E_MeV, double T_MeV) noexcept {
  const double f_eq = 1.0 / (std::exp(E_MeV / T_MeV) + 1.0);
  return c_cgs * 4.0 * detail::pi * (E_MeV * E_MeV * E_MeV) /
         (hc_MeV_cm * hc_MeV_cm * hc_MeV_cm) * f_eq;
}

// The two energy-independent legs of the analytic map, factored out so that
// their exact pass-through is checkable at compile time.
constexpr double analytic_kappa_a(const AnalyticOpacityParams &p,
                                  Species s) noexcept {
  return p.kappa_a0[species_index(s)];
}
constexpr double analytic_kappa_s(const AnalyticOpacityParams &p,
                                  Species s) noexcept {
  return p.kappa_s0[species_index(s)];
}

// The analytic coefficient triple, verification-suite-design.md:213-217
// verbatim.
inline Coefficients analytic_coefficients(const AnalyticOpacityParams &p,
                                          Species s, double E_MeV,
                                          double T_MeV) noexcept {
  const double ka = analytic_kappa_a(p, s);
  return {ka, analytic_kappa_s(p, s),
          p.eta_scale[species_index(s)] * ka *
              kirchhoff_eta_over_kappa_a(E_MeV, T_MeV)};
}

// ---------------------------------------------------------------------------
// Source-agnostic dispatch — the hook the table source plugs into
// ---------------------------------------------------------------------------

// Returns the identical Coefficients type whichever source is selected.
// `table` is the table-source callable with the contract
//   (Species, double E_MeV, const FluidState&) -> Coefficients
// (the [MCNX-OPA-04] assembly, a later task); it is only invoked under
// CoefficientSource::Table. The analytic leg reads only st.T_MeV.
template <class TableEval>
inline Coefficients evaluate_coefficients(CoefficientSource src,
                                          const AnalyticOpacityParams &ap,
                                          TableEval &&table, Species s,
                                          double E_MeV, const FluidState &st) {
  if (src == CoefficientSource::Analytic)
    return analytic_coefficients(ap, s, E_MeV, st.T_MeV);
  return table(s, E_MeV, st);
}

// ---------------------------------------------------------------------------
// Parameter glue (defined in mcnux_coefficients.cxx; reads
// MCNuX::opacity_source, kappa_a0[], kappa_s0[], eta_scale[] of param.ccl).
// ---------------------------------------------------------------------------

CoefficientSource selected_coefficient_source();
AnalyticOpacityParams analytic_params_from_parameters();

// ---------------------------------------------------------------------------
// Compile-time verification  [MCNX-VER-05] — only what is constexpr-checkable:
// the exact pass-through of kappa_a/kappa_s, the enumerator distinctness of
// the source switch, and the device-capturable layout of the parameter
// aggregate (the AxisBounds precedent of mcnux_opacity.hxx). The Kirchhoff
// leg (std::exp) is verified at runtime by the `unit-selftest` battery.
// ---------------------------------------------------------------------------

namespace detail {

// Fixture: exact binary fractions, distinct per species, kappa_a0(NuX) > 0.
inline constexpr AnalyticOpacityParams coeffix_params = {
    {0.5, 0.25, 0.125}, {0.0625, 0.03125, 0.015625}, {1.0, 0.5, 1.0}};

} // namespace detail

static_assert(analytic_kappa_a(detail::coeffix_params, Species::NuE) == 0.5 &&
                  analytic_kappa_a(detail::coeffix_params, Species::NuEBar) ==
                      0.25 &&
                  analytic_kappa_a(detail::coeffix_params, Species::NuX) ==
                      0.125,
              "[MCNX-VER-05] kappa_a(s, E) = kappa_a0(s) must be an exact "
              "per-species pass-through");
static_assert(analytic_kappa_s(detail::coeffix_params, Species::NuE) ==
                      0.0625 &&
                  analytic_kappa_s(detail::coeffix_params, Species::NuEBar) ==
                      0.03125 &&
                  analytic_kappa_s(detail::coeffix_params, Species::NuX) ==
                      0.015625,
              "[MCNX-VER-05] kappa_s(s, E) = kappa_s0(s) must be an exact "
              "per-species pass-through");

static_assert(CoefficientSource::Table != CoefficientSource::Analytic,
              "[MCNX-VER-05] the two coefficient sources must be distinct");

static_assert(sizeof(AnalyticOpacityParams) ==
                  3 * NUM_SPECIES * sizeof(double),
              "AnalyticOpacityParams must be exactly 3 x NUM_SPECIES doubles");
static_assert(std::is_trivially_copyable<AnalyticOpacityParams>::value,
              "AnalyticOpacityParams must be trivially copyable (captured by "
              "value into device kernels)");
static_assert(std::is_trivially_copyable<Coefficients>::value &&
                  std::is_trivially_copyable<FluidState>::value,
              "Coefficients and FluidState must be trivially copyable");

} // namespace MCNuX

#endif // MCNUX_COEFFICIENTS_HXX
