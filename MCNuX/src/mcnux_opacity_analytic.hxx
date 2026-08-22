#ifndef MCNUX_OPACITY_ANALYTIC_HXX
#define MCNUX_OPACITY_ANALYTIC_HXX

// The analytic (gray) opacity mode — a first-class physics capability of
// specs/verification-suite-design.md [MCNX-VER-05] (formulas normative),
// feeding the identical per-species coefficient interface of
// specs/opacity-eos-evaluation.md:107-126 (kappa_a, kappa_s, eta per single
// species, g = 1). Per species s, with E and T in MeV:
//
//   kappa_a(s, E) = kappa_a0(s)
//   kappa_s(s, E) = kappa_s0(s)
//   eta(s, E)     = eta_scale(s) * c kappa_a0(s) * 4 pi E^3/(hc)^3
//                   * 1/(exp(E/T) + 1)
//
// Units are those of the coefficient interface (microphysics system): kappa
// in cm^-1, eta in MeV cm^-3 s^-1 MeV^-1. No code-unit conversion happens
// here — opacity_code_to_cgs / emissivity_code_to_cgs of mcnux_units.hxx are
// the downstream consumer's concern. The constants c and hc are
// MCNuX::c_cgs and MCNuX::hc_MeV_cm (the single-literal site), 4 pi reuses
// detail::pi of mcnux_tetrad.hxx; no literal is restated here.
//
// eta_scale = 1 gives exact Kirchhoff equilibrium at mu_nu = 0: the
// equilibrium spectral energy density J_eq(E) = eta/kappa_a * 1/c
// = 4 pi E^3/(hc)^3 * f_eq(E) (the analytic target of the equilibration
// benchmarks); eta_scale = 0 gives an absorption-only medium (beam
// benchmarks). The mode is the nu_x coverage mechanism: kappa_a0(nu_x) > 0
// with eta_scale(nu_x) = 1 exercises nu_x emission and absorption, which the
// weaklib tables cannot express. No degeneracy factor enters anywhere here
// (g = 4 enters exactly once, at packet creation — conventions-and-units.md,
// [MCNX-CNV-06]).
//
// Runtime selection between this mode and the production tables lives in
// mcnux_coefficients.hxx (the source-agnostic dispatch); this header is the
// pure formula.
//
// Deliberately free of cctk.h, AMReX, CarpetX, and WeakLibInterp includes:
// plain portable constexpr C++ (the mcnux_units.hxx / mcnux_trp.hxx
// pattern) so host code, device code, and the self-test battery can all
// include it. The formula core is templated on the exponential, exactly as
// build_tetrad_with(..., sqrt_fn) in mcnux_tetrad.hxx: production passes
// std::exp (host entry point below, or the device functor of
// mcnux_coefficients.hxx), the compile-time self-tests pass a constexpr
// stand-in — std::exp is never called in a constant expression (clang
// rejects it; GCC's builtin extension would only mask the portability
// defect).

#include "mcnux_tetrad.hxx" // detail::pi
#include "mcnux_units.hxx"

#include <cmath>
#include <type_traits>

namespace MCNuX {

// THE coefficient interface triple (specs/opacity-eos-evaluation.md:107-116):
// kappa_a [cm^-1], kappa_s [cm^-1], eta [MeV cm^-3 s^-1 MeV^-1], per single
// species, g = 1. Source-agnostic: the table leg (T10) and the analytic leg
// both produce exactly this.
struct Coefficients {
  double kappa_a, kappa_s, eta;
};

// The per-species analytic-mode parameters of [MCNX-VER-05] (the parameter
// table of specs/verification-suite-design.md:112-118), indexed by
// species_index(Species) — the binding nu_e / nu_e_bar / nu_x order of
// mcnux_units.hxx, never a separate ordering. kappa_a0, kappa_s0 in cm^-1
// (>= 0); eta_scale dimensionless (>= 0). Range enforcement is param.ccl's.
struct AnalyticOpacityParams {
  double kappa_a0[NUM_SPECIES];
  double kappa_s0[NUM_SPECIES];
  double eta_scale[NUM_SPECIES];
};

static_assert(std::is_trivially_copyable_v<AnalyticOpacityParams>,
              "AnalyticOpacityParams must be trivially copyable "
              "(captured by value into device kernels)");
static_assert(std::is_trivially_copyable_v<Coefficients>,
              "Coefficients must be trivially copyable (device-returnable)");

// Planck phase-space factor 4 pi E^3/(hc)^3 [MeV^3 / (MeV cm)^3 = cm^-3],
// E in MeV. The spectral energy density of a unit-occupancy single-species
// gas; multiplied by an occupancy f it is J(E), and c kappa_a J(E) is the
// Kirchhoff emissivity.
constexpr double planck_phase_space_factor(double E_MeV) noexcept {
  return 4.0 * detail::pi * E_MeV * E_MeV * E_MeV /
         (hc_MeV_cm * hc_MeV_cm * hc_MeV_cm);
}

// The formula core of [MCNX-VER-05], templated on the exponential.
// The occupancy 1/(exp(E/T) + 1) is formed FIRST so that an overflowing
// exp(E/T) (E >> T) yields an occupancy of exactly 0 and hence eta = 0 —
// never the 0 * inf = NaN a product-ordered evaluation would produce.
template <class ExpF>
constexpr Coefficients
analytic_coefficients_with(const AnalyticOpacityParams &p, Species s,
                           double E_MeV, double T_MeV, ExpF exp_fn) noexcept {
  const int i = species_index(s);
  const double occupancy = 1.0 / (exp_fn(E_MeV / T_MeV) + 1.0);
  return Coefficients{p.kappa_a0[i], p.kappa_s0[i],
                      p.eta_scale[i] * c_cgs * p.kappa_a0[i] *
                          planck_phase_space_factor(E_MeV) * occupancy};
}

// Production host entry point: the core with std::exp (the build_tetrad
// precedent). Pure function of its arguments.
inline Coefficients analytic_coefficients(const AnalyticOpacityParams &p,
                                          Species s, double E_MeV,
                                          double T_MeV) noexcept {
  return analytic_coefficients_with(
      p, s, E_MeV, T_MeV, [](double x) { return std::exp(x); });
}

// ---------------------------------------------------------------------------
// Compile-time verification — [MCNX-VER-05] on exact binary-fraction
// fixtures, with a constexpr stand-in for the exponential (returns 1.0, so
// the occupancy is exactly 0.5 and no constexpr exp is needed).
// ---------------------------------------------------------------------------

namespace detail {

// exp stand-in for constant evaluation only: exp_fn(x) = 1 for every x,
// hence occupancy = 1/(1 + 1) = 0.5 exactly.
struct HalfOccupancyExp {
  constexpr double operator()(double) const noexcept { return 1.0; }
};

// Fixture: per-species opacities that are exact binary fractions, all
// distinct so a species-index slip is visible; eta_scale = 1 (Kirchhoff).
inline constexpr AnalyticOpacityParams anafix_kirchhoff = {
    {0.5, 0.25, 0.125}, {2.0, 1.0, 0.75}, {1.0, 1.0, 1.0}};

// Same opacities, eta_scale = 0 (absorption-only medium).
inline constexpr AnalyticOpacityParams anafix_absorb_only = {
    {0.5, 0.25, 0.125}, {2.0, 1.0, 0.75}, {0.0, 0.0, 0.0}};

inline constexpr double anafix_E = 4.0; // MeV
inline constexpr double anafix_T = 2.0; // MeV

inline constexpr Coefficients anafix_nue = analytic_coefficients_with(
    anafix_kirchhoff, Species::NuE, anafix_E, anafix_T, HalfOccupancyExp{});
inline constexpr Coefficients anafix_nuebar = analytic_coefficients_with(
    anafix_kirchhoff, Species::NuEBar, anafix_E, anafix_T,
    HalfOccupancyExp{});
inline constexpr Coefficients anafix_nux = analytic_coefficients_with(
    anafix_kirchhoff, Species::NuX, anafix_E, anafix_T, HalfOccupancyExp{});
inline constexpr Coefficients anafix_nux_absorb_only =
    analytic_coefficients_with(anafix_absorb_only, Species::NuX, anafix_E,
                               anafix_T, HalfOccupancyExp{});

} // namespace detail

// kappa_a(s, E) = kappa_a0(s) and kappa_s(s, E) = kappa_s0(s): bitwise
// pass-through, per species, in the binding index order.
static_assert(detail::anafix_nue.kappa_a == 0.5 &&
                  detail::anafix_nuebar.kappa_a == 0.25 &&
                  detail::anafix_nux.kappa_a == 0.125,
              "[MCNX-VER-05] kappa_a must be the per-species kappa_a0, exactly");
static_assert(detail::anafix_nue.kappa_s == 2.0 &&
                  detail::anafix_nuebar.kappa_s == 1.0 &&
                  detail::anafix_nux.kappa_s == 0.75,
              "[MCNX-VER-05] kappa_s must be the per-species kappa_s0, exactly");

// eta_scale = 0 gives eta = 0 exactly (absorption-only medium), with the
// opacities untouched.
static_assert(detail::anafix_nux_absorb_only.eta == 0.0 &&
                  detail::anafix_nux_absorb_only.kappa_a == 0.125 &&
                  detail::anafix_nux_absorb_only.kappa_s == 0.75,
              "[MCNX-VER-05] eta_scale = 0 must give eta = 0 exactly and "
              "leave the opacities unchanged");

// Kirchhoff / J_eq identity at eta_scale = 1: eta/(c kappa_a) =
// 4 pi E^3/(hc)^3 * f_eq, here with f_eq = 0.5 by the stand-in. Machine tier.
static_assert(detail::approx_eq(detail::anafix_nue.eta /
                                    (c_cgs * detail::anafix_nue.kappa_a),
                                planck_phase_space_factor(detail::anafix_E) *
                                    0.5,
                                detail::rtol_machine) &&
                  detail::approx_eq(
                      detail::anafix_nux.eta / (c_cgs * detail::anafix_nux.kappa_a),
                      planck_phase_space_factor(detail::anafix_E) * 0.5,
                      detail::rtol_machine),
              "[MCNX-VER-05] eta/(c kappa_a) must equal 4 pi E^3/(hc)^3 f_eq "
              "(the J_eq Kirchhoff identity)");

// The phase-space factor scales as E^3 (ratio 8 at 2E vs E), machine tier.
static_assert(detail::approx_eq(planck_phase_space_factor(2.0 * detail::anafix_E) /
                                    planck_phase_space_factor(detail::anafix_E),
                                8.0, detail::rtol_machine),
              "4 pi E^3/(hc)^3 must scale as E^3");

// nu_x coverage mechanism: kappa_a0(nu_x) > 0 with eta_scale(nu_x) = 1 yields
// a strictly positive nu_x emissivity — what the tables can never produce.
static_assert(detail::anafix_nux.eta > 0.0,
              "[MCNX-VER-05] nu_x with kappa_a0 > 0 and eta_scale = 1 must "
              "emit (eta > 0)");

// The factor itself is positive and finite at the fixture (guards an hc
// mis-derivation that would send it to 0 or inf).
static_assert(planck_phase_space_factor(detail::anafix_E) > 0.0 &&
                  planck_phase_space_factor(detail::anafix_E) < 1.0e300,
              "4 pi E^3/(hc)^3 must be positive and finite");

} // namespace MCNuX

#endif // MCNUX_OPACITY_ANALYTIC_HXX
