#ifndef MCNUX_COEFFICIENTS_HXX
#define MCNUX_COEFFICIENTS_HXX

// Source-agnostic transport-coefficient dispatch — the runtime-parameter
// selection of specs/opacity-eos-evaluation.md:118-126 ("the interface is
// source-agnostic: a runtime parameter selects the production weaklib tables
// or the analytic per-species formulas") and the parameter-selected opacity
// source of specs/verification-suite-design.md [MCNX-VER-05].
//
// Everything downstream of coefficients() (sampling, events, tallies, RNG,
// the trapped-regime relabeling of mcnux_trp.hxx) sees only the
// Coefficients triple and is identical in both modes. The selection is made
// ONCE, from MCNuX::opacity_source (param.ccl), into a trivially copyable
// CoefficientSource that consumers capture by value into device kernels;
// the per-call dispatch is a branch on an int.
//
// Table leg status (LOUD, for T10): the Table source is a PLACEHOLDER that
// returns Coefficients{0, 0, 0} — a transparent, non-emitting medium. This is
// acceptable only because no runtime consumer of coefficients() exists yet
// (emission sampling and event sampling are unbuilt); T10 replaces the
// placeholder with the [MCNX-OPA-04/05] table assembly (EmAb/Iso lookups in
// Kelvin/log10 through mcnux_opacity.hxx, Kirchhoff closure with
// planck_phase_space_factor, nu_x mapping) by adding its table views as
// further members of CoefficientSource WITHOUT changing the coefficients()
// call below. Until then, rho and Ye are unused.
//
// AMReX-qualified (unlike the pure headers it includes) because the table leg
// is not constexpr (WeakLibInterp kernels, mcnux_opacity.hxx) and because the
// entry point must be callable from device kernels in either mode; the
// analytic leg's constexpr template is device-callable under AMReX's
// --expt-relaxed-constexpr (see mcnux_rng.hxx), through the device-capable
// exp functor below. Never includes cctk_Parameters.h: the param.ccl ->
// CoefficientSource plumbing is make_coefficient_source(), called from a
// scheduled routine's DECLARE_CCTK_PARAMETERS scope (the mcnux_paramcheck.cxx
// idiom), passing the KEYWORD string and the three CCTK_REAL arrays.

#include "mcnux_opacity.hxx"
#include "mcnux_opacity_analytic.hxx"
#include "mcnux_units.hxx"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>
#include <cstring>
#include <type_traits>

namespace MCNuX {

// The runtime-selected source. Values are pinned (a CoefficientSource is
// captured by value; the int is what crosses into device code).
enum class OpacitySource : int {
  Table = 0,   // production weaklib tables (specs/opacity-eos-evaluation.md)
  Analytic = 1 // per-species gray Kirchhoff-form coefficients [MCNX-VER-05]
};

static_assert(static_cast<int>(OpacitySource::Table) == 0 &&
                  static_cast<int>(OpacitySource::Analytic) == 1,
              "OpacitySource enum values are pinned (0 = table, 1 = analytic)");

// Everything coefficients() needs, built once on the host from the
// parameters and captured by value. T10 appends its table views here.
struct CoefficientSource {
  OpacitySource source;
  AnalyticOpacityParams analytic;
};

static_assert(std::is_trivially_copyable_v<CoefficientSource>,
              "CoefficientSource must be trivially copyable "
              "(captured by value into device kernels)");

namespace detail {

// Device-capable exponential for the analytic leg (std::exp is supported in
// device code by nvcc/hip; the host build is the plain <cmath> call).
struct DeviceExp {
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
  operator()(double x) const noexcept {
    return std::exp(x);
  }
};

} // namespace detail

// THE consumer entry point: the per-species coefficient triple at a fluid
// state (rho [g/cm^3], T [MeV], Ye) and fluid-frame energy E [MeV], for
// species s — specs/opacity-eos-evaluation.md:107-126. Microphysics units
// (cm^-1, cm^-1, MeV cm^-3 s^-1 MeV^-1); g = 1 per single species.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Coefficients
coefficients(const CoefficientSource &src, Species s, double rho_cgs,
             double T_MeV, double Ye, double E_MeV) noexcept {
  switch (src.source) {
  case OpacitySource::Analytic:
    return analytic_coefficients_with(src.analytic, s, E_MeV, T_MeV,
                                      detail::DeviceExp{});
  case OpacitySource::Table:
  default:
    // T10 PLACEHOLDER — see the header comment. Transparent, non-emitting.
    (void)rho_cgs;
    (void)Ye;
    return Coefficients{0.0, 0.0, 0.0};
  }
}

// ---------------------------------------------------------------------------
// Host-only plumbing from MCNuX/param.ccl
// ---------------------------------------------------------------------------

namespace detail {

// Constexpr string equality (so the keyword map can be static_asserted).
constexpr bool cstr_eq(const char *a, const char *b) noexcept {
  while (*a != '\0' && *a == *b) {
    ++a;
    ++b;
  }
  return *a == *b;
}

} // namespace detail

// The KEYWORD values of MCNuX::opacity_source (param.ccl). The flesh
// validates a KEYWORD parameter against its declared value set before any
// scheduled routine runs, so no other string can reach this function from a
// parameter; the fall-through to Table exists only to keep the function total.
constexpr OpacitySource
opacity_source_from_keyword(const char *keyword) noexcept {
  if (detail::cstr_eq(keyword, "analytic"))
    return OpacitySource::Analytic;
  return OpacitySource::Table;
}

static_assert(opacity_source_from_keyword("table") == OpacitySource::Table,
              "param.ccl keyword \"table\" must select the table source");
static_assert(opacity_source_from_keyword("analytic") ==
                  OpacitySource::Analytic,
              "param.ccl keyword \"analytic\" must select the analytic source");

// Build the captured source from the parameters as they arrive after
// DECLARE_CCTK_PARAMETERS: the KEYWORD as `const char *`, each CCTK_REAL[3]
// array as `const CCTK_REAL *` (NUM_SPECIES entries, species-index order).
inline CoefficientSource
make_coefficient_source(const char *keyword, const double *kappa_a0,
                        const double *kappa_s0,
                        const double *eta_scale) noexcept {
  CoefficientSource src{};
  src.source = opacity_source_from_keyword(keyword);
  for (int i = 0; i < NUM_SPECIES; ++i) {
    src.analytic.kappa_a0[i] = kappa_a0[i];
    src.analytic.kappa_s0[i] = kappa_s0[i];
    src.analytic.eta_scale[i] = eta_scale[i];
  }
  return src;
}

} // namespace MCNuX

#endif // MCNUX_COEFFICIENTS_HXX
