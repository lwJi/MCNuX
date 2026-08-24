#include "mcnux_coefficients.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

namespace MCNuX {

// Parameter glue for the source-agnostic coefficient interface of
// mcnux_coefficients.hxx (specs/verification-suite-design.md [MCNX-VER-05]):
// the only place the MCNuX::opacity_source keyword and the per-species
// analytic arrays kappa_a0[], kappa_s0[], eta_scale[] of param.ccl are read.
// Both are STEERABLE=NEVER, so consumers capture the returned values once at
// startup and carry them (AnalyticOpacityParams is trivially copyable) into
// their kernels.

CoefficientSource selected_coefficient_source() {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_EQUALS(opacity_source, "analytic"))
    return CoefficientSource::Analytic;
  return CoefficientSource::Table;
}

AnalyticOpacityParams analytic_params_from_parameters() {
  DECLARE_CCTK_PARAMETERS;
  AnalyticOpacityParams p{};
  for (int i = 0; i < NUM_SPECIES; ++i) {
    p.kappa_a0[i] = kappa_a0[i];
    p.kappa_s0[i] = kappa_s0[i];
    p.eta_scale[i] = eta_scale[i];
  }
  return p;
}

} // namespace MCNuX
