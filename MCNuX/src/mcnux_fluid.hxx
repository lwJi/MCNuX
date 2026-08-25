#ifndef MCNUX_FLUID_HXX
#define MCNUX_FLUID_HXX

// Fluid-state / four-velocity gather at a packet's cell.
//
// Governing spec: specs/neutrino-matter-interactions.md.
//   * [MCNX-INT-05] coefficient provenance — the opacities of [MCNX-INT-01]
//     are evaluated "at the cell-centered fluid state of the packet's current
//     cell — never re-derived, re-interpolated sub-cell, or modified here".
//     The gather below is therefore a DIRECT ONE-CELL READ of the
//     cell-centered (ccc) HydroBaseX fields: containing-cell anchor, then the
//     stored value, bitwise — no trilinear (or any) sub-cell interpolation.
//     (Contrast VertexMetricGather in mcnux_geodesic.hxx, whose fields are
//     vertex-centered and whose spec mandates interpolation.)
//   * The inputs table of the same spec pins u^mu as the CELL-CENTERED fluid
//     four-velocity with g_munu u^mu u^nu = -1; the metric entering the lift
//     is therefore evaluated at the containing cell's CENTER coordinate, not
//     at the packet position.
//
// Units ([MCNX-HYD-04], specs/hydro-coupling-source-terms.md:165-175, and the
// TestMCNuX writer contract): rho in g/cm^3, temperature in MeV, Ye
// dimensionless — passed through VERBATIM into the coefficient-interface
// FluidState (mcnux_coefficients.hxx). No unit conversion happens here; the
// MeV -> Kelvin conversion of mcnux_units.hxx belongs to the WeakLibInterp
// call boundary only.
//
// v^i -> u^mu lift (the Valencia 3-velocity of HydroBaseX::vel,
// CarpetX/HydroBaseX/interface.ccl "velocity v^i"):
//
//   W   = 1 / sqrt(1 - gamma_ij v^i v^j)     (Lorentz factor)
//   u^t = W / alpha
//   u^i = W (v^i - beta^i / alpha)
//
// NO SPEC PINS THIS FORMULA — it is the standard Valencia-formulation lift
// (v^i = (u^i/u^t + beta^i)/alpha, W = alpha u^t), recorded here as this
// task's one unpinned physics decision. Its consumers are the tetrad's
// timelike leg ([MCNX-PKT-04], build_tetrad in mcnux_tetrad.hxx) and the
// fluid-frame energy nu = -p_mu u^mu ([MCNX-INT-01], fluid_frame_energy in
// mcnux_interactions.hxx), both of which take u^mu as a plain double[4] with
// the normalization contract g_munu u^mu u^nu = -1 (verified in the selftest
// battery via spacetime_metric_from_adm/metric_dot). gamma_ij v^i v^j is a
// direct contraction against the six stored components — never a re-derived
// 3x3 inverse (the single-inverse rule of mcnux_tetrad.hxx).
//
// The only AMReX dependence is the amrex::Array4 views and the host/device
// annotations of the gather functor (the VertexMetricGather precedent).

#include "mcnux_coefficients.hxx" // FluidState
#include "mcnux_tetrad.hxx"       // SpatialMetric, detail::csqrt

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// One-cell fluid sample  [MCNX-INT-05]
// ---------------------------------------------------------------------------

// The gather's output: the coefficient-interface FluidState (units verbatim)
// plus the Valencia 3-velocity v^i of the same cell. u^mu is deliberately NOT
// a member: the FluidState consumers (mcnux_table_coeffs.hxx,
// mcnux_table_range.hxx) never see a velocity, and the u^mu consumers take a
// plain double[4] (see the lift below).
struct FluidSample {
  FluidState state; // rho [g/cm^3], T [MeV], Ye — verbatim
  double vel[3];    // Valencia v^i, dimensionless
};

// Cell-centered one-cell gather functor over the four HydroBaseX groups
// (component order per CarpetX/HydroBaseX/interface.ccl: vel = velx vely
// velz; rho, temperature, Ye single-component). Trivially copyable by-value
// capture into packet kernels ([MCNX-GPU-02]). The containing-cell anchor is
// the identical arithmetic as VertexMetricGather (mcnux_geodesic.hxx): for
// a ccc Array4 the global index i covers the cell
// [prob_lo + i dx, prob_lo + (i+1) dx), value stored at its center, so
// i0 = floor((x - prob_lo)/dx) IS the containing cell — one contains() guard,
// no i0+1 corner (nothing beyond the single cell is read).
struct CellFluidGather {
  amrex::Array4<const amrex::Real> rho;         // 1 component, g/cm^3
  amrex::Array4<const amrex::Real> vel;         // 3 components, v^i
  amrex::Array4<const amrex::Real> temperature; // 1 component, MeV
  amrex::Array4<const amrex::Real> ye;          // 1 component
  double prob_lo[3];
  double dx[3];

  // Containing-cell index of a position (the shared anchor arithmetic).
  AMREX_GPU_HOST_DEVICE void containing_cell(double x, double y, double z,
                                             int i0[3]) const noexcept {
    using std::floor;
    const double pos[3] = {x, y, z};
    for (int d = 0; d < 3; ++d)
      i0[d] = static_cast<int>(floor((pos[d] - prob_lo[d]) / dx[d]));
  }

  // Center coordinate of cell index i along direction d — the metric
  // evaluation point of the cell-centered u^mu (spec inputs table).
  AMREX_GPU_HOST_DEVICE double cell_center(int i, int d) const noexcept {
    return prob_lo[d] + (double(i) + 0.5) * dx[d];
  }

  AMREX_GPU_HOST_DEVICE FluidSample operator()(double x, double y,
                                               double z) const noexcept {
    int i0[3];
    containing_cell(x, y, z, i0);
    return at_cell(i0);
  }

  // The direct one-cell read, bitwise ([MCNX-INT-05]: no sub-cell
  // re-interpolation).
  AMREX_GPU_HOST_DEVICE FluidSample at_cell(const int i0[3]) const noexcept {
    AMREX_ASSERT(rho.contains(i0[0], i0[1], i0[2]) &&
                 vel.contains(i0[0], i0[1], i0[2]) &&
                 temperature.contains(i0[0], i0[1], i0[2]) &&
                 ye.contains(i0[0], i0[1], i0[2]));
    FluidSample s{};
    s.state.rho_cgs = rho(i0[0], i0[1], i0[2]);
    s.state.T_MeV = temperature(i0[0], i0[1], i0[2]);
    s.state.Ye = ye(i0[0], i0[1], i0[2]);
    for (int k = 0; k < 3; ++k)
      s.vel[k] = vel(i0[0], i0[1], i0[2], k);
    return s;
  }
};

static_assert(std::is_trivially_copyable_v<CellFluidGather> &&
                  std::is_trivially_copyable_v<FluidSample>,
              "[MCNX-GPU-02] the fluid gather functor and its sample are "
              "captured by value into packet kernels and must be trivially "
              "copyable");

// ---------------------------------------------------------------------------
// Valencia lift v^i -> u^mu (formula documented in the header comment)
// ---------------------------------------------------------------------------

// Core, templated on the sqrt so the compile-time self-tests can run it with
// detail::csqrt (the build_tetrad_with pattern); production callers use
// valencia_four_velocity below. gamma_ij v^i v^j by direct contraction of the
// six stored components — no 3x3 inverse is derived here.
template <class SqrtF>
constexpr void valencia_four_velocity_with(double alpha,
                                           const double beta_up[3],
                                           const SpatialMetric &gamma,
                                           const double v_up[3], double u[4],
                                           SqrtF sqrt_fn) noexcept {
  const double v2 =
      gamma.xx * v_up[0] * v_up[0] + gamma.yy * v_up[1] * v_up[1] +
      gamma.zz * v_up[2] * v_up[2] +
      2.0 * (gamma.xy * v_up[0] * v_up[1] + gamma.xz * v_up[0] * v_up[2] +
             gamma.yz * v_up[1] * v_up[2]);
  const double W = 1.0 / sqrt_fn(1.0 - v2);
  u[0] = W / alpha;
  for (int i = 0; i < 3; ++i)
    u[i + 1] = W * (v_up[i] - beta_up[i] / alpha);
}

AMREX_GPU_HOST_DEVICE inline void
valencia_four_velocity(double alpha, const double beta_up[3],
                       const SpatialMetric &gamma, const double v_up[3],
                       double u[4]) noexcept {
  valencia_four_velocity_with(alpha, beta_up, gamma, v_up, u,
                              [](double x) { return std::sqrt(x); });
}

// ---------------------------------------------------------------------------
// Compile-time verification (static limit; the measured rows live in the
// runtime `unit-selftest` battery, mcnux_selftest_fluidgather.cxx)
// ---------------------------------------------------------------------------

namespace detail {

// v = 0 gives W = 1 exactly (csqrt(1.0) == 1.0), so on any (alpha, beta) the
// lift is u = (1/alpha, -beta^i/alpha) exactly; on flat static data that is
// the normalized static observer u = (1, 0, 0, 0).
constexpr bool static_valencia_lift_holds() noexcept {
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double v0[3] = {0.0, 0.0, 0.0};
  double u[4] = {0.0, 0.0, 0.0, 0.0};
  valencia_four_velocity_with(1.0, beta0, SpatialMetric{1, 0, 0, 1, 0, 1}, v0,
                              u, [](double x) { return csqrt(x); });
  return u[0] == 1.0 && u[1] == 0.0 && u[2] == 0.0 && u[3] == 0.0;
}

} // namespace detail

static_assert(detail::static_valencia_lift_holds(),
              "the v = 0 Valencia lift on flat static data must give the "
              "static observer u = (1, 0, 0, 0) exactly");

} // namespace MCNuX

#endif // MCNUX_FLUID_HXX
