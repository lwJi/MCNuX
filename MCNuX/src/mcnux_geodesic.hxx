#ifndef MCNUX_GEODESIC_HXX
#define MCNUX_GEODESIC_HXX

// Geodesic push of the packet population: the metric gather, the null-geodesic
// right-hand side, and the time integrator.
//
// Governing spec: specs/geodesic-propagation.md.
//   * [MCNX-GEO-01] equations of motion (restated there, normative):
//         dx^i/dt = gamma^{ij} p_j/p^t - beta^i
//         dp_i/dt = -alpha (d_i alpha) p^t + (d_i beta^k) p_k
//                   - 1/2 (d_i gamma^{jk}) p_j p_k/p^t
//     with d_i gamma^{jk} = -gamma^{ja} (d_i gamma_ab) gamma^{bk}.
//   * [MCNX-GEO-02] null closure: p^t = sqrt(gamma^{ij} p_i p_j)/alpha,
//     recomputed inside every right-hand-side evaluation (never stored) —
//     through the single closure definition of mcnux_tetrad.hxx.
//   * [MCNX-GEO-03] gather contract: trilinear interpolation of the
//     vertex-centered ADMBaseX lapse/shift/metric, derivatives by
//     differentiating the interpolant; exact (machine tier) on constant and
//     linear fields; ONE scheme feeds every term of one RHS evaluation.
//   * [MCNX-GEO-04] integrator: classical RK4 (order 4 >= the required 2),
//     re-gathering at every stage position; in flat spacetime dp_i/dt
//     vanishes identically and the push is a straight line at |dx/dt| = 1.
//
// Layering: MetricSnapshot and geodesic_rhs are pure functions of their
// arguments; geodesic_step is templated on the gather so that a higher-order
// gather (T16's Schwarzschild legs) can be swapped in without touching the
// RHS or the integrator, and so that the self-test battery can drive the
// integrator with analytic metrics. The only AMReX dependence is the
// amrex::Array4 views and the host/device annotations of the runtime gather
// functor (the mcnux_particles.hxx precedent); CarpetX's interpolator
// (interp.hxx/interpolate.cxx) is deliberately NOT reused ([MCNX-GPU-07]).
//
// Conventions: signature (-,+,+,+); alpha is the lapse; geometrized units
// G = c = M_sun = 1; Latin indices spatial; the vertex-centered Array4 global
// index i sits at x = prob_lo + i dx (the plain AMReX NODE convention, which
// is also what CarpetX's TSV writer and GridDesc use for vertex-centered
// groups: CarpetX/CarpetX/src/io_tsv.cxx, schedule.cxx).

#include "mcnux_tetrad.hxx" // SpatialMetric, spatial_metric_inverse, p_t_closure_with
#include "mcnux_units.hxx"  // detail::approx_eq, detail::rtol_machine

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// The gathered metric snapshot at one point  [MCNX-GEO-03]
// ---------------------------------------------------------------------------
// Every field the right-hand side consumes, gathered by ONE scheme at ONE
// position. Index conventions: beta[k] = beta^k, dalpha[i] = d_i alpha,
// dbeta[i][k] = d_i beta^k, dg[i] = d_i gamma_jk (six components each, in
// the gxx, gxy, gxz, gyy, gyz, gzz order of SpatialMetric).
struct MetricSnapshot {
  double alpha;
  double beta[3];
  SpatialMetric g;
  double dalpha[3];
  double dbeta[3][3];
  SpatialMetric dg[3];
};

// The exact flat snapshot (alpha = 1, beta = 0, gamma = delta, all
// derivatives 0): the [MCNX-GEO-04] flat limit and the self-test anchor.
constexpr MetricSnapshot flat_metric_snapshot() noexcept {
  MetricSnapshot s{};
  s.alpha = 1.0;
  s.g = SpatialMetric{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  return s;
}

// ---------------------------------------------------------------------------
// Right-hand side  [MCNX-GEO-01], [MCNX-GEO-02]
// ---------------------------------------------------------------------------

struct GeodesicRhs {
  double dxdt[3]; // dx^i/dt
  double dpdt[3]; // dp_i/dt
  double pt;      // the p^t of this evaluation (closure; never state)
};

namespace detail {

// The production square root as a host/device callable object (a lambda
// would do on a CPU toolchain, but the explicit functor keeps the device
// annotation visible).
struct SqrtStd {
  AMREX_GPU_HOST_DEVICE double operator()(double v) const noexcept {
    return std::sqrt(v);
  }
};

// No-op stage observer for geodesic_step_with (production path).
struct NoStageObserver {
  AMREX_GPU_HOST_DEVICE constexpr void
  operator()(const MetricSnapshot &, const double *, const GeodesicRhs &) const
      noexcept {}
};

} // namespace detail

// The equations of motion of [MCNX-GEO-01], verbatim, with p^t closed by
// [MCNX-GEO-02] inside the call (never stored). Templated on the sqrt so the
// compile-time self-tests below can run it with detail::csqrt; production
// callers use geodesic_rhs. The inverse-metric derivative is eliminated
// through the spec's own identity d_i gamma^{jk} = -gamma^{ja} (d_i gamma_ab)
// gamma^{bk}: contracting it with p_j p_k gives -p^a (d_i gamma_ab) p^b with
// p^a = gamma^{aj} p_j, so the third term of dp_i/dt is
//   -1/2 (d_i gamma^{jk}) p_j p_k / p^t = +1/2 p^a p^b (d_i gamma_ab) / p^t .
template <class SqrtF>
constexpr GeodesicRhs geodesic_rhs_with(const MetricSnapshot &m,
                                        const double p[3],
                                        SqrtF sqrt_fn) noexcept {
  const InverseSpatialMetric gu = spatial_metric_inverse(m.g);
  const double pt = p_t_closure_with(p[0], p[1], p[2], gu, m.alpha, sqrt_fn);

  // Raised momentum p^a = gamma^{aj} p_j.
  const double pu[3] = {gu.xx * p[0] + gu.xy * p[1] + gu.xz * p[2],
                        gu.xy * p[0] + gu.yy * p[1] + gu.yz * p[2],
                        gu.xz * p[0] + gu.yz * p[1] + gu.zz * p[2]};

  GeodesicRhs r{};
  r.pt = pt;
  for (int i = 0; i < 3; ++i)
    r.dxdt[i] = pu[i] / pt - m.beta[i];

  for (int i = 0; i < 3; ++i) {
    const SpatialMetric &dg = m.dg[i];
    // p^a p^b d_i gamma_ab (symmetric: off-diagonal terms twice).
    const double pdgp =
        dg.xx * pu[0] * pu[0] + dg.yy * pu[1] * pu[1] + dg.zz * pu[2] * pu[2] +
        2.0 * (dg.xy * pu[0] * pu[1] + dg.xz * pu[0] * pu[2] +
               dg.yz * pu[1] * pu[2]);
    const double dbeta_p =
        m.dbeta[i][0] * p[0] + m.dbeta[i][1] * p[1] + m.dbeta[i][2] * p[2];
    r.dpdt[i] = -m.alpha * m.dalpha[i] * pt + dbeta_p + 0.5 * pdgp / pt;
  }
  return r;
}

AMREX_GPU_HOST_DEVICE inline GeodesicRhs
geodesic_rhs(const MetricSnapshot &m, const double p[3]) noexcept {
  return geodesic_rhs_with(m, p, detail::SqrtStd{});
}

// ---------------------------------------------------------------------------
// Integrator  [MCNX-GEO-04]
// ---------------------------------------------------------------------------
// Classical RK4 on the coupled (x^i, p_i) system over one push interval dt.
// `gather(x, y, z) -> MetricSnapshot` is re-evaluated at every stage
// position, so each RHS evaluation sees one consistent gathered snapshot
// ([MCNX-GEO-03]: no mixing of schemes within a step). `observe(snapshot, p,
// rhs)` is called once per stage (verification hook; a no-op in production).
//
// Flat-spacetime exactness mechanism: with all gathered derivatives exactly
// zero every stage has dp/dt == 0 exactly, so the stage momenta and the
// final momentum equal the input bitwise, every stage velocity is the same
// p^i/|p|, and the position update is x + dt * ((k + 2(k + k) + k)/6) — a
// straight line to roundoff. Keep the stage arithmetic in this simple form.
template <class Gather, class SqrtF, class Observer>
constexpr void geodesic_step_with(const Gather &gather, double x[3],
                                  double p[3], double dt, SqrtF sqrt_fn,
                                  Observer &&observe) noexcept {
  GeodesicRhs k[4]{};
  double xs[3], ps[3];

  const MetricSnapshot m1 = gather(x[0], x[1], x[2]);
  k[0] = geodesic_rhs_with(m1, p, sqrt_fn);
  observe(m1, p, k[0]);

  for (int i = 0; i < 3; ++i) {
    xs[i] = x[i] + 0.5 * dt * k[0].dxdt[i];
    ps[i] = p[i] + 0.5 * dt * k[0].dpdt[i];
  }
  const MetricSnapshot m2 = gather(xs[0], xs[1], xs[2]);
  k[1] = geodesic_rhs_with(m2, ps, sqrt_fn);
  observe(m2, ps, k[1]);

  for (int i = 0; i < 3; ++i) {
    xs[i] = x[i] + 0.5 * dt * k[1].dxdt[i];
    ps[i] = p[i] + 0.5 * dt * k[1].dpdt[i];
  }
  const MetricSnapshot m3 = gather(xs[0], xs[1], xs[2]);
  k[2] = geodesic_rhs_with(m3, ps, sqrt_fn);
  observe(m3, ps, k[2]);

  for (int i = 0; i < 3; ++i) {
    xs[i] = x[i] + dt * k[2].dxdt[i];
    ps[i] = p[i] + dt * k[2].dpdt[i];
  }
  const MetricSnapshot m4 = gather(xs[0], xs[1], xs[2]);
  k[3] = geodesic_rhs_with(m4, ps, sqrt_fn);
  observe(m4, ps, k[3]);

  for (int i = 0; i < 3; ++i) {
    x[i] += dt * ((k[0].dxdt[i] + 2.0 * (k[1].dxdt[i] + k[2].dxdt[i]) +
                   k[3].dxdt[i]) /
                  6.0);
    p[i] += dt * ((k[0].dpdt[i] + 2.0 * (k[1].dpdt[i] + k[2].dpdt[i]) +
                   k[3].dpdt[i]) /
                  6.0);
  }
}

// Production entry point: std::sqrt, no observer. Host/device.
template <class Gather>
AMREX_GPU_HOST_DEVICE inline void geodesic_step(const Gather &gather,
                                                double x[3], double p[3],
                                                double dt) noexcept {
  geodesic_step_with(gather, x, p, dt, detail::SqrtStd{},
                     detail::NoStageObserver{});
}

// ---------------------------------------------------------------------------
// Vertex-centered trilinear gather  [MCNX-GEO-03]
// ---------------------------------------------------------------------------
// A trivially-copyable device functor over the three ADMBaseX groups, held as
// amrex::Array4<const Real> views (component order = declaration order in
// CarpetX/ADMBaseX/interface.ccl: metric = gxx gxy gxz gyy gyz gzz, lapse =
// alp, shift = betax betay betaz) plus the level geometry. The views are
// those of the packet's owning box (mfab->array(box index)), so the 8
// corners of the containing cell are always inside the nodal fab for a
// redistributed packet; an AMREX_ASSERT guards the view bounds.
//
// Scheme (documented per the spec): with s_d = (x_d - prob_lo_d)/dx_d,
// anchor i0_d = floor(s_d), weight f_d = s_d - i0_d, the value is the nested
// linear interpolation lerp(a, b, f) = a + f (b - a) along x, then y, then z
// of the 8 corner values, and d_x is the nested lerp along y, then z of the
// 4 corner differences (F(i0+1) - F(i0))/dx_x (likewise d_y, d_z). Writing
// it as differences makes the derivatives of a constant field identically
// zero and its value bitwise exact (lerp(C, C, f) == C), which is what
// makes flat spacetime exact rather than merely machine-close; linear
// fields are reproduced to roundoff. Trilinear derivatives are first-order
// in h on curved metrics — admissible per geodesic-propagation.md:102-105.
struct VertexMetricGather {
  amrex::Array4<const amrex::Real> metric; // 6 components
  amrex::Array4<const amrex::Real> lapse;  // 1 component
  amrex::Array4<const amrex::Real> shift;  // 3 components
  double prob_lo[3];
  double dx[3];

  AMREX_GPU_HOST_DEVICE MetricSnapshot operator()(double x, double y,
                                                  double z) const noexcept {
    using std::floor;
    const double pos[3] = {x, y, z};
    int i0[3];
    double f[3];
    for (int d = 0; d < 3; ++d) {
      const double s = (pos[d] - prob_lo[d]) / dx[d];
      const double fl = floor(s);
      i0[d] = static_cast<int>(fl);
      f[d] = s - fl;
    }
    AMREX_ASSERT(metric.contains(i0[0], i0[1], i0[2]) &&
                 metric.contains(i0[0] + 1, i0[1] + 1, i0[2] + 1) &&
                 lapse.contains(i0[0], i0[1], i0[2]) &&
                 lapse.contains(i0[0] + 1, i0[1] + 1, i0[2] + 1) &&
                 shift.contains(i0[0], i0[1], i0[2]) &&
                 shift.contains(i0[0] + 1, i0[1] + 1, i0[2] + 1));

    MetricSnapshot m{};
    double dv[3];
    m.alpha = sample(lapse, 0, i0, f, dv);
    for (int d = 0; d < 3; ++d)
      m.dalpha[d] = dv[d];
    for (int k = 0; k < 3; ++k) {
      m.beta[k] = sample(shift, k, i0, f, dv);
      for (int d = 0; d < 3; ++d)
        m.dbeta[d][k] = dv[d];
    }
    double *const g[6] = {&m.g.xx, &m.g.xy, &m.g.xz,
                          &m.g.yy, &m.g.yz, &m.g.zz};
    for (int c = 0; c < 6; ++c) {
      *g[c] = sample(metric, c, i0, f, dv);
      double *const dg[3][6] = {
          {&m.dg[0].xx, &m.dg[0].xy, &m.dg[0].xz, &m.dg[0].yy, &m.dg[0].yz,
           &m.dg[0].zz},
          {&m.dg[1].xx, &m.dg[1].xy, &m.dg[1].xz, &m.dg[1].yy, &m.dg[1].yz,
           &m.dg[1].zz},
          {&m.dg[2].xx, &m.dg[2].xy, &m.dg[2].xz, &m.dg[2].yy, &m.dg[2].yz,
           &m.dg[2].zz}};
      for (int d = 0; d < 3; ++d)
        *dg[d][c] = dv[d];
    }
    return m;
  }

private:
  AMREX_GPU_HOST_DEVICE static constexpr double lerp(double a, double b,
                                                     double f) noexcept {
    return a + f * (b - a);
  }

  // Value and gradient of component `comp` of one view at the anchored
  // position; the nested-lerp scheme documented above.
  AMREX_GPU_HOST_DEVICE double sample(const amrex::Array4<const amrex::Real> &a,
                                      int comp, const int i0[3],
                                      const double f[3],
                                      double dv[3]) const noexcept {
    double c[2][2][2];
    for (int cz = 0; cz < 2; ++cz)
      for (int cy = 0; cy < 2; ++cy)
        for (int cx = 0; cx < 2; ++cx)
          c[cx][cy][cz] = a(i0[0] + cx, i0[1] + cy, i0[2] + cz, comp);

    // Value: x, then y, then z.
    const double v = lerp(
        lerp(lerp(c[0][0][0], c[1][0][0], f[0]),
             lerp(c[0][1][0], c[1][1][0], f[0]), f[1]),
        lerp(lerp(c[0][0][1], c[1][0][1], f[0]),
             lerp(c[0][1][1], c[1][1][1], f[0]), f[1]),
        f[2]);

    // d/dx: differences along x, lerped along y then z.
    dv[0] = lerp(lerp((c[1][0][0] - c[0][0][0]) / dx[0],
                      (c[1][1][0] - c[0][1][0]) / dx[0], f[1]),
                 lerp((c[1][0][1] - c[0][0][1]) / dx[0],
                      (c[1][1][1] - c[0][1][1]) / dx[0], f[1]),
                 f[2]);
    // d/dy: differences along y, lerped along x then z.
    dv[1] = lerp(lerp((c[0][1][0] - c[0][0][0]) / dx[1],
                      (c[1][1][0] - c[1][0][0]) / dx[1], f[0]),
                 lerp((c[0][1][1] - c[0][0][1]) / dx[1],
                      (c[1][1][1] - c[1][0][1]) / dx[1], f[0]),
                 f[2]);
    // d/dz: differences along z, lerped along x then y.
    dv[2] = lerp(lerp((c[0][0][1] - c[0][0][0]) / dx[2],
                      (c[1][0][1] - c[1][0][0]) / dx[2], f[0]),
                 lerp((c[0][1][1] - c[0][1][0]) / dx[2],
                      (c[1][1][1] - c[1][1][0]) / dx[2], f[0]),
                 f[1]);
    return v;
  }
};

static_assert(std::is_trivially_copyable_v<VertexMetricGather>,
              "[MCNX-GPU-02] the gather functor is captured by value into "
              "packet kernels and must be trivially copyable");
static_assert(std::is_trivially_copyable_v<MetricSnapshot> &&
                  std::is_trivially_copyable_v<GeodesicRhs>,
              "snapshot and RHS PODs must be trivially copyable");

// ---------------------------------------------------------------------------
// Compile-time verification (flat limit of [MCNX-GEO-01/04])
// ---------------------------------------------------------------------------

namespace detail {

struct FlatGather {
  constexpr MetricSnapshot operator()(double, double, double) const noexcept {
    return flat_metric_snapshot();
  }
};

struct CSqrt {
  constexpr double operator()(double v) const noexcept { return csqrt(v); }
};

// On the flat snapshot dp_i/dt is identically zero and |dx/dt| = 1
// (geodesic-propagation.md:113-115). p = (0, 0, 2) has an exact norm, so the
// velocity is exactly (0, 0, 1).
constexpr bool flat_rhs_holds() noexcept {
  const double p[3] = {0.0, 0.0, 2.0};
  const GeodesicRhs r = geodesic_rhs_with(flat_metric_snapshot(), p, CSqrt{});
  if (r.dpdt[0] != 0.0 || r.dpdt[1] != 0.0 || r.dpdt[2] != 0.0)
    return false;
  if (r.dxdt[0] != 0.0 || r.dxdt[1] != 0.0 || r.dxdt[2] != 1.0)
    return false;
  return r.pt == 2.0;
}

// One RK4 step on flat data with an axis-aligned exact-norm momentum: the
// momentum is bitwise unchanged and the position advances by exactly dt.
constexpr bool flat_step_is_exact() noexcept {
  double x[3] = {0.125, -0.25, 0.0625};
  double p[3] = {0.0, 0.0, 2.0};
  geodesic_step_with(FlatGather{}, x, p, 0.125, CSqrt{}, NoStageObserver{});
  return x[0] == 0.125 && x[1] == -0.25 && x[2] == 0.0625 + 0.125 &&
         p[0] == 0.0 && p[1] == 0.0 && p[2] == 2.0;
}

} // namespace detail

static_assert(detail::flat_rhs_holds(),
              "flat-spacetime right-hand side must give dp_i/dt == 0 and "
              "|dx/dt| == 1 exactly ([MCNX-GEO-04])");
static_assert(detail::flat_step_is_exact(),
              "one RK4 step on flat data must be an exact straight line with "
              "bitwise-constant momentum ([MCNX-GEO-04])");

} // namespace MCNuX

#endif // MCNUX_GEODESIC_HXX
