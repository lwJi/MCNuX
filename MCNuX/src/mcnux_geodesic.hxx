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
#include <limits>
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
  // The domain's upper corner (geom.ProbHi), carried purely for the
  // [MCNX-GPU-05] escape predicate of mcnux_escape.hxx at the push/episode
  // sites; no gather arithmetic reads it.
  double prob_hi[3];
  double dx[3];

  // Un-pinned entry point: floor the anchor from the position, then sample.
  // f[d] = s - fl with fl = floor(s) equals s - double(i0[d]) bitwise (the
  // int round-trips exactly), so this path is arithmetic-identical to the
  // pre-at_anchor single-function form.
  AMREX_GPU_HOST_DEVICE MetricSnapshot operator()(double x, double y,
                                                  double z) const noexcept {
    using std::floor;
    const double pos[3] = {x, y, z};
    int i0[3];
    for (int d = 0; d < 3; ++d)
      i0[d] = static_cast<int>(floor((pos[d] - prob_lo[d]) / dx[d]));
    return at_anchor(i0, x, y, z);
  }

  // Anchored evaluation: sample with a CALLER-pinned anchor cell i0. The
  // weight f[d] = (x_d - prob_lo_d)/dx_d - i0_d is NOT clamped to [0, 1]:
  // for positions slightly outside the anchored cell the trilinear
  // polynomial of that cell extrapolates smoothly (each sub-step of the
  // cellwise integrator below then integrates ONE C-infinity polynomial,
  // satisfying the one-scheme-per-RHS contract of [MCNX-GEO-03]), and the
  // difference-form lerp keeps constant fields bitwise exact for ANY f —
  // the flat-exactness mechanism documented above is anchor-independent.
  AMREX_GPU_HOST_DEVICE MetricSnapshot at_anchor(const int i0[3], double x,
                                                 double y,
                                                 double z) const noexcept {
    const double pos[3] = {x, y, z};
    double f[3];
    for (int d = 0; d < 3; ++d)
      f[d] = (pos[d] - prob_lo[d]) / dx[d] - i0[d];
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

// Anchor-pinning adapter: presents one trilinear piece of the vertex gather
// as a position-only gather, so the geodesic_step_with template above is
// reused UNMODIFIED as the per-sub-step integrator of the cellwise driver.
struct AnchoredVertexGather {
  VertexMetricGather g;
  int i0[3];

  AMREX_GPU_HOST_DEVICE MetricSnapshot operator()(double x, double y,
                                                  double z) const noexcept {
    return g.at_anchor(i0, x, y, z);
  }
};

namespace detail {

// Minimum time along the FROZEN coordinate velocity dxdt to leave the box
// [face_lo, face_hi]: per axis (face - x)/dxdt toward the faced boundary,
// +inf for a zero component (branch BEFORE any division: 0/0 must never
// reach the FPU). This mirrors dt_to_cell_exit_faces of
// mcnux_interactions.hxx verbatim; it is duplicated here rather than
// included because mcnux_geodesic.hxx must not depend on the interactions
// layer (layering inversion) — a future-consolidation candidate for a
// shared kinematics header.
constexpr double dt_to_face_exit(const double x[3], const double dxdt[3],
                                 const double face_lo[3],
                                 const double face_hi[3]) noexcept {
  double dt_min = std::numeric_limits<double>::infinity();
  for (int d = 0; d < 3; ++d) {
    double cand = std::numeric_limits<double>::infinity();
    if (dxdt[d] > 0.0)
      cand = (face_hi[d] - x[d]) / dxdt[d];
    else if (dxdt[d] < 0.0)
      cand = (face_lo[d] - x[d]) / dxdt[d];
    if (cand < dt_min)
      dt_min = cand;
  }
  return dt_min;
}

} // namespace detail

// ---------------------------------------------------------------------------
// Face-aware substepped push  [MCNX-GEO-03]/[MCNX-GEO-04]
// ---------------------------------------------------------------------------
// The production push entry point on the trilinear gather. The trilinear
// interpolant is only C0 across the vertex planes x_d = prob_lo_d + i dx_d
// (the NODE convention above): its derivative fields jump there, so an RK4
// step integrated blindly across a face injects an O(dt * [f]) local error
// at a pseudo-random phase and the accumulated p_t drift degrades to a
// first-order random walk (measured on the schwarzschild-pt benchmark).
// This driver therefore splits the push interval at face crossings so that
// every RK4 (sub-)step integrates a SINGLE trilinear piece — one smooth
// polynomial, one gather scheme per RHS evaluation ([MCNX-GEO-03]) — which
// restores the integrator's design order on curved metrics.
//
// Per sub-step: (1) anchor i0 = floor((x - prob_lo)/dx) and one frozen RHS
// at the sub-step start; a packet AT-or-below its cell's lower face moving
// toward -d is re-anchored to the LOWER neighbor (the direction-aware entry
// adjustment of mcnux_interactions.cxx — without it an exact face landing
// with inward motion yields an infinite dt_exit = 0 loop). (2) The exit
// bound is the frozen-velocity time to leave [face_lo, face_hi], inflated
// by 1 + 4 eps so a crossing sub-step lands in (or bitwise on) the entered
// cell instead of roundoff-short of the face (the boundary-crawl mechanism
// documented at mcnux_interactions.cxx). (3) If the bound covers the
// remainder, one anchored RK4 step finishes the interval — a step that
// crosses no face floors to the same i0 and computes the same f as the
// un-pinned gather, so no-crossing steps are arithmetic-identical to the
// plain geodesic_step (the flat benchmark stays within roundoff: on
// constant fields the anchored values are bitwise anchor-independent).
// Otherwise the crossing time is located precisely before the sub-step is
// taken, because the accuracy of the located crossing bounds the accuracy
// of the whole push: the derivative fields JUMP by O(1) across the face,
// so a sub-step that ends a time delta-t past the true crossing integrates
// the wrong polynomial's dp_i/dt for that long — an O(delta-t) error with
// a pseudo-random sign at every crossing. The frozen-velocity bound alone
// locates the crossing only to O(dt^2) and the accumulated random walk
// caps the push at observed order ~2 (measured on the schwarzschild-pt
// legs; a one-shot Heun-averaged velocity still left order-breaking
// realizations at ~1e-11). The production sequence is therefore
// (a) Heun refinement — re-evaluate the RHS at the frozen-estimate
// crossing state on the SAME anchored polynomial and redo the bound with
// the averaged velocity (also fixing the crossing DIRECTION) — then
// (b) Newton polish on the actual integrated trajectory: trial anchored
// RK4 steps over the estimate, updating it by the landing's face residual
// over the averaged velocity until the landing sits within a few ulps of
// the face. The polished sub-step ends on the face to machine precision,
// so the wrong-polynomial interval shrinks to the 4-eps inflation tier
// (~1e-16 relative over a whole leg) and the push recovers the clean
// smooth-RK4 O(dt^4) drift (measured: ~1e-11 at Delta t = M/16 on the
// schwarzschild-pt base leg). Then one anchored RK4 sub-step of dt_exit,
// and loop. RK4 stage positions may still overshoot the face by O(dt^2)
// mid-step; the anchored gather extrapolates its polynomial smoothly
// across it (a consistent one-scheme RK4 on the extended polynomial). A
// safety cap bounds the sub-step count; on hitting it the whole remainder
// is taken in one anchored step — residual time is NEVER dropped (dropped
// time would silently corrupt the convergence observable).
AMREX_GPU_HOST_DEVICE inline void
geodesic_step_cellwise(const VertexMetricGather &gather, double x[3],
                       double p[3], double dt) noexcept {
  using std::floor;
  constexpr int max_substeps = 64;
  constexpr double exit_inflate =
      1.0 + 4.0 * std::numeric_limits<double>::epsilon();

  double dt_rem = dt;
  for (int sub = 0; sub < max_substeps; ++sub) {
    AnchoredVertexGather ag{gather, {0, 0, 0}};
    double face_lo[3], face_hi[3];
    for (int d = 0; d < 3; ++d) {
      ag.i0[d] =
          static_cast<int>(floor((x[d] - gather.prob_lo[d]) / gather.dx[d]));
      face_lo[d] = gather.prob_lo[d] + ag.i0[d] * gather.dx[d];
    }

    // Frozen sub-step-start velocity (the episode-driver convention): one
    // RHS at x on the floored anchor — identical to the un-pinned gather.
    const MetricSnapshot m = ag(x[0], x[1], x[2]);
    const GeodesicRhs rhs = geodesic_rhs(m, p);

    // Direction-aware entry adjustment (mcnux_interactions.cxx precedent).
    for (int d = 0; d < 3; ++d) {
      if (x[d] <= face_lo[d] && rhs.dxdt[d] < 0.0) {
        --ag.i0[d];
        face_lo[d] -= gather.dx[d];
      }
      face_hi[d] = face_lo[d] + gather.dx[d];
    }

    const double dt0 = detail::dt_to_face_exit(x, rhs.dxdt, face_lo, face_hi);
    if (dt0 * exit_inflate >= dt_rem || sub == max_substeps - 1) {
      geodesic_step(ag, x, p, dt_rem);
      return;
    }

    // (a) Heun refinement of the crossing time (see the header comment):
    // one Euler predictor to the frozen-estimate crossing state, an RHS on
    // the same anchored polynomial there, and the bound redone with the
    // averaged velocity — tracking the crossing direction and target face
    // for the Newton polish below.
    double xc[3], pc[3];
    for (int d = 0; d < 3; ++d) {
      xc[d] = x[d] + dt0 * rhs.dxdt[d];
      pc[d] = p[d] + dt0 * rhs.dpdt[d];
    }
    const GeodesicRhs rhc = geodesic_rhs(ag(xc[0], xc[1], xc[2]), pc);
    double vavg[3];
    for (int d = 0; d < 3; ++d)
      vavg[d] = 0.5 * (rhs.dxdt[d] + rhc.dxdt[d]);

    double dt_est = std::numeric_limits<double>::infinity();
    int dstar = -1;
    for (int d = 0; d < 3; ++d) {
      double cand = std::numeric_limits<double>::infinity();
      if (vavg[d] > 0.0)
        cand = (face_hi[d] - x[d]) / vavg[d];
      else if (vavg[d] < 0.0)
        cand = (face_lo[d] - x[d]) / vavg[d];
      if (cand < dt_est) {
        dt_est = cand;
        dstar = d;
      }
    }
    // Anti-stall guard: the entry adjustment guarantees a strictly positive
    // bound only for the FROZEN velocity; if an averaged component flips
    // sign against an exactly-landed face the refined bound can degenerate
    // to 0 (or lose every finite candidate) — fall back to the frozen
    // bound and its direction rather than loop in place.
    double slope;
    if (dt_est > 0.0 && dstar >= 0) {
      slope = vavg[dstar];
    } else {
      dt_est = dt0;
      dstar = -1;
      slope = 0.0;
      for (int d = 0; d < 3; ++d) {
        const double cand =
            rhs.dxdt[d] > 0.0   ? (face_hi[d] - x[d]) / rhs.dxdt[d]
            : rhs.dxdt[d] < 0.0 ? (face_lo[d] - x[d]) / rhs.dxdt[d]
                                : std::numeric_limits<double>::infinity();
        if (cand == dt0 && dstar < 0) {
          dstar = d;
          slope = rhs.dxdt[d];
        }
      }
    }
    if (dt_est * exit_inflate >= dt_rem) {
      geodesic_step(ag, x, p, dt_rem);
      return;
    }

    // (b) Newton polish on the actual integrated trajectory: trial
    // anchored RK4 steps, residual = landing minus target face, slope the
    // averaged velocity component (the true trajectory slope to O(dt) —
    // ample for a few near-quadratic corrections). Converged when the
    // landing is within a few ulps of the face; a correction leaving
    // (0, dt_rem) stops the polish with the last in-range estimate.
    {
      const double face_target = slope > 0.0 ? face_hi[dstar] : face_lo[dstar];
      const double res_tol = 4.0 * std::numeric_limits<double>::epsilon() *
                             (std::fabs(face_target) + gather.dx[dstar]);
      for (int nit = 0; nit < 4; ++nit) {
        double xt[3] = {x[0], x[1], x[2]};
        double pt[3] = {p[0], p[1], p[2]};
        geodesic_step(ag, xt, pt, dt_est);
        const double resid = xt[dstar] - face_target;
        if (std::fabs(resid) <= res_tol)
          break;
        const double dt_new = dt_est - resid / slope;
        if (!(dt_new > 0.0) || dt_new >= dt_rem)
          break;
        dt_est = dt_new;
      }
    }

    const double dt_exit = dt_est * exit_inflate;
    if (dt_exit >= dt_rem) {
      geodesic_step(ag, x, p, dt_rem);
      return;
    }
    geodesic_step(ag, x, p, dt_exit);
    dt_rem -= dt_exit;
  }
}

static_assert(std::is_trivially_copyable_v<VertexMetricGather> &&
                  std::is_trivially_copyable_v<AnchoredVertexGather>,
              "[MCNX-GPU-02] the gather functors are captured by value into "
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
