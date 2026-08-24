// ---------------------------------------------------------------------------
// Rows 79..85 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the geodesic push of specs/geodesic-propagation.md
// (T15) through mcnux_geodesic.hxx: the vertex-centered trilinear gather
// [MCNX-GEO-03] on synthetic in-memory vertex data (host buffers wrapped as
// amrex::Array4 views over a small vertex box — the in-memory synthetic-table
// discipline of the opacity rows; no grid, no driver), the null-closure
// identity [MCNX-GEO-02] at every RK stage of a curved analytic metric, the
// flat-spacetime exactness of [MCNX-GEO-04] (straight line to rtol_machine,
// momentum bitwise constant, |dx/dt| = 1), and the convergence order of the
// integrator on a smooth analytic metric (RK4 is used; the row asserts the
// spec's order >= 2 bound, ratio >= 3.9 on step halving). Booleans only.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_geodesic.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace MCNuX {

namespace {

// Synthetic vertex box: 5 x 5 x 5 vertices (indices 0..4), anisotropic
// spacing and an off-origin prob_lo so that index/coordinate mix-ups cannot
// pass. Ten fields in ADMBaseX component order (metric 6, lapse 1, shift 3).
struct SyntheticVertexData {
  static constexpr int n = 5;
  static constexpr double prob_lo[3] = {-1.0, -0.5, 0.25};
  static constexpr double dx[3] = {0.25, 0.5, 0.125};
  std::vector<double> metric, lapse, shift; // component-major buffers

  SyntheticVertexData()
      : metric(std::size_t(n) * n * n * 6), lapse(std::size_t(n) * n * n),
        shift(std::size_t(n) * n * n * 3) {}

  static std::size_t at(int i, int j, int k, int c) {
    return std::size_t(i) + std::size_t(n) * (std::size_t(j) + std::size_t(n) * (std::size_t(k) + std::size_t(n) * c));
  }

  // Fill every field as the affine function a + b . (x, y, z) of the vertex
  // coordinates (a, b per field from `coef`, ten rows of four).
  void fill_affine(const double coef[10][4]) {
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          const double x = prob_lo[0] + i * dx[0];
          const double y = prob_lo[1] + j * dx[1];
          const double z = prob_lo[2] + k * dx[2];
          for (int c = 0; c < 10; ++c) {
            const double v = coef[c][0] + coef[c][1] * x + coef[c][2] * y +
                             coef[c][3] * z;
            if (c < 6)
              metric[at(i, j, k, c)] = v;
            else if (c == 6)
              lapse[at(i, j, k, 0)] = v;
            else
              shift[at(i, j, k, c - 7)] = v;
          }
        }
  }

  VertexMetricGather gather() const {
    const amrex::Dim3 lo{0, 0, 0}, hi{n, n, n}; // hi exclusive
    VertexMetricGather g{};
    g.metric = amrex::Array4<const amrex::Real>(metric.data(), lo, hi, 6);
    g.lapse = amrex::Array4<const amrex::Real>(lapse.data(), lo, hi, 1);
    g.shift = amrex::Array4<const amrex::Real>(shift.data(), lo, hi, 3);
    for (int d = 0; d < 3; ++d) {
      g.prob_lo[d] = prob_lo[d];
      g.dx[d] = dx[d];
    }
    return g;
  }
};

// A smooth curved analytic metric with exact derivatives (polynomial, SPD
// and diagonally dominant on |x^i| <~ 1.5): the fixture for the closure
// identity and convergence rows. Host-only test scaffolding.
struct CurvedAnalyticMetric {
  MetricSnapshot operator()(double x, double y, double z) const noexcept {
    MetricSnapshot m{};
    m.alpha = 1.0 + 0.2 * x + 0.1 * y * z;
    m.dalpha[0] = 0.2;
    m.dalpha[1] = 0.1 * z;
    m.dalpha[2] = 0.1 * y;
    m.beta[0] = 0.1 * y;
    m.beta[1] = 0.05 * z;
    m.beta[2] = 0.1 * x * y;
    m.dbeta[1][0] = 0.1;      // d_y beta^x
    m.dbeta[2][1] = 0.05;     // d_z beta^y
    m.dbeta[0][2] = 0.1 * y;  // d_x beta^z
    m.dbeta[1][2] = 0.1 * x;  // d_y beta^z
    m.g = SpatialMetric{1.0 + 0.1 * x * x, 0.1 * z,  0.05 * y,
                        1.0 + 0.2 * y * y, 0.1 * x,  1.0 + 0.1 * z * z};
    m.dg[0] = SpatialMetric{0.2 * x, 0.0, 0.0, 0.0, 0.1, 0.0};
    m.dg[1] = SpatialMetric{0.0, 0.0, 0.05, 0.4 * y, 0.0, 0.0};
    m.dg[2] = SpatialMetric{0.0, 0.1, 0.0, 0.0, 0.0, 0.2 * z};
    return m;
  }
};

} // namespace

void append_geodesic_rows(Battery &b) {
  // --- Rows 79..80: affine fields, values and derivatives at off-vertex
  // points (machine tier). Distinct, exact-binary coefficients per field. ---
  const double coef[10][4] = {
      {1.25, 0.125, -0.0625, 0.03125},  // gxx
      {0.0625, -0.03125, 0.125, 0.25},  // gxy
      {-0.125, 0.25, 0.0625, -0.125},   // gxz
      {1.5, -0.25, 0.125, 0.0625},      // gyy
      {0.03125, 0.0625, -0.125, 0.25},  // gyz
      {1.125, 0.0625, 0.25, -0.03125},  // gzz
      {0.875, 0.25, -0.125, 0.0625},    // alp
      {0.125, -0.0625, 0.25, 0.125},    // betax
      {-0.25, 0.125, 0.0625, -0.25},    // betay
      {0.0625, 0.03125, -0.25, 0.125},  // betaz
  };
  SyntheticVertexData data;
  data.fill_affine(coef);
  const VertexMetricGather gather = data.gather();

  // Query points strictly inside the 4 x 4 x 4 cell box, at generic
  // (non-binary) fractions of a cell, plus one exactly on a vertex.
  const double queries[][3] = {{-0.63, 0.37, 0.41},
                               {-0.07, 1.23, 0.66},
                               {-0.9, -0.45, 0.27},
                               {-0.5, 0.5, 0.5},
                               {-0.26, 1.49, 0.74}};
  bool values_ok = true, derivs_ok = true;
  for (const auto &q : queries) {
    const MetricSnapshot m = gather(q[0], q[1], q[2]);
    double want[10], got[10], dgot[10][3];
    for (int c = 0; c < 10; ++c)
      want[c] = coef[c][0] + coef[c][1] * q[0] + coef[c][2] * q[1] +
                coef[c][3] * q[2];
    const double *const gv[6] = {&m.g.xx, &m.g.xy, &m.g.xz,
                                 &m.g.yy, &m.g.yz, &m.g.zz};
    for (int c = 0; c < 6; ++c) {
      got[c] = *gv[c];
      for (int d = 0; d < 3; ++d) {
        const double *const dgv[6] = {&m.dg[d].xx, &m.dg[d].xy, &m.dg[d].xz,
                                      &m.dg[d].yy, &m.dg[d].yz, &m.dg[d].zz};
        dgot[c][d] = *dgv[c];
      }
    }
    got[6] = m.alpha;
    for (int d = 0; d < 3; ++d)
      dgot[6][d] = m.dalpha[d];
    for (int k = 0; k < 3; ++k) {
      got[7 + k] = m.beta[k];
      for (int d = 0; d < 3; ++d)
        dgot[7 + k][d] = m.dbeta[d][k];
    }
    for (int c = 0; c < 10; ++c) {
      values_ok =
          values_ok && detail::approx_eq(got[c], want[c], detail::rtol_machine);
      for (int d = 0; d < 3; ++d)
        derivs_ok = derivs_ok && detail::approx_eq(dgot[c][d], coef[c][1 + d],
                                                   detail::rtol_machine);
    }
  }
  b.add_boolean("geo.gather.linear_values", values_ok);
  b.add_boolean("geo.gather.linear_derivs", derivs_ok);

  // --- Row 81: constant fields — derivatives identically zero and values
  // reproduced bitwise (exact tier; the mechanism behind flat exactness). ---
  double const_coef[10][4] = {};
  for (int c = 0; c < 10; ++c)
    const_coef[c][0] = coef[c][0];
  SyntheticVertexData cdata;
  cdata.fill_affine(const_coef);
  const VertexMetricGather cgather = cdata.gather();
  bool const_ok = true;
  for (const auto &q : queries) {
    const MetricSnapshot m = cgather(q[0], q[1], q[2]);
    const double *const gv[6] = {&m.g.xx, &m.g.xy, &m.g.xz,
                                 &m.g.yy, &m.g.yz, &m.g.zz};
    for (int c = 0; c < 6; ++c)
      const_ok = const_ok && *gv[c] == coef[c][0];
    const_ok = const_ok && m.alpha == coef[6][0];
    for (int k = 0; k < 3; ++k)
      const_ok = const_ok && m.beta[k] == coef[7 + k][0];
    for (int d = 0; d < 3; ++d) {
      const_ok = const_ok && m.dalpha[d] == 0.0 && m.dg[d].xx == 0.0 &&
                 m.dg[d].xy == 0.0 && m.dg[d].xz == 0.0 && m.dg[d].yy == 0.0 &&
                 m.dg[d].yz == 0.0 && m.dg[d].zz == 0.0;
      for (int k = 0; k < 3; ++k)
        const_ok = const_ok && m.dbeta[d][k] == 0.0;
    }
  }
  b.add_boolean("geo.gather.constant_derivs_zero", const_ok);

  // --- Row 82: null-closure identity alpha^2 (p^t)^2 = gamma^{ij} p_i p_j
  // at EVERY RHS evaluation (all four RK stages of every step) along a
  // curved trajectory, machine tier. ---
  {
    const CurvedAnalyticMetric curved;
    bool closure_ok = true;
    int nstages = 0;
    const auto observe = [&](const MetricSnapshot &m, const double *p,
                             const GeodesicRhs &r) {
      const InverseSpatialMetric gu = spatial_metric_inverse(m.g);
      const double lhs = (m.alpha * r.pt) * (m.alpha * r.pt);
      const double rhs = momentum_norm2(gu, p[0], p[1], p[2]);
      closure_ok = closure_ok && detail::approx_eq(lhs, rhs, detail::rtol_machine);
      ++nstages;
    };
    double x[3] = {0.1, -0.2, 0.15};
    double p[3] = {0.6, 0.3, -0.4};
    for (int step = 0; step < 8; ++step)
      geodesic_step_with(curved, x, p, 0.125, detail::SqrtStd{}, observe);
    b.add_boolean("geo.rhs.null_closure_identity", closure_ok && nstages == 32);
  }

  // --- Rows 83..84: flat spacetime. Straight line x(t) = x0 + t p/|p| to
  // rtol_machine over the benchmark's four steps with p bitwise unchanged;
  // |dx/dt| = 1 to rtol_machine for several momenta. ---
  {
    const detail::FlatGather flat;
    const double x0[3] = {0.125, -0.25, 0.0625};
    const double p0[3] = {0.5, -0.25, 0.125};
    double x[3] = {x0[0], x0[1], x0[2]};
    double p[3] = {p0[0], p0[1], p0[2]};
    const double dt = 0.125;
    const int nsteps = 4;
    bool line_ok = true;
    const double pnorm = std::sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]);
    for (int step = 1; step <= nsteps; ++step) {
      geodesic_step(flat, x, p, dt);
      const double t = step * dt;
      for (int i = 0; i < 3; ++i)
        line_ok = line_ok &&
                  detail::approx_eq(x[i], x0[i] + t * p0[i] / pnorm,
                                    detail::rtol_machine) &&
                  p[i] == p0[i];
    }
    b.add_boolean("geo.push.flat_straight_line", line_ok);

    const double momenta[][3] = {{1.0, 0.0, 0.0},     {0.0, -2.0, 0.0},
                                 {0.0, 0.0, 0.25},    {0.5, -0.25, 0.125},
                                 {-0.375, 0.75, 0.5}, {3.0, 4.0, 12.0}};
    bool speed_ok = true;
    for (const auto &pm : momenta) {
      const GeodesicRhs r = geodesic_rhs(flat_metric_snapshot(), pm);
      const double v2 = r.dxdt[0] * r.dxdt[0] + r.dxdt[1] * r.dxdt[1] +
                        r.dxdt[2] * r.dxdt[2];
      speed_ok = speed_ok && detail::approx_eq(std::sqrt(v2), 1.0, detail::rtol_machine) &&
                 r.dpdt[0] == 0.0 && r.dpdt[1] == 0.0 && r.dpdt[2] == 0.0;
    }
    b.add_boolean("geo.push.flat_unit_speed", speed_ok);
  }

  // --- Row 85: convergence order >= 2 on the curved analytic metric:
  // errors against a 32x finer reference at dt and dt/2 must shrink by at
  // least 3.9 (RK4 gives ~16). ---
  {
    const CurvedAnalyticMetric curved;
    const double x0[3] = {0.1, -0.2, 0.15};
    const double p0[3] = {0.6, 0.3, -0.4};
    const double T = 1.0;
    const auto integrate = [&](int nsteps, double out[6]) {
      double x[3] = {x0[0], x0[1], x0[2]};
      double p[3] = {p0[0], p0[1], p0[2]};
      for (int step = 0; step < nsteps; ++step)
        geodesic_step(curved, x, p, T / nsteps);
      for (int i = 0; i < 3; ++i) {
        out[i] = x[i];
        out[3 + i] = p[i];
      }
    };
    double ref[6], coarse[6], fine[6];
    integrate(128, ref);
    integrate(4, coarse);
    integrate(8, fine);
    double e_coarse = 0.0, e_fine = 0.0;
    for (int i = 0; i < 6; ++i) {
      e_coarse = std::max(e_coarse, detail::cabs(coarse[i] - ref[i]));
      e_fine = std::max(e_fine, detail::cabs(fine[i] - ref[i]));
    }
    const bool order_ok = e_fine > 0.0 && e_coarse / e_fine >= 3.9;
    CCTK_VINFO("MCNuX selftest geo.push.order_ge2: error(dt) = %.3e, "
               "error(dt/2) = %.3e, ratio = %.3f",
               e_coarse, e_fine, e_fine > 0.0 ? e_coarse / e_fine : 0.0);
    b.add_boolean("geo.push.order_ge2", order_ok);
  }
}

} // namespace MCNuX
