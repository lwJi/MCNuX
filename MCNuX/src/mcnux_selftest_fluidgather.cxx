// ---------------------------------------------------------------------------
// Rows 114..117 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the fluid-state/four-velocity gather of
// specs/neutrino-matter-interactions.md through mcnux_fluid.hxx: the
// cell-centered ONE-CELL gather of [MCNX-INT-05] ("never re-interpolated
// sub-cell") on synthetic in-memory ccc data (host buffers wrapped as
// amrex::Array4 views over a small cell box — the in-memory synthetic
// discipline of the geodesic rows; no grid, no driver), and the Valencia
// v^i -> u^mu lift: the v = 0 static limit exactly, and the normalization
// identity g_munu u^mu u^nu = -1 on a non-flat gamma_ij at machine tier —
// the battery's ONLY coverage of the W != 1 branch (every harness fixture
// has vel = 0). Booleans only.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_fluid.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace MCNuX {

namespace {

// Synthetic cell-centered box: 4 x 4 x 4 cells (indices 0..3), anisotropic
// spacing and an off-origin prob_lo (the SyntheticVertexData values) so that
// index/coordinate mix-ups cannot pass. Six fields in the gather's component
// order: rho, velx, vely, velz, temperature, Ye.
struct SyntheticCellData {
  static constexpr int n = 4;
  static constexpr double prob_lo[3] = {-1.0, -0.5, 0.25};
  static constexpr double dx[3] = {0.25, 0.5, 0.125};
  std::vector<double> rho, vel, temperature, ye; // component-major buffers

  SyntheticCellData()
      : rho(std::size_t(n) * n * n), vel(std::size_t(n) * n * n * 3),
        temperature(std::size_t(n) * n * n), ye(std::size_t(n) * n * n) {}

  static std::size_t at(int i, int j, int k, int c) {
    return std::size_t(i) +
           std::size_t(n) * (std::size_t(j) +
                             std::size_t(n) * (std::size_t(k) +
                                               std::size_t(n) * c));
  }

  static double center(int i, int d) {
    return prob_lo[d] + (double(i) + 0.5) * dx[d];
  }

  // Fill every field as the affine function a + b . (xc, yc, zc) of the CELL
  // CENTER coordinates (a, b per field from `coef`, six rows of four).
  void fill_affine(const double coef[6][4]) {
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          const double x = center(i, 0);
          const double y = center(j, 1);
          const double z = center(k, 2);
          for (int c = 0; c < 6; ++c) {
            const double v = coef[c][0] + coef[c][1] * x + coef[c][2] * y +
                             coef[c][3] * z;
            if (c == 0)
              rho[at(i, j, k, 0)] = v;
            else if (c < 4)
              vel[at(i, j, k, c - 1)] = v;
            else if (c == 4)
              temperature[at(i, j, k, 0)] = v;
            else
              ye[at(i, j, k, 0)] = v;
          }
        }
  }

  CellFluidGather gather() const {
    const amrex::Dim3 lo{0, 0, 0}, hi{n, n, n}; // hi exclusive
    CellFluidGather g{};
    g.rho = amrex::Array4<const amrex::Real>(rho.data(), lo, hi, 1);
    g.vel = amrex::Array4<const amrex::Real>(vel.data(), lo, hi, 3);
    g.temperature =
        amrex::Array4<const amrex::Real>(temperature.data(), lo, hi, 1);
    g.ye = amrex::Array4<const amrex::Real>(ye.data(), lo, hi, 1);
    for (int d = 0; d < 3; ++d) {
      g.prob_lo[d] = prob_lo[d];
      g.dx[d] = dx[d];
    }
    return g;
  }
};

// One gathered sample against the stored buffer values of cell (i, j, k),
// bitwise (the six fields in coefficient order).
bool sample_matches_cell(const SyntheticCellData &data, const FluidSample &s,
                         int i, int j, int k) {
  return s.state.rho_cgs == data.rho[SyntheticCellData::at(i, j, k, 0)] &&
         s.vel[0] == data.vel[SyntheticCellData::at(i, j, k, 0)] &&
         s.vel[1] == data.vel[SyntheticCellData::at(i, j, k, 1)] &&
         s.vel[2] == data.vel[SyntheticCellData::at(i, j, k, 2)] &&
         s.state.T_MeV ==
             data.temperature[SyntheticCellData::at(i, j, k, 0)] &&
         s.state.Ye == data.ye[SyntheticCellData::at(i, j, k, 0)];
}

} // namespace

void append_fluidgather_rows(Battery &b) {
  // Distinct, exact-binary affine coefficients per field (a, b_x, b_y, b_z).
  const double coef[6][4] = {
      {2.5, 0.5, -0.25, 0.125},          // rho
      {0.125, 0.0625, -0.03125, 0.25},   // velx
      {-0.0625, 0.125, 0.25, -0.125},    // vely
      {0.25, -0.125, 0.0625, 0.03125},   // velz
      {3.0, -0.5, 0.25, 0.0625},         // temperature
      {0.25, 0.03125, -0.0625, 0.125},   // Ye
  };

  // Query points strictly inside the box, each OFF-CENTER inside its cell
  // (generic non-binary fractions), with the containing cell worked out by
  // hand; plus one exactly at a cell center.
  struct Query {
    double x, y, z;
    int i, j, k;
  };
  const Query queries[] = {
      {-0.93, -0.37, 0.29, 0, 0, 0},  // fractions (0.28, 0.26, 0.32)
      {-0.51, 0.63, 0.66, 1, 2, 3},   // fractions (0.96, 0.26, 0.28)
      {-0.26, -0.45, 0.41, 2, 0, 1},  // fractions (0.96, 0.10, 0.28)
      {-0.07, 1.23, 0.52, 3, 3, 2},   // fractions (0.72, 0.46, 0.16)
      {-0.375, 0.75, 0.5625, 2, 2, 2} // exactly the cell (2,2,2) center
  };

  // --- Row 114: constant fields — every gathered value is the stored value
  // bitwise (exact tier; lerp cannot be blamed for passing this one, but it
  // anchors the read path before the affine tripwire below). ---
  double const_coef[6][4] = {};
  for (int c = 0; c < 6; ++c)
    const_coef[c][0] = coef[c][0];
  SyntheticCellData cdata;
  cdata.fill_affine(const_coef);
  const CellFluidGather cgather = cdata.gather();
  bool const_ok = true;
  for (const Query &q : queries) {
    const FluidSample s = cgather(q.x, q.y, q.z);
    const_ok = const_ok && sample_matches_cell(cdata, s, q.i, q.j, q.k) &&
               s.state.rho_cgs == coef[0][0] && s.vel[0] == coef[1][0] &&
               s.vel[1] == coef[2][0] && s.vel[2] == coef[3][0] &&
               s.state.T_MeV == coef[4][0] && s.state.Ye == coef[5][0];
  }
  b.add_boolean("fluid.gather.constant_exact", const_ok);

  // --- Row 115: affine fields — an off-center query returns EXACTLY the
  // stored value of its containing cell (bitwise), never the affine value at
  // the query position (which is what any interpolating gather would
  // reproduce exactly on affine data): [MCNX-INT-05] "never re-interpolated
  // sub-cell". The inequality legs give the row its teeth. ---
  SyntheticCellData data;
  data.fill_affine(coef);
  const CellFluidGather gather = data.gather();
  bool cell_ok = true;
  for (const Query &q : queries) {
    const FluidSample s = gather(q.x, q.y, q.z);
    int i0[3];
    gather.containing_cell(q.x, q.y, q.z, i0);
    cell_ok = cell_ok && i0[0] == q.i && i0[1] == q.j && i0[2] == q.k &&
              sample_matches_cell(data, s, q.i, q.j, q.k);
    // Tripwire: at the four off-center queries the affine value at the query
    // position differs from the cell-center value (the gradients are
    // nonzero), so an interpolating gather could not pass.
    const bool at_center = q.x == SyntheticCellData::center(q.i, 0) &&
                           q.y == SyntheticCellData::center(q.j, 1) &&
                           q.z == SyntheticCellData::center(q.k, 2);
    if (!at_center) {
      const double rho_at_query =
          coef[0][0] + coef[0][1] * q.x + coef[0][2] * q.y + coef[0][3] * q.z;
      cell_ok = cell_ok && s.state.rho_cgs != rho_at_query;
    }
  }
  b.add_boolean("fluid.gather.cell_value_no_interp", cell_ok);

  // --- Row 116: v = 0 Valencia lift on nontrivial (alpha, beta, gamma):
  // W = 1 exactly (std::sqrt(1.0) == 1.0), so u = (1/alpha, -beta^i/alpha)
  // exactly (exact tier). ---
  {
    const double alpha = 1.25;
    const double beta[3] = {0.5, -0.25, 0.125};
    const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
    const double v0[3] = {0.0, 0.0, 0.0};
    double u[4];
    valencia_four_velocity(alpha, beta, gamma, v0, u);
    const bool static_ok = u[0] == 1.0 / alpha && u[1] == -(beta[0] / alpha) &&
                           u[2] == -(beta[1] / alpha) &&
                           u[3] == -(beta[2] / alpha);
    b.add_boolean("fluid.u.static_lift_exact", static_ok);
  }

  // --- Row 117: nonzero-v lift on a non-flat gamma_ij (the only W != 1
  // coverage in the corpus — all harness fixtures have vel = 0):
  // g_munu u^mu u^nu = -1 via the independent
  // spacetime_metric_from_adm/metric_dot route, and alpha u^t = W with W
  // recomputed by an independent 3x3-loop contraction; machine tier. ---
  {
    const double alpha = 1.2;
    const double beta[3] = {0.1, -0.05, 0.2};
    const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
    const double v[3] = {0.2, -0.1, 0.15};
    double u[4];
    valencia_four_velocity(alpha, beta, gamma, v, u);

    const SpacetimeMetric g4 = spacetime_metric_from_adm(alpha, beta, gamma);
    const double norm = metric_dot(g4, u, u);

    // Independent W: gamma_ij v^i v^j by the full 3x3 loop.
    const double gij[3][3] = {{gamma.xx, gamma.xy, gamma.xz},
                              {gamma.xy, gamma.yy, gamma.yz},
                              {gamma.xz, gamma.yz, gamma.zz}};
    double v2 = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        v2 += gij[i][j] * v[i] * v[j];
    const double W_indep = 1.0 / std::sqrt(1.0 - v2);

    const bool norm_ok =
        detail::approx_eq(norm, -1.0, detail::rtol_machine) && W_indep > 1.0 &&
        detail::approx_eq(alpha * u[0], W_indep, detail::rtol_machine);
    b.add_boolean("fluid.u.normalization", norm_ok);
  }
}

} // namespace MCNuX
