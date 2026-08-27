// ---------------------------------------------------------------------------
// Rows 140..145 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the cell-exit-time helper of mcnux_interactions.hxx
// feeding compete_episode's dt_cell_exit ([MCNX-INT-02] of
// specs/neutrino-matter-interactions.md; the algorithm itself is
// implementation freedom, so these rows pin the DOCUMENTED choice: half-open
// [i dx, (i+1) dx) cell membership via the shared floor anchor, per-axis
// (face - x)/dxdt candidates with an explicit +inf branch on zero velocity
// components, min over axes). Exact-tier rows re-assert the axis-aligned flat
// cases through the production dt_to_cell_exit wrapper (the std::floor face
// derivation the header's static_asserts cannot cover) on the off-origin
// anisotropic fixture geometry of the geodesic rows; a boolean row covers the
// +inf sentinels; a boolean row pins the on-lower-face -> dt = 0 -> immediate
// CellExit convention; and a machine-tier row checks agreement with an
// RK4-stepped flat boundary crossing: geodesic_rhs supplies the frozen dxdt,
// one geodesic_step of exactly dt_to_cell_exit lands the exit coordinate on
// the independently predicted face to rtol_machine. Deterministic fixtures
// only — no scheduled transport, no containers, no grid.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_geodesic.hxx"
#include "mcnux_interactions.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <limits>

namespace MCNuX {

namespace {

// Off-origin, anisotropic level geometry — the SyntheticVertexData fixture of
// mcnux_selftest_geodesic.cxx, so index/coordinate mix-ups cannot pass. All
// values exact binary.
constexpr double fix_prob_lo[3] = {-1.0, -0.5, 0.25};
constexpr double fix_dx[3] = {0.25, 0.5, 0.125};

} // namespace

void append_cell_exit_rows(Battery &b) {
  constexpr double inf = std::numeric_limits<double>::infinity();

  // The pinned in-cell position, exact binary per axis. Containing cells and
  // faces (half-open convention, i0 = floor((x - prob_lo)/dx)):
  //   x: i0 = 1, faces [-0.75, -0.5);  y: i0 = 1, faces [0.0, 0.5);
  //   z: i0 = 0, faces [0.25, 0.375).
  const double x0[3] = {-0.625, 0.125, 0.3125};

  // --- Row 140: +x toward the upper face through the production wrapper:
  // dt = (-0.5 - (-0.625))/0.5 = 0.25 exactly. ---
  {
    const double v[3] = {0.5, 0.0, 0.0};
    b.add_exact("int.cellexit.pos_axis_exact",
                dt_to_cell_exit(x0, v, fix_prob_lo, fix_dx), 0.25);
  }

  // --- Row 141: -y toward the lower face: dt = (0.0 - 0.125)/(-0.25) = 0.5
  // exactly. ---
  {
    const double v[3] = {0.0, -0.25, 0.0};
    b.add_exact("int.cellexit.neg_axis_exact",
                dt_to_cell_exit(x0, v, fix_prob_lo, fix_dx), 0.5);
  }

  // --- Row 142: mixed velocity, min over axes: candidates 0.25 (+x),
  // 0.5 (-y), (0.375 - 0.3125)/1.0 = 0.0625 (+z) -> the z crossing wins,
  // exactly. ---
  {
    const double v[3] = {0.5, -0.25, 1.0};
    b.add_exact("int.cellexit.min_axis_exact",
                dt_to_cell_exit(x0, v, fix_prob_lo, fix_dx), 0.0625);
  }

  // --- Row 143: +inf sentinels, never NaN: an all-zero velocity returns
  // +inf (a legal compete_episode input), and zero components contribute
  // +inf so a single finite axis still wins: (0.5 - 0.125)/0.25 = 1.5. ---
  {
    const double v0[3] = {0.0, 0.0, 0.0};
    const double dt0 = dt_to_cell_exit(x0, v0, fix_prob_lo, fix_dx);
    const double v1[3] = {0.0, 0.25, 0.0};
    const double dt1 = dt_to_cell_exit(x0, v1, fix_prob_lo, fix_dx);
    const bool ok = std::isinf(dt0) && dt0 > 0.0 && !std::isnan(dt0) &&
                    dt1 == 1.5;
    b.add_boolean("int.cellexit.zero_velocity_inf", ok);
  }

  // --- Row 144: the documented half-open on-face convention: a packet
  // exactly on its lower x face (x = -0.75, i0 = floor(1.0) = 1) moving in
  // -x yields dt == 0 (a signed zero), and compete_episode with finite
  // channel times sees an immediate CellExit (0 strictly precedes them). ---
  {
    const double xf[3] = {-0.75, 0.125, 0.3125};
    const double v[3] = {-1.0, 0.0, 0.0};
    const double dt = dt_to_cell_exit(xf, v, fix_prob_lo, fix_dx);
    const EpisodeOutcome o = compete_episode(2.0, 3.0, dt, 4.0);
    const bool ok = dt == 0.0 && o.kind == EpisodeEnd::CellExit &&
                    o.dt_event == dt;
    b.add_boolean("int.cellexit.on_face_zero", ok);
  }

  // --- Row 145: RK4 agreement, machine tier. For two generic flat fixtures:
  // dxdt from the episode-start geodesic_rhs (the frozen kinematics of the
  // header), dt_exit from the production wrapper, then ONE geodesic_step of
  // exactly dt_exit; the coordinate of the independently predicted winning
  // axis must land on its predicted face to rtol_machine (flat RK4 is a
  // straight line to roundoff, [MCNX-GEO-04]). The prediction below is
  // test-local arithmetic: per-axis faces from the floor anchor, candidate
  // times, argmin. ---
  {
    const detail::FlatGather flat;
    const double fixtures[2][6] = {
        {-0.6, 0.3, 0.34, 0.5, -0.25, 0.125},
        {-0.2, -0.3, 0.53, -0.3, 0.7, -0.2}};
    bool ok = true;
    for (const auto &fx : fixtures) {
      double x[3] = {fx[0], fx[1], fx[2]};
      double p[3] = {fx[3], fx[4], fx[5]};

      const GeodesicRhs r = geodesic_rhs(flat_metric_snapshot(), p);
      const double dt_exit = dt_to_cell_exit(x, r.dxdt, fix_prob_lo, fix_dx);

      // Independent prediction of the winning axis and its face coordinate.
      int d_star = -1;
      double dt_star = inf, face_star = 0.0;
      for (int d = 0; d < 3; ++d) {
        if (r.dxdt[d] == 0.0)
          continue;
        const double i0 = std::floor((x[d] - fix_prob_lo[d]) / fix_dx[d]);
        const double face =
            fix_prob_lo[d] + (r.dxdt[d] > 0.0 ? i0 + 1.0 : i0) * fix_dx[d];
        const double cand = (face - x[d]) / r.dxdt[d];
        if (cand < dt_star) {
          dt_star = cand;
          d_star = d;
          face_star = face;
        }
      }

      geodesic_step(flat, x, p, dt_exit);
      ok = ok && d_star >= 0 && std::isfinite(dt_exit) && dt_exit > 0.0 &&
           dt_exit == dt_star &&
           detail::approx_eq(x[d_star], face_star, detail::rtol_machine);
    }
    b.add_boolean("int.cellexit.rk4_face_crossing", ok);
  }
}

} // namespace MCNuX
