// ---------------------------------------------------------------------------
// Rows 146..149 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the elastic scatter redraw of mcnux_interactions.hxx
// ([MCNX-INT-04] of specs/neutrino-matter-interactions.md: a new isotropic
// fluid-frame direction at FIXED fluid-frame energy nu, through the same
// pinned tetrad transform as creation). Deterministic fixtures only — no
// scheduled transport, no containers, no grid. The curved/boosted fixture is
// that of mcnux_selftest_tetrad.cxx (SPD gamma with off-diagonals,
// alpha != 1, beta != 0, u^t solved from the normalization), so the machine-
// tier identities are probed away from every flat degeneracy; the flat-
// static-limit row is exact (identity tetrad, linear transform); the draw-
// composition row pins scatter_uniforms -> scatter_redraw against the direct
// u(S, q, e, k) pair bitwise ([MCNX-INT-06] k = 2, 3).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_interactions.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <cstdint>

namespace MCNuX {

void append_scatter_rows(Battery &b) {
  // --- Curved/boosted fixture (the mcnux_selftest_tetrad.cxx inputs). ---
  const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
  const double alp = 0.83;
  const double betau[3] = {0.06, -0.11, 0.04};
  const SpacetimeMetric m = spacetime_metric_from_adm(alp, betau, gamma);
  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);

  // Boosted u^mu: spatial components chosen, u^t fixed by g u u = -1
  // (future-pointing root).
  double uvec[4] = {0.0, 0.15, -0.05, 0.10};
  {
    const double a = m.g[0][0];
    double bb = 0.0, cc = 1.0;
    for (int i = 1; i < 4; ++i) {
      bb += m.g[0][i] * uvec[i];
      for (int j = 1; j < 4; ++j)
        cc += m.g[i][j] * uvec[i] * uvec[j];
    }
    uvec[0] = (bb + std::sqrt(bb * bb - a * cc)) / (-a);
  }

  const Tetrad tet = build_tetrad(alp, betau, gamma, uvec);
  const double nu = 3.25;

  // One redraw with fixed literal uniforms (the FROZEN episode kinematics:
  // tetrad, metric, nu all from the fixture above).
  const CoordinateMomentum pm = scatter_redraw(tet, m, nu, 0.35, 0.62);

  // --- Row 146: [MCNX-INT-04] nu preservation — the fluid-frame energy
  // recomputed from the redraw output with the SAME u^mu/metric equals the
  // input nu at machine tier (the elasticity identity of
  // neutrino-matter-interactions.md:208-211). ---
  {
    const double p_cov[3] = {pm.px, pm.py, pm.pz};
    const double nu_back = fluid_frame_energy(alp, betau, p_cov, pm.pt, uvec);
    b.add_boolean("int.scatter.nu_preserved",
                  detail::approx_eq(nu_back, nu, detail::rtol_machine));
  }

  // --- Row 147: [MCNX-GEO-02] null closure of the redrawn momentum:
  // p^t = sqrt(gamma^{ij} p_i p_j)/alpha recovers the transform's p^t and
  // alpha^2 (p^t)^2 = gamma^{ij} p_i p_j, machine tier. ---
  {
    const double pt = p_t_closure(pm.px, pm.py, pm.pz, gu, alp);
    const double n2 = momentum_norm2(gu, pm.px, pm.py, pm.pz);
    b.add_boolean("int.scatter.null_closure",
                  detail::approx_eq(pt, pm.pt, detail::rtol_machine) &&
                      detail::approx_eq(alp * alp * pt * pt, n2,
                                        detail::rtol_machine));
  }

  // --- Row 148: flat-static-limit isotropy, exact: with Minkowski and
  // u = (1, 0, 0, 0) the tetrad is the identity and the transform linear,
  // so the redraw output is bitwise nu * (sin(theta) cos(phi),
  // sin(theta) sin(phi), cos(theta)) with cos(theta) = 2 u1 - 1 and
  // phi = 2 pi u2 (the [MCNX-PKT-04] draw), and p^t = nu exactly. Exact
  // binary uniforms so cos(theta) is exact. ---
  {
    const SpatialMetric flat{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    const double beta0[3] = {0.0, 0.0, 0.0};
    const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
    const SpacetimeMetric mf = spacetime_metric_from_adm(1.0, beta0, flat);
    const Tetrad tf = build_tetrad(1.0, beta0, flat, ustat);

    const double u1 = 0.375, u2 = 0.25; // exact binary
    const CoordinateMomentum pf = scatter_redraw(tf, mf, nu, u1, u2);

    // Independent composition of the pinned draw formulas. The direction
    // components are formed FIRST (the [MCNX-PKT-04] n_hat), then scaled by
    // nu — the same association as the fluid-frame momentum nu * n_hat, so
    // the comparison is bitwise.
    const double costheta = 2.0 * u1 - 1.0;
    const double sintheta = std::sqrt(1.0 - costheta * costheta);
    const double phi = 2.0 * detail::pi * u2;
    const bool ok = pf.pt == nu && pf.px == nu * (sintheta * std::cos(phi)) &&
                    pf.py == nu * (sintheta * std::sin(phi)) &&
                    pf.pz == nu * costheta;
    b.add_boolean("int.scatter.flat_isotropy", ok);
  }

  // --- Row 149: [MCNX-INT-06] draw composition — the scatter_uniforms
  // wrapper fed through scatter_redraw equals the redraw on the direct
  // u(S, q, e, 2/3) pair bitwise (the row-113 DrawRecord discipline, on the
  // curved fixture so nothing degenerates). ---
  {
    const std::uint64_t S = 1296518744ull;
    const std::uint64_t q = 42ull;
    const std::uint32_t e = 1u;
    const ScatterUniforms su = scatter_uniforms(S, q, e);
    const CoordinateMomentum via_wrap = scatter_redraw(tet, m, nu, su.u1, su.u2);
    const CoordinateMomentum via_direct = scatter_redraw(
        tet, m, nu, u(S, q, e, draw_k_cos_theta), u(S, q, e, draw_k_phi));
    const bool ok = su.u1 == u(S, q, e, 2u) && su.u2 == u(S, q, e, 3u) &&
                    via_wrap.pt == via_direct.pt &&
                    via_wrap.px == via_direct.px &&
                    via_wrap.py == via_direct.py &&
                    via_wrap.pz == via_direct.pz;
    b.add_boolean("int.scatter.draw_composition", ok);
  }
}

} // namespace MCNuX
