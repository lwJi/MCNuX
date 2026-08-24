// ---------------------------------------------------------------------------
// Rows 53..59 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — tetrad construction + null-closure helpers of
// specs/geodesic-propagation.md [MCNX-GEO-02] and
// specs/packet-representation-and-sampling.md [MCNX-PKT-04] (T13).
//
// Inputs are synthetic and self-contained: a flat-static pair (Minkowski,
// u^mu = (1,0,0,0)) probing the closed-form limit of
// packet-representation-and-sampling.md:187-190, and a curved pair (SPD
// gamma_ij with off-diagonal terms, alpha != 1, beta^i != 0, and a boosted
// u^mu whose u^t is solved at runtime from g_munu u^mu u^nu = -1) probing the
// generic identities. All rows are machine-tier booleans (error 0 or 1), so
// the golden data is bitwise-stable across libm ulp differences. The
// direction draw takes fixed literal uniforms — the helper deliberately
// performs no u(S, q, e, k) call ([MCNX-PKT-05] draw-map bookkeeping belongs
// to the consuming code, not to these helpers).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <cmath>

namespace MCNuX {

void append_tetrad_rows(Battery &b) {
  // Absolute machine-tier closeness for O(1) quantities whose target may be
  // exactly 0 (orthonormality off-diagonals, the null norm), where a relative
  // comparison is undefined.
  const auto near_abs = [](double got, double want) {
    return detail::cabs(got - want) <= detail::rtol_machine;
  };

  // --- Synthetic curved input. ---
  const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
  const double alp = 0.83;
  const double betau[3] = {0.06, -0.11, 0.04};
  const SpacetimeMetric m = spacetime_metric_from_adm(alp, betau, gamma);

  // Boosted u^mu: spatial components chosen, u^t fixed by the normalization
  // g_munu u^mu u^nu = -1 (future-pointing root of the quadratic
  // g_tt (u^t)^2 + 2 g_ti u^i u^t + (gamma_ij u^i u^j + 1) = 0).
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

  // Row 53 — [MCNX-GEO-02] analytic inverse: gamma^{ij} gamma_jk = delta^i_k
  // (machine) and det(gamma_ij) > 0 on the SPD test metric.
  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
  {
    const double lo[3][3] = {{gamma.xx, gamma.xy, gamma.xz},
                             {gamma.xy, gamma.yy, gamma.yz},
                             {gamma.xz, gamma.yz, gamma.zz}};
    const double up[3][3] = {
        {gu.xx, gu.xy, gu.xz}, {gu.xy, gu.yy, gu.yz}, {gu.xz, gu.yz, gu.zz}};
    bool ok = gu.det > 0.0;
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) {
        double s = 0.0;
        for (int j = 0; j < 3; ++j)
          s += up[i][j] * lo[j][k];
        ok = ok && near_abs(s, i == k ? 1.0 : 0.0);
      }
    b.add_boolean("geo.inverse_metric.identity", ok);
  }

  // One drawn packet on the curved background: fixed literal uniforms into
  // the isotropic draw, then the pinned x->y->z Gram-Schmidt tetrad and the
  // [MCNX-PKT-04] transform.
  const double nu = 3.25;
  const UnitVector3 nhat = fluid_frame_direction(0.35, 0.62);
  const Tetrad tet = build_tetrad(alp, betau, gamma, uvec);
  const FourVectorUp pup = fluid_momentum_to_coordinate(tet, nu, nhat);
  const CoordinateMomentum pm = transform_to_coordinate_frame(tet, m, nu, nhat);

  // Row 54 — [MCNX-GEO-02] null closure on the transformed momentum:
  // p^t = sqrt(gamma^{ij} p_i p_j)/alpha recovers the tetrad's p^t (the
  // physical nullness content) and alpha^2 (p^t)^2 = gamma^{ij} p_i p_j
  // (geodesic-propagation.md:96-97), both machine tier.
  {
    const double pt = p_t_closure(pm.px, pm.py, pm.pz, gu, alp);
    const double n2 = momentum_norm2(gu, pm.px, pm.py, pm.pz);
    b.add_boolean("geo.null_closure",
                  detail::approx_eq(pt, pm.pt, detail::rtol_machine) &&
                      detail::approx_eq(alp * alp * pt * pt, n2,
                                        detail::rtol_machine));
  }

  // Row 55 — tetrad orthonormality: g_munu e^mu_(a) e^nu_(b) = eta_ab for
  // all 10 independent pairs, machine tier.
  {
    bool ok = true;
    for (int a = 0; a < 4; ++a)
      for (int c = a; c < 4; ++c) {
        const double eta = a != c ? 0.0 : (a == 0 ? -1.0 : 1.0);
        ok = ok && near_abs(metric_dot(m, tet.e[a], tet.e[c]), eta);
      }
    b.add_boolean("tetrad.orthonormality", ok);
  }

  // Row 56 — transformed momentum is null: g_munu p^mu p^nu = 0 normalized
  // by (alpha p^t)^2 (packet-representation-and-sampling.md:138-140).
  {
    const double norm = metric_dot(m, pup.v, pup.v);
    const double scale = (alp * pm.pt) * (alp * pm.pt);
    b.add_boolean("tetrad.null_norm", near_abs(norm / scale, 0.0));
  }

  // Row 57 — nu recovery: -p_mu u^mu equals the drawn nu, machine tier.
  {
    double pu = 0.0;
    for (int mu = 0; mu < 4; ++mu)
      for (int nu4 = 0; nu4 < 4; ++nu4)
        pu += m.g[mu][nu4] * pup.v[nu4] * uvec[mu];
    b.add_boolean("tetrad.nu_recovery",
                  detail::approx_eq(-pu, nu, detail::rtol_machine));
  }

  // --- Flat-static limit: Minkowski, u^mu = (1, 0, 0, 0). ---
  const SpatialMetric flat{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
  const SpacetimeMetric mf = spacetime_metric_from_adm(1.0, beta0, flat);
  const Tetrad tf = build_tetrad(1.0, beta0, flat, ustat);
  const CoordinateMomentum pf = transform_to_coordinate_frame(tf, mf, nu, nhat);

  // Row 58 — flat-static momentum: p^t = nu and p_i = nu n_hat_i, machine
  // tier (packet-representation-and-sampling.md:187-190).
  b.add_boolean("tetrad.flat_static.momentum",
                detail::approx_eq(pf.pt, nu, detail::rtol_machine) &&
                    detail::approx_eq(pf.px, nu * nhat.x,
                                      detail::rtol_machine) &&
                    detail::approx_eq(pf.py, nu * nhat.y,
                                      detail::rtol_machine) &&
                    detail::approx_eq(pf.pz, nu * nhat.z,
                                      detail::rtol_machine));

  // Row 59 — flat-static tetrad legs equal the coordinate basis, exactly
  // (the x->y->z Gram-Schmidt touches only exact values on this input).
  {
    bool ok = true;
    for (int a = 0; a < 4; ++a)
      for (int mu = 0; mu < 4; ++mu)
        ok = ok && tf.e[a][mu] == (a == mu ? 1.0 : 0.0);
    b.add_boolean("tetrad.flat_static.identity", ok);
  }
}

} // namespace MCNuX
