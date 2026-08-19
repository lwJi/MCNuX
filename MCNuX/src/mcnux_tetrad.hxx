#ifndef MCNUX_TETRAD_HXX
#define MCNUX_TETRAD_HXX

// Tetrad construction and null-closure helpers.
//
// Governing specs:
//   * specs/geodesic-propagation.md [MCNX-GEO-02] — the evolved packet state
//     is (x^i, p_i); p^t is closed by the null condition,
//         p^t = sqrt(gamma^{ij} p_i p_j) / alpha ,
//     RECOMPUTED from the current (p_i, gamma^{ij}, alpha) at every
//     evaluation. Nothing in this header caches or memoizes: every function
//     is a pure function of its arguments (callers may cache the returned
//     PODs, which the spec's implementation-freedom section permits).
//   * specs/packet-representation-and-sampling.md [MCNX-PKT-04] — isotropic
//     fluid-frame direction draw and the orthonormal-tetrad transform to grid
//     coordinates (arXiv:1708.08452 Eqs. 19-21). The tetrad is pinned:
//     timelike leg e^mu_(0') = u^mu (an already-normalized input; its
//     provider is the hydro coupling), spatial legs by Gram-Schmidt of the
//     coordinate basis vectors against g_munu in FIXED coordinate order
//     x, then y, then z — determinism-critical, never reordered.
//
// RNG discipline: no function here calls u(S, q, e, k) or any RNG. The
// direction draw takes its two uniforms as plain arguments; the (e, k)
// draw-map bookkeeping is owned by the consuming specs (emission k = 4/5,
// scattering k = 2/3, diffusion k = 1/2).
//
// Conventions: signature (-,+,+,+); alpha is the LAPSE throughout this header
// (specs/conventions-and-units.md:46-58 — never the trapped-regime fraction
// or a spectral index); geometrized units G = c = M_sun = 1; Latin indices
// spatial; gamma = det(gamma_ij).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// C++ (constexpr where the math is rational; <cmath> sqrt/sin/cos otherwise,
// which nvcc supports in device code) so host code, device code, and
// standalone unit tests can all include it. The static_asserts at the bottom
// are the compile-time half of the T13 test surface; the runtime half lives
// in the `unit-selftest` battery (mcnux_selftest.cxx rows 53..59).

#include "mcnux_units.hxx" // detail::cabs, detail::approx_eq, detail::rtol_machine

#include <cmath>

namespace MCNuX {

namespace detail {

// pi as a mathematical (not physical) constant; the single-literal rule of
// conventions-and-units.md:132 covers physical constants only. Correctly
// rounded binary64 value of pi.
inline constexpr double pi = 3.14159265358979323846;

// Constexpr Newton square root, used ONLY inside compile-time test
// predicates; every production path uses std::sqrt. Monotone-decreasing
// iteration from above, so termination is guaranteed; exact on exact squares
// of interest here (csqrt(1.0) == 1.0 identically, which is what the
// flat-static compile-time tetrad check needs).
constexpr double csqrt(double x) noexcept {
  if (x <= 0.0)
    return 0.0;
  double y = x >= 1.0 ? x : 1.0;
  while (true) {
    const double next = 0.5 * (y + x / y);
    if (next >= y)
      return y;
    y = next;
  }
}

} // namespace detail

// ---------------------------------------------------------------------------
// 3x3 spatial metric and its analytic inverse  [MCNX-GEO-02]
// ---------------------------------------------------------------------------

// Symmetric spatial metric gamma_ij: the six independent components in the
// (gxx, gxy, gxz, gyy, gyz, gzz) order of the ADMBaseX grid functions.
struct SpatialMetric {
  double xx, xy, xz, yy, yz, zz;
};

// Symmetric inverse spatial metric gamma^{ij} plus det(gamma_ij) (the
// Conventions use gamma = det(gamma_ij), e.g. in sqrt(-g) = alpha sqrt(gamma);
// exposing it here saves consumers a second cofactor expansion).
struct InverseSpatialMetric {
  double xx, xy, xz, yy, yz, zz;
  double det; // det(gamma_ij), NOT det(gamma^{ij})
};

// Analytic 3x3 symmetric inverse by adjugate/determinant. Pure and
// constexpr: recomputed from the six gamma_ij at every call, never cached.
// This is the single 3x3 inverse of the thorn — [MCNX-GEO-01/03] consumers
// must call this function rather than duplicating the cofactors.
constexpr InverseSpatialMetric
spatial_metric_inverse(const SpatialMetric &g) noexcept {
  const double cxx = g.yy * g.zz - g.yz * g.yz;
  const double cxy = g.xz * g.yz - g.xy * g.zz;
  const double cxz = g.xy * g.yz - g.xz * g.yy;
  const double cyy = g.xx * g.zz - g.xz * g.xz;
  const double cyz = g.xy * g.xz - g.xx * g.yz;
  const double czz = g.xx * g.yy - g.xy * g.xy;
  const double det = g.xx * cxx + g.xy * cxy + g.xz * cxz;
  return {cxx / det, cxy / det, cxz / det,
          cyy / det, cyz / det, czz / det, det};
}

// gamma^{ij} p_i p_j — the square of alpha p^t under the null closure. Kept
// separate (and constexpr) so the identity alpha^2 (p^t)^2 = gamma^{ij} p_i p_j
// of geodesic-propagation.md:96-97 can be checked without a redundant sqrt.
constexpr double momentum_norm2(const InverseSpatialMetric &gu, double px,
                                double py, double pz) noexcept {
  return gu.xx * px * px + gu.yy * py * py + gu.zz * pz * pz +
         2.0 * (gu.xy * px * py + gu.xz * px * pz + gu.yz * py * pz);
}

// The null closure of [MCNX-GEO-02]: p^t = sqrt(gamma^{ij} p_i p_j)/alpha,
// exactly as written at geodesic-propagation.md:91. Pure — recomputed from
// the current (p_i, gamma^{ij}, alpha) at every call; p^t is never state.
inline double p_t_closure(double px, double py, double pz,
                          const InverseSpatialMetric &gu,
                          double alpha) noexcept {
  return std::sqrt(momentum_norm2(gu, px, py, pz)) / alpha;
}

// ---------------------------------------------------------------------------
// Spacetime metric from the ADM pieces
// ---------------------------------------------------------------------------

// Full g_munu, index order (t, x, y, z). Built from (alpha, beta^i, gamma_ij):
// g_tt = -alpha^2 + beta_k beta^k, g_ti = beta_i = gamma_ij beta^j,
// g_ij = gamma_ij.
struct SpacetimeMetric {
  double g[4][4];
};

constexpr SpacetimeMetric
spacetime_metric_from_adm(double alpha, const double beta_up[3],
                          const SpatialMetric &gamma) noexcept {
  const double gij[3][3] = {{gamma.xx, gamma.xy, gamma.xz},
                            {gamma.xy, gamma.yy, gamma.yz},
                            {gamma.xz, gamma.yz, gamma.zz}};
  double beta_lo[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      beta_lo[i] += gij[i][j] * beta_up[j];
  double beta2 = 0.0;
  for (int i = 0; i < 3; ++i)
    beta2 += beta_lo[i] * beta_up[i];

  SpacetimeMetric m{};
  m.g[0][0] = -alpha * alpha + beta2;
  for (int i = 0; i < 3; ++i) {
    m.g[0][i + 1] = beta_lo[i];
    m.g[i + 1][0] = beta_lo[i];
    for (int j = 0; j < 3; ++j)
      m.g[i + 1][j + 1] = gij[i][j];
  }
  return m;
}

// g_munu a^mu b^nu for two contravariant 4-vectors.
constexpr double metric_dot(const SpacetimeMetric &m, const double a[4],
                            const double b[4]) noexcept {
  double s = 0.0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu)
      s += m.g[mu][nu] * a[mu] * b[nu];
  return s;
}

// ---------------------------------------------------------------------------
// Orthonormal tetrad  [MCNX-PKT-04]
// ---------------------------------------------------------------------------

// The four tetrad legs e^mu_(a) as coordinate-frame contravariant 4-vectors:
// e[a][mu], tetrad index a (0 = the timelike leg u^mu), coordinate component
// mu in (t, x, y, z). A public POD by design: [MCNX-TRP-10] composes a
// caller-supplied polar axis from these legs and caches one Tetrad per
// scattering episode (caller-side caching is spec-permitted; this header
// never caches).
struct Tetrad {
  double e[4][4];
};

// Gram-Schmidt core, templated on the sqrt so the compile-time self-tests can
// run it with detail::csqrt; production entry point below uses std::sqrt.
// Leg order and projection order are FIXED: timelike leg first, then the
// coordinate basis vectors x, y, z in exactly that order
// (packet-representation-and-sampling.md:135-138 pins the order for
// golden-data determinism — any other order is physically fine but forbidden).
template <class SqrtF>
constexpr Tetrad build_tetrad_with(const SpacetimeMetric &m, const double u[4],
                                   SqrtF sqrt_fn) noexcept {
  // Normalization signs of the already-built legs: g(e_(0'), e_(0')) = -1
  // (u^mu is normalized by contract), g(e_(i'), e_(i')) = +1.
  const double sigma[4] = {-1.0, 1.0, 1.0, 1.0};

  Tetrad t{};
  for (int mu = 0; mu < 4; ++mu)
    t.e[0][mu] = u[mu]; // timelike leg: e^mu_(0') = u^mu, verbatim

  for (int a = 1; a < 4; ++a) { // coordinate basis x (a=1), y (a=2), z (a=3)
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    v[a] = 1.0;
    for (int b = 0; b < a; ++b) {
      const double d = metric_dot(m, v, t.e[b]);
      for (int mu = 0; mu < 4; ++mu)
        v[mu] -= sigma[b] * d * t.e[b][mu];
    }
    const double n = sqrt_fn(metric_dot(m, v, v));
    for (int mu = 0; mu < 4; ++mu)
      t.e[a][mu] = v[mu] / n;
  }
  return t;
}

// Production tetrad construction from the ADM pieces and the (normalized)
// fluid four-velocity. Pure: rebuilt from its arguments at every call.
inline Tetrad build_tetrad(double alpha, const double beta_up[3],
                           const SpatialMetric &gamma,
                           const double u[4]) noexcept {
  return build_tetrad_with(spacetime_metric_from_adm(alpha, beta_up, gamma), u,
                           [](double x) { return std::sqrt(x); });
}

// ---------------------------------------------------------------------------
// Isotropic fluid-frame direction draw + transform  [MCNX-PKT-04]
// ---------------------------------------------------------------------------

// Fluid-frame unit direction n_hat' (orthonormal-frame spatial components).
struct UnitVector3 {
  double x, y, z;
};

// Isotropic direction from two uniforms in [0, 1):
//   cos(theta) = 2 u1 - 1 in [-1, 1),  phi = 2 pi u2 in [0, 2 pi)
// (packet-representation-and-sampling.md:123-126). The uniforms are plain
// ARGUMENTS: this function performs no RNG call, so the per-consumer
// (e, k) draw maps stay auditable at the call sites.
inline UnitVector3 fluid_frame_direction(double u1, double u2) noexcept {
  const double costheta = 2.0 * u1 - 1.0;
  const double sintheta = std::sqrt(1.0 - costheta * costheta);
  const double phi = 2.0 * detail::pi * u2;
  return {sintheta * std::cos(phi), sintheta * std::sin(phi), costheta};
}

// Contravariant coordinate-frame 4-vector (t, x, y, z), the intermediate
// p^mu = e^mu_(lambda') p^{lambda'} — exposed so verification code can form
// full-metric contractions (null norm, nu recovery) without re-deriving it.
struct FourVectorUp {
  double v[4];
};

// p^mu = e^mu_(lambda') p^{lambda'} with the fluid-frame null momentum
// p^{lambda'} = nu (1, n_hat') of [MCNX-PKT-04]. Accepts an arbitrary
// fluid-frame unit direction (not just a (u1, u2) pair) so [MCNX-TRP-10]'s
// correlated-axis draw can reuse it.
constexpr FourVectorUp
fluid_momentum_to_coordinate(const Tetrad &t, double nu,
                             const UnitVector3 &nhat) noexcept {
  const double pf[4] = {nu, nu * nhat.x, nu * nhat.y, nu * nhat.z};
  FourVectorUp p{};
  for (int mu = 0; mu < 4; ++mu)
    for (int a = 0; a < 4; ++a)
      p.v[mu] += t.e[a][mu] * pf[a];
  return p;
}

// The transform's packet-facing output: p^t (contravariant) and the lower
// spatial momentum p_i that the packet state stores ([MCNX-PKT-01]).
struct CoordinateMomentum {
  double pt;           // p^t
  double px, py, pz;   // p_i (LOWER spatial components)
};

// p^t = p^t,  p_i = g_{i mu} p^mu.
constexpr CoordinateMomentum
lower_spatial_momentum(const SpacetimeMetric &m,
                       const FourVectorUp &p) noexcept {
  CoordinateMomentum r{p.v[0], 0.0, 0.0, 0.0};
  double lo[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int mu = 0; mu < 4; ++mu)
      lo[i] += m.g[i + 1][mu] * p.v[mu];
  r.px = lo[0];
  r.py = lo[1];
  r.pz = lo[2];
  return r;
}

// The full [MCNX-PKT-04] transform, exactly as written at
// packet-representation-and-sampling.md:131-133:
//   p^t = e^t_(lambda') p^{lambda'} ,   p_i = g_{i mu} e^mu_(lambda') p^{lambda'}
constexpr CoordinateMomentum
transform_to_coordinate_frame(const Tetrad &t, const SpacetimeMetric &m,
                              double nu, const UnitVector3 &nhat) noexcept {
  return lower_spatial_momentum(m, fluid_momentum_to_coordinate(t, nu, nhat));
}

// ---------------------------------------------------------------------------
// Compile-time verification (the structural/flat-limit half of T13; the
// measured-error rows live in the runtime `unit-selftest` battery)
// ---------------------------------------------------------------------------

namespace detail {

// Absolute machine-tier closeness for O(1) quantities whose target may be 0
// (a relative test is undefined there; approx_eq's 1e-30 floor is far below
// the ~1e-16 rounding of these contractions).
constexpr bool near_machine_abs(double got, double want) noexcept {
  return cabs(got - want) <= rtol_machine;
}

// The constant SPD test metric (off-diagonal terms so a diagonal-only inverse
// cannot pass), shared with nothing — a pure test fixture.
inline constexpr SpatialMetric test_gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};

// gamma^{ij} gamma_jk = delta^i_k on the constant test metric, machine tier.
constexpr bool inverse_metric_identity_holds() noexcept {
  const InverseSpatialMetric gu = spatial_metric_inverse(test_gamma);
  if (!(gu.det > 0.0))
    return false;
  const double lo[3][3] = {{test_gamma.xx, test_gamma.xy, test_gamma.xz},
                           {test_gamma.xy, test_gamma.yy, test_gamma.yz},
                           {test_gamma.xz, test_gamma.yz, test_gamma.zz}};
  const double up[3][3] = {
      {gu.xx, gu.xy, gu.xz}, {gu.xy, gu.yy, gu.yz}, {gu.xz, gu.yz, gu.zz}};
  for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 3; ++k) {
      double s = 0.0;
      for (int j = 0; j < 3; ++j)
        s += up[i][j] * lo[j][k];
      if (!near_machine_abs(s, i == k ? 1.0 : 0.0))
        return false;
    }
  return true;
}

// Minkowski from the ADM pieces is exactly diag(-1, 1, 1, 1).
constexpr bool flat_adm_metric_holds() noexcept {
  const double beta0[3] = {0.0, 0.0, 0.0};
  const SpacetimeMetric m =
      spacetime_metric_from_adm(1.0, beta0, SpatialMetric{1, 0, 0, 1, 0, 1});
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) {
      const double want = mu != nu ? 0.0 : (mu == 0 ? -1.0 : 1.0);
      if (m.g[mu][nu] != want)
        return false;
    }
  return true;
}

// Flat-static limit (packet-representation-and-sampling.md:187-190): with
// Minkowski and u^mu = (1, 0, 0, 0) the tetrad is the identity, exactly.
constexpr bool flat_static_tetrad_is_identity() noexcept {
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
  const SpacetimeMetric m =
      spacetime_metric_from_adm(1.0, beta0, SpatialMetric{1, 0, 0, 1, 0, 1});
  const Tetrad t =
      build_tetrad_with(m, ustat, [](double x) { return csqrt(x); });
  for (int a = 0; a < 4; ++a)
    for (int mu = 0; mu < 4; ++mu)
      if (t.e[a][mu] != (a == mu ? 1.0 : 0.0))
        return false;
  return true;
}

// Flat-static transform: p^t = nu and p_i = nu n_hat_i, exactly (the
// transform is linear; no sqrt is exercised on this path).
constexpr bool flat_static_momentum_holds() noexcept {
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
  const SpacetimeMetric m =
      spacetime_metric_from_adm(1.0, beta0, SpatialMetric{1, 0, 0, 1, 0, 1});
  const Tetrad t =
      build_tetrad_with(m, ustat, [](double x) { return csqrt(x); });
  const double nu = 2.5;
  const UnitVector3 nhat{0.6, 0.0, 0.8};
  const CoordinateMomentum p = transform_to_coordinate_frame(t, m, nu, nhat);
  return p.pt == nu && p.px == nu * nhat.x && p.py == nu * nhat.y &&
         p.pz == nu * nhat.z;
}

// Null-closure identity on the flat limit, sqrt-free: alpha^2 (p^t)^2 =
// gamma^{ij} p_i p_j with p^t = nu on a unit direction (exact components).
constexpr bool flat_null_closure_holds() noexcept {
  const InverseSpatialMetric gu =
      spatial_metric_inverse(SpatialMetric{1, 0, 0, 1, 0, 1});
  const double nu = 2.5;
  // (3/5, 0, 4/5) is unit only up to rounding in binary64; compare against
  // the recomputed sum rather than nu^2 to keep the check exact-tier.
  const double px = nu * 0.6, py = 0.0, pz = nu * 0.8;
  const double want = px * px + py * py + pz * pz;
  return momentum_norm2(gu, px, py, pz) == want;
}

} // namespace detail

static_assert(detail::inverse_metric_identity_holds(),
              "gamma^{ij} gamma_jk != delta^i_k on the constant test metric "
              "([MCNX-GEO-02] analytic 3x3 inverse)");
static_assert(detail::flat_adm_metric_holds(),
              "spacetime_metric_from_adm must reproduce Minkowski exactly on "
              "alpha = 1, beta = 0, gamma_ij = delta_ij");
static_assert(detail::flat_static_tetrad_is_identity(),
              "flat-static tetrad must be the identity "
              "(packet-representation-and-sampling.md:187-190)");
static_assert(detail::flat_static_momentum_holds(),
              "flat-static transform must give p^t = nu, p_i = nu n_hat_i "
              "([MCNX-PKT-04] flat limit)");
static_assert(detail::flat_null_closure_holds(),
              "flat-limit null closure alpha^2 (p^t)^2 = gamma^{ij} p_i p_j "
              "failed ([MCNX-GEO-02])");
static_assert(detail::csqrt(1.0) == 1.0 && detail::csqrt(4.0) == 2.0,
              "the compile-time test sqrt must be exact on exact squares");

} // namespace MCNuX

#endif // MCNUX_TETRAD_HXX
