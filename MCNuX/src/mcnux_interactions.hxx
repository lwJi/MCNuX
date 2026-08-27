#ifndef MCNUX_INTERACTIONS_HXX
#define MCNUX_INTERACTIONS_HXX

// Interaction-time sampling and sampling-episode competition — the pure
// functions of specs/neutrino-matter-interactions.md.
//
// Governing spec requirements:
//   * [MCNX-INT-01] — the interaction-time draw (provenance arXiv:1708.08452
//     Eq. 24, restated in the corpus's [0, 1) convention):
//         Dt_{s,a} = -ln(1 - u_{s,a}) p^t/(kappa_{s,a} nu)
//     with kappa the FLUID-FRAME opacity in geometrized code units, nu the
//     fluid-frame energy nu = -p_mu u^mu, and p^t/nu a dimensionless ratio
//     (both geometrized) converting the fluid-frame rate kappa*nu to
//     coordinate time. A zero-opacity channel is Dt = +inf, realized by an
//     explicit branch BEFORE the division: u(S, q, e, k) == 0.0 is a
//     reachable exact RNG output (mcnux_rng.hxx, the u64 -> [0, 1) closed
//     forms), so the naive formula would evaluate 0/0 = NaN there — the
//     branch IS the spec's "no NaN" requirement, not defensiveness.
//   * [MCNX-INT-02] — episode competition: the candidate event is the
//     shorter of Dt_s and Dt_a; it fires only if it STRICTLY precedes both
//     the packet's cell-exit time and the step remainder Dt_rem. Ties
//     (Dt_s == Dt_a bitwise, including both +inf) resolve to SCATTERING —
//     implemented structurally as `dt_s <= dt_a -> scatter`, so the
//     tie-break cannot drift. If no event fires: cell exit if it precedes
//     the step end, else step end. The Dt_cell_exit == Dt_rem tie resolving
//     to step end is implementation freedom (the spec pins only the episode
//     structure and draw order); golden data depends on this choice. Packet
//     weight N is NOT an input — competition never touches it (Discreteness:
//     no path-length attenuation anywhere in the baseline).
//   * [MCNX-INT-06] — RNG consumption map. All draws are evaluations of the
//     pure u(S, q, e, k); within an episode
//         k = 0 -> u_s (scattering time),  k = 1 -> u_a (absorption time),
//         k = 2 -> u1 (cos theta),         k = 3 -> u2 (phi)
//     in exactly this order and no other draws. The event counter e
//     increments by 1 at each episode START (creation is e = 0); every
//     episode consumes both k = 0, 1 unconditionally; k = 2, 3 are consumed
//     only if the episode ends in scattering; episodes ending without an
//     event (cell crossing or step end) consume exactly k = 0, 1, and the
//     next episode uses the next e. The wrappers below carry that (e, k)
//     bookkeeping; per the repo convention (mcnux_tetrad.hxx: RNG
//     discipline) the math functions themselves take already-drawn doubles.
//
// Coefficient provenance ([MCNX-INT-05]): the opacities enter as
// already-evaluated plain doubles (kappa_a, kappa_s) from the coefficient
// interface of specs/opacity-eos-evaluation.md (callers use
// evaluate_coefficients / RangedTableCoefficients); when the trapped-regime
// relabeling of specs/trapped-regime-treatment.md is active, the primed
// kappa_a', kappa_s' substitute TRANSPARENTLY at the call site — nothing in
// this header knows or cares which it received. This header therefore
// includes neither mcnux_coefficients.hxx nor mcnux_trp.hxx.
//
// Spec-normative deviations, preserved verbatim (neutrino-matter-
// interactions.md, Open questions — spec text normative per specs/README.md):
//   * Per-cell episode structure: unlike the source paper's
//     start-of-step draw, a new episode (with fresh draws) begins at every
//     cell entry; statistically exact for piecewise-constant coefficients by
//     memorylessness, and the undrawn remainder is discarded on crossing.
//   * Frozen episode kinematics: p^t, nu, and the opacities are evaluated
//     ONCE at episode start and held fixed for the episode, even on curved
//     metrics where p^t drifts within the cell — re-evaluating mid-episode
//     would be a spec change (it alters golden data), not an implementation
//     choice.
//   * Ties resolve to scattering: probability ~2^-53 per episode; the rule
//     exists so the branch is deterministic, not because it is
//     statistically consequential.
//
// Units at the call boundary: kappa arrives in code units via
//   kappa[code] = kappa[cm^-1] * length_code_to_cgs
// (conventions-and-units.md:88; opacity transforms INVERSELY to a length, so
// opacity_code_to_cgs is the wrong direction here) — opacity_cgs_to_code
// below wraps exactly that product. MeV is for the coefficient lookup only:
// nu[MeV] = nu[code] * energy_code_to_MeV (fluid_frame_energy_MeV below);
// the Dt ratio itself stays entirely in code units.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// C++ (constexpr where the math is rational; <cmath> log otherwise, which
// nvcc supports in device code) so host code, device code, and standalone
// unit tests can all include it. The static_asserts at the bottom are the
// compile-time half of the T18a-core test surface; the runtime half lives in
// the `unit-selftest` battery (mcnux_selftest_interactions.cxx, rows
// 106..113, and mcnux_selftest_cellexit.cxx, rows 140..145).

#include "mcnux_rng.hxx"    // u(S, q, e, k)
#include "mcnux_tetrad.hxx" // MCNUX_HOST_DEVICE, metric helpers, csqrt
#include "mcnux_units.hxx"  // conversion factors, approx_eq, rtol_machine

#include <cmath>
#include <cstdint>
#include <limits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Unit wrappers at the interaction call boundary  [MCNX-INT-01]
// ---------------------------------------------------------------------------

// kappa[code] = kappa[cm^-1] * length_code_to_cgs — the binding direction of
// conventions-and-units.md:88. Reuses the existing factor; no new literal.
constexpr double opacity_cgs_to_code(double kappa_cgs) noexcept {
  return kappa_cgs * length_code_to_cgs;
}

// nu[MeV] = nu[code] * energy_code_to_MeV — for the coefficient LOOKUP only
// ([MCNX-INT-01]: "converted to MeV for the coefficient lookup"); the draw
// law itself consumes nu in code units.
constexpr double fluid_frame_energy_MeV(double nu_code) noexcept {
  return nu_code * energy_code_to_MeV;
}

// ---------------------------------------------------------------------------
// Interaction-time draw  [MCNX-INT-01]
// ---------------------------------------------------------------------------

// Dt = -ln(1 - u) p^t/(kappa nu), all of (p^t, kappa, nu) in geometrized code
// units, u in [0, 1). The kappa == 0.0 branch (the zero-opacity channel:
// nu_x absorption in table mode, the transparency floor) returns +inf BEFORE
// the division: at u == 0 the naive formula is 0 * inf = NaN, and the spec
// pins "never fires, draw still consumed, no NaN". u == 0 with kappa > 0
// gives exactly Dt = 0 (-ln(1) == -0.0; the event fires immediately).
// Non-constexpr by design (std::log; the p_t_closure precedent) — covered by
// the runtime selftest rows, not the static_asserts below.
MCNUX_HOST_DEVICE inline double interaction_time(double u_draw, double p_t_up,
                                                 double kappa_code,
                                                 double nu_code) noexcept {
  if (kappa_code == 0.0)
    return std::numeric_limits<double>::infinity();
  return -std::log(1.0 - u_draw) * p_t_up / (kappa_code * nu_code);
}

// ---------------------------------------------------------------------------
// Fluid-frame energy  nu = -p_mu u^mu  [MCNX-INT-01]
// ---------------------------------------------------------------------------

// ADM shortcut, exact algebraically (machine tier numerically):
//   p_t = -alpha^2 p^t + beta^i p_i
//   nu  = -(p_t u^t + p_x u^x + p_y u^y + p_z u^z)
// The stored packet momentum p_i is COVARIANT and contracts directly against
// the contravariant u^i — no index raising. p^t is the null closure of
// [MCNX-GEO-02] (p_t_closure, mcnux_tetrad.hxx), passed in per the frozen
// episode kinematics above. Cross-checked in the selftest battery against
// the independent metric_dot route (raise p_i, contract with the full
// g_munu) on flat and pinned curved fixtures.
constexpr double fluid_frame_energy(double alpha, const double beta_up[3],
                                    const double p_cov[3], double p_t_up,
                                    const double u_up[4]) noexcept {
  const double p_t = -alpha * alpha * p_t_up + beta_up[0] * p_cov[0] +
                     beta_up[1] * p_cov[1] + beta_up[2] * p_cov[2];
  return -(p_t * u_up[0] + p_cov[0] * u_up[1] + p_cov[1] * u_up[2] +
           p_cov[2] * u_up[3]);
}

// ---------------------------------------------------------------------------
// Episode competition  [MCNX-INT-02]
// ---------------------------------------------------------------------------

// How a sampling episode ends. Scatter/Absorb: the channel fired (a new
// episode follows a scattering; an absorption removes the packet,
// [MCNX-INT-03]). CellExit: the boundary preceded the candidate — the packet
// crosses and a NEW episode begins (undrawn remainder discarded). StepEnd:
// the packet finishes the step with no event.
enum class EpisodeEnd : int { Scatter = 0, Absorb = 1, CellExit = 2, StepEnd = 3 };

struct EpisodeOutcome {
  EpisodeEnd kind;
  double dt_event; // coordinate-time interval from episode start, code units
};

// The competition of [MCNX-INT-02]. All four inputs are plain coordinate-time
// intervals measured from the episode start (+inf is legal for any of them;
// dt_s/dt_a are interaction_time outputs, dt_cell_exit and dt_rem arrive as
// synthetic/runtime inputs from the push). The candidate is chosen by
// `dt_s <= dt_a -> scatter`, which makes the bitwise-tie -> scattering rule
// structural; it fires only if it STRICTLY precedes both bounds.
constexpr EpisodeOutcome compete_episode(double dt_s, double dt_a,
                                         double dt_cell_exit,
                                         double dt_rem) noexcept {
  const bool scatter_wins = dt_s <= dt_a; // tie -> scattering, [MCNX-INT-02]
  const double dt_cand = scatter_wins ? dt_s : dt_a;
  if (dt_cand < dt_cell_exit && dt_cand < dt_rem)
    return {scatter_wins ? EpisodeEnd::Scatter : EpisodeEnd::Absorb, dt_cand};
  if (dt_cell_exit < dt_rem)
    return {EpisodeEnd::CellExit, dt_cell_exit};
  return {EpisodeEnd::StepEnd, dt_rem};
}

// ---------------------------------------------------------------------------
// Cell-exit time  (feeds compete_episode's dt_cell_exit, [MCNX-INT-02])
// ---------------------------------------------------------------------------
// Coordinate time from the current position, along the FROZEN coordinate-
// frame velocity dx^i/dt of the episode-start geodesic_rhs evaluation
// ([MCNX-GEO-01]: dxdt[i] = p^i/p^t - beta^i), to the boundary of the
// containing cell. NO SPEC PINS THIS ALGORITHM — cell-exit-time tracking is
// implementation freedom (neutrino-matter-interactions.md, Implementation
// freedom); only the [MCNX-INT-02] competition it feeds is normative. The
// frozen-kinematics assumption above applies: the caller evaluates
// geodesic_rhs ONCE at episode start and passes dxdt here — this helper
// never re-derives velocity or p^t, and takes no gather functor.
//
// Boundary convention (deliberate choice, recorded per the decision-recording
// idiom of mcnux_fluid.hxx): cell membership is HALF-OPEN, [i dx, (i+1) dx)
// per axis — i0 = floor((x - prob_lo)/dx) is the containing cell, the
// identical one-line anchor arithmetic as CellFluidGather::containing_cell
// (mcnux_fluid.hxx) and the emission creation cell (cell_lo = prob_lo + i dx,
// mcnux_emission.cxx). Consequence: a packet exactly on its lower face moving
// in -x yields dt = 0 (a signed zero, == 0.0) — an immediate CellExit in
// compete_episode, since finite candidates fire only if STRICTLY smaller.
//
// Per axis, the candidate time is (face - x)/dxdt toward the faced boundary;
// a zero velocity component contributes +inf via an explicit branch BEFORE
// any division (the kappa == 0.0 -> +inf idiom of interaction_time above:
// 0/0 must never reach the FPU), so the result is never NaN and never
// negative for in-cell positions; an all-zero velocity returns +inf, a legal
// compete_episode input (the channel/boundary sentinels already compete).

// Core on explicit faces: constexpr (rational arithmetic only), the
// p_t_closure_with-style testable layer. face_lo/face_hi are the containing
// cell's boundary coordinates per axis.
constexpr double dt_to_cell_exit_faces(const double x[3], const double dxdt[3],
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

// Production wrapper on the level geometry: derives the containing cell's
// faces from (prob_lo, dx) via the shared floor anchor (non-constexpr because
// of std::floor — the interaction_time/std::log precedent; covered by the
// runtime selftest rows, not the static_asserts below).
MCNUX_HOST_DEVICE inline double dt_to_cell_exit(const double x[3],
                                                const double dxdt[3],
                                                const double prob_lo[3],
                                                const double dx[3]) noexcept {
  using std::floor;
  double face_lo[3], face_hi[3];
  for (int d = 0; d < 3; ++d) {
    // i0 = floor((x - prob_lo)/dx): the CellFluidGather::containing_cell
    // anchor (mcnux_fluid.hxx), kept as a double so face_lo is computed in
    // one rounding step per term.
    const double i0 = floor((x[d] - prob_lo[d]) / dx[d]);
    face_lo[d] = prob_lo[d] + i0 * dx[d];
    face_hi[d] = face_lo[d] + dx[d];
  }
  return dt_to_cell_exit_faces(x, dxdt, face_lo, face_hi);
}

// ---------------------------------------------------------------------------
// RNG draw map  [MCNX-INT-06]
// ---------------------------------------------------------------------------

// The pinned draw indices within an episode, in exactly this order.
inline constexpr std::uint32_t draw_k_scatter_time = 0u; // u_s
inline constexpr std::uint32_t draw_k_absorb_time = 1u;  // u_a
inline constexpr std::uint32_t draw_k_cos_theta = 2u;    // u1 (scatter only)
inline constexpr std::uint32_t draw_k_phi = 3u;          // u2 (scatter only)

// The two channel-time uniforms, drawn unconditionally at every episode
// start (both are consumed regardless of the coefficient values). e is a
// plain argument: the caller owns the counter and increments it by 1 at each
// episode start (creation is e = 0).
struct EpisodeUniforms {
  double u_s, u_a;
};

constexpr EpisodeUniforms episode_uniforms(std::uint64_t S, std::uint64_t q,
                                           std::uint32_t e) noexcept {
  return {u(S, q, e, draw_k_scatter_time), u(S, q, e, draw_k_absorb_time)};
}

// The two direction uniforms of the elastic redraw ([MCNX-INT-04] via
// fluid_frame_direction(u1, u2)), consumed ONLY if the episode ends in
// scattering. Eventless episodes consume exactly k = 0, 1 and never call
// this.
struct ScatterUniforms {
  double u1, u2;
};

constexpr ScatterUniforms scatter_uniforms(std::uint64_t S, std::uint64_t q,
                                           std::uint32_t e) noexcept {
  return {u(S, q, e, draw_k_cos_theta), u(S, q, e, draw_k_phi)};
}

// ---------------------------------------------------------------------------
// Elastic scatter redraw  [MCNX-INT-04]
// ---------------------------------------------------------------------------

// The outgoing coordinate-frame momentum of an elastic scattering event: a
// NEW isotropic fluid-frame direction from the two already-drawn uniforms
// (u1, u2) — k = 2, 3 of the episode draw map above, via scatter_uniforms —
// at the FIXED fluid-frame energy nu, through the SAME pinned tetrad
// transform as packet creation ([MCNX-INT-04]: "the same pinned tetrad
// transform"; fluid_frame_direction then transform_to_coordinate_frame of
// mcnux_tetrad.hxx, the exact chain of the emission creation path). The
// tetrad, metric, and nu are the FROZEN episode-start kinematics (spec
// deviation above): the caller evaluates them once per episode and passes
// them here — this function re-derives nothing. Non-constexpr
// (fluid_frame_direction's sqrt/sin/cos); the explicit-direction composition
// is asserted below, the drawn path by the runtime selftest rows 146..149.
MCNUX_HOST_DEVICE inline CoordinateMomentum
scatter_redraw(const Tetrad &t, const SpacetimeMetric &m, double nu, double u1,
               double u2) noexcept {
  return transform_to_coordinate_frame(t, m, nu, fluid_frame_direction(u1, u2));
}

// ---------------------------------------------------------------------------
// Compile-time verification (the constexpr half of the T18a-core test
// surface; the std::log draw path is runtime-only and lives in the
// `unit-selftest` battery rows 106..113)
// ---------------------------------------------------------------------------

namespace detail {

inline constexpr double inf = std::numeric_limits<double>::infinity();

// Flat fixture: Minkowski, static fluid u = (1, 0, 0, 0) => nu == p^t
// exactly (p_t = -p^t with alpha = 1, beta = 0; no rounding on this path).
constexpr bool flat_fluid_frame_energy_holds() noexcept {
  const double beta0[3] = {0.0, 0.0, 0.0};
  const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
  const double p_cov[3] = {0.3, -0.4, 0.12};
  const double p_t_up = 2.5; // arbitrary: the flat-static identity is linear
  return fluid_frame_energy(1.0, beta0, p_cov, p_t_up, ustat) == p_t_up;
}

// Pinned curved fixture (nonzero shift, non-diagonal gamma): the ADM
// shortcut must agree with the independent metric_dot route
// nu = -g_munu p^mu u^nu at machine tier. Entirely constexpr via csqrt.
constexpr bool curved_fluid_frame_energy_holds() noexcept {
  const double alpha = 1.1;
  const double beta_up[3] = {0.1, -0.05, 0.2};
  const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
  const double p_cov[3] = {0.4, -0.25, 0.3};

  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
  const double p_t_up = p_t_closure_with(p_cov[0], p_cov[1], p_cov[2], gu,
                                         alpha, [](double v) { return csqrt(v); });

  // Normalized fluid four-velocity from an arbitrary timelike seed.
  const SpacetimeMetric m = spacetime_metric_from_adm(alpha, beta_up, gamma);
  const double w[4] = {1.0, 0.2, -0.1, 0.15};
  const double n = csqrt(-metric_dot(m, w, w));
  const double u_up[4] = {w[0] / n, w[1] / n, w[2] / n, w[3] / n};

  // Independent route: raise p_i (p^i = gamma^{ij} p_j - beta^i p^t) and
  // contract two CONTRAVARIANT vectors against the full metric.
  const double p_up[4] = {
      p_t_up,
      gu.xx * p_cov[0] + gu.xy * p_cov[1] + gu.xz * p_cov[2] -
          beta_up[0] * p_t_up,
      gu.xy * p_cov[0] + gu.yy * p_cov[1] + gu.yz * p_cov[2] -
          beta_up[1] * p_t_up,
      gu.xz * p_cov[0] + gu.yz * p_cov[1] + gu.zz * p_cov[2] -
          beta_up[2] * p_t_up};

  const double nu_shortcut =
      fluid_frame_energy(alpha, beta_up, p_cov, p_t_up, u_up);
  const double nu_metric = -metric_dot(m, p_up, u_up);
  return near_machine_abs(nu_shortcut, nu_metric) &&
         near_machine_abs(metric_dot(m, u_up, u_up), -1.0) && nu_shortcut > 0.0;
}

// The full competition case table of [MCNX-INT-02], including the strict
// precedence at equality with a bound and the +inf sentinel competitions.
constexpr bool competition_table_holds() noexcept {
  // Event fires: shorter channel wins.
  {
    const EpisodeOutcome o = compete_episode(1.0, 2.0, 3.0, 4.0);
    if (!(o.kind == EpisodeEnd::Scatter && o.dt_event == 1.0))
      return false;
  }
  {
    const EpisodeOutcome o = compete_episode(2.0, 1.0, 3.0, 4.0);
    if (!(o.kind == EpisodeEnd::Absorb && o.dt_event == 1.0))
      return false;
  }
  // Bitwise tie -> scattering (finite and both-+inf-under-finite-bounds).
  {
    const EpisodeOutcome o = compete_episode(1.5, 1.5, 10.0, 20.0);
    if (!(o.kind == EpisodeEnd::Scatter && o.dt_event == 1.5))
      return false;
  }
  // Cell boundary first: cross, new episode (remainder discarded upstream).
  {
    const EpisodeOutcome o = compete_episode(3.0, 4.0, 1.0, 2.0);
    if (!(o.kind == EpisodeEnd::CellExit && o.dt_event == 1.0))
      return false;
  }
  // Step end first: finish the step with no event.
  {
    const EpisodeOutcome o = compete_episode(3.0, 4.0, 2.0, 1.0);
    if (!(o.kind == EpisodeEnd::StepEnd && o.dt_event == 1.0))
      return false;
  }
  // Strict precedence: candidate == cell exit does NOT fire.
  {
    const EpisodeOutcome o = compete_episode(2.0, 3.0, 2.0, 4.0);
    if (!(o.kind == EpisodeEnd::CellExit && o.dt_event == 2.0))
      return false;
  }
  // Zero-opacity sentinels: +inf channels never fire.
  {
    const EpisodeOutcome o = compete_episode(inf, inf, 2.0, 3.0);
    if (!(o.kind == EpisodeEnd::CellExit && o.dt_event == 2.0))
      return false;
  }
  {
    const EpisodeOutcome o = compete_episode(inf, 5.0, 7.0, 6.0);
    if (!(o.kind == EpisodeEnd::Absorb && o.dt_event == 5.0))
      return false;
  }
  {
    const EpisodeOutcome o = compete_episode(inf, inf, inf, 3.0);
    if (!(o.kind == EpisodeEnd::StepEnd && o.dt_event == 3.0))
      return false;
  }
  return true;
}

// Exact axis-aligned flat cases of the cell-exit core, on faces with exact
// binary coordinates, plus the sentinel and composed-competition checks.
constexpr bool cell_exit_faces_hold() noexcept {
  const double face_lo[3] = {0.0, 0.0, 0.0};
  const double face_hi[3] = {1.0, 2.0, 4.0};
  const double x[3] = {0.25, 1.0, 2.0};
  // +x toward the upper face: dt = (1 - 0.25)/0.5 = 1.5 exactly.
  {
    const double v[3] = {0.5, 0.0, 0.0};
    if (dt_to_cell_exit_faces(x, v, face_lo, face_hi) != 1.5)
      return false;
  }
  // -y toward the lower face: dt = (0 - 1)/(-0.5) = 2 exactly.
  {
    const double v[3] = {0.0, -0.5, 0.0};
    if (dt_to_cell_exit_faces(x, v, face_lo, face_hi) != 2.0)
      return false;
  }
  // Mixed velocity: per-axis candidates 3 (+x), 4 (-y), 2 (+z); min picks z.
  {
    const double v[3] = {0.25, -0.25, 1.0};
    if (dt_to_cell_exit_faces(x, v, face_lo, face_hi) != 2.0)
      return false;
  }
  // A zero component contributes +inf; the min still picks the finite axis.
  {
    const double v[3] = {0.0, 0.0, 1.0};
    if (dt_to_cell_exit_faces(x, v, face_lo, face_hi) != 2.0)
      return false;
  }
  // All-zero velocity: +inf (legal compete_episode input), never NaN.
  {
    const double v[3] = {0.0, 0.0, 0.0};
    if (dt_to_cell_exit_faces(x, v, face_lo, face_hi) != inf)
      return false;
  }
  // Composed with [MCNX-INT-02]: the produced dt_cell_exit precedes both
  // channels and dt_rem -> CellExit at exactly that time.
  {
    const double v[3] = {0.5, 0.0, 0.0};
    const EpisodeOutcome o = compete_episode(
        5.0, 6.0, dt_to_cell_exit_faces(x, v, face_lo, face_hi), 10.0);
    if (!(o.kind == EpisodeEnd::CellExit && o.dt_event == 1.5))
      return false;
  }
  return true;
}

// [MCNX-INT-04] the redraw composition preserves nu: on the pinned curved
// fixture (the curved_fluid_frame_energy_holds inputs), the fluid-frame
// energy recomputed from the transformed momentum — with the SAME u^mu and
// metric — equals the input nu at machine tier, and the transformed p^t
// satisfies the [MCNX-GEO-02] null closure. Constexpr via csqrt and an
// EXPLICIT exact unit direction (the drawn fluid_frame_direction path is
// runtime-only; scatter_redraw is exactly this composition preceded by it).
constexpr bool scatter_composition_preserves_nu() noexcept {
  const double alpha = 1.1;
  const double beta_up[3] = {0.1, -0.05, 0.2};
  const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
  const SpacetimeMetric m = spacetime_metric_from_adm(alpha, beta_up, gamma);

  // Normalized fluid four-velocity from an arbitrary timelike seed.
  const double w[4] = {1.0, 0.2, -0.1, 0.15};
  const double n = csqrt(-metric_dot(m, w, w));
  const double u_up[4] = {w[0] / n, w[1] / n, w[2] / n, w[3] / n};

  const Tetrad t =
      build_tetrad_with(m, u_up, [](double v) { return csqrt(v); });
  const double nu = 3.25;
  const UnitVector3 nhat{0.6, 0.0, 0.8}; // exactly unit: 0.36 + 0.64 = 1
  const CoordinateMomentum cm = transform_to_coordinate_frame(t, m, nu, nhat);

  const double p_cov[3] = {cm.px, cm.py, cm.pz};
  const double nu_back = fluid_frame_energy(alpha, beta_up, p_cov, cm.pt, u_up);

  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
  const double pt_closed = p_t_closure_with(cm.px, cm.py, cm.pz, gu, alpha,
                                            [](double v) { return csqrt(v); });
  return approx_eq(nu_back, nu, rtol_machine) &&
         approx_eq(pt_closed, cm.pt, rtol_machine);
}

// The wrappers are exactly the pinned direct u(S, q, e, k) calls.
constexpr bool draw_map_wrappers_hold() noexcept {
  const std::uint64_t S = 1296518744ull;
  const std::uint64_t q = 0x0123456789ABCDEFull;
  for (std::uint32_t e = 0u; e < 3u; ++e) {
    const EpisodeUniforms eu = episode_uniforms(S, q, e);
    if (!(eu.u_s == u(S, q, e, 0u) && eu.u_a == u(S, q, e, 1u)))
      return false;
    const ScatterUniforms su = scatter_uniforms(S, q, e);
    if (!(su.u1 == u(S, q, e, 2u) && su.u2 == u(S, q, e, 3u)))
      return false;
  }
  return true;
}

} // namespace detail

static_assert(detail::flat_fluid_frame_energy_holds(),
              "[MCNX-INT-01] flat-static fluid-frame energy must be exactly "
              "p^t (nu = -p_mu u^mu with u = (1,0,0,0) on Minkowski)");
static_assert(detail::curved_fluid_frame_energy_holds(),
              "[MCNX-INT-01] the ADM p_t shortcut must agree with the "
              "metric_dot route nu = -g_munu p^mu u^nu on the pinned curved "
              "fixture (machine tier)");
static_assert(detail::competition_table_holds(),
              "[MCNX-INT-02] episode competition: min rule, strict precedence "
              "over cell-exit/step-end, tie -> scattering, +inf sentinels");
static_assert(detail::cell_exit_faces_hold(),
              "cell-exit core: exact axis-aligned flat cases, +inf on zero "
              "velocity components (never NaN), min over axes, and CellExit "
              "composition with compete_episode ([MCNX-INT-02])");
static_assert(detail::scatter_composition_preserves_nu(),
              "[MCNX-INT-04] the elastic-redraw composition (tetrad "
              "transform at fixed nu) must preserve the fluid-frame energy "
              "and satisfy the null closure on the pinned curved fixture "
              "(machine tier)");
static_assert(detail::draw_map_wrappers_hold(),
              "[MCNX-INT-06] episode_uniforms/scatter_uniforms must be "
              "exactly the pinned u(S, q, e, k) evaluations");
static_assert(draw_k_scatter_time == 0u && draw_k_absorb_time == 1u &&
                  draw_k_cos_theta == 2u && draw_k_phi == 3u,
              "[MCNX-INT-06] draw-index assignment is pinned: k = 0 u_s, "
              "k = 1 u_a, k = 2 u1, k = 3 u2");
static_assert(opacity_cgs_to_code(1.0) == length_code_to_cgs,
              "[MCNX-INT-01] kappa[code] = kappa[cm^-1] * length_code_to_cgs "
              "(opacity transforms inversely to a length; the code->cgs "
              "factor is the WRONG direction here)");
static_assert(detail::approx_eq(opacity_cgs_to_code(3.25) *
                                    opacity_code_to_cgs,
                                3.25, detail::rtol_machine),
              "opacity cgs -> code -> cgs must round-trip at machine tier");
static_assert(fluid_frame_energy_MeV(1.0) == energy_code_to_MeV,
              "nu[MeV] = nu[code] * energy_code_to_MeV (lookup-only "
              "conversion)");

} // namespace MCNuX

#endif // MCNUX_INTERACTIONS_HXX
