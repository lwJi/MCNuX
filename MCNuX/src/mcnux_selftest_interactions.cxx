// ---------------------------------------------------------------------------
// Rows 106..113 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the interaction-time draw, episode competition, and
// RNG draw map of specs/neutrino-matter-interactions.md
// [MCNX-INT-01/02/06] (T18a-core), through the mcnux_interactions.hxx pure
// functions. Deterministic pinned fixtures only (the `unit-selftest` parfile
// is a single-shot itlast=0 run): the std::log draw path — which the header's
// static_asserts cannot cover — is exercised against independently composed
// closed forms at the machine tier; the sentinel, tie, competition-ordering,
// and draw-map rows are exact/boolean; the nu = -p_mu u^mu ADM shortcut is
// cross-checked against the independent metric_dot route on flat and pinned
// curved fixtures. No scheduled transport, no containers, no live tables.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_interactions.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <cmath>
#include <cstdint>
#include <limits>

namespace MCNuX {

namespace {

// Pinned draw-map fixture identity (the pinned production seed of
// mcnux_stats.hxx is not required here: any pinned (S, q) works and the
// direct-u comparisons are exact either way).
constexpr std::uint64_t fix_S = 1296518744ull;
constexpr std::uint64_t fix_q = 42ull;

// The independent flat/curved nu route: raise p_i against gamma^{ij}, build
// p^mu, contract two contravariant vectors with the full ADM metric.
double nu_via_metric_dot(double alpha, const double beta_up[3],
                         const SpatialMetric &gamma, const double p_cov[3],
                         double p_t_up, const double u_up[4]) {
  const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
  const double p_up[4] = {
      p_t_up,
      gu.xx * p_cov[0] + gu.xy * p_cov[1] + gu.xz * p_cov[2] -
          beta_up[0] * p_t_up,
      gu.xy * p_cov[0] + gu.yy * p_cov[1] + gu.yz * p_cov[2] -
          beta_up[1] * p_t_up,
      gu.xz * p_cov[0] + gu.yz * p_cov[1] + gu.zz * p_cov[2] -
          beta_up[2] * p_t_up};
  const SpacetimeMetric m = spacetime_metric_from_adm(alpha, beta_up, gamma);
  return -metric_dot(m, p_up, u_up);
}

} // namespace

void append_interaction_rows(Battery &b) {
  constexpr double inf = std::numeric_limits<double>::infinity();

  // --- Row 106: [MCNX-INT-01] draw law vs independently composed closed
  // forms, machine tier. u = 1/2 -> Dt = ln(2) p^t/(kappa nu); u = 3/4 ->
  // Dt = ln(4) p^t/(kappa nu). kappa enters through the cgs -> code wrapper
  // so a direction-inverted conversion cannot pass (the closed form uses
  // the length factor directly). ---
  {
    const double p_t_up = 2.0;
    const double nu_code = 4.0;
    const double kappa_cgs = 0.5; // cm^-1
    const double kappa_code = opacity_cgs_to_code(kappa_cgs);
    const double scale = p_t_up / (kappa_cgs * length_code_to_cgs * nu_code);
    const bool ok =
        detail::approx_eq(interaction_time(0.5, p_t_up, kappa_code, nu_code),
                          std::log(2.0) * scale, detail::rtol_machine) &&
        detail::approx_eq(interaction_time(0.75, p_t_up, kappa_code, nu_code),
                          std::log(4.0) * scale, detail::rtol_machine);
    b.add_boolean("int.draw.closed_form", ok);
  }

  // --- Row 107: [MCNX-INT-01] u = 0 with kappa > 0 -> Dt == 0 exactly
  // (-ln(1) is -0.0; the event fires immediately). ---
  b.add_exact("int.draw.u_zero_exact", interaction_time(0.0, 2.0, 0.5, 4.0),
              0.0);

  // --- Row 108: [MCNX-INT-02] zero-opacity sentinel: kappa = 0 gives +inf
  // and NEVER NaN, at both the reachable exact RNG output u = 0 (the naive
  // 0/0 path) and a generic u > 0. ---
  {
    const double dt0 = interaction_time(0.0, 2.0, 0.0, 4.0);
    const double dt1 = interaction_time(0.9, 2.0, 0.0, 4.0);
    const bool ok = std::isinf(dt0) && dt0 > 0.0 && !std::isnan(dt0) &&
                    std::isinf(dt1) && dt1 > 0.0 && !std::isnan(dt1);
    b.add_boolean("int.draw.zero_opacity_no_nan", ok);
  }

  // --- Row 109: [MCNX-INT-02] bitwise tie Dt_s == Dt_a resolves to
  // SCATTERING (exact; the structural dt_s <= dt_a rule), on a finite tie
  // and on a tie between two interaction_time evaluations with identical
  // inputs (bitwise-equal by purity). ---
  {
    const EpisodeOutcome tie = compete_episode(1.5, 1.5, 10.0, 20.0);
    const double dt_s = interaction_time(0.5, 2.0, 0.5, 4.0);
    const double dt_a = interaction_time(0.5, 2.0, 0.5, 4.0);
    const EpisodeOutcome tie2 = compete_episode(dt_s, dt_a, 1.0e3, 2.0e3);
    const bool ok = tie.kind == EpisodeEnd::Scatter && tie.dt_event == 1.5 &&
                    dt_s == dt_a && tie2.kind == EpisodeEnd::Scatter &&
                    tie2.dt_event == dt_s;
    b.add_boolean("int.compete.tie_scatter", ok);
  }

  // --- Row 110: [MCNX-INT-02] competition ordering sweep: each of the four
  // outcomes wins in turn; strict precedence at candidate == cell-exit; the
  // +inf sentinel channels never fire. ---
  {
    bool ok = true;
    {
      const EpisodeOutcome o = compete_episode(1.0, 2.0, 3.0, 4.0);
      ok = ok && o.kind == EpisodeEnd::Scatter && o.dt_event == 1.0;
    }
    {
      const EpisodeOutcome o = compete_episode(2.0, 1.0, 3.0, 4.0);
      ok = ok && o.kind == EpisodeEnd::Absorb && o.dt_event == 1.0;
    }
    {
      const EpisodeOutcome o = compete_episode(3.0, 4.0, 1.0, 2.0);
      ok = ok && o.kind == EpisodeEnd::CellExit && o.dt_event == 1.0;
    }
    {
      const EpisodeOutcome o = compete_episode(3.0, 4.0, 2.0, 1.0);
      ok = ok && o.kind == EpisodeEnd::StepEnd && o.dt_event == 1.0;
    }
    { // strict precedence: candidate == cell exit does not fire
      const EpisodeOutcome o = compete_episode(2.0, 3.0, 2.0, 4.0);
      ok = ok && o.kind == EpisodeEnd::CellExit && o.dt_event == 2.0;
    }
    { // both channels transparent, cell boundary inside the step
      const EpisodeOutcome o = compete_episode(inf, inf, 2.0, 3.0);
      ok = ok && o.kind == EpisodeEnd::CellExit && o.dt_event == 2.0;
    }
    { // one transparent channel: the finite one still competes
      const EpisodeOutcome o = compete_episode(inf, 5.0, 7.0, 6.0);
      ok = ok && o.kind == EpisodeEnd::Absorb && o.dt_event == 5.0;
    }
    { // everything transparent/unbounded but the step: finish with no event
      const EpisodeOutcome o = compete_episode(inf, inf, inf, 3.0);
      ok = ok && o.kind == EpisodeEnd::StepEnd && o.dt_event == 3.0;
    }
    b.add_boolean("int.compete.ordering", ok);
  }

  // --- Row 111: [MCNX-INT-01] nu = -p_mu u^mu ADM shortcut, flat fixtures,
  // vs the independent metric_dot route: static fluid (nu == p^t exactly)
  // and a boosted fluid u = (1.25, 0.75, 0, 0) (v = 0.6, gamma = 1.25,
  // exact binary), machine tier. ---
  {
    const double beta0[3] = {0.0, 0.0, 0.0};
    const SpatialMetric flat{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    const InverseSpatialMetric gu = spatial_metric_inverse(flat);
    const double p_cov[3] = {0.3, -0.4, 0.12};
    const double p_t_up = p_t_closure(p_cov[0], p_cov[1], p_cov[2], gu, 1.0);

    const double ustat[4] = {1.0, 0.0, 0.0, 0.0};
    const double nu_stat = fluid_frame_energy(1.0, beta0, p_cov, p_t_up, ustat);
    const double uboost[4] = {1.25, 0.75, 0.0, 0.0};
    const double nu_boost =
        fluid_frame_energy(1.0, beta0, p_cov, p_t_up, uboost);

    const bool ok =
        nu_stat == p_t_up &&
        detail::approx_eq(
            nu_stat, nu_via_metric_dot(1.0, beta0, flat, p_cov, p_t_up, ustat),
            detail::rtol_machine) &&
        detail::approx_eq(
            nu_boost,
            nu_via_metric_dot(1.0, beta0, flat, p_cov, p_t_up, uboost),
            detail::rtol_machine);
    b.add_boolean("int.nu.flat_identity", ok);
  }

  // --- Row 112: [MCNX-INT-01] nu identity on the pinned curved fixture
  // (nonzero shift, non-diagonal gamma; p^t from the production null
  // closure; u^mu normalized from a timelike seed so g u u = -1 to machine),
  // vs the metric_dot route, machine tier. ---
  {
    const double alpha = 1.1;
    const double beta_up[3] = {0.1, -0.05, 0.2};
    const SpatialMetric gamma{1.3, 0.12, -0.05, 1.7, 0.08, 1.1};
    const double p_cov[3] = {0.4, -0.25, 0.3};
    const InverseSpatialMetric gu = spatial_metric_inverse(gamma);
    const double p_t_up =
        p_t_closure(p_cov[0], p_cov[1], p_cov[2], gu, alpha);

    const SpacetimeMetric m = spacetime_metric_from_adm(alpha, beta_up, gamma);
    const double w[4] = {1.0, 0.2, -0.1, 0.15};
    const double n = std::sqrt(-metric_dot(m, w, w));
    const double u_up[4] = {w[0] / n, w[1] / n, w[2] / n, w[3] / n};

    const double nu_shortcut =
        fluid_frame_energy(alpha, beta_up, p_cov, p_t_up, u_up);
    const double nu_metric =
        nu_via_metric_dot(alpha, beta_up, gamma, p_cov, p_t_up, u_up);
    const bool ok =
        detail::approx_eq(nu_shortcut, nu_metric, detail::rtol_machine) &&
        detail::approx_eq(metric_dot(m, u_up, u_up), -1.0,
                          detail::rtol_machine) &&
        nu_shortcut > 0.0;
    b.add_boolean("int.nu.curved_identity", ok);
  }

  // --- Row 113: [MCNX-INT-06] k-order / e-increment policy on a synthetic
  // two-episode script, exact. Episode at e = 1 (creation was e = 0) ends in
  // scattering -> consumes exactly k = 0, 1, 2, 3; the next episode at
  // e = 2 is eventless -> exactly k = 0, 1. Every drawn value must equal the
  // pinned direct u(S, q, e, k) call bitwise, and the recorded (e, k) trace
  // must match the draw-map audit sequence exactly. ---
  {
    struct DrawRecord {
      std::uint32_t e, k;
      double v;
    };
    DrawRecord trace[8];
    int ntrace = 0;
    const auto record = [&](std::uint32_t e, std::uint32_t k, double v) {
      trace[ntrace].e = e;
      trace[ntrace].k = k;
      trace[ntrace].v = v;
      ++ntrace;
    };

    // Episode 1 (e incremented from creation's 0): unconditional channel
    // draws, then the episode ends in scattering, so the direction pair is
    // consumed too.
    const std::uint32_t e1 = 1u;
    const EpisodeUniforms eu1 = episode_uniforms(fix_S, fix_q, e1);
    record(e1, draw_k_scatter_time, eu1.u_s);
    record(e1, draw_k_absorb_time, eu1.u_a);
    const ScatterUniforms su1 = scatter_uniforms(fix_S, fix_q, e1);
    record(e1, draw_k_cos_theta, su1.u1);
    record(e1, draw_k_phi, su1.u2);

    // Episode 2 (post-scatter, e incremented by 1): eventless — exactly the
    // two channel draws, nothing else (the undrawn remainder is discarded).
    const std::uint32_t e2 = 2u;
    const EpisodeUniforms eu2 = episode_uniforms(fix_S, fix_q, e2);
    record(e2, draw_k_scatter_time, eu2.u_s);
    record(e2, draw_k_absorb_time, eu2.u_a);

    const std::uint32_t want_e[6] = {1u, 1u, 1u, 1u, 2u, 2u};
    const std::uint32_t want_k[6] = {0u, 1u, 2u, 3u, 0u, 1u};
    bool ok = ntrace == 6;
    for (int i = 0; ok && i < 6; ++i)
      ok = trace[i].e == want_e[i] && trace[i].k == want_k[i] &&
           trace[i].v == u(fix_S, fix_q, want_e[i], want_k[i]);
    b.add_boolean("int.rng.episode_draw_map", ok);
  }
}

} // namespace MCNuX
