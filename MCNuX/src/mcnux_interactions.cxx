// The runtime episode driver of specs/neutrino-matter-interactions.md
// [MCNX-INT-01/02/04/05/06] (scattering-only increment): per packet and per
// transport step, a sequence of sampling episodes — interaction-time draws
// against the cell-exit time and the step remainder — with elastic
// scattering applied at fixed fluid-frame energy nu through the pinned
// tetrad transform, the scattering momentum exchange deposited into the
// episode cell ([MCNX-HYD-02] via the shared C2 helper), and the identical
// pre-negation deltas accumulated into the per-step event-side audit
// ([MCNX-HYD-05]).
//
// All episode math is the already-verified pure functions of
// mcnux_interactions.hxx / mcnux_tetrad.hxx / mcnux_fluid.hxx /
// mcnux_geodesic.hxx — this file is glue: it calls, never re-derives.
//
// Scheduling (MCNuX/schedule.ccl): MCNuX_EpisodeDriver runs IN
// MCNuX_TransportStep AFTER MCNuX_ZeroSourceTerms BEFORE
// MCNuX_AddToSourceTerms — the MCNuX_GeodesicPush slot; when
// enable_interactions is set the plain push is NOT scheduled (its condition
// carries && !enable_interactions), so the two never double-push the same
// packets. The driver subsumes the push: every packet advances by exactly
// cctk_delta_time of coordinate time across its episode segments
// (geodesic_step per segment, [MCNX-GEO-01/04]).
//
// Spec-normative deviations preserved (neutrino-matter-interactions.md,
// Open questions — spec text normative per specs/README.md):
//   * Per-cell episodes: a NEW episode (fresh draws, e incremented) begins
//     at every cell crossing; the undrawn remainder is discarded.
//   * Frozen episode kinematics: p^t, nu, the opacities, the tetrad, and the
//     coordinate velocity dx^i/dt are evaluated ONCE at episode start and
//     held fixed for the episode (the segment advance itself is the RK4
//     geodesic_step, which re-gathers per stage as [MCNX-GEO-03] requires).
//   * Ties resolve to scattering (structural inside compete_episode).
//
// Scattering-only placeholder: an Absorb outcome ends the packet's step with
// NO deposit and NO removal — the [MCNX-INT-03] removal machinery is a later
// task. The channel draws are consumed either way, so the [MCNX-INT-06]
// audit trace stays exact; benchmarks keep the absorption channel
// transparent (kappa_a = 0 -> Dt_a = +inf) until removal lands.

#include "mcnux_gather.hxx" // does the CarpetX driver.hxx relative include

#include "mcnux_coefficients.hxx"
#include "mcnux_deposit.hxx"
#include "mcnux_emission.hxx" // packet_rng_key
#include "mcnux_fluid.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_interactions.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_table_range.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <AMReX_AmrParGDB.H> // completes AmrParGDB (mcnux_geodesic.cxx note)
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particle.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdint>
#include <limits>
#include <vector>

namespace MCNuX {

namespace {

// The per-step event-side audit accumulator behind interaction_step_audit()
// (the g_emission_audit pattern of mcnux_emission.cxx; the two contributors
// own their accumulators independently).
LedgerAudit g_interaction_audit{};

// Defensive per-packet episode cap. Unreachable safety net: dt_rem strictly
// decreases across episodes except on a dt = 0 CellExit (a packet bitwise ON
// a cell face moving toward it — measure-zero off contrived fixtures), and
// the cap turns even that into a plain break instead of a device-side hang.
constexpr int max_episodes_per_step = 10000;

} // namespace

LedgerAudit &interaction_step_audit() { return g_interaction_audit; }

extern "C" void MCNuX_EpisodeDriver(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_EpisodeDriver;
  DECLARE_CCTK_PARAMETERS;

  require_driver();

  // Captured once per step, host side (the mcnux_emission.cxx shape). The
  // table slot is a RangedTableCoefficients with default (null) views —
  // safe because MCNuX_ParamCheck forbids enable_interactions with
  // opacity_source = "table", so the Analytic branch of
  // evaluate_coefficients never invokes it; its clamp counters stay null
  // (host counters are not device-safe). [MCNX-INT-05]: this single narrow
  // evaluate_coefficients expression is where the trapped-regime relabeling
  // of mcnux_trp.hxx substitutes in a later task.
  const CoefficientSource src = selected_coefficient_source();
  const AnalyticOpacityParams ap = analytic_params_from_parameters();
  const RangedTableCoefficients table_eval{};

  const std::uint64_t S = static_cast<std::uint64_t>(rng_seed);
  const double dt_step = cctk_delta_time; // coordinate, code units

  const MetricGroups mgroups = metric_groups();
  const HydroGroups hgroups = hydro_groups();
  const SourceGroups sgroups = source_groups();

  // Per-step event-side audit ([MCNX-HYD-05]): 10 device accumulator slots
  // (net dPt..dL; gross |dPt|..|dL|) plus three counters for the log line
  // (scattering events, episodes, max episodes of any one packet — the
  // near-cap observable). Step-reset here.
  interaction_step_audit() = LedgerAudit{};
  amrex::Gpu::DeviceVector<double> audit_dev(13, 0.0);
  double *const audit_ptr = audit_dev.data();

  for_each_packet_tile_raw([&](const CarpetX::GHExt::PatchData &patchdata,
                               const CarpetX::GHExt::PatchData::LevelData
                                   &leveldata,
                               PacketContainer &, const PacketIter &pti) {
    const long np = pti.numParticles();
    if (np == 0)
      return;
    const auto ptd = pti.GetParticleTile().getParticleTileData();

    // Both gathers of the episode: the trilinear vertex metric ([MCNX-GEO-03])
    // and the one-cell fluid state ([MCNX-INT-05]), built for this tile's box.
    const VertexMetricGather mgather =
        make_gather(patchdata, leveldata, mgroups, pti);
    const CellFluidGather fgather =
        make_fluid_gather(patchdata, leveldata, hgroups, pti);

    // Source views of the SAME box: CarpetX builds the group MultiFabs on
    // the AmrCore's BoxArray/DistributionMapping — the ones the ParGDB
    // tracks — so the particle iterator's box index addresses the matching
    // grid fab (the box-index reasoning of mcnux_gather.hxx; first deposit
    // from a packet-tile walk).
    constexpr int tl = 0;
    amrex::MultiFab &rad_mf =
        *leveldata.groupdata.at(sgroups.rad_force)->mfab.at(tl);
    amrex::MultiFab &lep_mf =
        *leveldata.groupdata.at(sgroups.lepton_source)->mfab.at(tl);
    const int box = pti.index();
    const SourceViews views{rad_mf.array(box, 0), rad_mf.array(box, 1),
                            rad_mf.array(box, 2), rad_mf.array(box, 3),
                            lep_mf.array(box, 0)};

    // Coordinate cell volume (code units) — the deposit/audit normalization
    // pair with dt_step ([MCNX-HYD-02]/[MCNX-HYD-05]; NEVER the cgs pair).
    const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
    const double dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      double x[3] = {ptd.rdata(PIdx::x)[ip], ptd.rdata(PIdx::y)[ip],
                     ptd.rdata(PIdx::z)[ip]};
      double p[3] = {ptd.rdata(PIdx::px)[ip], ptd.rdata(PIdx::py)[ip],
                     ptd.rdata(PIdx::pz)[ip]};
      const double N = ptd.rdata(PIdx::w)[ip]; // weight: NEVER changed here
      const Species s =
          static_cast<Species>(ptd.idata(IntIdx::species)[ip]);
      // The per-packet draw key: the REPACKED logical key of packet_rng_key
      // (bit 63 provably clear), NEVER the raw idcpu word (the
      // mcnux_emission.hxx rationale).
      const std::uint64_t q = packet_rng_key(
          std::uint64_t(amrex::particle_impl::unpack_id(ptd.m_idcpu[ip])),
          std::uint64_t(amrex::particle_impl::unpack_cpu(ptd.m_idcpu[ip])));

      double dt_rem = dt_step;
      int it = 0;
      for (; it < max_episodes_per_step && dt_rem > 0.0; ++it) {
        // (1) [MCNX-INT-06]: e increments by 1 at EVERY episode start
        // (step start, post-scatter, cell crossing; creation left e = 0,
        // so the first episode uses e = 1).
        const std::uint32_t e = static_cast<std::uint32_t>(
            ++ptd.idata(IntIdx::event_counter)[ip]);

        // (2) Freeze the episode kinematics at the current position:
        // metric snapshot, null-closure p^t, coordinate velocity,
        // containing cell, fluid state, Valencia u^mu, fluid-frame energy
        // nu (all held fixed below). The geodesic RHS supplies the FROZEN
        // dx^i/dt of the cell-exit bound ([MCNX-GEO-01]; its pt member is
        // the same closure as pt_up).
        const MetricSnapshot ms = mgather(x[0], x[1], x[2]);
        const InverseSpatialMetric gu = spatial_metric_inverse(ms.g);
        const double pt_up = p_t_closure(p[0], p[1], p[2], gu, ms.alpha);
        const GeodesicRhs rhs = geodesic_rhs(ms, p);

        // Episode cell: the half-open anchor of CellFluidGather, with a
        // direction-aware entry adjustment. Two roundoff configurations of
        // the bare anchor cannot make progress or even move backwards:
        // (a) a packet BITWISE on its cell's lower face moving toward -d is
        // ENTERING the lower neighbor, but the anchor re-assigns it to the
        // cell it just left, whose exit time along -d is exactly 0 — an
        // infinite dt = 0 CellExit loop (exact face landings are routine on
        // flat metrics with binary-exact grids); (b) a packet a few ulps
        // BELOW a face (the floor ratio rounds the tiny offset away, so the
        // anchor still names the upper cell) gets a NEGATIVE exit candidate
        // and ping-pongs across the face with shrinking amplitude — dozens
        // of wasted episodes per crossing (observed at the x = 0 faces).
        // Cell-exit tracking is implementation freedom
        // (neutrino-matter-interactions.md, Implementation freedom); the
        // documented choice here is that a packet at-or-below its anchored
        // cell's lower face while moving toward -d belongs to the LOWER
        // neighbor — its coefficients, deposit target, and exit faces are
        // that cell's, and the exit candidate is strictly positive. The
        // pure dt_to_cell_exit convention (selftest row 144) is untouched —
        // the driver simply never feeds it a stuck configuration.
        int ic[3];
        fgather.containing_cell(x[0], x[1], x[2], ic);
        double face_lo[3], face_hi[3];
        for (int d = 0; d < 3; ++d) {
          face_lo[d] = fgather.prob_lo[d] + ic[d] * fgather.dx[d];
          if (x[d] <= face_lo[d] && rhs.dxdt[d] < 0.0) {
            --ic[d];
            face_lo[d] -= fgather.dx[d];
          }
          face_hi[d] = face_lo[d] + fgather.dx[d];
        }
        const FluidSample fs = fgather.at_cell(ic);
        double u4[4];
        valencia_four_velocity(ms.alpha, ms.beta, ms.g, fs.vel, u4);
        const double nu_code = fluid_frame_energy(ms.alpha, ms.beta, p,
                                                  pt_up, u4);

        // (3) Coefficients ONLY via the source-agnostic interface at the
        // cell fluid state ([MCNX-INT-05]), lookup energy in MeV, opacities
        // cm^-1 -> code via the pinned wrapper (never a retyped literal).
        const Coefficients co =
            evaluate_coefficients(src, ap, table_eval, s,
                                  fluid_frame_energy_MeV(nu_code), fs.state);
        const double kappa_s_code = opacity_cgs_to_code(co.kappa_s);
        const double kappa_a_code = opacity_cgs_to_code(co.kappa_a);

        // (4) [MCNX-INT-06]: both channel draws ALWAYS consumed; the
        // cell-exit bound uses the FROZEN episode-start coordinate velocity
        // and the episode cell's faces from above.
        const EpisodeUniforms eu = episode_uniforms(S, q, e);
        const double dt_s =
            interaction_time(eu.u_s, pt_up, kappa_s_code, nu_code);
        const double dt_a =
            interaction_time(eu.u_a, pt_up, kappa_a_code, nu_code);
        // The exit bound is inflated by a few ulps (an implementation-
        // freedom choice like the entry adjustment above): the RK4 segment
        // of an exactly-computed dt_exit can land roundoff-SHORT of the
        // face, and the follow-up micro-episodes then crawl the boundary
        // one ulp at a time (dozens of wasted draws per crossing; a
        // half-ulp round-to-even tie would even stall until the episode
        // cap). Inflating by 4 eps makes a crossing episode land in (or
        // bitwise on) the entered cell in one go; the perturbation is far
        // below every physical tolerance and the competition stays exact.
        constexpr double exit_inflate =
            1.0 + 4.0 * std::numeric_limits<double>::epsilon();
        const double dt_exit =
            dt_to_cell_exit_faces(x, rhs.dxdt, face_lo, face_hi) *
            exit_inflate;

        // (5) The [MCNX-INT-02] competition.
        const EpisodeOutcome out =
            compete_episode(dt_s, dt_a, dt_exit, dt_rem);

        // (6) Advance the geodesic over the episode segment.
        geodesic_step(mgather, x, p, out.dt_event);
        dt_rem -= out.dt_event;
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[11], 1.0); // episode count

        // (7) Branch on the outcome.
        if (out.kind == EpisodeEnd::StepEnd)
          break;
        if (out.kind == EpisodeEnd::CellExit)
          continue; // new episode, undrawn remainder discarded
        if (out.kind == EpisodeEnd::Absorb) {
          // Scattering-only placeholder for [MCNX-INT-03] (a later task):
          // no deposit, no removal — the packet simply ends its step here.
          // The channel draws were already consumed, so the [MCNX-INT-06]
          // trace stays exact when removal lands.
          break;
        }

        // Scatter ([MCNX-INT-04]): the direction pair (k = 2, 3) is
        // consumed ONLY on this branch; the redraw uses the FROZEN
        // episode-start tetrad/metric/nu. p_in is the packet's p_i at
        // event time (post-segment) with the frozen episode p^t.
        const ScatterUniforms su = scatter_uniforms(S, q, e);
        const Tetrad tet = build_tetrad(ms.alpha, ms.beta, ms.g, u4);
        const SpacetimeMetric gm =
            spacetime_metric_from_adm(ms.alpha, ms.beta, ms.g);
        const CoordinateMomentum p_out =
            scatter_redraw(tet, gm, nu_code, su.u1, su.u2);

        // Ledger: the pre-negation scattering delta, deposited into the
        // EPISODE cell with the code-unit (dV, dt) pair; the audit uses the
        // IDENTICAL delta and pair ([MCNX-HYD-05] pairing).
        const LedgerDelta d =
            scattering_delta(N, pt_up, p[0], p[1], p[2], p_out.pt, p_out.px,
                             p_out.py, p_out.pz);
        p[0] = p_out.px;
        p[1] = p_out.py;
        p[2] = p_out.pz;
        deposit_delta(views, ic[0], ic[1], ic[2], d, dV, dt_step);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[0], d.dPt);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[1], d.dPx);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[2], d.dPy);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[3], d.dPz);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[4], d.dL);
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[5], detail::cabs(d.dPt));
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[6], detail::cabs(d.dPx));
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[7], detail::cabs(d.dPy));
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[8], detail::cabs(d.dPz));
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[9], detail::cabs(d.dL));
        amrex::Gpu::Atomic::AddNoRet(&audit_ptr[10], 1.0); // scatter count
        // Continue: a new episode follows a scattering ([MCNX-INT-02]).
      }

      amrex::Gpu::Atomic::Max(&audit_ptr[12], double(it + 1));

      ptd.rdata(PIdx::x)[ip] = x[0];
      ptd.rdata(PIdx::y)[ip] = x[1];
      ptd.rdata(PIdx::z)[ip] = x[2];
      ptd.rdata(PIdx::px)[ip] = p[0];
      ptd.rdata(PIdx::py)[ip] = p[1];
      ptd.rdata(PIdx::pz)[ip] = p[2];
    });
  });
  amrex::Gpu::streamSynchronize();

  // Re-establish ownership after the step's motion (the push precedent).
  for (int patch = 0; patch < num_packet_patches(); ++patch)
    packet_population(patch).Redistribute();

  // Fold the device audit into the step accumulator.
  std::vector<double> audit_host(13, 0.0);
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, audit_dev.begin(), audit_dev.end(),
                   audit_host.begin());
  LedgerAudit &audit = interaction_step_audit();
  audit.net = LedgerDelta{audit_host[0], audit_host[1], audit_host[2],
                          audit_host[3], audit_host[4]};
  audit.gross = LedgerDelta{audit_host[5], audit_host[6], audit_host[7],
                            audit_host[8], audit_host[9]};

  CCTK_VINFO("MCNuX episode driver: %ld episodes, %ld scattering events, "
             "max %ld episodes/packet (cap %d)",
             long(audit_host[11]), long(audit_host[10]), long(audit_host[12]),
             max_episodes_per_step);

  // Refresh the per-packet golden table when the synthetic fixture is the
  // population (the `interactions-fixedseed` benchmark; rows keyed
  // id - 1, SIZE cross-checked inside). Production populations get their
  // own diagnostics elsewhere (MCNuX_EmissionDiag).
  if (test_synthetic_packets)
    fill_packet_diag(mgroups, {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz, pk_pt},
                     packet_diag_size());
}

} // namespace MCNuX
