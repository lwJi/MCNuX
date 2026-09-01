// The runtime episode driver of specs/neutrino-matter-interactions.md
// [MCNX-INT-01/02/03/04/05/06]: per packet and per transport step, a
// sequence of sampling episodes — interaction-time draws against the
// cell-exit time and the step remainder — with elastic scattering applied at
// fixed fluid-frame energy nu through the pinned tetrad transform,
// absorption removing the packet whole, every event's momentum/lepton
// exchange deposited into the episode cell ([MCNX-HYD-02] via the shared C2
// helper), and the identical pre-negation deltas accumulated into the
// per-step event-side audit ([MCNX-HYD-05]).
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
// Absorption ([MCNX-INT-03]): an Absorb outcome deposits the packet's FULL
// content (absorption_delta — no weight decay, the packet dies whole) and
// removes the packet: make_invalid clears the idcpu validity bit and the
// unconditional per-patch Redistribute() after the tile walk compacts the
// slot away. The id VALUE is retired per [MCNX-GPU-04] — AMReX's
// ParticleType::NextID() is monotonic and never rewinds, so a removed
// packet's id is never reissued (the storage slot is reused by compaction;
// the id value is what the contract covers).
//
// This file also owns MCNuX_IdContractCheck (parameter-gated by
// MCNuX::test_id_contract): the runtime [MCNX-GPU-04] id-contract guard over
// creation + removal — live logical keys pairwise distinct, ids in range,
// bit 63 of the logical key clear, no retired key ever reappears.

#include "mcnux_gather.hxx" // does the CarpetX driver.hxx relative include

#include "mcnux_coefficients.hxx"
#include "mcnux_deposit.hxx"
#include "mcnux_emission.hxx" // packet_rng_key
#include "mcnux_escape.hxx"
#include "mcnux_fluid.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_interactions.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_table_range.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_trp.hxx"
#include "mcnux_units.hxx"

#include <AMReX_AmrParGDB.H> // completes AmrParGDB (mcnux_geodesic.cxx note)
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particle.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <set>
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
  // (host counters are not device-safe). [MCNX-INT-05]: the single narrow
  // evaluate_coefficients expression below is where the trapped-regime
  // relabeling of mcnux_trp.hxx substitutes (event-sampling law only,
  // [MCNX-TRP-02]) when trapped_scheme = "relabeled"; with the default
  // explicit scheme the sampling expressions are byte-for-byte the baseline.
  const CoefficientSource src = selected_coefficient_source();
  const AnalyticOpacityParams ap = analytic_params_from_parameters();
  const TrpParams trp = trp_params_from_parameters();
  const RangedTableCoefficients table_eval{};

  const std::uint64_t S = static_cast<std::uint64_t>(rng_seed);
  const double dt_step = cctk_delta_time; // coordinate, code units

  const MetricGroups mgroups = metric_groups();
  const HydroGroups hgroups = hydro_groups();
  const SourceGroups sgroups = source_groups();

  // Per-step event-side audit ([MCNX-HYD-05]): 10 device accumulator slots
  // (net dPt..dL; gross |dPt|..|dL| — shared by the scattering AND
  // absorption channels, correct per the per-cell episode framing) plus
  // four counters for the log line (scattering events, episodes, max
  // episodes of any one packet — the near-cap observable — and absorption
  // events). Step-reset here.
  interaction_step_audit() = LedgerAudit{};
  amrex::Gpu::DeviceVector<double> audit_dev(14, 0.0);
  double *const audit_ptr = audit_dev.data();

  // ONE step-scoped escape buffer shared by every tile's kernel
  // ([MCNX-GPU-05], mcnux_escape.hxx; the push-site precedent in
  // mcnux_geodesic.cxx — per-tile buffers would lose counts).
  amrex::Gpu::DeviceVector<double> esc_dev(escape_num_slots, 0.0);
  double *const esc_ptr = esc_dev.data();

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

    // Cell light-crossing time of the [MCNX-TRP-04] alpha selection: code
    // units with c = 1, so Dt_c = min_d dx_d. Uniform over the tile (one
    // level's geometry), computed host-side and captured by value; consumed
    // only on the relabeled branch below.
    const double trp_dt_c =
        std::min(fgather.dx[0], std::min(fgather.dx[1], fgather.dx[2]));

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
      // Last-in-domain position anchor of the escape energy tally
      // (mcnux_escape.hxx convention): refreshed at every episode start
      // while the packet is still inside the domain, so an escaping
      // packet's p^t is evaluated on a metric gather that never leaves
      // the domain.
      double x_last_inside[3] = {x[0], x[1], x[2]};
      int it = 0;
      for (; it < max_episodes_per_step && dt_rem > 0.0; ++it) {
        if (!outside_domain(x, mgather.prob_lo, mgather.prob_hi)) {
          x_last_inside[0] = x[0];
          x_last_inside[1] = x[1];
          x_last_inside[2] = x[2];
        }
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

        // Trapped-regime relabeling ([MCNX-TRP-02], event-sampling law
        // substitution point): with trapped_scheme = "relabeled", alpha per
        // (cell, species, energy bin) is the [MCNX-TRP-04] selection rule on
        // the UNPRIMED code-unit kappa_a (or the fixed-alpha override), and
        // the primed pair feeds the two interaction_time draws below —
        // NOWHERE else: the RNG draw sites are scheme-independent (the
        // [MCNX-TRP-05] alpha = 1 bitwise obligation), and every deposit,
        // audit, and tally keeps the actual fired event's quantities. On the
        // default explicit branch these are bitwise copies (zero FP ops).
        double kappa_a_samp = kappa_a_code;
        double kappa_s_samp = kappa_s_code;
        if (trp.relabeled) {
          const double alpha =
              (trp.alpha_fixed > 0.0)
                  ? trp.alpha_fixed
                  : alpha_select(kappa_a_code, trp_dt_c, trp.xi);
          // Map by NAME (RelabeledCoefficients field order differs from
          // Coefficients); the eta slot is unused in the sampling law.
          const RelabeledCoefficients rc =
              relabel(0.0, kappa_a_code, kappa_s_code, alpha);
          kappa_a_samp = rc.kappa_a_p;
          kappa_s_samp = rc.kappa_s_p;
        }

        // (4) [MCNX-INT-06]: both channel draws ALWAYS consumed; the
        // cell-exit bound uses the FROZEN episode-start coordinate velocity
        // and the episode cell's faces from above.
        const EpisodeUniforms eu = episode_uniforms(S, q, e);
        const double dt_s =
            interaction_time(eu.u_s, pt_up, kappa_s_samp, nu_code);
        const double dt_a =
            interaction_time(eu.u_a, pt_up, kappa_a_samp, nu_code);
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
          // Absorption ([MCNX-INT-03]): deposit the packet's FULL content
          // (-N p^t, -N p_i, -l_s N via absorption_delta — no weight decay,
          // the packet dies whole) into the EPISODE cell with the identical
          // (dV, dt) pair, mirror the pre-negation delta into the audit,
          // then remove the packet: make_invalid clears the idcpu validity
          // bit and the unconditional Redistribute() below compacts the
          // slot away ([MCNX-GPU-04]: the id VALUE is retired — NextID()
          // never rewinds, so it is never reissued). p[0..2] are the
          // event-time momenta with the FROZEN episode p^t (the
          // frozen-kinematics convention). Both channel draws were already
          // consumed above, so the [MCNX-INT-06] trace stays exact.
          const LedgerDelta d =
              absorption_delta(s, N, pt_up, p[0], p[1], p[2]);
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
          amrex::Gpu::Atomic::AddNoRet(&audit_ptr[13], 1.0); // absorb count
          amrex::particle_impl::make_invalid(ptd.m_idcpu[ip]);
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

      // Escape check ([MCNX-GPU-05], mcnux_escape.hxx): a packet whose
      // final position left the domain is tallied per species and removed
      // via make_invalid (tally strictly BEFORE invalidation) — the same
      // sequence as the push site, using the last-in-domain metric anchor
      // for p^t. Skipped for packets the Absorb branch already invalidated
      // (their validity bit is clear); the unconditional Redistribute()
      // below compacts both removal routes alike.
      if ((ptd.m_idcpu[ip] >> 63) != 0 &&
          outside_domain(x, mgather.prob_lo, mgather.prob_hi)) {
        const MetricSnapshot me =
            mgather(x_last_inside[0], x_last_inside[1], x_last_inside[2]);
        const InverseSpatialMetric gue = spatial_metric_inverse(me.g);
        const double pt_esc = p_t_closure(p[0], p[1], p[2], gue, me.alpha);
        const int sidx = species_index(s);
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_count)], 1.0);
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_number)], N);
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_energy)], N * pt_esc);
        amrex::particle_impl::make_invalid(ptd.m_idcpu[ip]);
      }

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
  std::vector<double> audit_host(14, 0.0);
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, audit_dev.begin(), audit_dev.end(),
                   audit_host.begin());
  LedgerAudit &audit = interaction_step_audit();
  audit.net = LedgerDelta{audit_host[0], audit_host[1], audit_host[2],
                          audit_host[3], audit_host[4]};
  audit.gross = LedgerDelta{audit_host[5], audit_host[6], audit_host[7],
                            audit_host[8], audit_host[9]};

  // Fold the step's escapes into the run-cumulative tally ([MCNX-GPU-05];
  // NOT part of the ledger audit above — escapes deposit nothing and never
  // join the closure, per the mcnux_escape.hxx storage-surface decision).
  {
    std::vector<double> esc_host(escape_num_slots, 0.0);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, esc_dev.begin(), esc_dev.end(),
                     esc_host.begin());
    escape_run_tally().add_slots(esc_host.data());
  }

  CCTK_VINFO("MCNuX episode driver: %ld episodes, %ld scattering events, "
             "%ld absorption events, max %ld episodes/packet (cap %d)",
             long(audit_host[11]), long(audit_host[10]), long(audit_host[13]),
             long(audit_host[12]), max_episodes_per_step);

  // Refresh the per-packet golden table when the synthetic fixture is the
  // population (the `interactions-fixedseed`/`absorb-smoke` benchmarks; rows
  // keyed id - 1, SIZE cross-checked inside). It runs AFTER the
  // Redistribute() above (preserve that ordering): the table is zeroed then
  // filled from live packets only, so an absorbed packet's row reverting to
  // all-zero is the built-in discrete-removal observable. Production
  // populations get their own diagnostics elsewhere (MCNuX_EmissionDiag).
  if (test_synthetic_packets)
    fill_packet_diag(mgroups, {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz, pk_pt},
                     packet_diag_size());
}

// ---------------------------------------------------------------------------
// Runtime id-contract check  [MCNX-GPU-04]
// ---------------------------------------------------------------------------
// Parameter-gated verification scaffolding (MCNuX::test_id_contract; the
// `absorb-smoke` benchmark): once per transport step, after the episode
// driver (removal) — and generically covering creation too — sample every
// live packet's id and assert the [MCNX-GPU-04] contract. A hard-abort
// guard (CCTK_VERROR on any violation), no golden artifact.
//
// Bit-63 reading (decided here, the mcnux_emission.hxx precedent): the RAW
// AMReX idcpu word has bit 63 SET on every live particle by construction —
// it is AMReX's validity flag (pack_id, amrex/Src/Particle/AMReX_Particle.H)
// — so the Verification bullet "bit 63 clear" of
// specs/particle-container-and-gpu.md is testable only on the REPACKED
// logical key q = (cpu << 39) | id, exactly the packet_rng_key packing that
// the RNG already uses for the same reason. The check therefore samples ids
// (never storage-slot indices), repacks, and asserts on the logical key.
//
// Rank-local by construction (the harness runs NPROCS 1); the function-static
// sets persist across the run, which is what makes "no retired id ever
// reappears" checkable over the whole run.
extern "C" void MCNuX_IdContractCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_IdContractCheck;
  DECLARE_CCTK_PARAMETERS;

  require_driver();

  // Gather every live packet's raw idcpu word to the host (tile SoA vectors
  // live in device memory on a GPU build).
  std::vector<std::uint64_t> raw;
  for_each_packet_tile_raw([&](const CarpetX::GHExt::PatchData &,
                               const CarpetX::GHExt::PatchData::LevelData &,
                               PacketContainer &, const PacketIter &pti) {
    const auto &idcpu =
        pti.GetParticleTile().GetStructOfArrays().GetIdCPUData();
    const std::size_t n0 = raw.size();
    raw.resize(n0 + idcpu.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, idcpu.begin(), idcpu.end(),
                     raw.begin() + n0);
  });
  amrex::Gpu::streamSynchronize();

  static std::set<std::uint64_t> ever_live; // every logical key ever seen
  static std::set<std::uint64_t> retired;   // keys that have vanished
  static std::set<std::uint64_t> prev_live; // live keys at the last check

  std::set<std::uint64_t> live;
  for (const std::uint64_t word : raw) {
    // Raw-word validity: after the driver's Redistribute() (default
    // remove_negative = true) no invalidated slot may survive in the tiles.
    if (!(word >> 63))
      CCTK_VERROR("[MCNX-GPU-04] id contract: an INVALID packet slot "
                  "(idcpu = 0x%llx) survived Redistribute()",
                  static_cast<unsigned long long>(word));
    const long long id = amrex::particle_impl::unpack_id(word);
    const int cpu = amrex::particle_impl::unpack_cpu(word);
    if (id < 1 || id > (1ll << 39) - 3)
      CCTK_VERROR("[MCNX-GPU-04] id contract: packet id %lld outside "
                  "[1, 2^39 - 3]",
                  id);
    if (cpu < 0 || cpu >= (1 << 24))
      CCTK_VERROR("[MCNX-GPU-04] id contract: packet cpu %d outside "
                  "[0, 2^24)",
                  cpu);
    const std::uint64_t q = packet_rng_key(static_cast<std::uint64_t>(id),
                                           static_cast<std::uint64_t>(cpu));
    if (q >> 63)
      CCTK_VERROR("[MCNX-GPU-04] id contract: logical key 0x%llx has "
                  "bit 63 set",
                  static_cast<unsigned long long>(q));
    if (!live.insert(q).second)
      CCTK_VERROR("[MCNX-GPU-04] id contract: duplicate live packet id "
                  "%lld (cpu %d)",
                  id, cpu);
    if (retired.count(q))
      CCTK_VERROR("[MCNX-GPU-04] id contract: retired packet id %lld "
                  "(cpu %d) reappeared",
                  id, cpu);
    ever_live.insert(q);
  }

  // Keys live at the last check but absent now have been removed: their id
  // values are retired for the rest of the run.
  for (const std::uint64_t q : prev_live)
    if (!live.count(q))
      retired.insert(q);
  prev_live = std::move(live);

  CCTK_VINFO("MCNuX id-contract check: %zu live, %zu ever live, %zu retired "
             "(all [MCNX-GPU-04] assertions passed)",
             prev_live.size(), ever_live.size(), retired.size());
}

} // namespace MCNuX
