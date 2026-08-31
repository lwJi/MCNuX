// The runtime emission/creation loop of
// specs/packet-representation-and-sampling.md [MCNX-PKT-02/03/05]: per cell,
// species, and energy bin, the count law decides how many packets to create
// (one cell-keyed uniform draw per (cell, step, s, b)); created packets are
// appended to the per-patch PacketContainer with contiguous AMReX-reserved
// ids, initialized by the pinned creation draw map (six uniforms at event
// e = 0 of the packet's own logical key), and their emission LedgerDelta is
// deposited into the creation cell through the shared atomic helper of
// mcnux_deposit.hxx ([MCNX-HYD-02]/[MCNX-GPU-02]) while the identical deltas
// accumulate into the per-step event-side audit ([MCNX-HYD-05]).
//
// All sampling math is the already-verified pure functions of
// mcnux_emission.hxx / mcnux_tetrad.hxx / mcnux_fluid.hxx — this file is
// glue: it calls, never re-derives.
//
// Scheduling (MCNuX/schedule.ccl): MCNuX_Emission runs IN
// MCNuX_AddToSourceTerms (group membership gives the zero-then-add ordering
// of [MCNX-HYD-03] and automatic LedgerClosure coverage), which is ordered
// AFTER MCNuX_TransportObserver — so the observer has already incremented
// the counter and the 0-based transport-step index of the [MCNX-PKT-05]
// cell draws is n = *transport_step_count - 1. MCNuX_GeodesicPush (when
// active) runs BEFORE this group, so packets created at step n receive
// their first full-dt push at step n + 1 (their creation-time draw u_time
// is consumed by the pinned map but currently unused: the packet state
// stores no time).
//
// Append idiom (no in-repo precedent; WarpX exemplar
// warpx/Source/Particles/ParticleCreation/AddParticles.cpp,
// SmartUtils.H): count kernel -> exclusive scan -> host-side contiguous id
// reservation via ParticleType::NextID (monotonic from 1, never reused,
// unique per rank) -> tile resize -> fill kernel at slot old + offset + j.
// One Redistribute() per patch afterwards (creation writes into the owning
// cell's own tile, so this is a no-op guard).

#include "mcnux_gather.hxx" // does the CarpetX driver.hxx relative include

#include "mcnux_coefficients.hxx"
#include "mcnux_deposit.hxx"
#include "mcnux_emission.hxx"
#include "mcnux_fluid.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_srcterms.hxx"
#include "mcnux_stats.hxx"
#include "mcnux_stats_emission.hxx"
#include "mcnux_table_range.hxx"
#include "mcnux_tetrad.hxx"
#include "mcnux_units.hxx"

#include <AMReX_AmrParGDB.H> // completes AmrParGDB (mcnux_geodesic.cxx note)
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particle.H>
#include <AMReX_Scan.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdint>
#include <vector>

namespace MCNuX {

namespace {

// ---------------------------------------------------------------------------
// Analytic-mode energy-bin grid
// ---------------------------------------------------------------------------
// The analytic (gray) coefficient source has no tabulated energy grid, and
// no spec pins one for it — this array is the pinned ANALYTIC FIXTURE grid
// (the nodes-as-edges convention of mcnux_emission.hxx): linear-MeV nodes
// {5, 15, 25} -> 2 bins with centers 10 (= the default E_p) and 20 MeV.
// When table residency lands, the table path supersedes this per-table
// (log10 -> linear conversion of emab_LogEs — deliberately NOT done here);
// MCNuX_ParamCheck forbids enable_emission with opacity_source = "table"
// until then.
constexpr int emission_num_nodes = 3;
constexpr int emission_num_bins = emission_num_nodes - 1; // nodes are edges

// Device-capturable bin-edge table (function-local plain aggregate — the
// mcnux_srcterms.cxx CUDA note on namespace-scope struct constants).
struct BinTable {
  double lo[emission_num_bins];
  double hi[emission_num_bins];
};

// The per-step event-side audit accumulator behind emission_step_audit().
LedgerAudit g_emission_audit{};

// ---------------------------------------------------------------------------
// stats-emission accumulator ([MCNX-PKT-06], [MCNX-VER-07])
// ---------------------------------------------------------------------------
// Step-scoped species-keyed reduction of the emission pass, gated on
// test_stats_emission (the production kernels see only a null-pointer check
// when it is off). Per species s (base = stats_stride * s):
//   +0  Sum of the PRE-floor real-valued N_p over (cell, bin)   — the
//       expected count of the [MCNX-PKT-02] floor + Bernoulli law (using
//       the post-Bernoulli integer here would bias the count z; brief risk
//       (b)), accumulated in the count pass;
//   +1  Sum p(1 - p), p = frac(N_p)                             — the count
//       variance (independent per-cell Bernoulli remainders), count pass;
//   +2  realized created-packet count                           — fill pass;
//   +3  Sum cos(theta)                                          — fill pass
//       (cos(theta)/phi are local at draw time, never stored on the packet,
//       so they MUST be accumulated here);
//   +4  Sum phi                                                 — fill pass.
// Folded to g_stats_accum after the walk (the audit_dev idiom); consumed by
// MCNuX_StatsEmission below.
constexpr int stats_stride = 5;
constexpr int stats_num_slots = NUM_SPECIES * stats_stride;
double g_stats_accum[stats_num_slots] = {};

} // namespace

LedgerAudit &emission_step_audit() { return g_emission_audit; }

extern "C" void MCNuX_Emission(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_Emission;
  DECLARE_CCTK_PARAMETERS;

  require_driver();

  // Captured once per step, host side.
  const CoefficientSource src = selected_coefficient_source();
  const AnalyticOpacityParams ap = analytic_params_from_parameters();
  // The production-correct dispatch shape: the table slot is a
  // RangedTableCoefficients (never the bare assembly). No table-residency
  // layer exists yet, so its views are default (null) — safe because
  // MCNuX_ParamCheck forbids enable_emission with opacity_source = "table",
  // and the Analytic branch of evaluate_coefficients never invokes the slot.
  const RangedTableCoefficients table_eval{};

  const std::uint64_t S = static_cast<std::uint64_t>(rng_seed);
  // 0-based transport-step index of the [MCNX-PKT-05] cell draws: the
  // observer (this group is AFTER MCNuX_TransportObserver) has already
  // incremented the counter. This `e = n` step index is a DIFFERENT
  // namespace from the per-packet event_counter (creation draws are e = 0
  // of the packet's own key) — never conflate the two.
  const std::uint32_t n = static_cast<std::uint32_t>(*transport_step_count - 1);
  const double dt = cctk_delta_time; // coordinate, code units
  const double E_p_MeV = E_p;

  // The analytic fixture bin grid, through the sanctioned grid helpers.
  constexpr double nodes_MeV[emission_num_nodes] = {5.0, 15.0, 25.0};
  BinTable bins{};
  for (int b = 0; b < emission_num_bins; ++b) {
    const EnergyBin bin = energy_bin(nodes_MeV, emission_num_nodes, b);
    bins.lo[b] = bin.E_lo;
    bins.hi[b] = bin.E_hi;
  }

  const MetricGroups mgroups = metric_groups();
  const HydroGroups hgroups = hydro_groups();
  const SourceGroups sgroups = source_groups();

  const int myproc = amrex::ParallelDescriptor::MyProc();

  // Per-step event-side audit ([MCNX-HYD-05]): 10 device accumulator slots
  // (net dPt, dPx, dPy, dPz, dL; gross |dPt|, ..., |dL|), copied to the
  // host accumulator after the walk. Step-reset here.
  emission_step_audit() = LedgerAudit{};
  amrex::Gpu::DeviceVector<double> audit_dev(10, 0.0);
  double *const audit_ptr = audit_dev.data();

  // stats-emission accumulator ([MCNX-VER-07]): step-reset always (stale
  // values must never leak into MCNuX_StatsEmission), device slots only when
  // the gate is on.
  for (int i = 0; i < stats_num_slots; ++i)
    g_stats_accum[i] = 0.0;
  amrex::Gpu::DeviceVector<double> stats_dev(
      test_stats_emission ? stats_num_slots : 0, 0.0);
  double *const stats_ptr = test_stats_emission ? stats_dev.data() : nullptr;

  long total_created = 0;

  for (int patch = 0; patch < num_packet_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    PacketContainer &pc = packet_population(patch);
    for (const auto &leveldata : patchdata.leveldata) {
      const int level = leveldata.level;
      if (level > pc.finestLevel())
        break;
      const amrex::Geometry &geom = patchdata.amrcore->Geom(level);

      // Coordinate cell volume (code units) — the deposit/audit
      // normalization pair ([MCNX-HYD-02]: per unit coordinate volume per
      // unit coordinate time).
      const double dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
      // Count-law unit system (the one documented choice the packet_count
      // comment requires): dV -> cm^3 and dt -> s via the pinned
      // mcnux_units.hxx factors, so eta_b stays in its Kirchhoff
      // microphysics units (bin-integrated MeV cm^-3 s^-1) and E_p in MeV;
      // the product N_p is dimensionless.
      const double dV_cm3 = dV * length_code_to_cgs * length_code_to_cgs *
                            length_code_to_cgs;
      const double dt_s = dt * time_code_to_cgs;

      // The level's domain lower corner: the (i0, j0, k0) anchor of the
      // [MCNX-PKT-05] cell key.
      const amrex::IntVect dom_lo = geom.Domain().smallEnd();
      const int di0 = dom_lo[0], dj0 = dom_lo[1], dk0 = dom_lo[2];

      // Emission walks GRID boxes (tiling off) of the HydroBaseX rho group
      // MultiFab — the same BoxArray/DistributionMapping the ParGDB tracks,
      // so the MFIter box index addresses the matching particle tile (the
      // box-index reasoning of mcnux_gather.hxx).
      constexpr int tl = 0;
      amrex::MultiFab &rho_mf = *leveldata.groupdata.at(hgroups.rho)->mfab.at(tl);
      amrex::MultiFab &rad_mf =
          *leveldata.groupdata.at(sgroups.rad_force)->mfab.at(tl);
      amrex::MultiFab &lep_mf =
          *leveldata.groupdata.at(sgroups.lepton_source)->mfab.at(tl);

      for (amrex::MFIter mfi(rho_mf, false); mfi.isValid(); ++mfi) {
        const int box = mfi.index();
        const amrex::Box bx = mfi.validbox(); // valid cells only, no ghosts

        // Capacity guard of the cell-key packing, checked loudly on the box
        // corners before any packing ([MCNX-PKT-05]: overflow requires a
        // deliberate re-pin, never a silent wrap).
        if (!cell_key_capacity_ok(patch, level, bx.smallEnd(0), bx.smallEnd(1),
                                  bx.smallEnd(2), di0, dj0, dk0) ||
            !cell_key_capacity_ok(patch, level, bx.bigEnd(0), bx.bigEnd(1),
                                  bx.bigEnd(2), di0, dj0, dk0))
          CCTK_VERROR(
              "MCNuX emission cell-key capacity exceeded on patch %d level %d "
              "box %d (specs/packet-representation-and-sampling.md "
              "[MCNX-PKT-05]: patch < 128, level < 256, per-dimension offsets "
              "< 65536)",
              patch, level, box);

        const VertexMetricGather mgather =
            make_gather(patchdata, leveldata, mgroups, box);
        const CellFluidGather fgather =
            make_fluid_gather(patchdata, leveldata, hgroups, box);
        const SourceViews views{rad_mf.array(box, 0), rad_mf.array(box, 1),
                                rad_mf.array(box, 2), rad_mf.array(box, 3),
                                lep_mf.array(box, 0)};

        // Flattened (cell) x (species) x (bin) index space of both passes:
        // m = (icell * NUM_SPECIES + s) * nbins + b, icell in x-fastest
        // order over the valid box.
        const int ilo = bx.smallEnd(0), jlo = bx.smallEnd(1),
                  klo = bx.smallEnd(2);
        const int nx = bx.length(0), ny = bx.length(1);
        const amrex::Long M =
            bx.numPts() * amrex::Long(NUM_SPECIES) * emission_num_bins;

        // -------------------------------------------------------------
        // Count pass: one emission_count per (cell, s, b) — exactly one
        // cell_count_uniform draw and one bin_integrated_eta evaluation
        // per (species, bin, cell, step), no other RNG use.
        // -------------------------------------------------------------
        amrex::Gpu::DeviceVector<amrex::Long> counts(M);
        amrex::Long *const counts_ptr = counts.data();

        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE(amrex::Long m) noexcept {
          const int b = int(m % emission_num_bins);
          const amrex::Long m1 = m / emission_num_bins;
          const int sidx = int(m1 % NUM_SPECIES);
          const amrex::Long icell = m1 / NUM_SPECIES;
          const int i = ilo + int(icell % nx);
          const int j = jlo + int((icell / nx) % ny);
          const int k = klo + int(icell / (amrex::Long(nx) * ny));
          const Species s = static_cast<Species>(sidx);

          const int ic[3] = {i, j, k};
          const FluidSample fs = fgather.at_cell(ic);
          // sqrt(-g) = alpha sqrt(det gamma) at the cell CENTER (the
          // count law is a per-cell statement).
          const MetricSnapshot ms = mgather(fgather.cell_center(i, 0),
                                            fgather.cell_center(j, 1),
                                            fgather.cell_center(k, 2));
          const InverseSpatialMetric gu = spatial_metric_inverse(ms.g);
          const double sqrt_neg_g = ms.alpha * std::sqrt(gu.det);

          const EnergyBin bin{bins.lo[b], bins.hi[b]};
          const double eta_spectral =
              evaluate_coefficients(src, ap, table_eval, s,
                                    bin_center_energy(bin), fs.state)
                  .eta;
          const double eta_b = bin_integrated_eta(eta_spectral, bin_width(bin));
          const double N_p =
              packet_count(s, sqrt_neg_g, dV_cm3, dt_s, eta_b, E_p_MeV);
          const std::uint64_t K = cell_key(patch, level, i, j, k, di0, dj0, dk0);
          counts_ptr[m] = emission_count(N_p, cell_count_uniform(S, K, n, b, s));

          // stats-emission ([MCNX-VER-07]): expected count = the PRE-floor
          // real-valued N_p; count variance = p(1 - p) of the Bernoulli
          // remainder, p computed with emission_count's own truncation.
          if (stats_ptr != nullptr) {
            const double p =
                N_p - static_cast<double>(static_cast<std::int64_t>(N_p));
            amrex::Gpu::Atomic::AddNoRet(&stats_ptr[stats_stride * sidx + 0],
                                         N_p);
            amrex::Gpu::Atomic::AddNoRet(&stats_ptr[stats_stride * sidx + 1],
                                         p * (1.0 - p));
          }
        });

        // -------------------------------------------------------------
        // Scan -> offsets + total; id reservation; tile resize.
        // -------------------------------------------------------------
        amrex::Gpu::DeviceVector<amrex::Long> offsets(M);
        const amrex::Long total =
            amrex::Scan::ExclusiveSum(M, counts.data(), offsets.data());
        if (total == 0)
          continue;

        // Contiguous id range [base_id, base_id + total): NextID is
        // monotonic from 1, aborts past LastParticleID = 2^39 - 3, never
        // reused — per-rank unique for the whole run ([MCNX-GPU-04]).
        const amrex::Long base_id = PacketContainer::ParticleType::NextID();
        PacketContainer::ParticleType::NextID(base_id + total);

        auto &tile =
            pc.DefineAndReturnParticleTile(level, box, mfi.LocalTileIndex());
        const amrex::Long old_np = amrex::Long(tile.size());
        tile.resize(old_np + total);
        const auto ptd = tile.getParticleTileData();

        const amrex::Long *const counts_c = counts.data();
        const amrex::Long *const offs_c = offsets.data();

        // -------------------------------------------------------------
        // Fill pass: per created packet, the pinned six-draw creation map
        // and the [MCNX-PKT-03/04] initialization chain; deposit + audit.
        // -------------------------------------------------------------
        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE(amrex::Long m) noexcept {
          const amrex::Long cnt = counts_c[m];
          if (cnt == 0)
            return;
          const int b = int(m % emission_num_bins);
          const amrex::Long m1 = m / emission_num_bins;
          const int sidx = int(m1 % NUM_SPECIES);
          const amrex::Long icell = m1 / NUM_SPECIES;
          const int i = ilo + int(icell % nx);
          const int j = jlo + int((icell / nx) % ny);
          const int k = klo + int(icell / (amrex::Long(nx) * ny));
          const Species s = static_cast<Species>(sidx);

          const int ic[3] = {i, j, k};
          const FluidSample fs = fgather.at_cell(ic);

          const EnergyBin bin{bins.lo[b], bins.hi[b]};
          // [MCNX-PKT-03]: nu is the exact bin center, in MeV; the tetrad
          // transform and the stored p_i are geometrized code units, so the
          // fluid-frame energy fed to the transform is converted MeV ->
          // code HERE (the pinned mcnux_units.hxx factor; missing this
          // conversion corrupts every momentum by ~1e60).
          const double nu_MeV = bin_center_energy(bin);
          const double nu_code = nu_MeV / energy_code_to_MeV;
          // Weight N = E_p/nu, both in MeV ([MCNX-PKT-03]; the nu_x g = 4
          // entered ONCE, in packet_count — never here).
          const double N = packet_weight(E_p_MeV, nu_MeV);

          for (amrex::Long jp = 0; jp < cnt; ++jp) {
            const amrex::Long ord = offs_c[m] + jp; // packet ordinal in box
            const amrex::Long slot = old_np + ord;  // tile slot
            const std::uint64_t id = std::uint64_t(base_id + ord);
            // The per-packet draw key: the REPACKED logical key of
            // packet_rng_key (bit 63 provably clear), NOT the raw idcpu
            // word (see the mcnux_emission.hxx rationale). Computable
            // before any draw, as the e = 0 creation map requires.
            const std::uint64_t q = packet_rng_key(id, std::uint64_t(myproc));
            const CreationUniforms cu = creation_uniforms(S, q);

            // Position: uniform in the creation cell, per axis
            // (cell_lo = prob_lo + i dx, the CellFluidGather convention).
            const double x = creation_position(
                cu.u_pos_x, fgather.prob_lo[0] + i * fgather.dx[0],
                fgather.dx[0]);
            const double y = creation_position(
                cu.u_pos_y, fgather.prob_lo[1] + j * fgather.dx[1],
                fgather.dx[1]);
            const double z = creation_position(
                cu.u_pos_z, fgather.prob_lo[2] + k * fgather.dx[2],
                fgather.dx[2]);
            // cu.u_time is consumed by the pinned draw map (k = 3) but its
            // value is currently unused: the packet state stores no
            // creation time, and the full-dt push begins next step.

            // Direction/momentum chain of [MCNX-PKT-04], metric sampled at
            // the packet's creation position.
            const UnitVector3 nhat = fluid_frame_direction(cu.u1, cu.u2);

            // stats-emission ([MCNX-VER-07]): realized count and the angular
            // first moments, accumulated HERE because cos(theta)/phi exist
            // only at draw time. nhat.z IS cos(theta) = 2 u1 - 1 bitwise
            // (fluid_frame_direction returns it as the z component); phi is
            // recomputed by the identical 2 pi u2 expression of
            // fluid_frame_direction (it is never returned).
            if (stats_ptr != nullptr) {
              amrex::Gpu::Atomic::AddNoRet(&stats_ptr[stats_stride * sidx + 2],
                                           1.0);
              amrex::Gpu::Atomic::AddNoRet(&stats_ptr[stats_stride * sidx + 3],
                                           nhat.z);
              amrex::Gpu::Atomic::AddNoRet(&stats_ptr[stats_stride * sidx + 4],
                                           2.0 * detail::pi * cu.u2);
            }
            const MetricSnapshot mm = mgather(x, y, z);
            double u4[4];
            valencia_four_velocity(mm.alpha, mm.beta, mm.g, fs.vel, u4);
            const Tetrad tet = build_tetrad(mm.alpha, mm.beta, mm.g, u4);
            const SpacetimeMetric gm =
                spacetime_metric_from_adm(mm.alpha, mm.beta, mm.g);
            const CoordinateMomentum cm =
                transform_to_coordinate_frame(tet, gm, nu_code, nhat);

            // Packet state ([MCNX-GPU-03] schema): position, LOWER spatial
            // momentum (p^t is never stored — the null closure recomputes
            // it), weight, species, event counter 0. cpu = MyProc(), NOT 0:
            // cross-rank id uniqueness ([MCNX-GPU-04]).
            ptd.rdata(PIdx::x)[slot] = x;
            ptd.rdata(PIdx::y)[slot] = y;
            ptd.rdata(PIdx::z)[slot] = z;
            ptd.rdata(PIdx::px)[slot] = cm.px;
            ptd.rdata(PIdx::py)[slot] = cm.py;
            ptd.rdata(PIdx::pz)[slot] = cm.pz;
            ptd.rdata(PIdx::w)[slot] = N;
            ptd.idata(IntIdx::species)[slot] = sidx;
            ptd.idata(IntIdx::event_counter)[slot] = 0;
            ptd.m_idcpu[slot] =
                amrex::SetParticleIDandCPU(amrex::Long(id), myproc);

            // Deposit + audit from the SAME pre-negation LedgerDelta, with
            // the SAME code-unit (dV, dt) pair ([MCNX-HYD-05] pairing).
            const LedgerDelta d =
                emission_delta(s, N, cm.pt, cm.px, cm.py, cm.pz);
            deposit_delta(views, i, j, k, d, dV, dt);
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
          }
        });
        amrex::Gpu::streamSynchronize();

        total_created += long(total);
      }
    }
  }

  // Creation writes into the owning cell's own tile; Redistribute() is the
  // no-op ownership guard (the push precedent). Unlike the push/episode
  // Redistribute() sites this one is PROVABLY never an escape site
  // ([MCNX-GPU-05], mcnux_escape.hxx): creation_position maps u through
  // map_uniform_to_interval into the walked cell's own half-open
  // [cell_lo, cell_lo + dx) interval (mcnux_emission.hxx), so every created
  // packet starts strictly inside the domain — no escape detection here.
  for (int patch = 0; patch < num_packet_patches(); ++patch)
    packet_population(patch).Redistribute();

  // Fold the device audit into the step accumulator.
  std::vector<double> audit_host(10, 0.0);
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, audit_dev.begin(), audit_dev.end(),
                   audit_host.begin());
  LedgerAudit &audit = emission_step_audit();
  audit.net = LedgerDelta{audit_host[0], audit_host[1], audit_host[2],
                          audit_host[3], audit_host[4]};
  audit.gross = LedgerDelta{audit_host[5], audit_host[6], audit_host[7],
                            audit_host[8], audit_host[9]};

  // Fold the device stats accumulator into the step accumulator (same
  // sync-then-copy discipline as the audit fold above).
  if (test_stats_emission) {
    std::vector<double> stats_host(stats_num_slots, 0.0);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, stats_dev.begin(),
                     stats_dev.end(), stats_host.begin());
    for (int i = 0; i < stats_num_slots; ++i)
      g_stats_accum[i] = stats_host[i];
  }

  CCTK_VINFO("MCNuX emission step %u: created %ld packets", unsigned(n),
             total_created);
}

// ---------------------------------------------------------------------------
// Emission diagnostics  (MCNuX::mcnux_packet_diag, MCNuX::mcnux_emission_diag)
// ---------------------------------------------------------------------------
// Mirrors the created population into the per-packet golden tables (the
// shared fill_packet_diag of mcnux_gather.hxx / mcnux_geodesic.cxx: row =
// packet id - 1, hard CCTK_VERROR above the declared SIZE — the emit-smoke /
// emission-fixedseed benchmarks tune eta_scale so the total emitted count
// stays within it). The AT initial leg fills the (empty-population) all-zero
// tables so iteration-0 output is defined.
//
// The discrete identity half (id, weight, species, event counter — the
// "exact packet set (ids, states)" observable of the `emission-fixedseed`
// row, specs/verification-suite-design.md:243) goes into the sibling group
// mcnux_emission_diag via the local staging kernel below, modeled on
// fill_packet_diag: device buffer keyed by unpack_id - 1, zero-fill first,
// hard CCTK_VERROR past the declared SIZE, host copy-out. Single-rank
// assumption as in fill_packet_diag (every MCNuX test runs NPROCS 1).

namespace {

// Number of doubles per packet in the emission-identity staging buffer:
// id, weight, species, event counter.
constexpr int emdiag_width = 4;

struct EmissionDiagColumns {
  CCTK_REAL *id, *w, *species, *ec;
};

inline int emission_diag_size() {
  return diag_array_size("MCNuX::mcnux_emission_diag");
}

void fill_emission_diag(const EmissionDiagColumns &cols, const int nrows) {
  for (int r = 0; r < nrows; ++r) {
    cols.id[r] = 0.0;
    cols.w[r] = 0.0;
    cols.species[r] = 0.0;
    cols.ec[r] = 0.0;
  }

  for_each_packet_tile_raw([&](const CarpetX::GHExt::PatchData &,
                               const CarpetX::GHExt::PatchData::LevelData &,
                               PacketContainer &, const PacketIter &pti) {
    const long np = pti.numParticles();
    if (np == 0)
      return;
    const auto ptd = pti.GetParticleTile().getParticleTileData();

    amrex::Gpu::DeviceVector<double> staging(std::size_t(np) * emdiag_width);
    double *const out = staging.data();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      double *const row = out + ip * emdiag_width;
      row[0] = double(amrex::particle_impl::unpack_id(ptd.m_idcpu[ip]));
      row[1] = ptd.rdata(PIdx::w)[ip];
      row[2] = double(ptd.idata(IntIdx::species)[ip]);
      row[3] = double(ptd.idata(IntIdx::event_counter)[ip]);
    });
    amrex::Gpu::streamSynchronize();

    std::vector<double> host(std::size_t(np) * emdiag_width);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, staging.begin(), staging.end(),
                     host.begin());

    for (long ip = 0; ip < np; ++ip) {
      const double *const row = host.data() + ip * emdiag_width;
      const long id = long(row[0]);
      if (id < 1 || id > nrows)
        CCTK_VERROR("MCNuX packet id %ld is outside the emission diagnostic "
                    "table (MCNuX::mcnux_emission_diag has SIZE=%d rows)",
                    id, nrows);
      const int r = int(id - 1);
      cols.id[r] = row[0];
      cols.w[r] = row[1];
      cols.species[r] = row[2];
      cols.ec[r] = row[3];
    }
  });
}

} // namespace

extern "C" void MCNuX_EmissionDiag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_EmissionDiag;
  DECLARE_CCTK_PARAMETERS;

  fill_packet_diag(metric_groups(), {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz,
                                     pk_pt},
                   packet_diag_size());
  fill_emission_diag({em_id, em_w, em_species, em_ec}, emission_diag_size());
}

// ---------------------------------------------------------------------------
// Statistical-reduction writer  (MCNuX::mcnux_stats_diag)
// ---------------------------------------------------------------------------
// The first stats-tier writer ([MCNX-PKT-06], [MCNX-VER-07]; the
// `stats-emission` benchmark of specs/verification-suite-design.md:255-279):
// reduce the emission step into per-check (estimate, expected, sigma, z)
// rows through the frozen MCNuX::zscore() of mcnux_stats.hxx (called
// verbatim, never reimplemented) and mirror them into the
// MCNuX::mcnux_stats_diag table for golden output. Acceptance |z| <= 4 is
// judged once at golden capture; thereafter the archived z values are
// themselves golden numbers diffed at 1e-12.
//
// Row layout (SIZE=14), species blocks in the binding species_index order
// (0 nu_e, 1 nu_e_bar, 2 nu_x; mcnux_units.hxx):
//   4s + 0  realized created-packet count vs Sum N_p (pre-floor exact
//           counts), sigma = sqrt(Sum p(1 - p)) — the [MCNX-PKT-02] floor +
//           Bernoulli construction;
//   4s + 1  total emitted fluid-frame energy: estimate = E_p * realized
//           count (each packet carries w * nu = E_p exactly,
//           packet_weight of mcnux_emission.hxx), expected = E_p * Sum N_p
//           (the count-law identity E_p N_p = g_s sqrt(-g) dV dt eta_b —
//           the nu_x g = 4 leg of [MCNX-CNV-06] is carried through
//           degeneracy() inside the accumulated N_p), sigma = E_p *
//           sigma_count — so the energy and count z rows are exactly
//           proportional by construction;
//   4s + 2  cos(theta) sample mean vs 0, sigma = sqrt(1/(3 n_s)) with n_s
//           the realized count (documented sample-size choice);
//   4s + 3  phi sample mean vs pi, sigma = pi sqrt(1/(3 n_s));
//   12      RNG uniformity sample mean vs 1/2, sigma = 1/sqrt(12 N)
//           (spec-stated, rng-and-statistical-acceptance.md:118);
//   13      RNG uniformity second moment vs 1/3, sigma = sqrt(4/(45 N)).
// Sigma formulas are benchmark-owned (mcnux_stats_emission.hxx), never
// added to mcnux_stats.hxx (its binding non-goal).
//
// The RNG rows sweep the pure u(S, q, e, k) on the host over the documented
// tuple sweep q = 1..N at fixed e = 0, k = 0, with the pinned
// S = stats_seed_primary = 1296518744 and N = stats_default_num_packets =
// 2^20 — independent of the run's own rng_seed and packet draws, so they
// are computed identically on the AT initial leg. The packet rows are
// zero-filled at initial (no emission has run; deterministic iteration-0
// output) and computed from g_stats_accum in-step, where a species with
// zero realized packets is a hard error (its angular sigmas would be
// meaningless at sigma = 0; the benchmark is tuned so all three species
// emit).

namespace {

inline int stats_diag_size() {
  return diag_array_size("MCNuX::mcnux_stats_diag");
}

} // namespace

extern "C" void MCNuX_StatsEmission(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_StatsEmission;
  DECLARE_CCTK_PARAMETERS;

  constexpr int expected_rows = 14;
  const int nrows = stats_diag_size();
  if (nrows != expected_rows)
    CCTK_VERROR("MCNuX::mcnux_stats_diag has SIZE=%d but MCNuX_StatsEmission "
                "writes %d check rows; the two must agree "
                "(MCNuX/interface.ccl)",
                nrows, expected_rows);

  for (int r = 0; r < nrows; ++r) {
    stats_estimate[r] = 0.0;
    stats_expected[r] = 0.0;
    stats_sigma[r] = 0.0;
    stats_z[r] = 0.0;
  }

  const auto put = [&](int r, const ZScore &zs) {
    stats_estimate[r] = zs.estimate;
    stats_expected[r] = zs.expected;
    stats_sigma[r] = zs.sigma;
    stats_z[r] = zs.z;
    CCTK_VINFO("MCNuX stats-emission row %2d: estimate = %.17g, "
               "expected = %.17g, sigma = %.17g, z = %.17g  %s",
               r, zs.estimate, zs.expected, zs.sigma, zs.z,
               zs.pass ? "PASS" : "FAIL");
  };

  // Packet rows 0..11: in-step only (the AT initial leg keeps the
  // deterministic zero fill above — no emission has run).
  if (cctk_iteration > 0) {
    const double E_p_MeV = E_p;
    for (int s = 0; s < NUM_SPECIES; ++s) {
      const double exp_count = g_stats_accum[stats_stride * s + 0];
      const double var_count = g_stats_accum[stats_stride * s + 1];
      const double n_s = g_stats_accum[stats_stride * s + 2];
      const double sum_cos = g_stats_accum[stats_stride * s + 3];
      const double sum_phi = g_stats_accum[stats_stride * s + 4];

      if (n_s <= 0.0)
        CCTK_VERROR(
            "MCNuX stats-emission: species %d realized zero packets this "
            "step; its angular check rows are meaningless at sigma = 0 "
            "([MCNX-VER-07]). Tune the benchmark (E_p / eta_scale) so every "
            "species emits.",
            s);

      const double sigma_count = sigma_bernoulli_sum(var_count);
      put(4 * s + 0, zscore(n_s, exp_count, sigma_count));
      put(4 * s + 1, zscore(E_p_MeV * n_s, E_p_MeV * exp_count,
                            E_p_MeV * sigma_count));
      put(4 * s + 2, zscore(sum_cos / n_s, 0.0, sigma_isotropy_costheta(n_s)));
      put(4 * s + 3,
          zscore(sum_phi / n_s, detail::pi, sigma_isotropy_phi(n_s)));
    }
  }

  // RNG uniformity-at-the-pinned-bar rows 12..13
  // (rng-and-statistical-acceptance.md:118): host sweep q = 1..N at e = 0,
  // k = 0, pinned primary seed, N = 2^20.
  {
    const std::uint64_t S = stats_seed_primary;
    const std::int64_t N = stats_default_num_packets;
    double sum = 0.0, sum2 = 0.0;
    for (std::int64_t qi = 1; qi <= N; ++qi) {
      const double r = u(S, static_cast<std::uint64_t>(qi), 0u, 0u);
      sum += r;
      sum2 += r * r;
    }
    const double Nd = static_cast<double>(N);
    put(12, zscore(sum / Nd, 0.5, sigma_uniform_mean(Nd)));
    put(13, zscore(sum2 / Nd, 1.0 / 3.0, sigma_uniform_second_moment(Nd)));
  }
}

} // namespace MCNuX
