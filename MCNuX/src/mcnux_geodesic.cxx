// Runtime side of specs/geodesic-propagation.md: the packet population owner
// (the thorn's first live PacketContainer), the deterministic synthetic
// packet-seeding fixture, the scheduled geodesic push, and the per-packet
// diagnostic table that the `minkowski-freestream` benchmark golden-checks.
//
// Driver access follows the [MCNX-CTX-01] idiom of mcnux_cadence.cxx: the
// per-patch AmrCore and per-level group MultiFabs are reached through the
// CarpetX singleton `ghext`; grid data is read only through amrex::Array4
// views of the ADMBaseX groups declared in the routines' READS: clauses
// ([MCNX-CTX-04]). The particle operators here never use the grid-point loop
// machinery (loop_device.hxx): per-packet work is one amrex::ParallelFor
// over a particle tile capturing trivially-copyable values ([MCNX-GPU-02],
// the packet-kernel idiom documented in mcnux_particles.hxx).
//
// Scheduling ([MCNX-CTX-03]): both routines run in global/level mode —
// CarpetX invokes such a routine exactly once per traversal — and walk the
// patches, levels, and particle tiles themselves.

// The driver-access plumbing (require_driver, group lookups, gather
// constructors, tile walkers, the diagnostic-table declaration) lives in
// mcnux_gather.hxx (the shared C1 factoring; it performs the CarpetX
// relative include of driver.hxx itself). This file remains the owner of
// the packet population and of the fill_packet_diag definition.
#include "mcnux_gather.hxx"

#include "mcnux_escape.hxx"
#include "mcnux_fluid.hxx"
#include "mcnux_geodesic.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_tetrad.hxx"

#include <AMReX_GpuAtomic.H>

#include <AMReX_AmrParGDB.H> // completes AmrParGDB for the ParGDBBase* ctor
#include <AMReX_GpuContainers.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleTransformation.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdint>
#include <memory>
#include <vector>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Population owner
// ---------------------------------------------------------------------------

namespace {

// One container per patch, constructed lazily on the patch's ParGDB
// (amrex::ParticleContainer_impl(ParGDBBase*) is the only inherited
// constructor that tracks the driver's grid structure; the AmrCore*
// convenience constructor does not exist for pure-SoA containers).
std::vector<std::unique_ptr<PacketContainer> > populations;

} // namespace

int num_packet_patches() {
  require_driver();
  return int(CarpetX::ghext->patchdata.size());
}

PacketContainer &packet_population(const int patch) {
  require_driver();
  if (populations.size() != CarpetX::ghext->patchdata.size())
    populations.resize(CarpetX::ghext->patchdata.size());
  std::unique_ptr<PacketContainer> &slot = populations.at(patch);
  if (!slot) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    // A container built before the grids exist would size its particle
    // levels from finestLevel() == -1 (AMReX_ParticleContainer.H, the
    // resizeData() note); the first touch is AT initial, never earlier.
    if (!patchdata.amrcore || patchdata.amrcore->finestLevel() < 0)
      CCTK_VERROR("MCNuX packet population for patch %d requested before the "
                  "CarpetX grid hierarchy exists",
                  patch);
    slot = std::make_unique<PacketContainer>(patchdata.amrcore->GetParGDB());
  }
  return *slot;
}

// ---------------------------------------------------------------------------
// Diagnostic table  (MCNuX::mcnux_packet_diag)
// ---------------------------------------------------------------------------
// The group lookups, gather constructors, and tile walkers used below are
// the shared ones of mcnux_gather.hxx (C1 factoring — one copy, reused by
// the emission loop).

namespace {

// Number of doubles per packet in the device staging buffer:
// id, x, y, z, p_x, p_y, p_z, p_t.
constexpr int diag_width = 8;

} // namespace

// Fill the diagnostic table from the current population: row = packet id - 1
// (so the table is independent of tile order and iteration), columns the
// position, the lower momentum, and the covariant energy
// p_t = -alpha^2 p^t + beta^i p_i (geodesic-propagation.md:121-123) with p^t
// from the null closure on the gathered metric at the packet position. The
// per-tile values are produced by a device kernel into a staging buffer and
// copied to the host. Single-rank assumption: the host array is filled from
// the local tiles only (every MCNuX test runs NPROCS 1; a multi-rank table
// would need a reduction over ranks).
void fill_packet_diag(const MetricGroups &groups, const PacketDiagColumns &cols,
                      const int nrows) {
  for (int r = 0; r < nrows; ++r) {
    cols.x[r] = 0.0;
    cols.y[r] = 0.0;
    cols.z[r] = 0.0;
    cols.px[r] = 0.0;
    cols.py[r] = 0.0;
    cols.pz[r] = 0.0;
    cols.pt[r] = 0.0;
  }

  for_each_packet_tile(groups, [&](PacketContainer &, int, const PacketIter &pti,
                                   const VertexMetricGather gather) {
    const long np = pti.numParticles();
    if (np == 0)
      return;
    const auto ptd = pti.GetParticleTile().getParticleTileData();

    amrex::Gpu::DeviceVector<double> staging(std::size_t(np) * diag_width);
    double *const out = staging.data();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      const double x[3] = {ptd.rdata(PIdx::x)[ip], ptd.rdata(PIdx::y)[ip],
                           ptd.rdata(PIdx::z)[ip]};
      const double p[3] = {ptd.rdata(PIdx::px)[ip], ptd.rdata(PIdx::py)[ip],
                           ptd.rdata(PIdx::pz)[ip]};
      const MetricSnapshot m = gather(x[0], x[1], x[2]);
      const InverseSpatialMetric gu = spatial_metric_inverse(m.g);
      const double pt_up = p_t_closure(p[0], p[1], p[2], gu, m.alpha);
      const double pt_lo = -m.alpha * m.alpha * pt_up + m.beta[0] * p[0] +
                           m.beta[1] * p[1] + m.beta[2] * p[2];
      double *const row = out + ip * diag_width;
      row[0] = double(amrex::particle_impl::unpack_id(ptd.m_idcpu[ip]));
      row[1] = x[0];
      row[2] = x[1];
      row[3] = x[2];
      row[4] = p[0];
      row[5] = p[1];
      row[6] = p[2];
      row[7] = pt_lo;
    });
    amrex::Gpu::streamSynchronize();

    std::vector<double> host(std::size_t(np) * diag_width);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, staging.begin(), staging.end(),
                     host.begin());

    for (long ip = 0; ip < np; ++ip) {
      const double *const row = host.data() + ip * diag_width;
      const long id = long(row[0]);
      if (id < 1 || id > nrows)
        CCTK_VERROR("MCNuX packet id %ld is outside the diagnostic table "
                    "(MCNuX::mcnux_packet_diag has SIZE=%d rows)",
                    id, nrows);
      const int r = int(id - 1);
      cols.x[r] = row[1];
      cols.y[r] = row[2];
      cols.z[r] = row[3];
      cols.px[r] = row[4];
      cols.py[r] = row[5];
      cols.pz[r] = row[6];
      cols.pt[r] = row[7];
    }
  });
}

// ---------------------------------------------------------------------------
// Synthetic packet fixture
// ---------------------------------------------------------------------------

// Two deterministic seed sets, selected by MCNuX::synthetic_packet_fixture
// (exact binary fractions in both; row index in the diagnostic table =
// packet id - 1 = fixture index; function-local constexpr per the
// mcnux_srcterms.cxx CUDA note on namespace-scope struct constants):
//   * "minkowski" (the `minkowski-freestream` benchmark): eight packets
//     inside |x^i| <= 0.25 with momenta spanning the six axis directions
//     (at four different magnitudes — the coordinate speed must be 1
//     regardless) plus two oblique directions; over the benchmark's
//     4 x 0.125 coordinate time every packet stays well inside the
//     [-1, 1]^3 domain, so the escape policy is never exercised. The
//     `escape-freestream` benchmark reuses the SAME fixture over 8 steps,
//     where seven of the eight packets exit and the [MCNX-GPU-05] escape
//     tally IS the observable.
//   * "schwarzschild" (the `schwarzschild-pt` p_t-drift legs): eight
//     packets at isotropic r in [6.06, 9.53] M with dominant
//     outgoing-radial momenta; over the benchmark's T = 100 M none escapes
//     the [2, 130]^3 domain (final radii ~100-105 with every coordinate
//     inside the box, measured at capture).
namespace {

struct PacketFixture {
  double x, y, z;
  double px, py, pz;
};

} // namespace

extern "C" void MCNuX_SeedSyntheticPackets(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SeedSyntheticPackets;
  DECLARE_CCTK_PARAMETERS;

  constexpr int nfixture = 8;
  constexpr PacketFixture fixture_minkowski[nfixture] = {
      {0.125, 0.0, 0.0, 1.0, 0.0, 0.0},
      {-0.125, 0.0, 0.0, -1.0, 0.0, 0.0},
      {0.0, 0.125, 0.0, 0.0, 2.0, 0.0},
      {0.0, -0.125, 0.0, 0.0, -0.5, 0.0},
      {0.0, 0.0, 0.125, 0.0, 0.0, 0.25},
      {0.0, 0.0, -0.125, 0.0, 0.0, -4.0},
      {0.25, -0.25, 0.125, 0.5, -0.25, 0.125},
      {-0.0625, 0.1875, -0.25, -0.375, 0.75, 0.5},
  };
  constexpr PacketFixture fixture_schwarzschild[nfixture] = {
      {3.5, 3.5, 3.5, 1.0, 1.0, 1.0},
      {4.5, 4.5, 4.5, 0.5, 0.5, 0.5},
      {5.5, 5.5, 5.5, 2.0, 2.0, 2.0},
      {6.0, 4.0, 3.0, 1.5, 1.0, 0.75},
      {8.0, 3.0, 3.0, 1.0, 0.5, 0.25},
      {4.0, 6.0, 3.0, 0.5, 1.5, 0.5},
      {3.0, 4.0, 7.0, 0.25, 0.5, 1.0},
      {5.0, 5.0, 2.5, 1.0, 1.0, 0.5},
  };
  const PacketFixture *const fixture =
      CCTK_EQUALS(synthetic_packet_fixture, "schwarzschild")
          ? fixture_schwarzschild
          : fixture_minkowski;

  const int nrows = packet_diag_size();
  if (nrows != nfixture)
    CCTK_VERROR("MCNuX synthetic packet fixture has %d packets, but "
                "MCNuX::mcnux_packet_diag is declared with SIZE=%d in "
                "MCNuX/interface.ccl; the two must agree",
                nfixture, nrows);

  const MetricGroups groups = metric_groups();

  // Seed into patch 0, level 0, (grid 0, tile 0) on rank 0 through a
  // host-side (pinned) tile, then let Redistribute() assign owners. The
  // pinned-tile + copyParticles route is the portable form: on a GPU build
  // the container's SoA lives in device memory and cannot be written from
  // the host directly.
  PacketContainer &pc = packet_population(0);
  if (amrex::ParallelDescriptor::MyProc() == 0) {
    using PinnedTile =
        PacketContainer::ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;
    PinnedTile pinned;
    pinned.define(0, 0);
    pinned.resize(nfixture);
    const auto hd = pinned.getParticleTileData();
    for (int i = 0; i < nfixture; ++i) {
      const PacketFixture &f = fixture[i];
      hd.rdata(PIdx::x)[i] = f.x;
      hd.rdata(PIdx::y)[i] = f.y;
      hd.rdata(PIdx::z)[i] = f.z;
      hd.rdata(PIdx::px)[i] = f.px;
      hd.rdata(PIdx::py)[i] = f.py;
      hd.rdata(PIdx::pz)[i] = f.pz;
      hd.rdata(PIdx::w)[i] = 1.0;
      hd.idata(IntIdx::species)[i] = 0;
      hd.idata(IntIdx::event_counter)[i] = 0;
      hd.m_idcpu[i] = amrex::SetParticleIDandCPU(i + 1, 0);
    }
    auto &tile = pc.DefineAndReturnParticleTile(0, 0, 0);
    const auto old_np = tile.numParticles();
    tile.resize(old_np + nfixture);
    amrex::copyParticles(tile, pinned, 0, int(old_np), nfixture);
  }
  pc.Redistribute();

  CCTK_VINFO("MCNuX seeded %d synthetic packets (%ld in the population)",
             nfixture, long(pc.TotalNumberOfParticles()));

  fill_packet_diag(groups, {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz, pk_pt},
                   nrows);
}

// ---------------------------------------------------------------------------
// Beam seeding fixture  (the `stats-beam` benchmark)
// ---------------------------------------------------------------------------

// Deterministic collimated monoenergetic beam of the [MCNX-INT-01/03]
// beam-attenuation benchmark (specs/neutrino-matter-interactions.md:197-200;
// gated by MCNuX::test_stats_beam): beam_num_packets nu_e packets of weight
// 1.0 at the seed plane x = beam_x0, spread over a deterministic y-z lattice
// strictly inside the domain (index arithmetic only — NO RNG draw is
// consumed by seeding), all with the identical +x null momentum
// p_i = (E_code, 0, 0) (on the benchmark's Minkowski background the null
// closure gives p^t = E_code and coordinate speed 1). Ids are 1..N with
// cpu 0 (the synthetic-fixture idiom; MCNuX_ParamCheck forbids the beam
// together with enable_emission and test_synthetic_packets, so no id
// collision is reachable). Unlike the 8-row synthetic fixture this writes
// NO per-packet diagnostic table — the benchmark observable is the
// aggregate escape/absorption reduction of MCNuX_StatsBeam
// (mcnux_interactions.cxx), never a 2^20-row table.
extern "C" void MCNuX_SeedBeamPackets(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SeedBeamPackets;
  DECLARE_CCTK_PARAMETERS;

  require_driver();

  const long np = long(beam_num_packets);
  // Fluid-frame == coordinate-frame energy on the at-rest Minkowski
  // background: E in MeV -> code via the pinned factor (the emission loop's
  // nu_MeV / energy_code_to_MeV form; never a retyped literal).
  const double E_code = beam_energy_MeV / energy_code_to_MeV;

  // The seed plane and lattice live on the patch-0 level-0 geometry (the
  // benchmark is single-patch unigrid; Redistribute() would re-own packets
  // anyway on a multi-box layout).
  const auto &patchdata = CarpetX::ghext->patchdata.at(0);
  const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
  const double xlo = geom.ProbLo(0), xhi = geom.ProbHi(0);
  const double ylo = geom.ProbLo(1), yhi = geom.ProbHi(1);
  const double zlo = geom.ProbLo(2), zhi = geom.ProbHi(2);
  if (!(beam_x0 >= xlo) || !(beam_x0 < xhi))
    CCTK_VERROR("MCNuX::beam_x0 = %.17g is outside the domain [%.17g, %.17g) "
                "(the escape predicate is half-open: a seed ON the upper "
                "face would escape immediately)",
                double(beam_x0), xlo, xhi);

  // Smallest square lattice covering np: nyz = ceil(sqrt(np)) by integer
  // search (np = 2^20 gives nyz = 1024 exactly). Cell-center-style offsets
  // (j + 0.5)/nyz keep every packet strictly inside the y/z extent and off
  // every face for binary-exact domains.
  long nyz = 1;
  while (nyz * nyz < np)
    ++nyz;
  const double dy = (yhi - ylo) / double(nyz);
  const double dz = (zhi - zlo) / double(nyz);

  PacketContainer &pc = packet_population(0);
  if (amrex::ParallelDescriptor::MyProc() == 0) {
    using PinnedTile =
        PacketContainer::ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;
    PinnedTile pinned;
    pinned.define(0, 0);
    pinned.resize(np);
    const auto hd = pinned.getParticleTileData();
    for (long i = 0; i < np; ++i) {
      const long j = i % nyz;
      const long k = i / nyz;
      hd.rdata(PIdx::x)[i] = beam_x0;
      hd.rdata(PIdx::y)[i] = ylo + (double(j) + 0.5) * dy;
      hd.rdata(PIdx::z)[i] = zlo + (double(k) + 0.5) * dz;
      hd.rdata(PIdx::px)[i] = E_code;
      hd.rdata(PIdx::py)[i] = 0.0;
      hd.rdata(PIdx::pz)[i] = 0.0;
      hd.rdata(PIdx::w)[i] = 1.0;
      hd.idata(IntIdx::species)[i] = 0; // nu_e
      hd.idata(IntIdx::event_counter)[i] = 0;
      hd.m_idcpu[i] = amrex::SetParticleIDandCPU(i + 1, 0);
    }
    auto &tile = pc.DefineAndReturnParticleTile(0, 0, 0);
    const auto old_np = tile.numParticles();
    tile.resize(old_np + np);
    amrex::copyParticles(tile, pinned, 0, int(old_np), np);
  }
  pc.Redistribute();

  CCTK_VINFO("MCNuX seeded %ld beam packets at x = %.17g on a %ld x %ld "
             "y-z lattice, E = %.17g MeV (%ld in the population)",
             np, double(beam_x0), nyz, nyz, double(beam_energy_MeV),
             long(pc.TotalNumberOfParticles()));
}

// ---------------------------------------------------------------------------
// Scatterbox seeding fixture  (the `stats-scatterbox` benchmark)
// ---------------------------------------------------------------------------

// Deterministic monoenergetic packet lattice of the [MCNX-INT-02/04]
// collision-statistics benchmark (specs/neutrino-matter-interactions.md:
// 201-207; gated by MCNuX::test_stats_scatterbox), modeled line-for-line on
// MCNuX_SeedBeamPackets above: scatterbox_num_packets nu_e packets of
// weight 1.0, all with the identical +x null momentum p_i = (E_code, 0, 0)
// (on the benchmark's Minkowski background the null closure gives
// p^t = E_code and coordinate speed 1, so one transport step gives every
// packet path length exactly c dt = cctk_delta_time — the benchmark's
// pinned l_path). Ids are 1..N with cpu 0 (the synthetic-fixture idiom;
// MCNuX_ParamCheck forbids the scatterbox together with enable_emission,
// test_synthetic_packets, and test_stats_beam, so no id collision is
// reachable), and NO RNG draw is consumed by seeding.
//
// Geometry (the one change from the beam's y-z plane lattice): a 3D
// cell-center lattice n^3 with n = ceil(cbrt(np)) by integer search
// (np = 2^20 gives n = 102), filled by linear index
// i -> (i % n, (i/n) % n, i/n^2) over the inner cube [-h, h]^3 with
// h = scatterbox_half_width. The cube must be strictly inside the domain
// (hard error below); with h + c dt below the domain half-width no packet
// can reach a face within the single transport step (the maximum
// displacement is the path length c dt, direction changes only shorten the
// net displacement), so the run has ZERO escapes by construction — the
// writer MCNuX_StatsScatterbox hard-errors if any occur. No per-packet
// diagnostic table is written (the beam precedent: the observable is the
// aggregate statistical reduction, never a 2^20-row table).
extern "C" void MCNuX_SeedScatterboxPackets(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SeedScatterboxPackets;
  DECLARE_CCTK_PARAMETERS;

  require_driver();

  const long np = long(scatterbox_num_packets);
  // Fluid-frame == coordinate-frame energy on the at-rest Minkowski
  // background: E in MeV -> code via the pinned factor (the beam seeder's
  // form; never a retyped literal).
  const double E_code = scatterbox_energy_MeV / energy_code_to_MeV;

  // The lattice lives on the patch-0 level-0 geometry (the benchmark is
  // single-patch unigrid; Redistribute() would re-own packets anyway on a
  // multi-box layout).
  const auto &patchdata = CarpetX::ghext->patchdata.at(0);
  const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
  const double h = scatterbox_half_width;
  for (int d = 0; d < 3; ++d)
    if (!(-h > geom.ProbLo(d)) || !(h < geom.ProbHi(d)))
      CCTK_VERROR("MCNuX::scatterbox_half_width = %.17g puts the seeding "
                  "cube [-h, h]^3 outside the axis-%d domain extent "
                  "[%.17g, %.17g); the cube must be strictly inside the "
                  "domain (a seed on or beyond a face would escape "
                  "immediately — the escape predicate is half-open)",
                  h, d, geom.ProbLo(d), geom.ProbHi(d));

  // Smallest cubic lattice covering np: n = ceil(cbrt(np)) by integer
  // search (np = 2^20 gives n = 102: 101^3 = 1030301 < 2^20 <= 1061208 =
  // 102^3). Cell-center-style offsets (i + 0.5)/n keep every packet
  // strictly inside the cube and off its faces.
  long n = 1;
  while (n * n * n < np)
    ++n;
  const double cell = 2.0 * h / double(n);

  PacketContainer &pc = packet_population(0);
  if (amrex::ParallelDescriptor::MyProc() == 0) {
    using PinnedTile =
        PacketContainer::ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;
    PinnedTile pinned;
    pinned.define(0, 0);
    pinned.resize(np);
    const auto hd = pinned.getParticleTileData();
    for (long i = 0; i < np; ++i) {
      const long ix = i % n;
      const long iy = (i / n) % n;
      const long iz = i / (n * n);
      hd.rdata(PIdx::x)[i] = -h + (double(ix) + 0.5) * cell;
      hd.rdata(PIdx::y)[i] = -h + (double(iy) + 0.5) * cell;
      hd.rdata(PIdx::z)[i] = -h + (double(iz) + 0.5) * cell;
      hd.rdata(PIdx::px)[i] = E_code;
      hd.rdata(PIdx::py)[i] = 0.0;
      hd.rdata(PIdx::pz)[i] = 0.0;
      hd.rdata(PIdx::w)[i] = 1.0;
      hd.idata(IntIdx::species)[i] = 0; // nu_e
      hd.idata(IntIdx::event_counter)[i] = 0;
      hd.m_idcpu[i] = amrex::SetParticleIDandCPU(i + 1, 0);
    }
    auto &tile = pc.DefineAndReturnParticleTile(0, 0, 0);
    const auto old_np = tile.numParticles();
    tile.resize(old_np + np);
    amrex::copyParticles(tile, pinned, 0, int(old_np), np);
  }
  pc.Redistribute();

  CCTK_VINFO("MCNuX seeded %ld scatterbox packets on a %ld^3 lattice inside "
             "[-%.17g, %.17g]^3, E = %.17g MeV (%ld in the population)",
             np, n, h, h, double(scatterbox_energy_MeV),
             long(pc.TotalNumberOfParticles()));
}

// ---------------------------------------------------------------------------
// The geodesic push
// ---------------------------------------------------------------------------

// One transport step of free streaming for every packet: per tile, one
// device kernel applying geodesic_step_cellwise (RK4 on the [MCNX-GEO-01]
// RHS with the trilinear gather of [MCNX-GEO-03], sub-stepped at cell-face
// crossings so every RK4 piece integrates one smooth trilinear polynomial —
// see mcnux_geodesic.hxx) over the coarsest-level Delta t (cctk_delta_time
// in level/global mode, the cadence contract of mcnux_cadence.cxx).
//
// Escape handling ([MCNX-GPU-05], mcnux_escape.hxx): a packet whose
// post-step position satisfies the half-open outside_domain predicate is
// tallied per species (count, number Sum N, energy Sum N p^t — p^t the
// upper-index null closure with the FINAL momentum on the metric gathered
// at the PRE-step, last-in-domain position) into one step-scoped device
// buffer shared by every tile's kernel, then removed explicitly via
// make_invalid (the absorption idiom of mcnux_interactions.cxx; the id
// value stays retired, [MCNX-GPU-04]) — tally strictly BEFORE
// invalidation, so nothing is ever silently dropped inside Redistribute().
// The unconditional per-patch Redistribute() afterwards compacts the
// invalidated slots and re-establishes ownership ([MCNX-GPU-06] ordering;
// the bounded/local redistribution of specs/particle-container-and-gpu.md
// remains a later refinement of that call). The fill_packet_diag after it
// zeroes then refills the golden table from LIVE packets only, so an
// escaped packet's row flipping to all-zero is the built-in
// discrete-removal observable (the `escape-freestream` benchmark).
extern "C" void MCNuX_GeodesicPush(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_GeodesicPush;
  DECLARE_CCTK_PARAMETERS;

  const MetricGroups groups = metric_groups();
  const double dt = cctk_delta_time;

  // ONE step-scoped escape buffer, captured by raw pointer into every
  // tile's kernel (the mcnux_emission.cxx device-accumulator scope
  // precedent; per-tile buffers would lose counts across tiles).
  amrex::Gpu::DeviceVector<double> esc_dev(escape_num_slots, 0.0);
  double *const esc_ptr = esc_dev.data();

  for_each_packet_tile(groups, [&](PacketContainer &, int, const PacketIter &pti,
                                   const VertexMetricGather gather) {
    const long np = pti.numParticles();
    const auto ptd = pti.GetParticleTile().getParticleTileData();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      double x[3] = {ptd.rdata(PIdx::x)[ip], ptd.rdata(PIdx::y)[ip],
                     ptd.rdata(PIdx::z)[ip]};
      double p[3] = {ptd.rdata(PIdx::px)[ip], ptd.rdata(PIdx::py)[ip],
                     ptd.rdata(PIdx::pz)[ip]};
      // Pre-step position: the last-in-domain metric anchor of the escape
      // energy tally (mcnux_escape.hxx convention).
      const double x0[3] = {x[0], x[1], x[2]};
      geodesic_step_cellwise(gather, x, p, dt);
      if (outside_domain(x, gather.prob_lo, gather.prob_hi)) {
        const MetricSnapshot m = gather(x0[0], x0[1], x0[2]);
        const InverseSpatialMetric gu = spatial_metric_inverse(m.g);
        const double pt_up = p_t_closure(p[0], p[1], p[2], gu, m.alpha);
        const int sidx = ptd.idata(IntIdx::species)[ip];
        const double N = ptd.rdata(PIdx::w)[ip];
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_count)], 1.0);
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_number)], N);
        amrex::Gpu::Atomic::AddNoRet(
            &esc_ptr[escape_slot(sidx, escape_channel_energy)], N * pt_up);
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

  for (int patch = 0; patch < num_packet_patches(); ++patch)
    packet_population(patch).Redistribute();

  // Fold the step's escapes into the run-cumulative tally (per-step +=,
  // never reset; surfaced by MCNuX_EscapeDiag when gated on).
  {
    std::vector<double> esc_host(escape_num_slots, 0.0);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, esc_dev.begin(), esc_dev.end(),
                     esc_host.begin());
    escape_run_tally().add_slots(esc_host.data());
  }

  fill_packet_diag(groups, {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz, pk_pt},
                   packet_diag_size());
}

// ---------------------------------------------------------------------------
// Fluid-gather diagnostic  (MCNuX::mcnux_fluid_diag)
// ---------------------------------------------------------------------------

namespace {

// Number of doubles per packet in the device staging buffer:
// id, rho, T, Ye, u^0, u^1, u^2, u^3.
constexpr int fluid_diag_width = 8;

} // namespace

// Fill the per-packet fluid-gather diagnostic table
// (specs/neutrino-matter-interactions.md): for every packet, the direct
// one-cell read of the HydroBaseX (rho, T[MeV], Ye) fields at the packet's
// containing cell ([MCNX-INT-05], CellFluidGather — units verbatim, no
// interpolation) and the cell-centered fluid four-velocity u^mu from the
// Valencia lift of mcnux_fluid.hxx, with the metric gathered at the
// containing cell's CENTER coordinate (the spec's u^mu is cell-centered).
// Row keying, staging, and the single-rank assumption are those of
// fill_packet_diag above.
extern "C" void MCNuX_FluidGatherDiag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_FluidGatherDiag;
  DECLARE_CCTK_PARAMETERS;

  const MetricGroups mgroups = metric_groups();
  const HydroGroups hgroups = hydro_groups();
  const int nrows = diag_array_size("MCNuX::mcnux_fluid_diag");

  CCTK_REAL *const cols[7] = {fl_rho, fl_temp, fl_ye, fl_u0,
                              fl_u1,  fl_u2,   fl_u3};
  for (int c = 0; c < 7; ++c)
    for (int r = 0; r < nrows; ++r)
      cols[c][r] = 0.0;

  for_each_packet_tile_raw([&](const CarpetX::GHExt::PatchData &patchdata,
                               const CarpetX::GHExt::PatchData::LevelData
                                   &leveldata,
                               PacketContainer &, const PacketIter &pti) {
    const long np = pti.numParticles();
    if (np == 0)
      return;
    const auto ptd = pti.GetParticleTile().getParticleTileData();
    const VertexMetricGather mgather =
        make_gather(patchdata, leveldata, mgroups, pti);
    const CellFluidGather fgather =
        make_fluid_gather(patchdata, leveldata, hgroups, pti);

    amrex::Gpu::DeviceVector<double> staging(std::size_t(np) *
                                             fluid_diag_width);
    double *const out = staging.data();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      const double x[3] = {ptd.rdata(PIdx::x)[ip], ptd.rdata(PIdx::y)[ip],
                           ptd.rdata(PIdx::z)[ip]};
      int i0[3];
      fgather.containing_cell(x[0], x[1], x[2], i0);
      const FluidSample s = fgather.at_cell(i0);
      // Metric at the containing cell's center: u^mu is cell-centered.
      const MetricSnapshot m = mgather(fgather.cell_center(i0[0], 0),
                                       fgather.cell_center(i0[1], 1),
                                       fgather.cell_center(i0[2], 2));
      double u[4];
      valencia_four_velocity(m.alpha, m.beta, m.g, s.vel, u);
      double *const row = out + ip * fluid_diag_width;
      row[0] = double(amrex::particle_impl::unpack_id(ptd.m_idcpu[ip]));
      row[1] = s.state.rho_cgs;
      row[2] = s.state.T_MeV;
      row[3] = s.state.Ye;
      row[4] = u[0];
      row[5] = u[1];
      row[6] = u[2];
      row[7] = u[3];
    });
    amrex::Gpu::streamSynchronize();

    std::vector<double> host(std::size_t(np) * fluid_diag_width);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, staging.begin(), staging.end(),
                     host.begin());

    for (long ip = 0; ip < np; ++ip) {
      const double *const row = host.data() + ip * fluid_diag_width;
      const long id = long(row[0]);
      if (id < 1 || id > nrows)
        CCTK_VERROR("MCNuX packet id %ld is outside the fluid diagnostic "
                    "table (MCNuX::mcnux_fluid_diag has SIZE=%d rows)",
                    id, nrows);
      const int r = int(id - 1);
      for (int c = 0; c < 7; ++c)
        cols[c][r] = row[1 + c];
    }
  });
}

} // namespace MCNuX
