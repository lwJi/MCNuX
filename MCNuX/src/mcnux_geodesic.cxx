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
// the exemplar_packet_kernel pattern of mcnux_particles.hxx).
//
// Scheduling ([MCNX-CTX-03]): both routines run in global/level mode —
// CarpetX invokes such a routine exactly once per traversal — and walk the
// patches, levels, and particle tiles themselves.

// TODO: Don't include files from other thorns; create a proper interface
//
// The same relative include as mcnux_cadence.cxx (see the rationale there):
// CarpetX exposes no public capability header for `ghext`.
#include "../../../CarpetX/CarpetX/src/driver.hxx"

#include "mcnux_geodesic.hxx"
#include "mcnux_particles.hxx"
#include "mcnux_tetrad.hxx"

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

void require_driver() {
  if (!CarpetX::ghext)
    CCTK_VERROR("MCNuX reached a packet operator, but the CarpetX driver "
                "singleton `ghext` is null. MCNuX accesses grid data and the "
                "particle hierarchy through the driver's native AMReX objects "
                "(specs/carpetx-thorn-integration.md [MCNX-CTX-01]) and cannot "
                "run without them.");
}

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

namespace {

// ---------------------------------------------------------------------------
// Metric access
// ---------------------------------------------------------------------------

struct MetricGroups {
  int metric, lapse, shift;
};

MetricGroups metric_groups() {
  MetricGroups g{CCTK_GroupIndex("ADMBaseX::metric"),
                 CCTK_GroupIndex("ADMBaseX::lapse"),
                 CCTK_GroupIndex("ADMBaseX::shift")};
  if (g.metric < 0 || g.lapse < 0 || g.shift < 0)
    CCTK_VERROR("MCNuX geodesic push needs the ADMBaseX metric, lapse, and "
                "shift groups (specs/geodesic-propagation.md, Source of "
                "truth); is ADMBaseX active?");
  return g;
}

using PacketIter = amrex::ParIterSoA<PIdx::nattribs, IntIdx::nattribs>;

// The gather functor for one particle tile: the ADMBaseX views of the tile's
// box (CarpetX builds the group MultiFabs on the AmrCore's BoxArray and
// DistributionMapping — the ones the ParGDB tracks — so the particle
// iterator's box index addresses the matching fab) plus the level geometry.
VertexMetricGather
make_gather(const CarpetX::GHExt::PatchData &patchdata,
            const CarpetX::GHExt::PatchData::LevelData &leveldata,
            const MetricGroups &groups, const PacketIter &pti) {
  constexpr int tl = 0;
  const int box = pti.index();
  VertexMetricGather gather{};
  gather.metric = leveldata.groupdata.at(groups.metric)->mfab.at(tl)->array(box);
  gather.lapse = leveldata.groupdata.at(groups.lapse)->mfab.at(tl)->array(box);
  gather.shift = leveldata.groupdata.at(groups.shift)->mfab.at(tl)->array(box);
  const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
  for (int d = 0; d < 3; ++d) {
    gather.prob_lo[d] = geom.ProbLo(d);
    gather.dx[d] = geom.CellSize(d);
  }
  return gather;
}

// Walk every particle tile of the population: f(pc, level, pti, gather).
template <class F> void for_each_packet_tile(const MetricGroups &groups, F &&f) {
  for (int patch = 0; patch < num_packet_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    PacketContainer &pc = packet_population(patch);
    for (const auto &leveldata : patchdata.leveldata) {
      const int level = leveldata.level;
      if (level > pc.finestLevel())
        break;
      for (PacketIter pti(pc, level); pti.isValid(); ++pti)
        f(pc, level, pti, make_gather(patchdata, leveldata, groups, pti));
    }
  }
}

// ---------------------------------------------------------------------------
// Diagnostic table  (MCNuX::mcnux_packet_diag)
// ---------------------------------------------------------------------------

// Size of the MCNuX::mcnux_packet_diag array as declared in interface.ccl.
int packet_diag_size() {
  const int gi = CCTK_GroupIndex("MCNuX::mcnux_packet_diag");
  if (gi < 0)
    CCTK_VERROR("MCNuX::mcnux_packet_diag is not a known group");
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("MCNuX::mcnux_packet_diag must be a 1-dimensional grid array");
  return int(*sizes[0]);
}

struct PacketDiagColumns {
  CCTK_REAL *x, *y, *z, *px, *py, *pz, *pt;
};

// Number of doubles per packet in the device staging buffer:
// id, x, y, z, p_x, p_y, p_z, p_t.
constexpr int diag_width = 8;

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

// The deterministic seed set of the `minkowski-freestream` benchmark (and of
// T16's Schwarzschild legs): eight packets with exact binary-fraction
// positions inside |x^i| <= 0.25 and momenta spanning the six axis
// directions (at four different magnitudes — the coordinate speed must be 1
// regardless) plus two oblique directions. Over the benchmark's 4 x 0.125
// coordinate time every packet stays well inside the [-1, 1]^3 domain, so
// the escape policy is never exercised. Row index in the diagnostic table =
// packet id - 1 = fixture index. Function-local constexpr (the
// mcnux_srcterms.cxx CUDA note on namespace-scope struct constants).
struct PacketFixture {
  double x, y, z;
  double px, py, pz;
};

} // namespace

extern "C" void MCNuX_SeedSyntheticPackets(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SeedSyntheticPackets;
  DECLARE_CCTK_PARAMETERS;

  constexpr int nfixture = 8;
  constexpr PacketFixture fixture[nfixture] = {
      {0.125, 0.0, 0.0, 1.0, 0.0, 0.0},
      {-0.125, 0.0, 0.0, -1.0, 0.0, 0.0},
      {0.0, 0.125, 0.0, 0.0, 2.0, 0.0},
      {0.0, -0.125, 0.0, 0.0, -0.5, 0.0},
      {0.0, 0.0, 0.125, 0.0, 0.0, 0.25},
      {0.0, 0.0, -0.125, 0.0, 0.0, -4.0},
      {0.25, -0.25, 0.125, 0.5, -0.25, 0.125},
      {-0.0625, 0.1875, -0.25, -0.375, 0.75, 0.5},
  };

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
// The geodesic push
// ---------------------------------------------------------------------------

// One transport step of free streaming for every packet: per tile, one
// device kernel applying geodesic_step (RK4 on the [MCNX-GEO-01] RHS with
// the trilinear gather of [MCNX-GEO-03]) over the coarsest-level Delta t
// (cctk_delta_time in level/global mode, the cadence contract of
// mcnux_cadence.cxx). Afterwards a plain Redistribute() re-establishes
// ownership (packets that left the domain are dropped by AMReX; escape
// tallies are a later task's, as is the bounded/local redistribution of
// specs/particle-container-and-gpu.md — T23 refines this call).
extern "C" void MCNuX_GeodesicPush(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_GeodesicPush;
  DECLARE_CCTK_PARAMETERS;

  const MetricGroups groups = metric_groups();
  const double dt = cctk_delta_time;

  for_each_packet_tile(groups, [&](PacketContainer &, int, const PacketIter &pti,
                                   const VertexMetricGather gather) {
    const long np = pti.numParticles();
    const auto ptd = pti.GetParticleTile().getParticleTileData();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
      double x[3] = {ptd.rdata(PIdx::x)[ip], ptd.rdata(PIdx::y)[ip],
                     ptd.rdata(PIdx::z)[ip]};
      double p[3] = {ptd.rdata(PIdx::px)[ip], ptd.rdata(PIdx::py)[ip],
                     ptd.rdata(PIdx::pz)[ip]};
      geodesic_step(gather, x, p, dt);
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

  fill_packet_diag(groups, {pk_x, pk_y, pk_z, pk_px, pk_py, pk_pz, pk_pt},
                   packet_diag_size());
}

} // namespace MCNuX
