#ifndef MCNUX_GATHER_HXX
#define MCNUX_GATHER_HXX

// Shared driver-access plumbing for packet operators (the C1 factoring):
// the group-index lookups, the per-box gather constructors, the packet-tile
// walkers, and the per-packet diagnostic-table fill — hoisted verbatim from
// the anonymous namespace of mcnux_geodesic.cxx so the emission/creation
// loop (mcnux_emission.cxx) and later packet operators reuse ONE copy
// (never duplicated).
//
// Driver access follows the [MCNX-CTX-01] idiom of mcnux_cadence.cxx: the
// per-patch AmrCore and per-level group MultiFabs are reached through the
// CarpetX singleton `ghext`; grid data is read only through amrex::Array4
// views of groups declared in the calling routines' READS: clauses
// ([MCNX-CTX-04]).
//
// This header deliberately carries NO static_asserts: it is CarpetX-facing
// glue (the driver include below), not portable constexpr math. It is still
// #included from stub.cxx per the every-shared-header idiom so that it is
// compiled on every build.

// TODO: Don't include files from other thorns; create a proper interface
//
// The same relative include as mcnux_cadence.cxx (see the rationale there):
// CarpetX exposes no public capability header for `ghext`.
#include "../../../CarpetX/CarpetX/src/driver.hxx"

#include "mcnux_fluid.hxx"    // CellFluidGather
#include "mcnux_geodesic.hxx" // VertexMetricGather
#include "mcnux_particles.hxx"

#include <AMReX_ParIter.H>

#include <cctk.h>

namespace MCNuX {

inline void require_driver() {
  if (!CarpetX::ghext)
    CCTK_VERROR("MCNuX reached a packet operator, but the CarpetX driver "
                "singleton `ghext` is null. MCNuX accesses grid data and the "
                "particle hierarchy through the driver's native AMReX objects "
                "(specs/carpetx-thorn-integration.md [MCNX-CTX-01]) and cannot "
                "run without them.");
}

// ---------------------------------------------------------------------------
// Group-index lookups
// ---------------------------------------------------------------------------

struct MetricGroups {
  int metric, lapse, shift;
};

inline MetricGroups metric_groups() {
  MetricGroups g{CCTK_GroupIndex("ADMBaseX::metric"),
                 CCTK_GroupIndex("ADMBaseX::lapse"),
                 CCTK_GroupIndex("ADMBaseX::shift")};
  if (g.metric < 0 || g.lapse < 0 || g.shift < 0)
    CCTK_VERROR("MCNuX needs the ADMBaseX metric, lapse, and shift groups "
                "(specs/geodesic-propagation.md, Source of truth); is "
                "ADMBaseX active?");
  return g;
}

// The four HydroBaseX groups of the fluid-state gather
// (specs/neutrino-matter-interactions.md [MCNX-INT-05]; vel is the only
// multi-member group). Mirrors metric_groups().
struct HydroGroups {
  int rho, vel, temperature, ye;
};

inline HydroGroups hydro_groups() {
  HydroGroups g{CCTK_GroupIndex("HydroBaseX::rho"),
                CCTK_GroupIndex("HydroBaseX::vel"),
                CCTK_GroupIndex("HydroBaseX::temperature"),
                CCTK_GroupIndex("HydroBaseX::Ye")};
  if (g.rho < 0 || g.vel < 0 || g.temperature < 0 || g.ye < 0)
    CCTK_VERROR("MCNuX fluid-state gather needs the HydroBaseX rho, vel, "
                "temperature, and Ye groups "
                "(specs/neutrino-matter-interactions.md [MCNX-INT-05]); is "
                "HydroBaseX active?");
  return g;
}

using PacketIter = amrex::ParIterSoA<PIdx::nattribs, IntIdx::nattribs>;

// ---------------------------------------------------------------------------
// Per-box gather constructors
// ---------------------------------------------------------------------------

// The gather functor for one grid box: the ADMBaseX views of the box
// (CarpetX builds the group MultiFabs on the AmrCore's BoxArray and
// DistributionMapping — the ones the ParGDB tracks — so a particle
// iterator's box index and a grid MFIter's box index both address the
// matching fab) plus the level geometry.
inline VertexMetricGather
make_gather(const CarpetX::GHExt::PatchData &patchdata,
            const CarpetX::GHExt::PatchData::LevelData &leveldata,
            const MetricGroups &groups, const int box) {
  constexpr int tl = 0;
  VertexMetricGather gather{};
  gather.metric = leveldata.groupdata.at(groups.metric)->mfab.at(tl)->array(box);
  gather.lapse = leveldata.groupdata.at(groups.lapse)->mfab.at(tl)->array(box);
  gather.shift = leveldata.groupdata.at(groups.shift)->mfab.at(tl)->array(box);
  const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
  for (int d = 0; d < 3; ++d) {
    gather.prob_lo[d] = geom.ProbLo(d);
    gather.prob_hi[d] = geom.ProbHi(d);
    gather.dx[d] = geom.CellSize(d);
  }
  return gather;
}

inline VertexMetricGather
make_gather(const CarpetX::GHExt::PatchData &patchdata,
            const CarpetX::GHExt::PatchData::LevelData &leveldata,
            const MetricGroups &groups, const PacketIter &pti) {
  return make_gather(patchdata, leveldata, groups, pti.index());
}

// The fluid-gather functor for one grid box: the four HydroBaseX ccc views
// of the box plus the level geometry (same box-index reasoning as
// make_gather above).
inline CellFluidGather
make_fluid_gather(const CarpetX::GHExt::PatchData &patchdata,
                  const CarpetX::GHExt::PatchData::LevelData &leveldata,
                  const HydroGroups &groups, const int box) {
  constexpr int tl = 0;
  CellFluidGather gather{};
  gather.rho = leveldata.groupdata.at(groups.rho)->mfab.at(tl)->array(box);
  gather.vel = leveldata.groupdata.at(groups.vel)->mfab.at(tl)->array(box);
  gather.temperature =
      leveldata.groupdata.at(groups.temperature)->mfab.at(tl)->array(box);
  gather.ye = leveldata.groupdata.at(groups.ye)->mfab.at(tl)->array(box);
  const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
  for (int d = 0; d < 3; ++d) {
    gather.prob_lo[d] = geom.ProbLo(d);
    gather.dx[d] = geom.CellSize(d);
  }
  return gather;
}

inline CellFluidGather
make_fluid_gather(const CarpetX::GHExt::PatchData &patchdata,
                  const CarpetX::GHExt::PatchData::LevelData &leveldata,
                  const HydroGroups &groups, const PacketIter &pti) {
  return make_fluid_gather(patchdata, leveldata, groups, pti.index());
}

// ---------------------------------------------------------------------------
// Packet-tile walkers
// ---------------------------------------------------------------------------

// Walk every particle tile of the population, raw form:
// f(patchdata, leveldata, pc, pti) — callers build the gathers they need.
template <class F> void for_each_packet_tile_raw(F &&f) {
  for (int patch = 0; patch < num_packet_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    PacketContainer &pc = packet_population(patch);
    for (const auto &leveldata : patchdata.leveldata) {
      const int level = leveldata.level;
      if (level > pc.finestLevel())
        break;
      for (PacketIter pti(pc, level); pti.isValid(); ++pti)
        f(patchdata, leveldata, pc, pti);
    }
  }
}

// Walk every particle tile of the population: f(pc, level, pti, gather).
template <class F>
void for_each_packet_tile(const MetricGroups &groups, F &&f) {
  for_each_packet_tile_raw(
      [&](const CarpetX::GHExt::PatchData &patchdata,
          const CarpetX::GHExt::PatchData::LevelData &leveldata,
          PacketContainer &pc, const PacketIter &pti) {
        f(pc, leveldata.level, pti,
          make_gather(patchdata, leveldata, groups, pti));
      });
}

// ---------------------------------------------------------------------------
// Per-packet diagnostic table  (MCNuX::mcnux_packet_diag)
// ---------------------------------------------------------------------------

// Size of a 1d DISTRIB=CONSTANT diagnostic array as declared in
// interface.ccl (mcnux_packet_diag, mcnux_fluid_diag).
inline int diag_array_size(const char *const group) {
  const int gi = CCTK_GroupIndex(group);
  if (gi < 0)
    CCTK_VERROR("%s is not a known group", group);
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("%s must be a 1-dimensional grid array", group);
  return int(*sizes[0]);
}

inline int packet_diag_size() {
  return diag_array_size("MCNuX::mcnux_packet_diag");
}

struct PacketDiagColumns {
  CCTK_REAL *x, *y, *z, *px, *py, *pz, *pt;
};

// Fill the diagnostic table from the current population: row = packet
// id - 1, columns the position, the lower momentum, and the covariant
// energy p_t. Defined in mcnux_geodesic.cxx (the population owner); shared
// so the emission diagnostic (mcnux_emission.cxx) fills the same table.
void fill_packet_diag(const MetricGroups &groups, const PacketDiagColumns &cols,
                      int nrows);

} // namespace MCNuX

#endif // MCNUX_GATHER_HXX
