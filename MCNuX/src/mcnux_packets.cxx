// Construction of the MCNuX packet container on the CarpetX driver's AMR
// hierarchy, and the pinned tile-kernel exemplar.
//
// Governing spec: specs/particle-container-and-gpu.md [MCNX-GPU-01] (the
// container family, pinned in mcnux_packets.hxx) and [MCNX-GPU-02] (device
// kernels over particle tiles). The driver objects this file consumes are the
// ones specs/carpetx-thorn-integration.md [MCNX-CTX-01] makes available; the
// packet spec consumes them as given.

// TODO: Don't include files from other thorns; create a proper interface
//
// Same driver-access idiom, and the same recorded caveat, as
// src/mcnux_cadence.cxx: CarpetX exposes no public capability header for
// `ghext`, so its data model is reached through a relative include of the
// driver's own header. MCNuX lives one directory deeper than CarpetX's
// ODESolvers relative to the CarpetX arrangement -- hence the extra `../`.
#include "../../../CarpetX/CarpetX/src/driver.hxx"

#include "mcnux_packets.hxx"

#include <cctk.h>

#include <AMReX_GpuLaunch.H>
#include <AMReX_ParIter.H>

#include <memory>

namespace MCNuX {

std::unique_ptr<PacketContainer> make_packet_container() {
  if (!CarpetX::ghext)
    CCTK_VERROR("MCNuX tried to build its packet container, but the CarpetX "
                "driver singleton `ghext` is null. The packet population is "
                "defined on the driver's own AMR hierarchy "
                "(specs/particle-container-and-gpu.md [MCNX-GPU-01], "
                "specs/carpetx-thorn-integration.md [MCNX-CTX-01]) and cannot "
                "be built without it.");

  // Single-patch assumption, recorded: particle-container-and-gpu.md:303-307
  // organizes the packet population *per patch* (one container keyed to each
  // patch's AmrCore) and runs the verification baselines single-patch, because
  // the cross-patch transfer contract for packets crossing patch boundaries is
  // not yet specified. Until it is, a multipatch domain is an error rather
  // than a silently patch-0-only run.
  const std::size_t num_patches = CarpetX::ghext->patchdata.size();
  if (num_patches != 1)
    CCTK_VERROR("MCNuX found %d driver patches. The packet population is "
                "organized per patch, but the cross-patch transfer contract "
                "for packets is not yet specified "
                "(specs/particle-container-and-gpu.md, Open questions), so "
                "MCNuX currently requires a single-patch domain.",
                int(num_patches));

  amrex::AmrCore *const amrcore = CarpetX::ghext->patchdata.at(0).amrcore.get();
  if (!amrcore)
    CCTK_VERROR("MCNuX found driver patch 0 without an AmrCore; the packet "
                "container needs its Geometry, BoxArray, and "
                "DistributionMapping.");

  // Timing constraint, recorded: the level-0 BoxArray and DistributionMapping
  // only exist once the driver has run MakeNewGrids (just before
  // CCTK_BASEGRID, CarpetX/CarpetX/src/schedule.cxx:1239,1261), so a caller
  // must not invoke this before then. Nothing in MCNuX/schedule.ccl calls it
  // yet -- wiring the container into the schedule is a later task.
  return std::make_unique<PacketContainer>(amrcore);
}

long touch_packets(PacketContainer &pc) {
  long nvisited = 0;

  for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
    // Host-side tile walk; the per-packet work below is the device kernel.
    // This is emphatically *not* the CarpetX interpolator's model
    // (CarpetX/CarpetX/src/interpolate.cxx: `#pragma omp parallel for simd`
    // over particles, no amrex::ParallelFor anywhere) -- [MCNX-GPU-07].
    for (PacketContainer::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
      const long np = pti.numParticles();
      nvisited += np;

      // The by-value kernel view of [MCNX-GPU-02]: a trivially copyable
      // struct of raw SoA pointers (amrex/Src/Particle/AMReX_ParticleTile.H:32,
      // trivial-copyability asserted at :467 and re-asserted for this
      // container in mcnux_packets.hxx). It is the *only* thing the lambda
      // captures -- no `this`, no references, no host containers.
      auto ptd = pti.GetParticleTile().getParticleTileData();

      // One fused launch over the particle index, exactly the generic pure-SoA
      // structure of amrex/Tests/Particles/SOAParticle/main.cpp:99-112 and the
      // per-tile push loop of
      // warpx/Source/Particles/PhysicalParticleContainer.cpp. The later
      // gather / geodesic push / event sampling / deposition kernels replace
      // the body; the launch shape and capture rule stay as they are here.
      //
      // The body is deliberately an identity map on the packet state: it
      // reads every schema component of [MCNX-GPU-03] through the tile view
      // and writes each back unchanged, so the exemplar exercises the whole
      // schema without asserting any physics. Note what it does *not* do: no
      // amrex::Random() / ParallelForRNG (stateful device RNG is forbidden for
      // physics draws, specs/rng-and-statistical-acceptance.md; MCNuX draws
      // are evaluations of the pure u(S, q, e, k) of mcnux_rng.hxx), and no
      // per-packet host round-trip.
      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
        const amrex::ParticleReal xi[3] = {ptd.pos(RIdx::x, ip),
                                           ptd.pos(RIdx::y, ip),
                                           ptd.pos(RIdx::z, ip)};
        const amrex::ParticleReal pi[3] = {ptd.rdata(RIdx::p_x)[ip],
                                           ptd.rdata(RIdx::p_y)[ip],
                                           ptd.rdata(RIdx::p_z)[ip]};
        const amrex::ParticleReal weight = ptd.rdata(RIdx::w)[ip];
        const int species = ptd.idata(IIdx::species)[ip];
        const int event = ptd.idata(IIdx::event)[ip];

        ptd.pos(RIdx::x, ip) = xi[0];
        ptd.pos(RIdx::y, ip) = xi[1];
        ptd.pos(RIdx::z, ip) = xi[2];
        ptd.rdata(RIdx::p_x)[ip] = pi[0];
        ptd.rdata(RIdx::p_y)[ip] = pi[1];
        ptd.rdata(RIdx::p_z)[ip] = pi[2];
        ptd.rdata(RIdx::w)[ip] = weight;
        ptd.idata(IIdx::species)[ip] = species;
        ptd.idata(IIdx::event)[ip] = event;
        // The packet id q lives in the built-in 64-bit `idcpu` word, not in a
        // storage component ([MCNX-GPU-03]); it is immutable for the packet's
        // life ([MCNX-GPU-04]), so the exemplar only reads it, and only
        // through AMReX's sanctioned wrapper.
        const auto id = ptd.id(ip);
        amrex::ignore_unused(id);
      });
    }
  }

  return nvisited;
}

} // namespace MCNuX
