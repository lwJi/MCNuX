#ifndef MCNUX_PARTICLES_HXX
#define MCNUX_PARTICLES_HXX

// The MCNuX packet container: the pure-SoA storage of the packet population.
//
// Governing spec: specs/particle-container-and-gpu.md, requirements
// [MCNX-GPU-01] (pure-SoA container pin), [MCNX-GPU-02] (device-kernel
// execution over particle tiles), [MCNX-GPU-03] (component schema), and the
// storage half of [MCNX-GPU-04] (the packet id lives in the built-in 64-bit
// `idcpu` word). The AMReX legs of the precision gate
// specs/build-and-integration.md [MCNX-BLD-03] are asserted at the bottom of
// this header (this is the thorn's first AMReX-consuming header, so the gate
// lives here; stub.cxx keeps the plain double/CCTK_REAL legs).
//
// Storage mapping of the PKT state table (particle-container-and-gpu.md,
// "Packet component schema"): position x^i -> SoA reals 0..2 (the container's
// position components), momentum p_i -> SoA reals 3..5, weight N -> SoA real 6
// (an ordinary SoA component, exactly like WarpX's `w`), species s -> SoA int
// 0, event counter e -> SoA int 1 (consumed as uint32 by the RNG packing),
// packet id q -> the built-in `idcpu` (no component slot).
//
// In `amrex::ParticleContainerPureSoA<NArrayReal, NArrayInt>` the position is
// NOT a separate storage slot: it occupies real components
// 0..AMREX_SPACEDIM-1 of the NArrayReal compile-time reals
// (amrex/Src/Particle/AMReX_ParticleContainer.H, the
// "Pure SoA particle containers require at least AMREX_SPACEDIM real
// components for positions" static_assert). Hence NArrayReal = 3 position +
// 3 momentum + 1 weight = 7 and NArrayInt = 2.
//
// The index enums below follow the WarpX production container, the corpus'
// sole particle-operator execution reference
// (warpx/Source/Particles/WarpXParticleContainer.H, struct PIdx).
//
// [MCNX-GPU-07] anti-pattern warning: the CarpetX driver interpolator's
// legacy AoS `AmrParticleContainer` and host per-particle loops
// (CarpetX/CarpetX/src/interp.hxx, interpolate.cxx) are NOT the pattern here;
// only the driver-object *access* idiom (mcnux_cadence.cxx) is taken from it.
//
// [MCNX-RNG-05]: nothing in this header (or in any packet kernel) may use
// AMReX's stateful per-thread device RNG; physics draws are evaluations of
// the pure function u(S, q, e, k) of mcnux_rng.hxx.

#include <AMReX_ParticleContainer.H>

#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Compile-time SoA component indices  [MCNX-GPU-03]
// ---------------------------------------------------------------------------

// Real components. The first AMREX_SPACEDIM slots are the container's
// position components; MCNuX's own components start at index AMREX_SPACEDIM.
struct PIdx {
  enum { x = 0, y, z, px, py, pz, w, nattribs };
  static constexpr auto names = {"x", "y", "z", "px", "py", "pz", "w"};
  static_assert(names.size() == nattribs,
                "PIdx names must match the component count");
};

// Integer components: species s in {0, 1, 2} and the RNG event counter e.
struct IntIdx {
  enum { species = 0, event_counter, nattribs };
  static constexpr auto names = {"species", "event_counter"};
  static_assert(names.size() == nattribs,
                "IntIdx names must match the component count");
};

static_assert(PIdx::px == AMREX_SPACEDIM,
              "MCNuX components must start after the position slots");
static_assert(PIdx::nattribs == AMREX_SPACEDIM + 3 + 1 &&
                  IntIdx::nattribs == 2,
              "schema: 3 position + 3 momentum + 1 weight reals, 2 ints");

// ---------------------------------------------------------------------------
// The packet container  [MCNX-GPU-01]
// ---------------------------------------------------------------------------
// Pure SoA (zero AoS bytes per particle; per-component contiguous arrays);
// exemplar: WarpXParticleContainer. Legacy AoS containers are forbidden for
// the packet population. Constructors are inherited from the base: the live
// wiring to the driver's per-patch ParGDB
// (ghext->patchdata.at(p).amrcore->GetParGDB()) belongs to the first task
// that schedules a population owner, not to this header.
class PacketContainer
    : public amrex::ParticleContainerPureSoA<PIdx::nattribs, IntIdx::nattribs> {
public:
  using Base = amrex::ParticleContainerPureSoA<PIdx::nattribs, IntIdx::nattribs>;
  using Base::Base;
};

// The by-value tile view every packet kernel captures ([MCNX-GPU-02]).
using PacketTileData = PacketContainer::ParticleTileType::ParticleTileDataType;

// ---------------------------------------------------------------------------
// Exemplar tile kernel  [MCNX-GPU-02]
// ---------------------------------------------------------------------------
// The execution idiom every per-packet operator follows: one
// amrex::ParallelFor over the particle index of a tile, whose lambda captures
// only trivially-copyable values by value (here the ParticleTileData view).
// Pattern: amrex/Tests/Particles/SOAParticle/main.cpp
// (ptd = pti.GetParticleTile().getParticleTileData(); ParallelFor(np, ...)).
// This exemplar touches every schema component and bumps the event counter;
// it is compile-only in this increment (no scheduled caller launches it yet).
inline void exemplar_packet_kernel(PacketTileData ptd, long np) {
  amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
    const amrex::ParticleReal p2 =
        ptd.rdata(PIdx::px)[ip] * ptd.rdata(PIdx::px)[ip] +
        ptd.rdata(PIdx::py)[ip] * ptd.rdata(PIdx::py)[ip] +
        ptd.rdata(PIdx::pz)[ip] * ptd.rdata(PIdx::pz)[ip];
    amrex::ignore_unused(p2, ptd.rdata(PIdx::w)[ip],
                         ptd.idata(IntIdx::species)[ip], ptd.m_idcpu[ip]);
    ptd.idata(IntIdx::event_counter)[ip] += 1;
  });
}

// ---------------------------------------------------------------------------
// Compile-time layout pin  (particle-container-and-gpu.md, Verification)
// ---------------------------------------------------------------------------
static_assert(
    std::is_base_of_v<
        amrex::ParticleContainerPureSoA<PIdx::nattribs, IntIdx::nattribs>,
        PacketContainer>,
    "[MCNX-GPU-01] the packet container must derive from an "
    "amrex::ParticleContainerPureSoA instantiation");
static_assert(std::is_trivially_copyable_v<PacketTileData>,
              "[MCNX-GPU-02] the tile view must be trivially copyable "
              "(mirrors AMReX's own ParticleTileData guard)");
static_assert(PacketContainer::ParticleType::is_soa_particle,
              "[MCNX-GPU-01] zero AoS bytes per particle");
static_assert(PacketContainer::NArrayReal == PIdx::nattribs &&
                  PacketContainer::NArrayInt == IntIdx::nattribs,
              "[MCNX-GPU-03] one storage component per physical component");

// ---------------------------------------------------------------------------
// Precision gate, AMReX legs  [MCNX-BLD-03]
// ---------------------------------------------------------------------------
// All reals are IEEE-754 binary64; the build must pin both amrex::Real and
// amrex::ParticleReal to double (build-and-integration.md). A failure here is
// the gate catching a misconfigured AMReX build, not a bad assert.
static_assert(sizeof(amrex::Real) == 8,
              "[MCNX-BLD-03] AMReX must be built with double-precision Real");
static_assert(sizeof(amrex::ParticleReal) == 8,
              "[MCNX-BLD-03] AMReX must be built with double-precision "
              "ParticleReal");
static_assert(std::is_same_v<amrex::Real, double> &&
                  std::is_same_v<amrex::ParticleReal, double>,
              "[MCNX-BLD-03] amrex::Real and amrex::ParticleReal must be "
              "double (mirrors CarpetX's own driver.hxx pin)");

} // namespace MCNuX

#endif // MCNUX_PARTICLES_HXX
