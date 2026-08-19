#ifndef MCNUX_PACKETS_HXX
#define MCNUX_PACKETS_HXX

// The MCNuX packet container: the sole authoritative storage of the packet
// population.
//
// Governing spec: specs/particle-container-and-gpu.md, requirements
// [MCNX-GPU-01] (pure-SoA container family, binding), [MCNX-GPU-02]
// (device-kernel execution over particle tiles), [MCNX-GPU-03] (component
// schema: one storage component per physical component of the PKT state
// table, all reals binary64). The compile-time layout pin the spec's
// Verification section demands (particle-container-and-gpu.md:236-240) is the
// block of static_asserts at the bottom of this header, so every translation
// unit that includes it re-verifies the layout at build time.
//
// This header also carries the AMReX legs of the precision gate of
// specs/build-and-integration.md [MCNX-BLD-03]: it is the first (and, by
// convention, the permanent) home of the sizeof(amrex::Real) and
// sizeof(amrex::ParticleReal) assertions, because this is the translation-unit
// family that consumes AMReX headers. The language-level legs
// (sizeof(double), sizeof(CCTK_REAL)) stay in src/stub.cxx.
//
// DELIBERATE EXCEPTION to the AMReX-free convention of mcnux_units.hxx and
// mcnux_rng.hxx: those headers are plain portable C++ so that host code,
// device code, and standalone unit tests can all include them. This header
// cannot be — the container type *is* an AMReX type — so it is the thorn's
// first AMReX consumer and includes the AMReX particle headers itself.
// CarpetX's driver.hxx does not pull them in. No configuration.ccl change
// follows from this: `REQUIRES CarpetX` already puts AMReX on the include
// path, and AMREX_PARTICLES=ON is baked into the sandbox and CI toolchains.
// The driver singleton `ghext` is deliberately *not* reachable from here;
// everything that needs it lives in mcnux_packets.cxx.

#include "mcnux_units.hxx"

#include <AMReX_AmrCore.H>
// AMReX_AmrCore.H only forward-declares AmrParGDB, so the derived-to-base
// conversion used below needs its definition.
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParticleContainer.H>

#include <initializer_list>
#include <memory>
#include <type_traits>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Component schema  [MCNX-GPU-03]
// ---------------------------------------------------------------------------
// particle-container-and-gpu.md:127-140 maps each physical component of the
// PKT state table onto exactly one storage slot:
//
//   position x^i  -> SoA reals 0..2   (the container's built-in positions)
//   momentum p_i  -> SoA reals 3..5
//   weight   N    -> SoA real  6      (an ordinary SoA component, WarpX's `w`)
//   species  s    -> SoA int   0      (values in {0, 1, 2})
//   event    e    -> SoA int   1      (consumed as uint32 by the RNG packing)
//   packet id q   -> the built-in 64-bit `idcpu` -- NOT a storage component
//
// Component *names* are free (spec:138); the mapping is binding. The slot
// names below follow the WarpX `PIdx`/`IntIdx` idiom
// (warpx/Source/Particles/WarpXParticleContainer.H:64-99): an enum whose last
// enumerator is the attribute count, paired with a matching name list that is
// handed to AMReX via SetSoACompileTimeNames so that plotfiles and checkpoints
// carry the names too.
//
// AMReX requires the first AMREX_SPACEDIM real components of a pure-SoA
// container to *be* the positions (amrex/Src/Particle/AMReX_ParticleContainer.H
// :212-213), which is why x, y, z lead the enumeration.

struct RIdx {
  enum {
    x,   //!< position x^1
    y,   //!< position x^2
    z,   //!< position x^3
    p_x, //!< momentum p_1
    p_y, //!< momentum p_2
    p_z, //!< momentum p_3
    w,   //!< weight N (number of physical neutrinos the packet represents)
    nattribs
  };

  static constexpr auto names = {"x", "y", "z", "p_x", "p_y", "p_z", "w"};

  static_assert(names.size() == nattribs);
};

struct IIdx {
  enum {
    species, //!< species index s, in {0, 1, 2}: see MCNuX::Species
    event,   //!< event counter e, consumed as uint32 by u(S, q, e, k)
    nattribs
  };

  static constexpr auto names = {"species", "event"};

  static_assert(names.size() == nattribs);
};

// ---------------------------------------------------------------------------
// The container  [MCNX-GPU-01]
// ---------------------------------------------------------------------------
// Pure SoA, zero AoS bytes per particle, mirroring the production exemplar
// `class WarpXParticleContainer : public amrex::ParticleContainerPureSoA<...>`
// (warpx/Source/Particles/WarpXParticleContainer.H:193). Legacy AoS flavors
// (amrex::AmrParticleContainer, struct-particle ParticleContainer) are
// forbidden for the packet population -- the CarpetX driver interpolator's
// `amrex::AmrParticleContainer<3, 2>` (CarpetX/CarpetX/src/interp.hxx:38) is
// the anti-pattern of [MCNX-GPU-07], not the pattern to copy.
//
// The default allocator is AMReX's arena allocator, i.e. device memory in a
// GPU build; the container is move-only, inheriting the base's deleted copy
// operations (AMReX_ParticleContainer.H:301-305).
class PacketContainer
    : public amrex::ParticleContainerPureSoA<RIdx::nattribs, IIdx::nattribs> {
public:
  using Base =
      amrex::ParticleContainerPureSoA<RIdx::nattribs, IIdx::nattribs>;

  //! Construct the packet population on a driver patch's AMR hierarchy. The
  //! ParGDB carries the Geometry / BoxArray / DistributionMapping of that
  //! hierarchy and keeps the container tracking regrids automatically
  //! (AMReX_ParticleContainer.H:231-238); this is the WarpX construction path
  //! (warpx/Source/Particles/WarpXParticleContainer.cpp:94-95).
  explicit PacketContainer(amrex::AmrCore *amr_core)
      : Base(par_gdb_of(amr_core)) {
    SetSoACompileTimeNames({RIdx::names.begin(), RIdx::names.end()},
                           {IIdx::names.begin(), IIdx::names.end()});
  }

private:
  //! Null-checked accessor, so that the check happens *before* the base
  //! subobject is constructed from the pointer.
  static amrex::ParGDBBase *par_gdb_of(amrex::AmrCore *amr_core) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amr_core != nullptr,
        "MCNuX::PacketContainer requires a non-null amrex::AmrCore");
    return amr_core->GetParGDB();
  }
};

// ---------------------------------------------------------------------------
// Factory and exemplar kernel (defined in mcnux_packets.cxx)
// ---------------------------------------------------------------------------

//! Build a packet container on the CarpetX driver's AmrCore. Reaching the
//! driver singleton is what keeps this out of the header; see
//! mcnux_packets.cxx for the single-patch assumption and the schedule-timing
//! constraint (the level-0 grids only exist after MakeNewGrids).
std::unique_ptr<PacketContainer> make_packet_container();

//! The pinned [MCNX-GPU-02] tile-kernel exemplar: one amrex::ParallelFor per
//! (level, tile) over the particle index, capturing only the trivially
//! copyable ParticleTileData by value. Its body is deliberately an identity
//! map on the packet state -- it exists to establish the execution idiom that
//! the gather / push / sampling / deposition kernels of the later specs
//! follow, not to implement physics. Returns the number of packets visited.
long touch_packets(PacketContainer &pc);

// ---------------------------------------------------------------------------
// Compile-time layout pin  (particle-container-and-gpu.md:236-240)
// ---------------------------------------------------------------------------
// "The packet container type is (or derives from) an
// amrex::ParticleContainerPureSoA instantiation; a compile-time check asserts
// this and asserts the tile view is trivially copyable (mirroring AMReX's own
// static_assert(std::is_trivially_copyable<PTD>()))."
//
// ParticleContainerPureSoA is an alias template for ParticleContainer_impl
// (AMReX_ParticleContainer.H:1624-1625), so it names a real class type and
// std::is_base_of_v is well defined on it.

static_assert(
    std::is_base_of_v<
        amrex::ParticleContainerPureSoA<RIdx::nattribs, IIdx::nattribs>,
        PacketContainer>,
    "[MCNX-GPU-01] the packet population must live in a container of the "
    "amrex::ParticleContainerPureSoA family; legacy AoS containers are "
    "forbidden");

static_assert(std::is_trivially_copyable_v<PacketContainer::PTDType>,
              "[MCNX-GPU-02] the by-value kernel view "
              "(ParticleTileData) must be trivially copyable, so that device "
              "lambdas can capture it by value");

static_assert(std::is_trivially_copyable_v<PacketContainer::ConstPTDType>,
              "[MCNX-GPU-02] the const by-value kernel view must be trivially "
              "copyable too");

// The schema mapping itself, [MCNX-GPU-03]: one storage component per physical
// component of the PKT state table, and no spares.
static_assert(AMREX_SPACEDIM == 3,
              "MCNuX packets carry three spatial coordinates");
static_assert(RIdx::nattribs == AMREX_SPACEDIM + 3 + 1,
              "[MCNX-GPU-03] SoA reals are exactly 3 positions + 3 momenta + "
              "1 weight");
static_assert(IIdx::nattribs == 2,
              "[MCNX-GPU-03] SoA ints are exactly the species index and the "
              "event counter");
static_assert(PacketContainer::NArrayReal == RIdx::nattribs &&
                  PacketContainer::NArrayInt == IIdx::nattribs,
              "[MCNX-GPU-03] the container's component counts must be the "
              "schema's");
static_assert(PacketContainer::ParticleType::is_soa_particle,
              "[MCNX-GPU-01] zero AoS bytes per packet");
static_assert(RIdx::x == 0 && RIdx::y == 1 && RIdx::z == 2,
              "[MCNX-GPU-03] AMReX pins the first AMREX_SPACEDIM SoA reals to "
              "be the positions");
static_assert(IIdx::species == 0 && IIdx::event == 1,
              "[MCNX-GPU-03] species and event-counter slot indices");
// The species component holds MCNuX::Species values; the enumeration itself is
// owned by mcnux_units.hxx [MCNX-CNV-06] and is never re-authored here.
static_assert(NUM_SPECIES == 3 && species_index(Species::NuE) == 0 &&
                  species_index(Species::NuEBar) == 1 &&
                  species_index(Species::NuX) == 2,
              "[MCNX-GPU-03] the species component takes values in {0, 1, 2}");

// ---------------------------------------------------------------------------
// Precision gate, AMReX legs  (build-and-integration.md [MCNX-BLD-03])
// ---------------------------------------------------------------------------
// "All reals are IEEE-754 binary64" ([MCNX-GPU-03] restates it for the packet
// components). amrex::Real and amrex::ParticleReal are gated by *independent*
// build macros -- AMREX_SINGLE_PRECISION_PARTICLES demotes ParticleReal alone
// (amrex/Src/Base/AMReX_REAL.H:61-65,90) -- so neither leg implies the other.
static_assert(sizeof(amrex::Real) == 8,
              "[MCNX-BLD-03] MCNuX requires amrex::Real == binary64; build "
              "AMReX in double precision (AMREX_PRECISION=DOUBLE)");
static_assert(std::is_same_v<amrex::Real, double>,
              "[MCNX-BLD-03] amrex::Real must be `double`");
static_assert(sizeof(amrex::ParticleReal) == 8,
              "[MCNX-BLD-03] MCNuX requires amrex::ParticleReal == binary64; "
              "AMREX_SINGLE_PRECISION_PARTICLES must not be set");
static_assert(std::is_same_v<amrex::ParticleReal, double>,
              "[MCNX-BLD-03] amrex::ParticleReal must be `double`");

} // namespace MCNuX

#endif // MCNUX_PACKETS_HXX
