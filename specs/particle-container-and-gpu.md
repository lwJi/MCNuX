# Particle container and GPU execution model

> Technical leaf spec. Self-contained: an agent can implement the packet container
> layout, the device-kernel execution model, the packet-id contract, and the
> redistribution invariants from this file alone. It restates only the conventions it
> uses and references [conventions-and-units](./conventions-and-units.md) and
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) for the rest.
> `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope

MCNuX's stated purpose is GPU Monte Carlo transport, so — as a sanctioned exception to
the corpus rule against prescribing internals — this spec pins the two layout/execution
choices that correctness of that purpose demands, plus the invariants that make a
distributed AMR particle population well-defined:

- The **pure-SoA container pin**: the packet population lives in the
  `amrex::ParticleContainerPureSoA` container family, with the WarpX production
  container as exemplar and the rationale recorded.
- The **device-kernel execution requirement**: every per-packet operator runs as device
  kernels over particle tiles, with only trivially-copyable views captured by value and
  no per-step host round-trips — with the pinned kernel exemplar list (gather, push,
  deposition, tile view, generic structure).
- The **packet component schema** mapping (one storage component per physical component
  of [packet-representation-and-sampling](./packet-representation-and-sampling.md)'s
  state table) and the **packet-id uniqueness contract** that
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) and
  [packet-representation-and-sampling](./packet-representation-and-sampling.md) defer
  to this spec.
- The **after-each-step ownership invariant** and the **bounded-motion precondition**
  of local redistribution, including `Redistribute` timing under AMR.
- The **anti-pattern warning** against the nearest in-stack particle code.

Out of scope:

- The physical meaning, units, and valid ranges of the packet components, and how
  packets are created ([packet-representation-and-sampling](./packet-representation-and-sampling.md)).
- The equations the kernels implement: geodesic push
  ([geodesic-propagation](./geodesic-propagation.md)), event sampling and outcomes
  ([neutrino-matter-interactions](./neutrino-matter-interactions.md),
  [trapped-regime-treatment](./trapped-regime-treatment.md)), coefficient evaluation
  ([opacity-eos-evaluation](./opacity-eos-evaluation.md)).
- How the driver's grid objects (AmrCore, MultiFabs, `Array4` views, centering) are
  reached, the schedule placement of transport, and the READS/WRITES obligations —
  [carpetx-thorn-integration](./carpetx-thorn-integration.md); this spec consumes those
  objects as given.
- The source-term grid variables and the deposition *protocol* (zero-then-add,
  conservation ledger) — [hydro-coupling-source-terms](./hydro-coupling-source-terms.md);
  this spec pins only the execution idiom deposition kernels use.
- The RNG algorithm and the prohibition rationale for stateful device RNG
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)); this spec
  inherits the prohibition.
- Concrete benchmarks and golden data
  ([verification-suite-design](./verification-suite-design.md)).

## Source of truth

WarpX is the **sole particle-operator execution reference** in this corpus: it is the
reference stack's only production GPU particle code. AMReX supplies the container
family, the tile-view contract, and the redistribution guarantee. All cited paths
resolve under the validator's reference roots.

- `warpx/Source/Particles/WarpXParticleContainer.H` — the production pure-SoA container:
  `class WarpXParticleContainer : public amrex::ParticleContainerPureSoA<PIdx::nattribs, 0, …>`
  (~line 193), with the compile-time SoA real components enumerated in `struct PIdx`
  (position, weight `w` as an *ordinary SoA component*, momenta) and runtime components
  added by name.
- `warpx/Source/Particles/PhysicalParticleContainer.cpp` — the push loop (~lines
  1387–1568): per tile, POD position functors and raw
  `ParticleReal* AMREX_RESTRICT` SoA pointers (`attribs[…].dataPtr() + offset`) are
  built once, then one `amrex::ParallelFor` over `np_to_push` with a
  `[=] AMREX_GPU_DEVICE` lambda does gather → momentum push → position update; runtime
  components added by name via `AddRealComp("opticalDepthQSR")` etc. (~line 263).
- `warpx/Source/Particles/Gather/FieldGather.H` — `doGatherShapeN` (~line 348):
  `AMREX_GPU_HOST_DEVICE` device function reading `Array4<amrex::Real const>` field
  views (plus their `IndexType` staggerings) at a particle position — the gather
  exemplar.
- `warpx/Source/Particles/Pusher/UpdatePosition.H` — the `UpdatePosition` device
  function (~line 23), invoked from the same kernel body — the push exemplar.
- `warpx/Source/Particles/Pusher/GetAndSetPosition.H` — the trivially-copyable
  `GetParticlePosition`/`SetParticlePosition` functors holding raw
  `AMREX_RESTRICT` SoA pointers.
- `warpx/Source/Particles/Deposition/CurrentDeposition.H` — shape-stencil deposition:
  `amrex::Gpu::Atomic::AddNoRet(&jx_arr(i,j,k), …)` inside a `ParallelFor` over
  particles — the deposition exemplar.
- `warpx/Source/Particles/WarpXParticleContainer.cpp` — the compact self-contained
  deposition example (~lines 1988–2017): `auto ptd = tile.getParticleTileData()`
  captured by value, cell located from the position, `Gpu::Atomic::AddNoRet` into
  reduction-MultiFab `Array4`s.
- `warpx/Source/Evolve/WarpXEvolve.cpp` — per-step bounded redistribution (~lines
  723–776): `max_cells_travelled` computed from the step's known motion (window moves +
  `particle_max_grid_crossings`), then `RedistributeLocal(max_cells_travelled)`;
  `warpx/Source/Particles/MultiParticleContainer.cpp` forwards it as
  `Redistribute(0, 0, 0, /*local=*/max_cells_travelled)`.
- `amrex/Src/Particle/AMReX_ParticleTile.H` — `struct ParticleTileData` (~line 32): the
  by-value kernel view (`uint64_t* m_idcpu`, per-component real/int pointer arrays),
  guarded by `static_assert(std::is_trivially_copyable<PTD>())` (~line 467);
  `getParticleTileData()` fills it after refreshing runtime pointer caches.
- `amrex/Src/Particle/AMReX_ParticleContainer.H` — the `Redistribute` contract (doc
  comment, ~lines 476–519): assigns particles to the levels, grids, and tiles
  containing their current positions (default: finest level covering the position);
  a *local* redistribute assumes particles moved only a bounded distance since the last
  call; `OK()` (~line 589) verifies placement.
- `amrex/Tests/Particles/SOAParticle/main.cpp` — the generic pure-SoA structure:
  `ParticleContainerPureSoA<4, 2, …>`, `ptd` captured into `ParallelFor`,
  `ParticleType p(ptd, ip)` proxy access.
- Anti-pattern provenance (execution side; the object-access side is owned by
  [carpetx-thorn-integration](./carpetx-thorn-integration.md)):
  `CarpetX/CarpetX/src/interp.hxx` — legacy AoS `amrex::AmrParticleContainer<3, 2>`
  (~line 38); `CarpetX/CarpetX/src/interpolate.cxx` — per-particle
  `#pragma omp parallel for simd` host loops (~lines 241–280), a literal
  `// TODO: use OpenMP` above its tile walk (~line 693), and no `amrex::ParallelFor`
  anywhere in its particle path.

This spec cites no published papers; its sources of truth are the reference
repositories above, so its citation registry is empty by construction.

## Inputs & outputs

Inputs (supplied per [carpetx-thorn-integration](./carpetx-thorn-integration.md)): the
driver's per-patch `amrex::AmrCore` (geometry, `BoxArray`, `DistributionMapping` for
container construction), and the grid-function `MultiFab`s whose `Array4` views the
gather/deposition kernels consume. Output: the packet container itself — the sole
authoritative storage of the packet population — and its tile views, consumed by every
transport kernel.

### Packet component schema (storage mapping of the PKT state table)

| Physical component (PKT) | Storage | Notes |
|---|---|---|
| position `x^i` | the container's 3 SoA position components (`ParticleReal`) | geometrized coordinates |
| momentum `p_i` | 3 SoA real components | geometrized |
| weight `N` | 1 SoA real component | an ordinary SoA component, exactly like WarpX's `w` |
| species `s` | 1 SoA integer component | values in {0, 1, 2} |
| event counter `e` | 1 SoA integer component | consumed as uint32 by the RNG packing |
| packet id `q` | the built-in 64-bit `idcpu` | see [MCNX-GPU-04] |

Component *names* are free; the mapping (one storage component per physical component,
`double` precision reals) is binding. Runtime-added components (by name, WarpX-style)
and compile-time components are equally admissible.

## Correctness requirements

- **[MCNX-GPU-01] Pure-SoA container, binding.** The packet population is stored in a
  container of the `amrex::ParticleContainerPureSoA` family (zero AoS bytes per
  particle; per-component contiguous arrays; the WarpX production container is the
  exemplar). Legacy AoS containers (`amrex::AmrParticleContainer`,
  struct-particle `ParticleContainer` flavors) are **forbidden** for the packet
  population. Rationale (binding, recorded): "GPU-friendly" is a deliverable of the
  thorn, and coalesced per-component access is the performance-bearing layout choice;
  fixed-seed golden data ties observable behavior to the population's iteration
  semantics, so the container family cannot be left to drift between iterations; and
  the reference stack's only AoS particle user — the CarpetX driver interpolator — is
  exactly the counter-example ([MCNX-GPU-07]).

- **[MCNX-GPU-02] Device-kernel execution.** Every per-packet operator — metric/fluid
  gather, geodesic push, event sampling, source-term deposition, packet-creation
  initialization — executes as device kernels over particle tiles: an
  `amrex::ParallelFor`-family launch over the particle index, whose lambda captures
  **only trivially-copyable values** (`Array4` views, `ParticleTileData`, raw SoA
  `dataPtr()` pointers, POD functors and scalars). No per-packet virtual dispatch, no
  host-side per-particle loops, and **no per-step host round-trips** of particle data
  outside `Redistribute` internals and I/O. Deposition into grid `Array4`s uses
  unordered device atomics (`amrex::Gpu::Atomic::AddNoRet`) — which is exactly why grid
  tallies carry no bitwise cross-run/cross-rank/cross-device requirement
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)). The pinned
  exemplars: gather `doGatherShapeN`; push `UpdatePosition` invoked from the tile
  kernel of `PhysicalParticleContainer.cpp`; deposition `CurrentDeposition.H` and the
  compact `WarpXParticleContainer.cpp` example; tile view
  `pti.GetParticleTile().getParticleTileData()`; generic structure
  `amrex/Tests/Particles/SOAParticle/main.cpp` (all registered in Source of truth).

- **[MCNX-GPU-03] Component schema.** The container carries exactly the physical
  content of the [packet-representation-and-sampling](./packet-representation-and-sampling.md)
  state table, per the schema mapping above: each physical component maps to exactly
  one SoA component (or the built-in `idcpu`), all reals are IEEE-754 binary64, and no
  packet state is authoritative anywhere outside the container (host mirrors for I/O
  and diagnostics are copies, never masters).

- **[MCNX-GPU-04] Packet-id uniqueness (the contract PKT and RNG defer to).** The
  packet id `q` is the AMReX 64-bit `idcpu` word (40-bit id + 24-bit cpu). Binding:
  every packet receives at creation an `idcpu` that is **unique across the entire run
  and all ranks**, is immutable for the packet's life, and is **never reused** — ids of
  removed (absorbed/escaped) packets stay retired. Bit 63 of `idcpu` is never set
  (guaranteed for cpu < 2²³ ranks; see Open questions), preserving the flag bit that
  disambiguates the RNG cell keys of
  [packet-representation-and-sampling](./packet-representation-and-sampling.md) from
  packet ids.

- **[MCNX-GPU-05] After-step ownership invariant.** At the end of every transport step,
  every packet resides on the rank, level, grid, and tile that owns its current
  position, in the sense of the AMReX `Redistribute` contract (finest level covering
  the position; periodic shifts applied where periodicity is configured). Packets whose
  positions have left the outer domain boundary are removed and tallied as escapes
  (per the `README.md` boundary policy) rather than redistributed.

- **[MCNX-GPU-06] Redistribution: bounded-motion precondition and AMR timing.** A
  *local* redistribute (`Redistribute(…, local = N)` / `RedistributeLocal(N)`) is
  admissible **only** under the proven precondition that no packet has moved more than
  `N` cells in any direction since the last redistribution; `N` must be derived from an
  actual per-step motion bound (e.g. the step's maximum coordinate speed × Δt over the
  minimum cell width on the packet's level, plus a grid-crossing margin — the WarpX
  `max_cells_travelled` precedent), not assumed. Absent such a bound, the global
  `Redistribute()` must be used. Timing under AMR (single-rate stepping per
  [carpetx-thorn-integration](./carpetx-thorn-integration.md)): redistribution runs at
  least once per transport step, after the step's last position-changing operation and
  before any subsequent operation that relies on [MCNX-GPU-05]; and after any regrid or
  level change, packets are redistributed onto the new grid hierarchy before the next
  transport operation touches them.

- **[MCNX-GPU-07] Anti-pattern warning (binding, stated identically in both technical
  specs).** The CarpetX driver interpolator is cited for **object access only**; its
  legacy AoS particle container and host-side per-particle loops
  (`CarpetX/CarpetX/src/interp.hxx`, `CarpetX/CarpetX/src/interpolate.cxx`) are an
  anti-pattern for GPU particle operators and are not adopted. It is the nearest
  in-stack particle code, and a fresh iteration must not mistake it for the pattern to
  copy: no MCNuX per-packet operator may follow its execution model (host `ParConstIter`
  loops, `#pragma omp parallel for simd` per-particle kernels, AoS particle storage).
  What *is* validly taken from it — where the driver's grid objects live — is owned by
  [carpetx-thorn-integration](./carpetx-thorn-integration.md).

### Conventions restated (the subset this leaf uses)

- Packet physical components, units, and valid ranges are owned by
  [packet-representation-and-sampling](./packet-representation-and-sampling.md); this
  spec maps them to storage. Species indices {0 (νe), 1 (ν̄e), 2 (νx)}.
- All reals are IEEE-754 binary64 (`double`), per
  [conventions-and-units](./conventions-and-units.md).
- All physics draws are evaluations of the pure function `u(S, q, e, k)`; AMReX's
  stateful per-thread device RNG is forbidden for physics draws
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)). Kernel
  structure may batch or cache draws; observable values must equal the pure function.

## Verification

- **Layout pin (compile-time, exact).** The packet container type is (or derives from)
  an `amrex::ParticleContainerPureSoA` instantiation; a compile-time check
  (`static_assert` / type trait) in the test suite asserts this and asserts the tile
  view is trivially copyable (mirroring AMReX's own
  `static_assert(std::is_trivially_copyable<PTD>())`).
- **Device execution (build + inspection).** The suite builds and runs with a GPU
  backend, with every per-packet operator compiled as a device kernel; inspection
  confirms the by-value-POD capture rule and the absence of per-step host round-trips
  of particle data (the [MCNX-GPU-02] exceptions aside). This is the verification mode
  the execution requirement admits — there is no runtime oracle for "ran on device."
- **Ownership invariant (exact).** A pinned multi-box, multi-tile test (≥ 2 boxes per
  rank, packets crossing box boundaries during the step) asserts after each step that
  the container's placement check passes (`OK()`—style: every packet's position is
  contained in the box/tile it is stored in, on the finest level covering it).
- **Bounded-motion audit (exact).** In every deterministic benchmark, the
  implementation's derived per-step motion bound `N` is recorded and compared against
  the measured maximum per-step cell displacement; measured ≤ derived must hold in
  every step of every benchmark. A local redistribute configured with the derived bound
  passes the ownership check above.
- **Id-contract test (exact).** Over a pinned multi-step run with creation and removal,
  all live ids are pairwise distinct, no removed id reappears, ids are immutable, and
  bit 63 is clear on every `idcpu`.
- **Single-rank determinism (exact tier).** Two single-rank CPU runs at fixed seed
  produce bitwise-identical packet populations after every step (the
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) reproducibility
  contract exercised through this container); golden benchmark outputs are archived at
  the Cactus harness defaults `ABSTOL=RELTOL=1e-12`. Grid tallies deposited by
  unordered atomics are compared across launch configurations only with norm-based
  tolerances (owned by the comparing spec), never bitwise.

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the literal `ParticleContainerPureSoA` and the anti-pattern phrase
`are an anti-pattern for GPU particle operators and are not adopted`; that its cited
`warpx/`, `amrex/`, and `CarpetX/` paths resolve; and that its `[MCNX-GPU-NN]` ids are
declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Kernel fusion (one fused transport kernel vs separate gather/push/sample/deposit
  launches), tiling strategy, launch shape, block sizes, and streams.
- Component names, ordering, compile-time vs runtime-added components, and any extra
  *diagnostic* components (which must not carry authoritative physics state beyond the
  schema).
- The iteration machinery (`ParIter` flavor, direct tile loops), particle sorting or
  binning for locality, and buffer growth policy.
- Local vs global `Redistribute` per step — subject to [MCNX-GPU-06]'s precondition —
  and how the motion bound is computed or over-estimated.
- Guard-cell buffering strategy for deposition (direct atomics into the MultiFab vs
  per-tile buffers merged afterward), provided tallies land in the owning cells and the
  conservation ledger of
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) holds.
- Host-side mirrors for I/O, checkpointing, and debugging (non-authoritative copies).

## Open questions / assumptions

- **Assumption: rank-count headroom for bit 63.** The AMReX `idcpu` packing sets bit 63
  only when the creating cpu index reaches 2²³; MCNuX assumes < 8388608 ranks, so the
  RNG cell-key flag bit stays unambiguous. Exceeding this would require a deliberate,
  golden-data-invalidating re-pin of the RNG key packing (owned by
  [packet-representation-and-sampling](./packet-representation-and-sampling.md)).
- **Assumption: `double`-precision particle reals.** The build pins
  `amrex::ParticleReal` (and `amrex::Real`) to `double`, matching the corpus value
  type; the build contract is owned by
  [build-and-integration](./build-and-integration.md).
- **Assumption: per-patch container organization.** On multipatch domains the packet
  population is organized per patch (a container keyed to each patch's AmrCore, the
  driver interpolator's object-access precedent); verification baselines run
  single-patch. A cross-patch transfer contract (packets crossing patch boundaries) is
  not yet specified and would extend this spec, not silently appear.
- **Open question: redistribution under subcycling.** [MCNX-GPU-06]'s timing assumes
  single-rate stepping (`CarpetX::use_subcycling = no`, pinned in
  [carpetx-thorn-integration](./carpetx-thorn-integration.md)). Subcycled transport
  would need per-level redistribution cadence rules specified here before adoption.
