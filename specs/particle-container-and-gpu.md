# Particle container and GPU execution (technical leaf: where packets live and how they move across ranks, levels, and devices)

> Technical leaf spec. Self-contained: an agent can implement the packet-container layer — the observable component schema, the residency and redistribution invariants, the deposition contract, and the checkpoint/restart guarantees — from this file alone, referencing `conventions-and-units.md` for packet nomenclature and units, `packet-representation-and-sampling.md` for the physical state contract this schema stores, `geodesic-propagation.md` for the bounded-motion guarantee consumed here, and `rng-and-statistical-acceptance.md` for the reproducibility bounds and (q, e) semantics. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines the container layer that carries Monte Carlo packets on the AMReX grid hierarchy: the observable per-packet component schema (what is stored, at what precision, and what must survive every container operation), the residency invariant tying each packet to the rank/level/grid/tile owning its position, the bounded-motion precondition that licenses local redistribution, the deposition contract by which per-event tallies become cell-centered grid values (atomic scatter into guard cells reconciled by guard-cell folding), and the checkpoint/restart guarantees including dual-grid restart. It constrains *observable* container behavior only; container flavor and kernel structure are implementation freedom.

In scope:

- The component schema as an observable layout contract: components, storage precision, and the round-trip/communication/checkpoint completeness rules.
- The residency invariant and redistribution: the owner definition (finest covering level), the observable point at which the invariant holds, and the bounded-motion precondition any local redistribution relies on.
- The deposition contract: event-to-cell attribution, guard-cell folding, multi-level accounting, and the reproducibility bound tallies are held to.
- Checkpoint/restart: population survival, (q, e) preservation, single-rank CPU determinism, dual-grid restart.
- The GPU-execution constraints, which are deliberately few: physics-affecting draws obey the RNG contract, and tallies obey the reduction-tolerance bound — everything else about device execution is freedom.

Out of scope:

- The physical meaning, creation-side validity, and invariants of the packet components (`packet-representation-and-sampling.md`; `conventions-and-units.md` for nomenclature and units).
- *What values* are deposited — the source-term assembly S_μ, S_L is owned by `hydro-coupling-source-terms.md`; this spec pins that whatever is deposited arrives completely and correctly.
- The grid-access idiom, schedule placement, execution modes, and centering bridges (`carpetx-thorn-integration.md`).
- The random-number generator and the packet-id allocation policy (`rng-and-statistical-acceptance.md`; uniqueness of q is MCNX-RNG-05).
- Load balancing, sorting policy, and performance tuning (out of scope corpus-wide; the WarpX practices cited below are provenance, not requirements).

## Source of truth

- `amrex/Src/Particle/AMReX_ParticleContainer.H` — the container template (AoS, mixed, and pure-SoA layouts behind one interface; compile-time plus named runtime components) and the `Redistribute` semantics: default assignment moves every particle to the rank/level/grid/tile owning its position, the owning level being the finest level covering it; the `local`/`max_cells_moved` variant is documented as valid only when particles have moved no farther than the stated cell distance since the last redistribution.
- `amrex/Src/Particle/AMReX_ParticleContainerI.H`, `amrex/Src/Particle/AMReX_ParticleLocator.H`, `amrex/Src/Particle/AMReX_ParticleTile.H`, `amrex/Src/Particle/AMReX_ParIter.H` — the redistribution implementation, the bin-based position→grid locator (the owner definition made executable), tile storage/device views, and iteration.
- `amrex/Src/Particle/AMReX_ParticleMesh.H` — the particle-to-mesh scatter path: deposition into a (guard-cell-grown) FAB with atomic adds on GPU, then guard-cell contributions folded into valid cells by `SumBoundary` (same grids) or `ParallelAdd` (temporary MultiFab).
- `amrex/Src/Particle/AMReX_ParticleIO.H` and `amrex/Src/Particle/AMReX_WriteBinaryParticleData.H` — checkpoint/restart: component-complete particle data plus id bookkeeping written per level; `Restart` reads on possibly different grids and distribution mappings and ends with an unconditional `Redistribute()` — dual-grid restart is an explicitly supported guarantee (exercised by the `amrex/Tests/Particles/` dual-grid restart tests).
- `warpx/Source/Particles/WarpXParticleContainer.H` and `warpx/Source/Evolve/WarpXEvolve.cpp` — the at-scale exemplar (provenance for the shape of this contract, not normative): a pure-SoA container with weight as an ordinary component; a CFL-bounded `RedistributeLocal(max_cells_travelled)` every step at single level, falling back to a full `Redistribute()` when mesh-refinement levels exist; atomic deposition into guard cells followed by a guard-cell sum. `warpx/Source/Initialization/WarpXAMReXInit.cpp` — tiling disabled on GPU (an implementation choice, cited as freedom precedent).
- `packet-representation-and-sampling.md` — the packet state contract (components, types, units, valid ranges) that the schema below stores; `geodesic-propagation.md` — the per-substep cell-crossing cap (|Δx^i| ≤ one cell width per direction, MCNX-GEO-04) that makes per-step motion bounded; `rng-and-statistical-acceptance.md` — the grid-tally reproducibility bound (MCNX-RNG-04) and packet-id uniqueness (MCNX-RNG-05).

## Inputs & outputs

### Observable component schema (binding)

Every packet in the population carries exactly these components, restated from the state contract of `packet-representation-and-sampling.md` (which owns their physical validity) with the storage obligations added here:

| Component | Symbol | Stored as | Units | Schema obligation |
|---|---|---|---|---|
| position | x^i | 3 × `double` (full IEEE-754 binary64) | code length | the position the locator assigns ownership from |
| four-momentum | p_t, p_i | 4 × `double` | MeV | stored per packet; the p_t consistency identity is owned by the packet spec |
| weight | N | `double` | dimensionless | no precision reduction at any point |
| species | s | integer component | — | observable value ∈ {0, 1, 2}, exact |
| RNG identity | (q, e) | 2 × 64-bit unsigned integers | — | never truncated, re-derived, or recycled |

Schema rules (binding):

1. **Full-precision round trip.** Storing a packet and reading it back yields bitwise-identical component values; no component is stored at reduced precision, quantized, or reconstructed lossily from others.
2. **Communication completeness.** Every schema component travels with the packet through every redistribution (whatever mechanism moves packets between tiles, grids, levels, or ranks ships all components bitwise — for named runtime components this means their communicate flag is enabled).
3. **Checkpoint completeness.** Every schema component is written to and restored from checkpoints bitwise, including (q, e) (checkpoint/restart preserves the RNG identity per `rng-and-statistical-acceptance.md`).
4. Additional *internal* per-packet components (diagnostics, cached values) are permitted, provided they are reconstructible or explicitly non-normative and their presence never alters the observable components.

### Residency invariant and redistribution

- **Owner definition.** The owner of a packet is the (MPI rank, AMR level, grid, tile) that the container's position-based assignment gives its x^i, with the owning level the **finest level covering the position** (the `amrex/Src/Particle/AMReX_ParticleContainer.H` default assignment; the locator of `amrex/Src/Particle/AMReX_ParticleLocator.H` is the executable owner definition).
- **Invariant (binding, observable point).** At the end of each transport step — from the point MCNuX's per-step transport phase completes until the next step's transport phase begins — every packet in the population, including packets created by emission during that step, resides at its owner. During substeps within the transport phase the invariant may be transiently violated; no consumer observes the population mid-phase.
- **Bounded-motion precondition (binding).** The per-substep cell-crossing cap of `geodesic-propagation.md` (MCNX-GEO-04: |Δx^i| ≤ one cell width per direction per substep) bounds a packet's total per-step motion by n_sub cells per direction, n_sub the number of substeps it took. Any use of a *local* (bounded-distance) redistribution is licensed only by this bound: the distance parameter handed to it (the `max_cells_moved` sense of `amrex/Src/Particle/AMReX_ParticleContainer.H`) must be ≥ the actual maximum per-direction cell displacement since the previous redistribution. Local and global redistribution must be observably equivalent: same final assignments, same populations. (Provenance: WarpX's CFL-justified `RedistributeLocal` in `warpx/Source/Evolve/WarpXEvolve.cpp`, which itself falls back to a global `Redistribute()` when refinement levels exist — choosing local vs. global is freedom; the precondition and equivalence are not.)
- **Regrid obligation.** After any change of the grid hierarchy (regridding, load-balance redistribution of boxes), a global redistribution must restore the invariant before the next transport phase; the bounded-motion license does not survive a hierarchy change.

### Deposition contract

Per-event tallies (the deposits defined by `neutrino-matter-interactions.md` and assembled per `hydro-coupling-source-terms.md`) reach the cell-centered grid as follows (binding):

- **Attribution.** Each event contribution deposits into the cell containing the event position, on the packet's owning (finest covering) level — the same attribution rule MCNX-SRC-02 asserts on values.
- **Guard-cell folding.** Contributions scattered into guard cells of a tile/box (an event near a box boundary) are folded into the owning valid cells — the `SumBoundary` / `ParallelAdd` reconciliation of `amrex/Src/Particle/AMReX_ParticleMesh.H` — such that no contribution is lost or double-counted.
- **Multi-level accounting.** In a multi-level run, the sum of deposits over the hierarchy with fine-covered coarse cells excluded (fine-mask convention) equals the per-event ledger; where the consumer-visible level is coarser than the deposit level, the fine-level deposits reach it through the driver's conservative restriction of cell-centered groups (`CarpetX/CarpetX/src/sync_restrict.cxx`), preserving the volume-weighted total.
- **Tolerance.** Deposited totals are deterministic tally identities at the relaxed tier `rtol = 1e-10` (summation-order freedom: atomics, per-thread partials); repeated identical parallel runs agree in tally norms to relative `1e-10` per MCNX-RNG-04. Bitwise deposition order is explicitly *not* required.

### Checkpoint/restart

- A checkpoint captures the complete packet population (all schema components, bitwise) together with whatever id-allocation state is needed so that packet-id uniqueness (MCNX-RNG-05) holds across the restart — no post-restart packet may receive a q ever used before the checkpoint.
- **Same-grid restart, single-rank CPU:** bitwise state resumption — the continued run's per-packet event sequences and outputs are identical to an uninterrupted run (this is the checkpoint leg of MCNX-RNG-03 and MCNX-SRC-03, gated here as MCNX-PAR-05).
- **Dual-grid restart (binding, per the AMReX guarantee):** restarting onto a different BoxArray/DistributionMapping (different rank count, different max grid size) is supported; after restart the census content (the multiset of packets and their components) is unchanged and the residency invariant holds. Physical results after a dual-grid restart agree with the same-grid run to the tally-reproducibility bound (`1e-10` norms), not bitwise.

## Correctness requirements

- **[MCNX-PAR-01] Schema completeness and round trip (exact).** Every packet component of the schema table round-trips bitwise through container storage, through every redistribution (single- and multi-rank), and through checkpoint/write-read: for a synthetic population with distinguishable per-packet values, all components of all packets are bitwise-identical afterward, with s and (q, e) integer-exact. Reference: the schema rules above.
- **[MCNX-PAR-02] Residency invariant (exact).** In an instrumented multi-level fixture, at the end of every transport step an exhaustive owner sweep (recomputing each packet's owner from its position via the locator definition) finds zero packets away from their owning (rank, level, grid, tile), with the owning level the finest level covering the position; packets emitted during the step are included. Reference: the invariant above; provenance `amrex/Src/Particle/AMReX_ParticleContainer.H` default assignment.
- **[MCNX-PAR-03] Bounded motion and local/global equivalence (exact).** In instrumented fixed-seed runs: the recorded maximum per-direction cell displacement of every packet between redistributions never exceeds the bound licensed by the substep count (and never exceeds the distance parameter of any local redistribution used); and the same fixed-seed step finished with a local redistribution and, separately, with a global one yields identical final assignments and bitwise-identical populations. Reference: the bounded-motion precondition; provenance `warpx/Source/Evolve/WarpXEvolve.cpp`, MCNX-GEO-04.
- **[MCNX-PAR-04] Deposition integrity (relaxed `1e-10`).** For fixed-seed instrumented runs with independent per-event logging: (a) grid totals equal the per-event ledger at `rtol = 1e-10`, including events depositing through tile guard cells (verified by placing emitters within one cell of box boundaries) — zero lost or double-counted contributions beyond summation rounding; (b) in the two-level fixture, the fine-mask hierarchy sum equals the ledger at `rtol = 1e-10`, and restricted coarse-level values over refined regions equal the volume-weighted average of the fine deposits at `rtol = 1e-10`; (c) repeated identical parallel runs agree in tally norms to relative `1e-10` (MCNX-RNG-04). Reference: the deposition contract above.
- **[MCNX-PAR-05] Checkpoint/restart determinism and dual-grid restart (golden `1e-12` + exact).** A fixed-seed single-rank CPU run with a mid-run checkpoint/restart produces per-packet TSV output identical to the uninterrupted run — harness diff at `ABSTOL = RELTOL = 1e-12` with zero violations. A dual-grid restart leg (same checkpoint, different rank count and max grid size) preserves the census content exactly (multiset equality of all schema components), satisfies the residency invariant immediately after restart, and assigns no duplicate q for the remainder of the run. Reference: the checkpoint/restart contract above; provenance `amrex/Src/Particle/AMReX_ParticleIO.H`.

## Verification

- **MCNX-PAR-01**: family-B unit check: build a synthetic population (≥ `1e4` packets, per-packet distinguishable values in every component, spanning tiles and levels), push it through store/read-back, a 1-rank and a 2-rank redistribution, and a checkpoint write/read; assert bitwise equality per component per packet.
- **MCNX-PAR-02**: family-B check on a two-level fixture (base grid with a refined central box): ≥ 4 fixed-seed transport steps at 1 and 2 ranks; after each step run the owner sweep and assert zero misplacements, including for packets that crossed the fine/coarse boundary in either direction and packets created that step. (This is the corpus's redistribution-invariant check after a multi-level transport step.)
- **MCNX-PAR-03**: the same instrumented runs record per-packet displacement between redistributions and assert the bound; the equivalence leg reruns one step from an identical checkpointed state under the local and the global path and asserts identical assignments and bitwise-identical populations.
- **MCNX-PAR-04**: family-B check with per-event logs on (a) a single-level fixture whose emitting cells sit within one cell of box boundaries, (b) the two-level fixture, (c) ≥ 3 repeated threaded (and, where a device build exists, GPU) runs; assert the three `1e-10` identities and norm agreements. GPU legs report SKIP loudly where no device build exists.
- **MCNX-PAR-05**: family-A Cactus regression (fixed seed, single-rank CPU, ≥ 4 steps, checkpoint at the midpoint, golden per-packet TSV at `ABSTOL = RELTOL = 1e-12`) plus a family-B dual-grid leg asserting census multiset equality, the post-restart owner sweep, and q uniqueness via id enumeration (shared instrumentation with MCNX-RNG-05).
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `amrex/`, `warpx/`, and `CarpetX/` paths, and contains the required claim strings (the finest-covering-level owner rule, the bounded-motion license, the guard-cell folding names, the dual-grid guarantee, each tolerance).

## Implementation freedom

- **Container flavor** — AoS, mixed, or pure-SoA layout; compile-time vs. named runtime components; the concrete container type (an `amrex::AmrParticleContainer` specialization per `amrex/Src/AmrCore/AMReX_AmrParticles.H`, a pure-SoA container per the WarpX precedent, or another realization) — provided the observable schema rules hold.
- **Redistribution strategy** — local vs. global, once per step vs. per substep, and the distance parameter — provided the residency invariant holds at the observable point and every local use satisfies the bounded-motion license.
- **Kernel structure** — fusion of push/event/deposit work into one `ParallelFor` per tile (the WarpX idiom) or separate kernels; tiling on or off; particle sorting/binning and its cadence — provided per-packet observables are pure functions of (S, q) per the RNG contract and tallies meet MCNX-PAR-04.
- **Deposition mechanics** — GPU atomics into grown FABs, CPU thread-private tiles merged afterward, or serial accumulation — bounded by MCNX-PAR-04 and MCNX-RNG-04.
- **Checkpoint format and mechanics** — AMReX-native particle checkpoint files or driver-managed I/O — provided MCNX-PAR-05's observables (bitwise same-grid resumption, dual-grid census/residency/uniqueness) hold.
- Packet-id allocation (rank-prefixed counters, id-header bookkeeping, reserved ranges) — owned by `rng-and-statistical-acceptance.md`'s uniqueness requirement; any scheme satisfying it is acceptable here.

## Open questions / assumptions

- **Subcycled fine-level transport is out of scope (assumption, shared with the corpus).** The verified configurations step all levels at the coarsest Δt (`carpetx-thorn-integration.md` pins `use_subcycling = no`); the residency and deposition contracts above are stated for single-rate stepping. Subcycled transport would need per-level redistribution timing (the `nGrow` machinery of `amrex/Src/Particle/AMReX_ParticleContainer.H`) and a revised multi-level accounting — a spec change here and in the CarpetX-integration spec.
- **Load balancing is permitted but unpinned (assumption).** A weight-based DistributionMapping change (WarpX-style) is an allowed hierarchy change handled by the regrid obligation (global redistribution before the next transport phase); no cost model or interval is specified, and no requirement above depends on load balance.
- **Reflux-style corrections for source terms are not required (assumption, recorded).** MCNX-PAR-04 requires conservative restriction of deposited sources; no flux-register correction at fine/coarse boundaries is specified, matching the operator-split (non-flux) character of the deposits. If a consumer demonstrates a need, the multi-level accounting gains a reflux clause here.
- **Guard-region width (assumption).** The deposition contract assumes the source-term MultiFabs carry at least one guard cell so near-boundary events can scatter locally before folding; an implementation depositing directly into valid cells with its own reconciliation is acceptable provided MCNX-PAR-04 holds unchanged.
