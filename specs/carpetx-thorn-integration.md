# CarpetX thorn integration (technical leaf: how MCNuX lives inside the driver)

> Technical leaf spec. Self-contained: an agent can implement the thorn/driver integration surface — the thorn shape and capability declarations, the schedule placement and cadence, the execution-mode rules, the pinned AMReX-native grid-access idiom, and the metric/fluid read contracts with their centering bridge — from this file alone, referencing `conventions-and-units.md` for the 3+1 variable mapping and centering facts, `hydro-coupling-source-terms.md` for the source-term protocol whose schedule realization is pinned here, `geodesic-propagation.md` for the per-step transport semantics being placed, and `rng-and-statistical-acceptance.md` for the reproducibility bounds the regressions ride on. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines how MCNuX exists as a Cactus thorn on the CarpetX driver: the four-ccl thorn shape and its capability obligations (`REQUIRES CarpetX`), the mandatory READS/WRITES declarations that keep the driver's validity machinery sound, the pinned AMReX-native grid-access idiom (driver-owned MultiFabs reached through `ghext`, iterated with `MFIter`/`ParIter`, read/written through `Array4` views inside `amrex::ParallelFor`), the execution-mode rule for routines that walk the container or the grid themselves, the schedule placement of the transport step (own group at `CCTK_EVOL` after `ODESolvers_Solve`) with the cadence stated observably, the single-rate configuration pin (`use_subcycling = no`), and the metric/fluid data-access contracts against ADMBaseX (vertex) and HydroBaseX (cell) including the centering-bridge obligation.

In scope:

- Thorn shape: the four ccl files, the capability declarations, and the grid-variable declarations MCNuX owes the driver.
- The access idiom (binding): driver-native AMReX objects, not the Loop thorn's point-loop layer, for all transport and deposition data paths.
- Schedule placement and the observable cadence contract, including the schedule realization of the zero-then-add protocol and the optional `TmunuBaseX_AddToTmunu` routine.
- Execution modes for scheduled routines; READS/WRITES/SYNC soundness under the driver's enforcement machinery.
- The single-rate stepping pin and its parameter guard.
- Metric/fluid read contracts and the centering bridge, with its machine-tier linear-exactness requirement.

Out of scope:

- Build and delivery: how CarpetX, AMReX, and WeakLibInterp are obtained and linked, compiler handling of GPU code (the build-and-integration spec).
- The transport physics and the source-term values (the domain leaf specs); container mechanics, residency, and deposition tolerances (`particle-container-and-gpu.md`).
- The GRMHD consumer's own schedule; only MCNuX's side of the coupling surface is pinned (`hydro-coupling-source-terms.md` owns the variable contract).
- The Loop thorn's `GF3D2`/`loop_*_device` API — deliberately unused for transport data access (see the pinned idiom below); its existence is noted only to exclude it.

## Source of truth

- `CarpetX/CarpetX/src/driver.hxx` — the driver's data model this spec pins access to: the global `ghext`; per-patch `patchdata` with `amrcore` (a `CactusAmrCore final : public amrex::AmrCore`, supplying `Geom(lev)`, `boxArray(lev)`, `DistributionMap(lev)`); per-(patch, level) `leveldata` with per-group `groupdata` holding the grid-function storage as `amrex::MultiFab`s in `mfab` (indexed by time level, one component per group variable) plus the group's `indextype` (centering).
- `CarpetX/CarpetX/src/interpolate.cxx` — the driver's own interpolation service: the in-repo demonstration of exactly the pinned idiom (a particle container built on `patchdata.amrcore.get()`, a `ParConstIter` walk per level, `groupdata.mfab.at(tl)->array(pti)` yielding the `Array4` view consumed in the kernel, `groupdata.indextype` consulted for centering).
- `CarpetX/CarpetX/src/schedule.cxx` and `CarpetX/CarpetX/src/valid.hxx` — the `CallFunction` dispatcher enforcing declared READS/WRITES against per-variable validity state, with poisoning (`poison_undefined_values`), optional checksums, and presync; local-mode routines are dispatched per box/tile, level/global-mode routines get a single call.
- `CarpetX/ODESolvers/schedule.ccl` — `ODESolvers_Solve` scheduled `AT evol` with `OPTIONS: level`; the RK substeps and `ODESolvers_PostStep` run inside it, so "after `ODESolvers_Solve`" means after the completed fluid/spacetime step including post-step processing.
- `flesh/src/schedule.ccl` — the flesh time bins (`CCTK_EVOL`, `CCTK_PARAMCHECK`, …) MCNuX schedules against.
- `CarpetX/TestLoopX/configuration.ccl` — precedent for `REQUIRES CarpetX` (the Loop capability alone does not expose `ghext`); `CarpetX/Loop/src/loop_device.hxx` — the point-loop layer this spec excludes from the transport data path.
- `CarpetX/ADMBaseX/interface.ccl` and `CarpetX/HydroBaseX/interface.ccl` — the metric (vertex-centered) and fluid (cell-centered) variables MCNuX reads; `CarpetX/TmunuBaseX/schedule.ccl` — the `TmunuBaseX_AddToTmunu` hook group the optional stress-energy routine is scheduled `IN`.
- `CarpetX/CarpetX/param.ccl` — `use_subcycling` (the parameter the single-rate pin guards on) and the enforcement-machinery parameters the soundness regression enables.
- `hydro-coupling-source-terms.md` — the zero-then-add protocol and cadence contract whose schedule realization is pinned here; `geodesic-propagation.md` — the metric time-slice policy (one slice per transport step) the placement must honor.

## Inputs & outputs

### Thorn shape and capability declarations (binding)

- MCNuX is a Cactus thorn of the standard four ccl files: `interface.ccl` (grid-variable groups with explicit `CENTERING` — the source-term groups cell-centered `CCC` per `hydro-coupling-source-terms.md`), `schedule.ccl` (routines with complete READS/WRITES/SYNC), `param.ccl` (run parameters incl. seed S, N_tgt, `rad_tmunu`, scheme selectors), `configuration.ccl` (capabilities).
- `configuration.ccl` **must** declare `REQUIRES CarpetX`: the pinned access idiom reaches driver internals (`ghext`), which the Loop capability alone does not expose (`CarpetX/TestLoopX/configuration.ccl` is the precedent). Further capability requirements (WeakLibInterp, …) are owned by the build-and-integration spec.

### The pinned grid-access idiom (binding; design decision Q9)

All transport-phase grid access — metric/fluid gathers, source-term zeroing and deposition, container construction — uses the driver's native AMReX objects:

```text
per patch:                ghext->patchdata[patch].amrcore        # amrex::AmrCore: Geom, boxArray, DistributionMap
per (patch, level, group): ghext->patchdata[patch].leveldata[level].groupdata[gi]->mfab[tl]
                                                                  # amrex::MultiFab, one component per variable
per box/tile:             MFIter / ParIter walk; data as amrex::Array4 views (mfab[tl]->array(mfi))
kernels:                  amrex::ParallelFor with POD views captured by value; centering from groupdata.indextype
```

exactly as the driver's own interpolation service does (`CarpetX/CarpetX/src/interpolate.cxx`). The Loop thorn's `GF3D2`/`p.I`/`loop_*_device` point-loop layer is **not** used on any transport or deposition data path (it is designed for grid-point-per-thread stencil updates in local mode and does not fit particle-tile iteration); incidental grid-only diagnostics may use it, but nothing normative may depend on it.

Reaching data through `mfab` does not relax the declaration obligation: READS/WRITES remain mandatory (below) because the driver's validity machinery enforces them regardless of how the pointer is obtained.

### Schedule placement and cadence (binding; design decision Q6)

- **Placement.** MCNuX's transport phase is a thorn-owned schedule group placed `AT CCTK_EVOL` (bin `evol`) `AFTER ODESolvers_Solve`. Since the RK substeps and `ODESolvers_PostStep` (where ADM and hydro variables are finalized by their producers) run inside `ODESolvers_Solve` (`CarpetX/ODESolvers/schedule.ccl`), the transport phase sees the completed, post-step metric and fluid state of the new slice. That slice is the "slice t_n" of `geodesic-propagation.md`'s time-slice policy: one transport step, one slice, never updated mid-step.
- **Cadence (normative, stated observably).** Packet positions/momenta and the source-term grid variables advance exactly once per coarsest-level Δt — one transport phase per `CCTK_EVOL` traversal, never per RK substage. Combined with the zero-then-add protocol of `hydro-coupling-source-terms.md`, whose zero/add phases both live inside this transport group: **the source terms read by the fluid at step n are those tallied during step n−1's transport** (the fluid's RHS reads, executing inside the *next* `ODESolvers_Solve`, see the previous transport phase's completed values).
- **Optional Tmunu routine.** When `rad_tmunu` is enabled, the stress-energy contribution routine is scheduled `IN TmunuBaseX_AddToTmunu` (`CarpetX/TmunuBaseX/schedule.ccl`), inheriting that group's guarantees (prior zeroing, current ADM data); its accumulation semantics are owned by MCNX-SRC-05.
- **Ancillary placements.** Container/table initialization at startup/basegrid-time bins; the parameter guard at `CCTK_PARAMCHECK` (below); analysis/diagnostic output at the analysis bin. Exact routine decomposition within these placements is freedom.

### Execution modes (binding)

Every scheduled MCNuX routine that itself walks `MFIter`/`ParIter` (transport, deposition, zeroing, container maintenance) is scheduled level- or global-mode (`OPTIONS: level` or `OPTIONS: global`), never local-mode: the driver's local mode already dispatches per box/tile (`CarpetX/CarpetX/src/schedule.cxx`), and nesting an iterator walk inside it would multiply or misroute the work. Point-wise grid-only routines (if any) may be local-mode.

### READS/WRITES soundness (binding)

Every scheduled routine declares complete `READS:` and `WRITES:` clauses with region specifiers (and `SYNC:` where ghost validity is needed) for every grid variable it touches — ADMBaseX metric/lapse/shift reads, HydroBaseX rho/temperature/Ye/vel reads, MCNuX source-term writes, and (when enabled) TmunuBaseX accumulation reads+writes. The declarations must be sound under the driver's full enforcement: presync, NaN-poisoning of undeclared data (`poison_undefined_values`), and checksum detection of undeclared writes (`CarpetX/CarpetX/src/valid.hxx`).

### Single-rate configuration pin (binding; verified configuration)

The verified configuration steps all levels at the coarsest Δt: `CarpetX::use_subcycling = no`. MCNuX must refuse to run otherwise: a parameter check at `CCTK_PARAMCHECK` aborts with a diagnostic naming this spec when `use_subcycling` is enabled. Subcycled transport is an open question below, not a silently-degraded mode.

### Metric/fluid access and the centering bridge (binding)

- Metric data (α, β^i, γ_ij) is vertex-centered; fluid data (ρ, T, Ye, v^i) and MCNuX's source terms are cell-centered (`conventions-and-units.md` centering table). Two bridge situations arise: evaluating vertex-centered metric data at cell centers (emission's α√γ, Tmunu lowering) and at packet positions (geodesic gathers).
- The bridge is a deterministic, documented interpolation. **Linear-exactness requirement (machine tier):** the vertex-to-cell-center bridge and the gather-at-position bridge must reproduce any field that is linear in the coordinates exactly (to `rtol = 1e-14`) — the 8-point average of `hydro-coupling-source-terms.md` and standard trilinear interpolation both satisfy this; a lower-order or biased bridge does not. Accuracy beyond linear exactness is bounded observably by the Schwarzschild benchmark (MCNX-GEO-06), not prescribed here.
- Reads honor the producing thorns' validity regions: interior reads plus declared SYNCs, per the READS/WRITES soundness rule; MCNuX never writes any ADMBaseX or HydroBaseX variable.

## Correctness requirements

- **[MCNX-CTX-01] Thorn shape and capability (exact, structural).** The four ccl files exist; `configuration.ccl` contains `REQUIRES CarpetX`; every MCNuX grid group declares an explicit centering matching its spec (`CCC` for `rad_force`/`lep_source`); no transport or deposition code path depends on the Loop point-loop layer. Verified by structural audit (fixed-string and shape checks). Reference: the thorn-shape and access-idiom sections above.
- **[MCNX-CTX-02] Declaration soundness under enforcement (exact).** A full fixed-seed run (transport, deposition, and the optional Tmunu leg enabled) with the driver's enforcement active — presync mode, `poison_undefined_values = yes`, checksums enabled — completes with zero validity errors, zero reads of poisoned data, and zero checksum violations. Reference: the READS/WRITES section; provenance `CarpetX/CarpetX/src/schedule.cxx`.
- **[MCNX-CTX-03] Execution-mode conformance (exact).** Every scheduled routine that walks `MFIter`/`ParIter` carries `OPTIONS: level` or `OPTIONS: global`; an instrumented run shows each such routine invoked exactly once per (level or run) per traversal, never per tile. Verified by schedule audit plus invocation counters. Reference: the execution-mode rule above.
- **[MCNX-CTX-04] Cadence observability (exact + golden `1e-12`).** In an instrumented ≥ 4-step run: the transport phase executes exactly once per `CCTK_EVOL` traversal (counter equality with the iteration count); packet state and source-term variables change during transport phases only; and the one-step lag holds — the shared fixed-seed regression of MCNX-SRC-03 (source values at the read point vs. committed golden at `ABSTOL = RELTOL = 1e-12`) passes with the schedule realization pinned here. Reference: the placement/cadence section above.
- **[MCNX-CTX-05] Single-rate guard (exact).** With `CarpetX::use_subcycling = yes`, the run aborts during `CCTK_PARAMCHECK` with a diagnostic identifying the unsupported configuration; zero transport steps execute. With `use_subcycling = no` the guard is silent. Reference: the configuration pin above.
- **[MCNX-CTX-06] Centering-bridge linear exactness (machine).** For probe fields f(x) = a + b_i x^i installed in vertex-centered (and, for the gather leg, cell-centered) groups: the vertex-to-cell-center bridge reproduces f at every cell center, and the gather-at-position bridge reproduces f at ≥ `1e3` random probe positions (interior, away from boundaries by one cell), each at `rtol = 1e-14`; repeated evaluation is bitwise deterministic. Reference: the centering-bridge section above.

## Verification

- **MCNX-CTX-01**: a family-B structural check (script or unit test) over the thorn tree: file presence, `REQUIRES CarpetX` in `configuration.ccl`, centering strings in `interface.ccl`, and a source-tree scan asserting no Loop-layer accessor appears on transport/deposition paths.
- **MCNX-CTX-02**: a family-A run of the standard fixed-seed transport fixture with presync, poisoning, and checksums enabled in the parameter file; pass = clean completion (the driver hard-errors on violations, so completion is the assertion) plus the harness diff of its outputs.
- **MCNX-CTX-03**: family-B: parse the schedule declarations for the mode options and run with invocation counters per routine; assert the per-traversal counts.
- **MCNX-CTX-04**: family-B counters (transport executions vs. iterations; write-phase discipline shared with MCNX-SRC-01) plus the family-A MCNX-SRC-03 regression executed against this schedule realization.
- **MCNX-CTX-05**: family-B: launch with `use_subcycling = yes` and assert nonzero exit during parameter checking with the diagnostic present and zero transport-step side effects; a control launch with `no` proceeds.
- **MCNX-CTX-06**: family-B unit check on a small grid: install linear probe fields, evaluate both bridges, assert `rtol = 1e-14` agreement and bitwise determinism across repeated evaluation.
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` and `flesh/` paths, and contains the required claim strings (`REQUIRES CarpetX`, the `ghext`/`mfab` access idiom, the placement and cadence sentences, `use_subcycling`, the linear-exactness tolerance).

## Implementation freedom

- The decomposition of the transport phase into routines and sub-groups, their names, and the internal ordering — provided the group placement, cadence counters, and write-phase discipline of MCNX-CTX-04 / MCNX-SRC-01 hold.
- Level-mode vs. global-mode for each iterator-walking routine (both satisfy MCNX-CTX-03); how patches are looped when `ghext` carries more than one.
- The bridge realization — 8-point averaging, trilinear gather, or higher-order interpolation — provided MCNX-CTX-06's linear exactness and determinism hold (higher order is bounded by MCNX-GEO-06, not forbidden).
- Whether ghost validity is obtained via `SYNC:` declarations or presync-triggered driver syncs; use of `mfab->setVal(0)` vs. a kernel for zeroing.
- Parameter names/defaults beyond those pinned by other specs, diagnostic outputs, and startup-time initialization structure.
- Incidental use of the Loop layer for non-normative diagnostics.

## Open questions / assumptions

- **Subcycled transport (open question, guarded).** Stepping fine-level packets at fine Δt would change the cadence contract, the residency/redistribution timing (`particle-container-and-gpu.md`), and the source-term lag semantics. It is excluded by MCNX-CTX-05 until specified; the guard makes the exclusion observable rather than implicit.
- **Single-patch assumption (assumption, non-blocking).** The verified configurations run one patch (`ghext->num_patches() = 1`); the access idiom is written per patch and generalizes by looping `patchdata`, but multi-patch transport (inter-patch packet exchange) is unspecified and would land as a spec change here.
- **Regridding cadence (assumption).** Verified configurations use a fixed grid hierarchy. When regridding is enabled, the regrid obligation of `particle-container-and-gpu.md` (global redistribution before the next transport phase) applies; the schedule hook that triggers it (postregrid-time vs. start of the transport group) is implementation freedom until exercised, at which point the choice is recorded here.
- **Consumer read ordering (assumption, shared with `hydro-coupling-source-terms.md`).** The cadence contract guarantees stable source values at every point outside the transport phase; a GRMHD consumer reading them anywhere in its RHS therefore needs no ordering edge against MCNuX beyond the standard bin structure. If a consumer requires intra-bin ordering, it declares it on its side; MCNuX's schedule stays consumer-ignorant.
- **Checkpoint bins (assumption).** Packet-population checkpointing (`particle-container-and-gpu.md` MCNX-PAR-05) rides the driver's checkpoint machinery at its standard bins; no MCNuX-specific checkpoint schedule is pinned. If the driver's particle support requires thorn-side hooks, they are added as a spec change here.
