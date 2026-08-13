# CarpetX thorn integration

> Technical leaf spec. Self-contained: an agent can implement MCNuX's thorn wiring —
> configuration, grid access, scheduling, and the transport cadence against the CarpetX
> driver and base thorns — from this file alone. It restates only the conventions it
> uses and references [conventions-and-units](./conventions-and-units.md) for the rest.
> `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope

This spec pins how MCNuX exists as a Cactus thorn on the CarpetX driver:

- The **grid-access idiom**: MCNuX reaches grid data through the driver's native AMReX
  objects (`ghext` → per-patch `AmrCore` → per-(patch, level, group) `MultiFab`s and
  `Array4` views, centering from `indextype`), and the **three obligations** that
  travel with that idiom, stated as correctness requirements: `REQUIRES CarpetX`,
  level/global-mode scheduling for container-walking routines, and mandatory
  READS/WRITES declarations — including the **pinned declaration convention** for the
  foreign-thorn groups (`ADMBaseX::metric`, `HydroBaseX::rho`, …) MCNuX reaches through
  `ghext` rather than through `CCTK_ARGUMENTS`.
- The **transport cadence contract**, stated as observable behavior: transport is
  operator-split, once per fluid step, outside the RK substeps, consuming end-of-step
  fluid state and metric — with the **single-rate pin**
  (`CarpetX::use_subcycling = no`) and subcycling recorded as an open question.
- The **anti-pattern warning** against the driver interpolator's execution model
  (required verbatim in both technical specs).

Out of scope:

- The packet container layout, device-kernel execution idiom, and redistribution
  invariants — [particle-container-and-gpu](./particle-container-and-gpu.md); this
  spec supplies the objects and the schedule slots those kernels run in.
- The physics executed inside the transport step
  ([geodesic-propagation](./geodesic-propagation.md),
  [neutrino-matter-interactions](./neutrino-matter-interactions.md),
  [packet-representation-and-sampling](./packet-representation-and-sampling.md),
  [trapped-regime-treatment](./trapped-regime-treatment.md),
  [opacity-eos-evaluation](./opacity-eos-evaluation.md)).
- The source-term grid variables, the zero-then-add accumulation protocol, and every
  interface requirement on the GRMHD partner (including the MeV temperature
  requirement) — [hydro-coupling-source-terms](./hydro-coupling-source-terms.md); this
  spec pins only *when* in the step those variables advance.
- The `TestMCNuX` test thorn, benchmark parameter files, and golden data
  ([verification-suite-design](./verification-suite-design.md)); the ThornList and
  build wiring ([build-and-integration](./build-and-integration.md)).

## Source of truth

All contracts below are stated against the CarpetX driver and base thorns — no
production GRMHD code is cited anywhere in this corpus. Reference artifacts:

- `CarpetX/CarpetX/src/driver.hxx` — the driver data model: the singleton
  `extern std::unique_ptr<GHExt> ghext` (~line 596); per patch
  `std::unique_ptr<CactusAmrCore> amrcore` (~line 365, `final : public amrex::AmrCore`);
  per (patch, level, group) `std::vector<std::unique_ptr<amrex::MultiFab>> mfab`
  (~line 477, indexed by time level, `numvars` components each); group centering in
  `GroupData::indextype`.
- `CarpetX/CarpetX/src/interpolate.cxx` — provenance for *where the objects live*
  (container constructed on `patchdata.amrcore`, ~line 570; grid data as
  `groupdata.mfab.at(tl)->array(pti)` `Array4` views, ~lines 731–738; centering
  decoded from `groupdata.indextype`) — **object access only**;
  `CarpetX/CarpetX/src/interp.hxx` — the same service's AoS container (~line 38), the
  execution side this corpus rejects ([MCNX-CTX-07]).
- `CarpetX/CarpetX/src/schedule.cxx` — the declare-then-enforce validity machinery
  (`CallFunction`, ~lines 2178–2565): READS checked and presync'd, WRITES poisoned
  before the call (`poison_undefined_values`), CRC32 checksums over
  valid-but-undeclared data, validity updated after the call; and the schedule
  `OPTIONS` decode (`local` invokes a routine once per (patch, level, component) tile,
  `level`/`global` invoke it exactly once). Four further mechanics of the same file are
  binding on the declaration convention below: `check_storage_clauses` (~lines
  2115–2175) aborts via `CCTK_VERROR` when *any* READS/WRITES/INVALIDATES/SYNC clause
  names a group without active STORAGE; `decode_clauses` (~lines 1025–1060) resolves
  each clause entry to a *variable index* (`RDWR_entry::varindex`) and decodes the
  region tag into interior/boundary/ghosts validity bits; the poison pass (~lines
  2309–2352) poisons a written region **except** where the same routine also declares it
  read (`provided & ~need`); and `InvalidateTimelevels` is deliberately **not** called in
  the evolution loop (~lines 1948–1951, commented out "we want to use them in the next
  iteration"), so non-evolved groups keep their validity across iterations.
- The stack's live pattern for declaring clauses on another implementation's groups —
  always `inherits:` plus a fully qualified `Implementation::group` clause:
  `CarpetX/TestInterpolate/interface.ccl` (`INHERITS: CoordinatesX`, ~line 5) with
  `CarpetX/TestInterpolate/schedule.ccl` (`OPTIONS: global` routines declaring
  `READS: CoordinatesX::vertex_coords(everywhere)` / `cell_coords(everywhere)`,
  ~lines 8–24) whose bodies reach that data by name through a driver service
  (`CCTK_VarIndex("CoordinatesX::vcoordx")`, `Interpolate(...)`,
  `CarpetX/TestInterpolate/src/test_vertex_interpolation.cxx` ~lines 17–19, 89) rather
  than through the group's `CCTK_ARGUMENTS` pointers — the exact
  declare-in-CCL/access-outside-the-argument-machinery precedent MCNuX needs;
  `CarpetX/TestNorms/interface.ccl` (`inherits: Driver, CarpetXRegrid`, line 3) with
  `CarpetX/TestNorms/schedule.ccl` (`WRITES: CarpetXRegrid::regrid_error(interior)`,
  line 15); `CarpetX/PoissonX/interface.ccl` (~line 5) with
  `CarpetX/PoissonX/schedule.ccl` (`READS: PDESolvers::point_type(everywhere)`, ~line
  55); `CarpetX/BoxInBox/interface.ccl` with `CarpetX/BoxInBox/schedule.ccl` (~line 24);
  `CarpetX/MovingBoxToy/schedule.ccl` (~lines 6–7). Across the whole CarpetX tree every
  qualified `Implementation::group` clause occurs in a thorn that inherits that
  implementation — there is no counterexample.
- `CarpetX/ADMBaseX/schedule.ccl` — unconditional `STORAGE:` for `metric`/`lapse`/
  `shift` (lines 3–8) and initial-data/gauge writers declaring `WRITES: …(everywhere)`
  (lines 68–112); `CarpetX/HydroBaseX/schedule.ccl` — unconditional `STORAGE:` for the
  fluid groups (lines 3–13) and the default `initial_hydro = "vacuum"` writer declaring
  `WRITES: …(everywhere)` on all of them (lines 59–75). Both base thorns therefore
  supply storage and everywhere-validity whenever they are active.
- `CarpetX/CarpetX/configuration.ccl` — `PROVIDES CarpetX { }` is an **empty** block: a
  build-order capability, no headers or libraries;
  `CarpetX/TestProlongate/configuration.ccl` and
  `CarpetX/ODESolvers/configuration.ccl` — consumer thorns declaring
  `REQUIRES CarpetX`.
- `CarpetX/ODESolvers/src/solve.hxx` — the stack's live idiom for reaching driver
  internals: a relative include of the driver's own `driver.hxx` (two directories up,
  into the CarpetX thorn's `src/`) behind an in-source
  `// TODO: Don't include files from other thorns; create a proper interface`.
- `CarpetX/ODESolvers/schedule.ccl` — `ODESolvers_Solve AT evol` with
  `OPTIONS: level`; the RK-substep groups `ODESolvers_RHS` / `ODESolvers_PostStep` are
  invoked per stage from *inside* the solver (via the `CallScheduleGroup` alias,
  `CarpetX/ODESolvers/interface.ccl`), with the `ODESolvers_CalledBySolve AT EVOL
  AFTER ODESolvers_Solve` group present for schedule display only.
- `CarpetX/CarpetX/param.ccl` — `BOOLEAN use_subcycling` (default `no`);
  `CarpetX/ODESolvers/param.ccl` — the live parameter-sharing idiom
  (`SHARES: Driver` + `USES BOOLEAN use_subcycling`).
- `CarpetX/ADMBaseX/interface.ccl` — the metric/lapse/shift groups MCNuX reads (no
  `CENTERING` declaration ⇒ driver-default vertex-centered);
  `CarpetX/HydroBaseX/interface.ccl` — the cell-centered (`CENTERING={ccc}`) fluid
  groups `rho`, `vel`, `eps`, `press`, `temperature`, `Ye` MCNuX reads.

This spec cites no published papers; its sources of truth are the reference
repositories above, so its citation registry is empty by construction.

## Inputs & outputs

### What MCNuX reads and writes on the grid

| Variables | Owner | Centering | Declared as ([MCNX-CTX-04]) | MCNuX role |
|---|---|---|---|---|
| `alp`, `betax…betaz`, `gxx…gzz` | ADMBaseX | vertex (default) | `READS: ADMBaseX::lapse(everywhere)`, `ADMBaseX::shift(everywhere)`, `ADMBaseX::metric(everywhere)` | READS — metric gather ([geodesic-propagation](./geodesic-propagation.md)) |
| `rho`, `temperature`, `Ye`, `velx…velz`, `eps`, `press` | HydroBaseX | `ccc` | `READS: HydroBaseX::rho(everywhere)`, `…::vel(everywhere)`, `…::eps(everywhere)`, `…::press(everywhere)`, `…::temperature(everywhere)`, `…::Ye(everywhere)` | READS — fluid state for opacities and tetrads (temperature in MeV, an interface requirement owned by [hydro-coupling-source-terms](./hydro-coupling-source-terms.md)) |
| MCNuX source-term variables (G^μ, lepton source) | MCNuX | `ccc` | `WRITES:` (and, for accumulating routines, also `READS:`) on MCNuX's own groups | WRITES — defined in [hydro-coupling-source-terms](./hydro-coupling-source-terms.md); advanced once per transport step per [MCNX-CTX-05] |

`ADMBaseX::curv`/`dtlapse`/`dtshift` and `HydroBaseX::entropy`/`Bvec`/`Avecx…Avecz` are
**not** touched by MCNuX and are therefore not declared (see [MCNX-CTX-04](c)).

### Declared CCL surface

`MCNuX/MCNuX/interface.ccl` declares `IMPLEMENTS: MCNuX` plus
`INHERITS: ADMBaseX HydroBaseX` — the visibility that makes the qualified clauses above
legal. Inheritance is a declaration-scope statement only: MCNuX still reaches the data
through `ghext` per [MCNX-CTX-01], exactly as `TestInterpolate` inherits `CoordinatesX`
yet reads that data by variable index through a driver service (Source of truth).
Consequence: every MCNuX parameter file activates a thorn implementing ADMBaseX and one
implementing HydroBaseX (the pinned ThornList of
[build-and-integration](./build-and-integration.md) contains both). A build that also
compiles the optional TmunuBaseX contributor
([hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-06]) adds
`INHERITS: TmunuBaseX` and the corresponding activation obligation — see
[MCNX-CTX-04](e).

### Driver objects consumed (the grid-access surface)

| Object | Reached as | Used for |
|---|---|---|
| per-patch AMR core | `ghext->patchdata[patch].amrcore` | `Geom(lev)`, `boxArray(lev)`, `DistributionMap(lev)` — packet-container construction and deposition targets |
| grid-function data | `ghext->patchdata[patch].leveldata[level].groupdata[gi]->mfab[tl]` | raw `amrex::MultiFab` per (group, time level) |
| field views | `groupdata.mfab.at(tl)->array(mfi)` | `Array4` views captured by value into device kernels |
| centering | `groupdata.indextype` | vertex/cell stencil bookkeeping in gather and deposition |

Parameters: MCNuX's own parameters are free (names, defaults) except where sibling
specs pin them; MCNuX additionally *reads* `CarpetX::use_subcycling` (via the
`SHARES:`/`USES` idiom or `CCTK_ParameterGet`) to enforce [MCNX-CTX-06].

## Correctness requirements

- **[MCNX-CTX-01] Driver-native grid access.** MCNuX reaches grid data exclusively
  through the driver's native AMReX objects, per the access table above: the per-patch
  `AmrCore` for grid structure, group `MultiFab`s for data, `Array4` views for kernel
  access, `indextype` for centering. The Loop thorn's `GF3D2`/point-loop layer is
  **not used** for particle operators — MC transport is particle-centric (gather at
  packet positions, push, deposit), not grid-point-centric. Provenance for where these
  objects live is the driver's own interpolation service (Source of truth) — for
  object access only, per [MCNX-CTX-07].

- **[MCNX-CTX-02] Hard driver dependency (obligation 1).** MCNuX's
  `configuration.ccl` declares `REQUIRES CarpetX`. The Loop capability alone does not
  expose `ghext` or any driver internals, and CarpetX's `PROVIDES CarpetX { }` block is
  empty — the requirement expresses the build-order dependency that makes the
  driver-internal include path of [MCNX-CTX-01] well-defined.

- **[MCNX-CTX-03] Level/global-mode scheduling for container walkers (obligation 2).**
  Every scheduled MCNuX routine that itself iterates `ParIter`/`MFIter` (or otherwise
  walks tiles/boxes) is scheduled with `OPTIONS: level` (or `global`), never
  local-mode: the driver invokes a local-mode routine once per (patch, level,
  component) tile, so a self-walking routine scheduled local would have its walk nested
  inside the driver's own tile loop, visiting tiles multiply. (`ODESolvers_Solve`'s
  `OPTIONS: level` is the stack precedent.)

- **[MCNX-CTX-04] READS/WRITES declarations mandatory, and how foreign groups are named
  (obligation 3).** Every scheduled MCNuX routine declares complete `READS:`/`WRITES:`
  clauses for every grid variable it touches — regardless of the fact that
  [MCNX-CTX-01] reaches the data pointers outside the flesh's argument machinery. The
  driver's validity machinery (poisoning of written-but-undeclared-read regions, CRC32
  checksums catching undeclared writes, presync of invalid ghosts before declared reads)
  enforces these declarations at runtime; an undeclared access is a hard error, not a
  silent hazard. Because the groups MCNuX touches are owned by other thorns, the
  declaration mechanism is **pinned here**, not left to the implementer:

  - **(a) Inheritance-qualified naming (the convention).** `MCNuX/MCNuX/interface.ccl`
    declares `INHERITS: ADMBaseX HydroBaseX`, and every clause on a group MCNuX does not
    own uses the fully qualified `Implementation::group(region)` form —
    `READS: ADMBaseX::metric(everywhere)`, `READS: HydroBaseX::rho(everywhere)`, and so
    on per the I/O table. This is the stack's only live pattern for cross-implementation
    clauses and it is universal there (Source of truth: TestInterpolate/CoordinatesX,
    TestNorms and BoxInBox and TestProlongate/CarpetXRegrid, PoissonX/PDESolvers,
    MovingBoxToy/BoxInBox — every qualified clause sits in a thorn that inherits that
    implementation, with no counterexample). `INHERITS:` is a *declaration-visibility*
    statement: it neither implies nor requires that MCNuX obtain the data through
    `CCTK_ARGUMENTS`, and it does not weaken [MCNX-CTX-01] or [MCNX-CTX-07].

  - **(b) One name, two uses.** The qualified string in a clause is character-for-
    character the string the routine passes to `CCTK_GroupIndex`/`CCTK_VarIndex` when it
    resolves the `ghext` lookup (`…leveldata[level].groupdata[gi]`) for that group.
    Declaration and access therefore cannot drift, and the correspondence is
    mechanically auditable (Verification). A group index obtained by any other route —
    a hard-coded index, or a name assembled at run time from a parameter — is
    non-compliant even if the numbers happen to agree.

  - **(c) Regions are exactly the regions touched.** Fields gathered at packet positions
    are declared `(everywhere)`: a gather stencil around a packet near a tile edge reads
    that tile's ghost points, so interior-only validity is insufficient and the driver
    must presync before the routine runs. Symmetrically, a group MCNuX does not touch is
    not declared: an over-declared read makes the driver demand validity nothing in the
    configuration provides and converts a correct run into a spurious validity abort.

  - **(d) Accumulate-in-place requires both clauses.** With `poison_undefined_values`
    the driver poisons every written region that the same routine does not also declare
    as read. Any routine that accumulates (`+=`) into MCNuX's source-term variables
    therefore declares them in **both** `READS:` and `WRITES:` with the same region; only
    the zeroing routine of [hydro-coupling-source-terms](./hydro-coupling-source-terms.md)
    declares `WRITES:` alone.

  - **(e) Foreign groups are read-only, with one sanctioned exception.** The base-thorn
    groups of the I/O table appear in `READS:` clauses and never in
    `WRITES:`/`INVALIDATES:`. The single exception is the parameter-selected,
    default-off TmunuBaseX contributor of
    [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-06], which
    accumulates into TmunuBaseX's own variables through TmunuBaseX's own extension
    group; a build that compiles it declares `INHERITS: TmunuBaseX` as well and names
    those groups by the same qualified form, under the same rules (b)–(d) — in
    particular the accumulate-in-place pair of (d). Whether that contributor exists at
    all, and its default, stay owned by [MCNX-HYD-06]; only the naming and inheritance
    mechanics are pinned here. Because `INHERITS:` cannot be parameter-conditional, a
    build carrying the contributor requires TmunuBaseX active in every run of that
    build.

  - **(f) Storage and activation.** Every group named in any clause must have active
    STORAGE at the moment the routine is called, or the driver aborts with its
    storage-check error naming the group. ADMBaseX and HydroBaseX allocate their groups
    unconditionally while active, so the obligation reduces to: both are active in every
    MCNuX parameter file, and MCNuX declares no clause on a group whose owning thorn it
    does not require.

  - **(g) Fallback ladder, in binding order.** If the flesh in use rejects the primary
    form, the implementation descends this ladder and records which rung it used:
    (1) enumerate the variables instead of the group (`READS: ADMBaseX::gxx(everywhere)`
    … `gzz(everywhere)`), which produces the identical per-variable clause entries;
    (2) if a qualified foreign clause cannot be expressed at all, the affected routine is
    **not scheduled** and the build or run fails loudly with a message naming the group.
    Deleting a clause while keeping the access is never admissible: under every rung,
    a routine touches no group it has not declared, and re-enabling an access that no
    rung supports requires amending this spec.

- **[MCNX-CTX-05] Transport cadence contract (observable behavior).** Transport is
  operator-split from the fluid/spacetime evolution, running **outside** the RK
  substeps, after the fluid step (including its `ODESolvers_PostStep` con2prim/boundary
  pass) completes. The observable contract: packet positions/momenta and the MCNuX
  source-term grid variables advance exactly once per coarsest-level Δt; each transport
  step propagates packets over the full step Δt using the end-of-step fluid state and
  metric of that step; and the source terms read by the fluid at step n are those
  tallied during step n−1's transport. The canonical realization is an MCNuX-owned
  schedule group `AT evol AFTER ODESolvers_Solve`; any schedule text producing exactly
  this observable behavior is admissible (`IN ODESolvers_PostStep` and
  `AT CCTK_POSTSTEP` do not, and were considered and rejected — see Open questions).

- **[MCNX-CTX-06] Single-rate pin.** The verified configuration is
  `CarpetX::use_subcycling = no`: all refinement levels advance with the same Δt, so
  "once per coarsest-level Δt" is once per `CCTK_EVOL` traversal for every packet on
  every level. MCNuX checks the parameter at `CCTK_PARAMCHECK` and **aborts** with an
  explanatory error when subcycling is enabled — transport behavior under subcycling
  is unspecified (open question below), and running unspecified is not permitted.

- **[MCNX-CTX-07] Anti-pattern warning (binding, stated identically in both technical
  specs).** The CarpetX driver interpolator is cited for **object access only**; its
  legacy AoS particle container and host-side per-particle loops
  (`CarpetX/CarpetX/src/interp.hxx`, `CarpetX/CarpetX/src/interpolate.cxx`) are an
  anti-pattern for GPU particle operators and are not adopted. From that service this
  corpus takes exactly one thing — where the driver's objects live ([MCNX-CTX-01]) —
  and none of its execution model; the binding execution idiom is
  [particle-container-and-gpu](./particle-container-and-gpu.md)'s.

### Conventions restated (the subset this leaf uses)

- Geometrized units `G = c = M_sun = 1` for all grid coordinates, Δt, metric
  quantities, and source terms; the transport step is the coarsest-level Δt.
- Metric fields read are lapse `α` (`alp`), shift `β^i` (`betax…betaz`), spatial
  metric `γ_ij` (`gxx…gzz`), vertex-centered; fluid fields are cell-centered;
  `HydroBaseX::temperature` is consumed as MeV (requirement on the partner owned by
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md)).
- Golden benchmark outputs run at the Cactus harness defaults `ABSTOL=RELTOL=1e-12`
  (`flesh/lib/sbin/RunTestUtils.pl`).

## Verification

- **Configuration dependency (build-time, exact).** A configuration without CarpetX
  fails at configure time on MCNuX's `REQUIRES CarpetX`; the delivered ThornList
  ([build-and-integration](./build-and-integration.md)) builds clean.
- **Paramcheck guard (exact).** A test parameter file with
  `CarpetX::use_subcycling = yes` aborts at `CCTK_PARAMCHECK` with MCNuX's error; the
  identical file with `no` runs.
- **Cadence observables (golden parity, 1e-12 harness).** A deterministic fixed-seed
  benchmark outputs per-iteration diagnostics (a transport-step counter and a
  source-term tally checksum as grid scalars/norms): the counter increments exactly
  once per `CCTK_EVOL` iteration (never per RK stage — with RK4 active, four RHS
  evaluations per step must show one transport advance); and a step-function change
  injected into the fluid state at iteration n first affects the source terms the
  fluid reads at iteration n+1 (the one-step lag of [MCNX-CTX-05]), all archived as
  golden TSV at `ABSTOL=RELTOL=1e-12`.
- **Validity-machinery run (exact).** The benchmark suite runs with
  `poison_undefined_values` enabled (and a presync-mode test run): completing without
  driver validity errors demonstrates the READS/WRITES declarations of [MCNX-CTX-04]
  are complete; a deliberately under-declared routine in a meta-test triggers the
  driver's checksum/poison error, proving the enforcement is live.
- **Declaration/access name parity (exact, source sweep, non-harness gate).** Every
  fully qualified group/variable name string passed to `CCTK_GroupIndex`/`CCTK_VarIndex`
  on an MCNuX `ghext` lookup path occurs verbatim in a `READS:`/`WRITES:` clause of the
  routine performing that lookup, and every foreign group declared in a clause is
  actually looked up by the declaring routine — [MCNX-CTX-04](b)(c). The sweep exits
  non-zero on either direction of mismatch, so both under-declaration and
  over-declaration fail the build gate.
- **Accumulate-in-place declaration (exact, meta-test).** A contributor routine that
  declares only `WRITES:` on a source-term group it accumulates into must fail under
  `poison_undefined_values` (the driver poisons the region before the call, so the `+=`
  consumes poison); the compliant `READS:` + `WRITES:` form of [MCNX-CTX-04](d) completes
  and reproduces the golden tally at `ABSTOL=RELTOL=1e-12`.
- **Storage/activation gate (exact).** A parameter file that omits ADMBaseX or
  HydroBaseX from `ActiveThorns` fails — at inheritance resolution or with the driver's
  storage-check abort naming the group — while the pinned benchmark parameter files,
  which activate both, run; this proves [MCNX-CTX-04](f) is enforced rather than assumed.
- **Scheduling-mode check (exact).** On a pinned multi-box run (≥ 2 boxes per rank), a
  per-invocation counter in each container-walking routine shows exactly one invocation
  per (patch, level) traversal — not one per tile — confirming [MCNX-CTX-03].

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the literals `REQUIRES CarpetX` and `CarpetX::use_subcycling = no`,
the cadence phrase `advance exactly once per coarsest-level Δt`, and the anti-pattern
phrase `are an anti-pattern for GPU particle operators and are not adopted`; that its
cited `CarpetX/` and `flesh/` paths resolve; and that its `[MCNX-CTX-NN]` ids are
declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix. The declaration
convention of [MCNX-CTX-04] adds **no** validator needle: it is a claim about CCL text
and driver behavior, verified by the source sweep and the runtime gates above, not by
the spec-text linter.

## Implementation freedom

- Schedule group and routine names, how many routines the transport step is split
  into, and their ordering within the group — subject to the observable contract of
  [MCNX-CTX-05].
- How driver internals are included (the stack's relative-include idiom, a vendored
  declaration, or a future proper interface) — the accessed objects and obligations
  are the contract, not the include mechanism.
- How `CarpetX::use_subcycling` is read (`SHARES:`/`USES` vs `CCTK_ParameterGet`) and
  the wording of the paramcheck error.
- Caching of object lookups (`amrcore`, group indices, `MultiFab` pointers) across the
  step; per-patch vs fused iteration order over patches and levels — provided the cached
  group index was obtained from the qualified name of [MCNX-CTX-04](b).
- How the clauses are distributed over routines (one gather routine declaring all reads
  vs several), and the names of MCNuX's own groups (owned by
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md)).
- MCNuX's own parameter names and defaults, diagnostics, and output cadence — except
  where sibling specs pin them.

Explicitly **not** free: how `READS:`/`WRITES:` name foreign-thorn groups reached
through `ghext`. That was previously an unstated gap and is now pinned by
[MCNX-CTX-04](a)–(g) — `INHERITS: ADMBaseX HydroBaseX` plus fully qualified
`Implementation::group(region)` clauses, with the fallback ladder as the only permitted
deviation.

## Open questions / assumptions

- **Open question: subcycled transport.** With `CarpetX::use_subcycling = yes`,
  refined levels advance in per-level batches with different Δt, and "once per
  coarsest-level Δt" no longer determines when packets on refined levels advance, when
  redistribution runs ([particle-container-and-gpu](./particle-container-and-gpu.md)),
  or which fluid state a mid-batch transport step should consume. This is deliberately
  unspecified; [MCNX-CTX-06]'s paramcheck guard keeps the unspecified region
  unreachable until a subcycling extension is specified here.
- **Assumption: rejected schedule placements (recorded).** `IN ODESolvers_PostStep`
  runs once per RK stage (the solver invokes the group per stage), violating the
  once-per-step observable; `AT CCTK_POSTSTEP` runs after the evol traversal in a bin
  whose interleaving with regridding and analysis makes the "source terms ready for
  step n+1's RHS reads" ordering fragile. Both were considered and rejected in design
  review; the contract binds the observable, and the `AT evol AFTER ODESolvers_Solve`
  realization is the one verified.
- **Assumption: driver-internal access remains include-based.** CarpetX exposes no
  public capability header for `ghext` (empty `PROVIDES CarpetX { }`); the stack's own
  thorns use the TODO-flagged relative include (`CarpetX/ODESolvers/src/solve.hxx`).
  MCNuX assumes this idiom remains available; a driver-side interface change would be
  an implementation migration, not a contract change.
- **Assumption: the foreign-clause precedent is real but indirect (recorded).** An
  earlier audit reported "no in-stack precedent" for [MCNX-CTX-04]'s convention; that
  reading came from grepping only for `READS: ADMBaseX::` / `READS: HydroBaseX::`. The
  pattern does exist, on other implementations and always paired with `inherits:`
  (Source of truth), including one case — `TestInterpolate` — that declares a foreign
  `(everywhere)` read on a non-local-mode routine and then reaches the data by name
  through a driver service rather than through `CCTK_ARGUMENTS`. Residual risk: no
  in-stack example names an ADMBaseX or HydroBaseX group, and the flesh's CST (which
  decides whether a qualified clause is in scope) is not part of this checkout, so the
  scope rule is inferred from the CarpetX consumer thorns rather than read from the
  parser. [MCNX-CTX-04](g)'s ladder exists exactly for that residual risk, and the first
  implementation to compile these clauses records which rung it used.
- **Assumption: base-thorn validity survives into `CCTK_EVOL`.** MCNuX's transport group
  runs `AT evol AFTER ODESolvers_Solve` and reads groups nothing evolves. This is
  satisfiable because CarpetX does not invalidate non-evolved groups inside the
  evolution loop (the `InvalidateTimelevels` call is commented out with the rationale
  that `ODESolvers_PostStep` results must survive into the next iteration) and both base
  thorns write their groups `(everywhere)` at initial data — so analytic-background
  benchmarks with no live GRMHD partner satisfy MCNuX's `(everywhere)` reads. If a
  driver revision re-enabled that invalidation, MCNuX's reads would have to be re-based
  on a partner that rewrites those groups each step: an assumption change, not a change
  to [MCNX-CTX-04].
- **Assumption: ADMBaseX default centering.** The metric groups declare no
  `CENTERING`, and the driver default is vertex-centered — the same recorded
  assumption as [geodesic-propagation](./geodesic-propagation.md); a change would
  alter gather bookkeeping, not this spec's obligations.
