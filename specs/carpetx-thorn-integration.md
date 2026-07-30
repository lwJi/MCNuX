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
  READS/WRITES declarations.
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
  `level`/`global` invoke it exactly once).
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

| Variables | Owner | Centering | MCNuX role |
|---|---|---|---|
| `alp`, `betax…betaz`, `gxx…gzz` | ADMBaseX | vertex (default) | READS — metric gather ([geodesic-propagation](./geodesic-propagation.md)) |
| `rho`, `temperature`, `Ye`, `velx…velz`, `eps`, `press` | HydroBaseX | `ccc` | READS — fluid state for opacities and tetrads (temperature in MeV, an interface requirement owned by [hydro-coupling-source-terms](./hydro-coupling-source-terms.md)) |
| MCNuX source-term variables (G^μ, lepton source) | MCNuX | `ccc` | WRITES — defined in [hydro-coupling-source-terms](./hydro-coupling-source-terms.md); advanced once per transport step per [MCNX-CTX-05] |

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

- **[MCNX-CTX-04] READS/WRITES declarations mandatory (obligation 3).** Every
  scheduled MCNuX routine declares complete `READS:`/`WRITES:` clauses for every grid
  variable it touches — regardless of the fact that [MCNX-CTX-01] reaches the data
  pointers outside the flesh's argument machinery. The driver's validity machinery
  (poisoning of written-but-undeclared-read regions, CRC32 checksums catching
  undeclared writes, presync of invalid ghosts before declared reads) enforces these
  declarations at runtime; an undeclared access is a hard error, not a silent hazard.

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
[verification-suite-design](./verification-suite-design.md) matrix.

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
  step; per-patch vs fused iteration order over patches and levels.
- MCNuX's own parameter names and defaults, diagnostics, and output cadence — except
  where sibling specs pin them.

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
- **Assumption: ADMBaseX default centering.** The metric groups declare no
  `CENTERING`, and the driver default is vertex-centered — the same recorded
  assumption as [geodesic-propagation](./geodesic-propagation.md); a change would
  alter gather bookkeeping, not this spec's obligations.
