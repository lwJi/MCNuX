# Build and integration

> Delivery spec. Self-contained: an agent can assemble, configure, and build the MCNuX
> delivery — the thorn layout, the dependency declarations, the pinned test ThornList,
> and the test/production activation split — from this file alone. It restates only the
> conventions it uses and references
> [carpetx-thorn-integration](./carpetx-thorn-integration.md) (thorn wiring semantics)
> and [verification-suite-design](./verification-suite-design.md) (what the tests are)
> for the rest. `README.md` is canonical if any restated convention here conflicts
> with it.

## Purpose & scope

This spec pins how MCNuX is delivered and built as part of a Cactus/CarpetX
configuration:

- The **repository/thorn assembly**: the MCNuX repository as a Cactus arrangement
  containing exactly the `MCNuX` transport thorn and the `TestMCNuX` sibling test
  thorn.
- The **dependency contract**: what MCNuX's `configuration.ccl` REQUIRES, and the
  capability chain those requirements pull in.
- The **build-mode contract**: the pinned `double` precision and the CPU/GPU build
  obligations.
- The **pinned test ThornList** — the verified configuration the regression suite
  runs in (flesh + CarpetX stack + WeakLibInterp + MCNuX + TestMCNuX).
- The **golden-data location** (`MCNuX/test/`) and the rule that production
  configurations never activate `TestMCNuX`.

Out of scope:

- The benchmarks' content, golden data, tolerances, and the `TestMCNuX` thorn's
  *design* — [verification-suite-design](./verification-suite-design.md); this spec
  pins where they live and how the configuration that runs them is assembled.
- Thorn wiring semantics — the grid-access idiom, scheduling obligations, and the
  meaning of `REQUIRES CarpetX` —
  [carpetx-thorn-integration](./carpetx-thorn-integration.md); this spec pins that the
  declaration exists in the delivered ccl.
- The device execution model and the `ParticleReal` consumption
  ([particle-container-and-gpu](./particle-container-and-gpu.md)); this spec pins the
  build-time precision configuration it assumes.
- WeakLibInterp's own build internals (its provider thorn is consumed as-is, never
  modified).

## Source of truth

The dependency chain and ThornList shape are pinned from the reference stack's own
configuration files:

- `CarpetX/CarpetX/configuration.ccl` — the driver's requirements
  (`REQUIRES AMReX IOUtil MPI yaml_cpp zlib` plus `Arith CarpetXRegrid Loop`) and its
  empty `PROVIDES CarpetX { }` block (a build-order capability, no headers).
- `CarpetX/ODESolvers/configuration.ccl` — `REQUIRES AMReX CarpetX Subcycling`: the
  evolution driver the cadence contract hooks pulls in `CarpetX/Subcycling`.
- `CarpetX/Loop/configuration.ccl`, `CarpetX/TestProlongate/configuration.ccl` — the
  capability idiom precedents (`PROVIDES Loop`; a consumer's bare `REQUIRES CarpetX`).
- `WeakLibInterp/cactus/thorns/WeakLibInterp/configuration.ccl` — the WeakLibInterp
  provider thorn: an ExternalLibraries-style wrapper (`PROVIDES WeakLibInterp` with a
  detect script; `REQUIRES AMReX HDF5 MPI`), so a consumer's whole obligation is
  `REQUIRES WeakLibInterp`.
- `WeakLibInterp/cactus/wli.th` — the in-stack precedent for a minimal pinned compile
  ThornList (a compile list over an existing Cactus checkout, not a fetch list).
- `flesh/lib/sbin/RunTestUtils.pl` — test discovery is keyed to ThornList membership:
  only thorns in the configuration's ThornList contribute `test/*.par` suites, which
  is why the ThornList below is part of the verification contract.

This spec cites no published papers; its sources of truth are the reference
repositories above, so its citation registry is empty by construction.

## Inputs & outputs

### Repository layout (the delivered arrangement)

The MCNuX repository is a Cactus arrangement named `MCNuX` containing exactly two
thorns:

| Thorn | Role | Test data |
|---|---|---|
| `MCNuX/MCNuX` | the transport thorn — all MCNuX physics, parameters, and grid variables | `MCNuX/test/` holds **every** benchmark parfile and golden directory of [verification-suite-design](./verification-suite-design.md), plus the suite's `test.ccl` |
| `MCNuX/TestMCNuX` | the sibling test thorn — analytic initial data only (design pinned in [verification-suite-design](./verification-suite-design.md)) | none — `TestMCNuX` ships **no** `test/` directory; keying all tests to `MCNuX` keeps discovery in one place |

The `specs/` directory (this corpus) and its validator live at the repository root,
outside both thorns; nothing under `specs/` participates in the build.

### Configuration declarations (binding)

| File | Must declare |
|---|---|
| `MCNuX/MCNuX/configuration.ccl` | `REQUIRES CarpetX` (the hard driver dependency of [carpetx-thorn-integration](./carpetx-thorn-integration.md) [MCNX-CTX-02]) and `REQUIRES WeakLibInterp` (the opacity/EOS library of [opacity-eos-evaluation](./opacity-eos-evaluation.md)) |
| `MCNuX/TestMCNuX/configuration.ccl` | `REQUIRES CarpetX` (it writes CarpetX-driver grid functions; no WeakLibInterp dependency — analytic initial data needs no tables) |

## Correctness requirements

- **[MCNX-BLD-01] Thorn assembly.** The delivered repository is the two-thorn
  arrangement of the layout table above — no third thorn, no `specs/fixtures/`
  directory, and no build products committed. All MCNuX physics lives in
  `MCNuX/MCNuX`; `TestMCNuX` contains initial data only.

- **[MCNX-BLD-02] Dependency contract.** The configuration declarations table is
  binding: `MCNuX/MCNuX/configuration.ccl` declares `REQUIRES CarpetX` and
  `REQUIRES WeakLibInterp`; a configuration lacking either capability must fail at
  configure time, before any compilation. MCNuX declares no dependency on any
  GRMHD or spacetime-evolution thorn (GRMHD-agnosticism,
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-07]);
  its couplings are grid-variable contracts against the base thorns only.

- **[MCNX-BLD-03] Precision and build modes.** The verified build configures AMReX
  with `amrex::Real` and `amrex::ParticleReal` both IEEE-754 binary64 (`double`) —
  the corpus value type every spec assumes; a configuration with single-precision
  particle reals is out of contract and must be rejected (a compile-time
  `static_assert` on both sizes is the enforcement; the `unit-selftest` battery
  re-checks `sizeof(amrex::Real) == 8` at run time). The delivered thorn must build
  in at least two modes: single-rank-capable CPU/MPI (the golden-data mode) and a
  GPU backend (the device-kernel verification mode of
  [particle-container-and-gpu](./particle-container-and-gpu.md) [MCNX-GPU-02]).

- **[MCNX-BLD-04] The pinned test ThornList.** The regression suite of
  [verification-suite-design](./verification-suite-design.md) runs in a configuration
  whose ThornList is exactly (flesh = the Cactus checkout itself; one thorn per line,
  `arrangement/thorn`):

  ```text
  CactusBase/IOUtil
  ExternalLibraries/AMReX
  ExternalLibraries/HDF5
  ExternalLibraries/MPI
  ExternalLibraries/yaml_cpp
  ExternalLibraries/zlib
  CarpetX/ADMBaseX
  CarpetX/Arith
  CarpetX/CarpetX
  CarpetX/CarpetXRegrid
  CarpetX/HydroBaseX
  CarpetX/Loop
  CarpetX/ODESolvers
  CarpetX/Subcycling
  CarpetX/TmunuBaseX
  WeakLibInterp/WeakLibInterp
  MCNuX/MCNuX
  MCNuX/TestMCNuX
  ```

  Membership rationale (each line is load-bearing): `CactusBase/IOUtil` and the five
  ExternalLibraries providers are required by the CarpetX driver
  (`AMReX IOUtil MPI yaml_cpp zlib`) and by WeakLibInterp (`HDF5`); `Arith`,
  `CarpetXRegrid`, `Loop` are the driver's capability chain; `ODESolvers` (with its
  required `Subcycling`) is the evolution loop the cadence contract of
  [carpetx-thorn-integration](./carpetx-thorn-integration.md) hooks; `ADMBaseX`,
  `HydroBaseX`, `TmunuBaseX` are the base-thorn interfaces MCNuX reads and
  contributes to; `WeakLibInterp/WeakLibInterp` provides the opacity/EOS capability;
  `MCNuX/MCNuX` and `MCNuX/TestMCNuX` are the delivery. This list is the *verified*
  configuration: it must configure, build, and run the full harness suite clean.
  Additions for a user's production configuration are free ([MCNX-BLD-06] aside);
  removals from the test configuration are not.

- **[MCNX-BLD-05] Golden-data location.** Every benchmark parfile, every golden
  directory, and the suite `test.ccl` live under `MCNuX/test/` — the harness
  discovers tests per ThornList thorn, and keying the whole suite to `MCNuX` puts the
  golden data where the tested code lives. `TestMCNuX` ships no `test/` directory;
  test parfiles activate `TestMCNuX` via their ActiveThorns lines where a non-trivial
  background is needed.

- **[MCNX-BLD-06] Test/production activation split.** `TestMCNuX` is compiled into
  the test configuration (it is in the ThornList) but a production configuration
  never activates `TestMCNuX`: no non-test parameter file lists it in ActiveThorns,
  and no MCNuX functionality may depend on its presence at run time (MCNuX must run
  with `TestMCNuX` absent from the ThornList entirely). Production parameter files
  supply real initial data through whatever spacetime/hydro thorns the user couples;
  MCNuX's contracts with them are exactly the interface requirements of
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-04].

### Conventions restated (the subset this leaf uses)

- All reals are IEEE-754 binary64 (`double`) throughout MCNuX; the build pins
  `amrex::Real` and `amrex::ParticleReal` accordingly
  ([conventions-and-units](./conventions-and-units.md) owns the value-type
  convention).
- Golden benchmark data is diffed at the Cactus harness defaults
  `ABSTOL=RELTOL=1e-12` (`flesh/lib/sbin/RunTestUtils.pl`); golden generation is
  single-rank CPU at the primary seed
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)).
- No production GRMHD code is cited or depended on anywhere in the corpus or the
  build (GRMHD-agnosticism).

## Verification

- **Configure-time gates (exact).** The pinned ThornList configures without error.
  Removing `CarpetX/CarpetX` (or the WeakLibInterp provider) from the list makes
  configuration fail on MCNuX's `REQUIRES` declarations — proving [MCNX-BLD-02] is
  live, not decorative.
- **Build gates (exact).** The pinned ThornList builds clean in the CPU/MPI mode; the
  GPU-backend build of the same list compiles every per-packet operator as device
  kernels ([particle-container-and-gpu](./particle-container-and-gpu.md) owns the
  inspection criteria). The precision `static_assert`s of [MCNX-BLD-03] compile
  (and fail the build under a single-precision AMReX configuration).
- **Harness run (exact / golden parity 1e-12).** The flesh harness, pointed at the
  built test configuration, discovers exactly the benchmark set of
  [verification-suite-design](./verification-suite-design.md) under `MCNuX/test/` and
  passes every test at the 1e-12 defaults.
- **Activation-split audit (exact).** A sweep of all committed parameter files:
  `TestMCNuX` appears in ActiveThorns only within `MCNuX/test/*.par`; and a smoke
  configuration whose ThornList omits both `MCNuX/TestMCNuX` and the test run still
  configures, builds, and starts MCNuX (no hidden dependency on the test thorn).
- **Layout audit (exact).** The repository contains exactly the two thorns of
  [MCNX-BLD-01]; `TestMCNuX` has no `test/` directory; `specs/` contains exactly the
  13 corpus files plus `tools/validate_specs.sh` (the corpus's own validator
  enforces the spec-side count).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the ThornList member literals `CarpetX/ODESolvers`,
`WeakLibInterp/WeakLibInterp`, `MCNuX/MCNuX`, and `MCNuX/TestMCNuX`, the dependency
literal `REQUIRES WeakLibInterp`, and the activation-split phrase
`never activates \`TestMCNuX\``; that its cited `CarpetX/`, `WeakLibInterp/`, and
`flesh/` paths resolve; and that its `[MCNX-BLD-NN]` ids are declared here and covered
by the [verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Everything inside the two thorns' source trees: file layout, languages, internal
  build fragments (`make.code.defn` structure), and header organization.
- How the build is driven (Simfactory, hand-run `make <config>`, CI scripts), option
  lists, compilers, and optimization flags — the contract is the configure/build/run
  outcomes, not the tooling.
- Production ThornLists (supersets of the capability chain with the user's GRMHD and
  initial-data thorns), subject to [MCNX-BLD-06].
- Adding OPTIONAL capabilities to `configuration.ccl` (e.g. `CUDA`) as the GPU build
  mode requires — the binding REQUIRES set of [MCNX-BLD-02] is a floor, not a
  ceiling.
- Where the non-harness CI gates of
  [verification-suite-design](./verification-suite-design.md) live and how they are
  invoked.

## Open questions / assumptions

- **Assumption: the surrounding checkout supplies the non-reference arrangements.**
  `CactusBase/IOUtil` and the `ExternalLibraries/*` providers are standard Einstein
  Toolkit components assumed present in the Cactus checkout the flesh root points at;
  they are outside this corpus's five validator-checked reference roots, so the
  ThornList names them without file-level citation (the CarpetX and WeakLibInterp
  configuration.ccl files cited above are the evidence that exactly these providers
  are required). If a checkout lacks them, configuration fails loudly at the
  ThornList stage.
- **Assumption: WeakLibInterp is consumed through its provider thorn.** The
  `WeakLibInterp/WeakLibInterp` wrapper (detect script; `WEAKLIBINTERP_DIR` /
  `WEAKLIBINTERP_INSTALL_DIR` options) is the delivery mechanism; MCNuX never links
  the library by hand-written paths. How the wrapper finds or builds the library is
  that thorn's concern, not MCNuX's.
- **Assumption: single AMReX.** The configuration's one `ExternalLibraries/AMReX`
  serves the driver, the particle containers, and WeakLibInterp alike; mixing AMReX
  versions between them is out of contract. The AMReX version floor is whatever the
  pinned CarpetX and WeakLibInterp checkouts require — MCNuX adds no floor of its
  own; recording one would become necessary only if MCNuX adopted an AMReX feature
  newer than theirs (none is anticipated by this corpus).
- **Open question: GetComponents fetch list.** The pinned ThornList is a *compile*
  list over an existing checkout (the `WeakLibInterp/cactus/wli.th` precedent). A
  public `GetComponents`-style fetch list (with repository URLs and revisions) is a
  distribution concern deferred until MCNuX is distributed beyond this development
  layout; the compile list above is the binding verified configuration either way.
