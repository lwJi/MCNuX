# Build and integration (delivery: assembling, building, and invoking the suite)

> Delivery spec. Self-contained: an agent can assemble a buildable Cactus/CarpetX configuration containing MCNuX, satisfy every capability and link-time obligation, and invoke the verification suite from this file alone, referencing `carpetx-thorn-integration.md` for the thorn's ccl shape and schedule surface (this spec owns only the capability line and the build/delivery consequences), `opacity-eos-evaluation.md` for the WeakLibInterp API contract whose link obligation is realized here, and `verification-suite-design.md` for what the invoked suite must contain. The tolerance tiers named here (golden `1e-12`, exact) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines the delivery contract: how MCNuX declares its capabilities (`REQUIRES CarpetX WeakLibInterp`), how a configuration obtains its dependencies (exactly one AMReX per executable, reached through the Cactus-provided install prefix; WeakLibInterp through its ExternalLibraries-style provider thorn hosted in the WeakLibInterp repository), how GPU code is compiled (inline, by the thorn's own configured compiler — no `.cu` files, the CarpetX convention), how the verification suite is invoked (family-A regression tests discovered and diffed by the Cactus RunTest harness; family-B unit checks as standalone host-runnable targets; statistical checks riding family-A runs as analysis-routine norm output), and the minimal-consumer smoke expectation — modeled on the `TestWeakLibInterp` scratch consumer — that every built executable must satisfy.

In scope:

- The capability declarations MCNuX owes its `configuration.ccl` and the transitive chain that makes them sufficient.
- The reference thornlist shape: the minimum thorn set that configures and builds MCNuX, and the option-list variables MCNuX reads but never redefines.
- The one-AMReX-per-executable rule and its symbol-level observability.
- The source/compilation convention for GPU code.
- The suite-invocation contract binding `verification-suite-design.md`'s two check families to concrete runners.
- The build/link/startup smoke expectation.

Out of scope:

- The thorn's schedule placement, grid-access idiom, and non-capability ccl content — owned by `carpetx-thorn-integration.md`.
- WeakLibInterp argument conventions, units, range policy, and table residency — owned by `opacity-eos-evaluation.md`; this spec only delivers the headers and library those contracts consume.
- The suite's rows, benchmarks, and tolerances — owned by `verification-suite-design.md`; this spec only pins how the suite is run.
- Modifying WeakLibInterp, CarpetX, AMReX, or the flesh; official Einstein Toolkit distribution packaging; cluster option lists (simfactory machine files); performance/optimization build variants.

## Source of truth

- `WeakLibInterp/cactus/thorns/WeakLibInterp/configuration.ccl` — the ExternalLibraries-style provider thorn this configuration includes: `REQUIRES AMReX HDF5 MPI` plus a `PROVIDES WeakLibInterp` block with a detect script and the option variables `WEAKLIBINTERP_DIR` / `WEAKLIBINTERP_INSTALL_DIR`. A consumer's whole obligation is `REQUIRES WeakLibInterp`.
- `WeakLibInterp/cactus/thorns/WeakLibInterp/src/detect.sh` and `WeakLibInterp/cactus/thorns/WeakLibInterp/src/build.sh` — find-or-build keyed on `WEAKLIBINTERP_DIR` (unset or `BUILD` → build the library from its repository checkout; a path → verified pre-installed prefix, loud failure if unusable); the build drives the library's own CMake against the ET AMReX thorn's install (`AMREX_DIR`), translating `AMREX_ENABLE_CUDA`/`AMREX_ENABLE_HIP` into the library's GPU backend and never forwarding an MPI guess (MPI is derived from what the AMReX prefix records).
- `WeakLibInterp/cactus/thorns/TestWeakLibInterp/configuration.ccl`, `WeakLibInterp/cactus/thorns/TestWeakLibInterp/src/test_wli.cxx`, and `WeakLibInterp/cactus/thorns/TestWeakLibInterp/schedule.ccl` — the minimal-consumer proof this spec's smoke expectation is modeled on: the consumer declares only `REQUIRES WeakLibInterp`, odr-uses one `_Point` entry point per family behind a never-true guard (compiling the header-inline device surface with the consumer's own compiler, no tables needed), and asserts `wli_value_type_size() == 8` at startup — the one out-of-line symbol (`WeakLibInterp/src/core/wli.cpp`) that forces real archive members into the link.
- `WeakLibInterp/specs/cactus-integration.md` — the library-side delivery contract MCNuX consumes: exactly one AMReX per executable, the configure-time configuration-consistency guard, consumers compile the device headers themselves, HDF5/MPI from the ET provider thorns. Nothing in it is normative for MCNuX beyond what is restated here.
- `CarpetX/CarpetX/configuration.ccl` — the driver capability: `PROVIDES CarpetX`, itself `REQUIRES AMReX IOUtil MPI yaml_cpp zlib` and `REQUIRES Arith CarpetXRegrid Loop`, so requiring CarpetX transitively pins the driver's whole dependency chain including the one AMReX install.
- `CarpetX/WaveToyX/src/wavetoyx.cxx` and its sibling `make.code.defn` (`SRCS = wavetoyx.cxx`) — the CarpetX source convention: thorn sources are plain `.cxx` files; no thorn in the CarpetX repository lists `.cu` files; device code is inline `amrex::ParallelFor` lambdas compiled once with whatever compiler the AMReX capability configures.
- `flesh/lib/sbin/RunTest.pl` and `flesh/lib/sbin/RunTestUtils.pl` — the Cactus test harness that invokes family-A tests: discovery requires `test/<name>.par` plus a same-named sibling golden directory; comparison is a numeric line-by-line diff at default `ABSTOL = RELTOL = 1e-12`, overridable per-thorn/per-test via `test.ccl`.
- `WeakLibInterp/cactus/wli.th` — precedent for a minimal compile thornlist carrying exactly one provider chain plus one consumer; the reference thornlist below follows the same discipline.

## Inputs & outputs

This spec defines no callable surface. Its inputs are a Cactus checkout plus a thornlist and option list; its outputs are a built executable and a runnable verification suite.

### Reference configuration (thornlist + option list)

The **reference thornlist** is the minimum set that configures and builds MCNuX (arrangement paths are checkout-layout freedom; the thorn set is binding):

| Ingredient | Thorns | Why |
|---|---|---|
| External libraries | ET `AMReX`, `HDF5` (+ its `zlib`), `MPI` providers | the one AMReX install everything links; HDF5 for WeakLibInterp's table reader; MPI for driver and library |
| Driver chain | `CarpetX` + its required capabilities (`Loop`, `Arith`, `CarpetXRegrid`, `IOUtil`, `yaml_cpp`, `zlib` per `CarpetX/CarpetX/configuration.ccl`) + `ODESolvers` | the driver MCNuX schedules against |
| Base variable thorns | `ADMBaseX`, `HydroBaseX`, `TmunuBaseX` | the metric/fluid/stress-energy groups MCNuX reads and (optionally) accumulates into |
| WeakLibInterp delivery | the `WeakLibInterp` provider thorn (hosted in the WeakLibInterp repository) | delivers headers + `libwli_lib.a` behind `REQUIRES WeakLibInterp` |
| MCNuX | the MCNuX thorn | this corpus's subject |

Option-list variables MCNuX **reads but never redefines** (each is owned by its provider thorn): `AMREX_DIR`, `AMREX_ENABLE_CUDA`, `AMREX_ENABLE_HIP` and the architecture variables (ET AMReX thorn — they select the single AMReX install and, under CUDA, switch consumer thorns to the CUDA compiler); `WEAKLIBINTERP_DIR`, `WEAKLIBINTERP_INSTALL_DIR` (WeakLibInterp provider thorn — find-or-build). MCNuX introduces **no** build-time option-list variables of its own; everything MCNuX-specific is a runtime Cactus parameter (`param.ccl`, owned by the other specs).

### Capability declarations (binding)

MCNuX's `configuration.ccl` declares exactly

```text
REQUIRES CarpetX WeakLibInterp
```

- `CarpetX` because the pinned access idiom reaches driver internals (`ghext`) — the obligation stated by `carpetx-thorn-integration.md` (MCNX-CTX-01); its capability chain transitively supplies AMReX, MPI, and the Loop/Arith/regrid capabilities.
- `WeakLibInterp` because MCNuX includes the library's headers and links its archive; its capability chain transitively supplies HDF5.
- MCNuX declares **no** direct `REQUIRES` on AMReX, HDF5, or MPI and **no** `PROVIDES` block of any kind: everything it needs beyond the two lines above arrives transitively, and duplicating a requirement would create a second place for version/backend skew to hide.

### Build outputs

- A Cactus executable in which MCNuX's objects (a) reference AMReX symbols from **exactly one** AMReX — the install the ET AMReX thorn provides, the same one CarpetX links and the same one the WeakLibInterp provider built the library against — and (b) contain pulled members of the WeakLibInterp archive (observable via the `wli_value_type_size` symbol).
- MCNuX contributes only plain `.cxx` sources listed in its `src/make.code.defn`; the build emits no MCNuX-owned shared library, no separate device-code objects, and no generated second build description.

### Suite invocation contract

How the two check families of `verification-suite-design.md` are invoked:

- **Family A (Cactus regression tests).** Each family-A matrix row is realized as `test/<name>.par` plus a same-named committed golden directory in the MCNuX thorn, and is discovered, run, and numerically diffed by the Cactus RunTest harness (`flesh/lib/sbin/RunTest.pl`) at the default `ABSTOL = RELTOL = 1e-12` (a per-test override via `test.ccl` is legal only where the matrix row records the rationale). Running the harness over the configuration executes every MCNuX family-A test; a diff violation or a runtime abort (including a statistical margin gate firing, `m < 0`) fails the test.
- **Family B (standalone unit checks).** Family-B checks build and run **without** a Cactus executable (host-only targets against the WeakLibInterp/AMReX libraries and MCNuX's host-testable kernels), each exiting nonzero on violation; an aggregate runner (make target, script, or CTest — freedom) propagates any failure to a nonzero exit. A check that cannot run in an environment (no device build, no production tables) reports a distinct, loud SKIP — never a silent pass.
- **Statistical checks** need no third runner: per `verification-suite-design.md`, they execute inside family-A runs (or family-B fixtures) as analysis-routine norm output, runtime-gated on the 4σ margin and committed as golden norm data.

## Correctness requirements

- **[MCNX-BLD-01] Capability closure (exact, structural).** `configuration.ccl` contains exactly the capability line `REQUIRES CarpetX WeakLibInterp`; MCNuX declares no direct requirement on AMReX, HDF5, or MPI and no `PROVIDES` block; the reference thornlist configures with zero unresolved-capability errors (the CST completes). Reference: the capability section above; provenance `CarpetX/CarpetX/configuration.ccl`, `WeakLibInterp/cactus/thorns/WeakLibInterp/configuration.ccl`.
- **[MCNX-BLD-02] One AMReX per executable (exact).** MCNuX never bundles, vendors, or builds AMReX or WeakLibInterp: its source tree contains no copied AMReX or WeakLibInterp sources and no nested build of either, and the built executable links exactly one AMReX per executable — the ET AMReX thorn's install, which the WeakLibInterp provider was also built against (`AMREX_DIR`, read-only). Two AMReX copies in one executable is a build error, not a degraded mode. Reference: the build-outputs section; provenance `WeakLibInterp/specs/cactus-integration.md`, `WeakLibInterp/cactus/thorns/WeakLibInterp/src/build.sh`.
- **[MCNX-BLD-03] Inline GPU compilation (exact, structural).** Every compiled MCNuX source is a `.cxx` file listed in `src/make.code.defn`; the thorn contains no `.cu` files, no `.hip` files, and no device-specific build description. All device code is inline (`amrex::ParallelFor` lambdas and header-inline WeakLibInterp `_Point` functions) and is compiled by the thorn's own configured compiler — under `AMREX_ENABLE_CUDA = yes` the ET AMReX thorn switches consumer thorns to the CUDA compiler, and MCNuX must compile unmodified under that switch. Reference: the CarpetX convention (`CarpetX/WaveToyX/src/wavetoyx.cxx` and its `make.code.defn`) and the consumer-compiles-device-headers rule of `WeakLibInterp/specs/cactus-integration.md`.
- **[MCNX-BLD-04] Suite invocation conformance (exact + golden `1e-12`).** Every family-A row of the coverage matrix is realized in the harness shape (`test/<name>.par` + same-named golden directory) and is discovered by `RunTest.pl`; a full harness run over the reference configuration executes all of them and passes at `ABSTOL = RELTOL = 1e-12` (or a row-recorded override). Every family-B check builds and runs without a Cactus executable and gates via exit status; the aggregate family-B run exits 0 with zero silent skips (environment-blocked checks report SKIP loudly). Reference: the suite-invocation section above; `verification-suite-design.md` owns the rows.
- **[MCNX-BLD-05] Minimal-consumer smoke (exact).** The reference configuration builds to completion in CPU mode, and the built executable passes a startup smoke assertion modeled on `TestWeakLibInterp` (`WeakLibInterp/cactus/thorns/TestWeakLibInterp/src/test_wli.cxx`): a scheduled startup routine asserts `wli_value_type_size() == 8` (proving real WeakLibInterp archive members were linked and the value type is double-sized), with one `_Point` entry point per consumed family odr-used behind a never-true guard (proving the device-header surface compiles with MCNuX's compiler, no tables needed). A zero-iteration smoke parameter file (driver up, MCNuX scheduled, no packets, no tables) runs to a clean exit. Reference: the smoke section of Inputs & outputs.

## Verification

- **MCNX-BLD-01**: a family-B structural check: fixed-string audit of `configuration.ccl` (the exact capability line present; no `REQUIRES AMReX`/`HDF5`/`MPI`; no `PROVIDES`) plus a configure of the reference thornlist asserting the CST completes with zero unresolved capabilities.
- **MCNX-BLD-02**: a family-B check in two parts: a source-tree scan asserting no vendored AMReX/WeakLibInterp sources or nested builds exist under the MCNuX thorn; and a symbol-level scan of the built executable (nm-style) asserting WeakLibInterp archive members are present and AMReX symbols resolve to the single Cactus-provided install (no duplicate-definition links).
- **MCNX-BLD-03**: a family-B structural check: `src/make.code.defn` lists only `.cxx` sources; a tree scan finds zero `.cu`/`.hip` files and zero secondary build descriptions. The CUDA-switch compile leg is environment-gated (loud SKIP where no CUDA toolchain exists) and asserts the unmodified tree compiles under `AMREX_ENABLE_CUDA = yes`.
- **MCNX-BLD-04**: run the Cactus harness (`RunTest.pl`) over the reference configuration and assert every MCNuX family-A test is discovered and passes; run the aggregate family-B runner and assert exit 0; cross-check discovered test names against the family-A rows of `verification-suite-design.md`'s matrix (every A-row realized, no orphan tests).
- **MCNX-BLD-05**: build the reference configuration in CPU mode; run the zero-iteration smoke parameter file; assert clean exit and the startup smoke assertion's success message (its failure is a `CCTK_ERROR` abort, so completion is the assertion).
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `WeakLibInterp/`, `CarpetX/`, and `flesh/` paths, and contains the required claim strings (the capability line, the one-AMReX rule, the `.cu` prohibition, the harness and smoke anchors).

## Implementation freedom

- Checkout layout: arrangement names/paths, the thornlist file's name and location, the Cactus configuration name, and whether GetComponents or a manual checkout assembles the tree — provided the thorn set and option-variable ownership above hold.
- The build environment: compilers, optimization flags, debug variants, option-list contents beyond the variables named above, and any simfactory machinery.
- How family-B checks are built and aggregated (Cactus make targets, a standalone CMake tree, a driver script) and how SKIP is reported — provided the exit-status gating and loud-SKIP rules hold.
- How the structural and symbol-level audits of MCNX-BLD-01/-02/-03 are scripted, and where the smoke parameter file lives.
- The name, wording, and schedule bin of the startup smoke routine — provided it runs before any transport work and aborts on failure.
- Richer configurations (a GRMHD evolution thorn, spacetime evolution, multi-patch) on top of the reference thornlist — out of scope for this spec's correctness, constrained only by the other specs.

## Open questions / assumptions

- **Provider-thorn checkout path (assumption, non-blocking).** The WeakLibInterp provider thorn is assumed reachable in the checkout as an arrangement entry into the WeakLibInterp repository (the `WeakLibInterp/cactus/wli.th` precedent lists `WeakLibInterp/WeakLibInterp` + `WeakLibInterp/TestWeakLibInterp`); the arrangement path is layout freedom, the thorn's contract is not.
- **GPU legs are compile-checked, not executed (assumption, mirrors WeakLibInterp).** The CUDA/HIP legs of MCNX-BLD-03 assert compilation only; GPU execution correctness is owned by the runtime specs (`particle-container-and-gpu.md`, `rng-and-statistical-acceptance.md`) and gated as SKIP where no device exists. The first real device build may surface flag-forwarding details; the configure-time consistency guard in the WeakLibInterp delivery makes a backend mismatch fail at configure time rather than at link.
- **Production tables are environment-gated (assumption, shared with `opacity-eos-evaluation.md`).** The reference build and the smoke run need no opacity/EOS tables; suite legs that need production tables follow the loud-SKIP discipline when the table root is absent.
- **Official ET distribution is deferred (assumption, non-blocking).** This contract is the group-use, build-from-checkout mode; ET release packaging (tarballs, licensing/doc review, thornlist registration) would land as a change to this spec without touching any other.
- **A single reference configuration is assumed sufficient (assumption).** The suite is defined against one reference thornlist/configuration; if a second configuration (e.g. a device build) becomes a standing gate, its thorn set and option pins are recorded here.
