# MCNuX

## Mission

MCNuX is a GPU Monte Carlo neutrino transport code implemented as a Cactus thorn on the CarpetX driver, with transport packets carried by the AMReX particle module and weak-interaction opacities/EOS quantities evaluated through the WeakLibInterp library. The canonical mission and correctness contract is `specs/README.md`; the 12 coded specs under `specs/` are the acceptance source of truth, and `specs/README.md` governs on any conflict. The build is driven by a Ralph loop reading `TODO.md`.

## Layout

- `specs/` — the spec corpus: `README.md` (canonical) + 12 coded leaf specs; `specs/tools/validate_specs.sh` is the corpus linter.
- `MCNuX/` — the production thorn: `interface.ccl`, `param.ccl`, `schedule.ccl`, `configuration.ccl` (`REQUIRES CarpetX WeakLibInterp Loop`), `src/` (see Source map), `test/` (see Testing).
- `TestMCNuX/` — the sibling test thorn (initial-data writers for benchmarks). `INHERITS: ADMBaseX HydroBaseX`, so every parfile activating it must activate both. Its `param.ccl` EXTENDS ADMBaseX/HydroBaseX keyword parameters cross-repo (`"Schwarzschild"`, `"uniform sphere"`, `"equilibration box"`). Ships no `test/` dir ([MCNX-BLD-05]); its parfiles live under `MCNuX/test/` ([MCNX-BLD-06]).
- `agent_scripts/` — `build.sh`, `test.sh`, `gate_paramcheck.sh` + `gates/*.par`, `sandbox/` (`setup.sh`, `mcnux.th` pinned ThornList, `README.md`), `ci/` (GitHub Actions helpers consumed by `.github/workflows/ci.yml`, never run locally).
- `TODO.md` — the durable build plan (written by the plan orchestrator only).
- `.research/` — transient plan-loop findings; regenerated every plan run.
- Reference checkouts are siblings: `../amrex`, `../warpx`, `../CarpetX`, `../WeakLibInterp` (flesh only inside the sandbox).

## Source map

One line per file; read the file itself for signatures and details. Shared headers live flat in `MCNuX/src/`, named `mcnux_*.hxx`, contents in `namespace MCNuX`.

Headers:

- `mcnux_units.hxx` — pinned constants, unit conversions, species enumeration.
- `mcnux_rng.hxx` — Philox4x32-10 stateless counter-based RNG.
- `mcnux_particles.hxx` — pure-SoA `PacketContainer`; declares the packet-population accessors owned by `mcnux_geodesic.cxx`.
- `mcnux_opacity.hxx` — WeakLibInterp call-boundary wrappers.
- `mcnux_tetrad.hxx` — inverse metric, null closure, tetrad construction.
- `mcnux_trp.hxx` — trapped-regime opacity relabeling map.
- `mcnux_srcterms.hxx` — exchange-ledger arithmetic.
- `mcnux_coefficients.hxx` — source-agnostic coefficient interface, analytic Kirchhoff-form formulas, table-vs-analytic dispatch.
- `mcnux_table_coeffs.hxx` — baseline table-source coefficient assembly + νx dataset mapping.
- `mcnux_table_range.hxx` — table-range policy: clamping with counters, table-derived transparency floor, `RangedTableCoefficients` wrapper.
- `mcnux_geodesic.hxx` — geodesic push: metric snapshot/gather, RK4 step.
- `mcnux_interactions.hxx` — interaction-time draw, fluid-frame energy ν, episode competition, episode draw-map wrappers.
- `mcnux_fluid.hxx` — cell-centered HydroBaseX fluid-state gather (`CellFluidGather`/`FluidSample`, containing-cell, no sub-cell interpolation) + Valencia v^i→u^μ lift.
- `mcnux_stats.hxx` — 4σ statistical-acceptance reduction; pinned seed and packet-count constants.

Compiled sources:

- `mcnux_paramcheck.cxx` — subcycling guard.
- `mcnux_cadence.cxx` — transport cadence group.
- `mcnux_coefficients.cxx` — parameter glue: coefficient-source selection and analytic params from `param.ccl` arrays.
- `mcnux_srcterms.cxx` — source-term zeroing + synthetic-deposit contributor.
- `mcnux_geodesic.cxx` — runtime packet-population owner (lazy per-patch containers), synthetic-packet seeding fixture, geodesic push scheduled `IN MCNuX_TransportStep`, fluid-gather diagnostic routine (`MCNuX_FluidGatherDiag`, gated `test_fluid_gather`).
- `mcnux_selftest.cxx` / `mcnux_selftest.hxx` / `mcnux_selftest_<domain>.cxx` — runtime selftest battery: core row table and ordered appender sequence in the first, shared machinery and appender declarations in the header, check code in the per-domain files.
- `stub.cxx` — compile-time selftest aggregation point.

`TestMCNuX/src/`:

- `schwarzschild.cxx` — isotropic Schwarzschild metric/lapse writers (vertex-centered).
- `hydro_profiles.cxx` — piecewise (ρ, T[MeV], Yₑ) hydro-profile writer (cell-centered, plus staggered `Avec` loops).

## Idioms & invariants

- Packet operators walk tiles via `for_each_packet_tile` + `amrex::ParallelFor` — never `loop_all_device` for packets. Grid functions use `loop_all_device` from `loop_device.hxx` with the field's centering (`<1,1,1>` ccc, `<0,0,0>` vertex).
- Table consumers use `RangedTableCoefficients`, never the bare `BaselineTableCoefficients` assembly.
- The selftest index→name row table in `mcnux_selftest.cxx` is APPEND-ONLY; new checks go in a per-domain `mcnux_selftest_<domain>.cxx` with its appender declared in `mcnux_selftest.hxx`.
- Every shared header is `#include`d from `stub.cxx` so its static_asserts run on every build.
- CarpetX's `ghext` is reached via the relative-include idiom (`mcnux_cadence.cxx` is the precedent).
- The `mcnux_packet_diag` grid array is the per-packet golden observable for transport tests; `mcnux_fluid_diag` is its analogue for the fluid-state gather.

## Testing

- `MCNuX/test/` is the Cactus test home: `test.ccl` (`NPROCS 1`), `<name>.par` parfiles, `<name>/` golden dirs. Current tests: `unit-selftest`, `source-zero-add`, `minkowski-freestream`, `schwarzschild-id`, `hydro-sphere-id`, `hydro-sphere-gather`.
- Golden TSVs are copied verbatim from sandbox harness runs, never hand-authored. Harness discovery: a `test/<name>.par` is only runnable once its golden dir `test/<name>/` exists; the first run against an empty golden dir is the capture step (vacuous "0 files identical").
- Parfile requirements: every MCNuX parfile needs `Cactus::presync_mode = "mixed-error"` and both `ADMBaseX` and `HydroBaseX` in `ActiveThorns` (MCNuX `INHERITS: ADMBaseX HydroBaseX`); runs that read the metric also need `ODESolvers` active (ADMBaseX initial data is scheduled only `IN ODESolvers_Initial`).
- CarpetX/AMReX rejects domains < 8 cells per direction (blocking_factor); off-origin test boxes must set `IO::out_{x,y,z}line_*` inside the domain.
- Benchmark parfiles set `IO::out_dir = $parfile`, `IO::out_fileinfo = "axis labels"`, `IO::parfile_write = no`, `CarpetX::out_metadata = no`.

## Build & run

- Build: `./agent_scripts/build.sh` from the repo root, inside the sandbox only (needs `$CACTUSX`/`$ETKCFG`, set by `agent_scripts/sandbox/setup.sh`; arrangement symlinks must exist). Incremental rebuilds take seconds. After a reverted/discarded task, stale objects can linger in `$CACTUSX/configs/mcnux/build/MCNuX/` (e.g. `.cxx/.cxx.d/.cxx.o` for deleted sources, or `.cxx.d` dependency files still naming a deleted/renamed *header* — "No rule to make target") and break the build — remove them (or clean-rebuild the thorn) before diagnosing further.
- Tests: `./agent_scripts/test.sh` (Cactus regression harness, sandbox only). Compile-time selftests (static_asserts in shared headers) run on every build via `stub.cxx`.
- Non-harness gates: `./agent_scripts/gate_paramcheck.sh` (sandbox only, after `build.sh`) — exact-criterion check for deliberately-aborting runs the harness cannot express; parfiles in `agent_scripts/gates/*.par`, transient output under `$CACTUSX/GATE/mcnux/`, never in the repo.
- Spec-corpus linter: `bash specs/tools/validate_specs.sh` from the repo root, sandbox only (no build needed) — the `flesh/` reference root resolves solely through `MCNX_FLESH_ROOT`, persisted by `agent_scripts/sandbox/setup.sh`; outside the sandbox every `flesh/...` citation check fails spuriously (environment defect, never a corpus defect).
