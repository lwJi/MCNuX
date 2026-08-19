# MCNuX

## Mission

MCNuX is a GPU Monte Carlo neutrino transport code implemented as a Cactus thorn on the CarpetX driver, with transport packets carried by the AMReX particle module and weak-interaction opacities/EOS quantities evaluated through the WeakLibInterp library. The canonical mission and correctness contract is `specs/README.md`; the 12 coded specs under `specs/` are the acceptance source of truth, and `specs/README.md` governs on any conflict. The build is driven by a Ralph loop reading `TODO.md`.

## Layout

- `specs/` — the spec corpus: `README.md` (canonical) + 12 coded leaf specs; `specs/tools/validate_specs.sh` is the corpus linter.
- `MCNuX/` — the production thorn (source root): `interface.ccl`, `param.ccl`, `configuration.ccl` (`REQUIRES CarpetX`, `REQUIRES WeakLibInterp`), `schedule.ccl`, `src/` (`make.code.defn`). Shared headers live flat in `MCNuX/src/`, named `mcnux_*.hxx`, contents in `namespace MCNuX` (`mcnux_units.hxx` constants/units/species, `mcnux_rng.hxx` Philox4x32-10 `u(S,q,e,k)`); compiled sources: `mcnux_paramcheck.cxx` (subcycling guard), `mcnux_cadence.cxx` (transport cadence group + `ghext` relative-include idiom), `mcnux_selftest.cxx` (append-only runtime battery), `stub.cxx` (compile-time selftest aggregation point — every shared header is `#include`d there). `MCNuX/test/` is the Cactus test home: `test.ccl` (`NPROCS 1`), `<name>.par` parfiles, `<name>/` golden dirs (golden TSVs copied verbatim from sandbox harness runs, never hand-authored); first test: `unit-selftest`. Benchmark parfiles set `IO::out_dir = $parfile`, `IO::out_fileinfo = "axis labels"`, `IO::parfile_write = no`, `CarpetX::out_metadata = no`.
- `TestMCNuX/` — the sibling test thorn (initial-data writers for benchmarks): `interface.ccl`, `param.ccl`, `configuration.ccl`, `schedule.ccl`, `src/make.code.defn` (no sources yet).
- `agent_scripts/` — `build.sh`, `test.sh`, `gate_paramcheck.sh` + `gates/*.par` (non-harness exact-criterion gates), `sandbox/` (`setup.sh`, `mcnux.th` pinned ThornList, `README.md`), and `ci/` (GitHub Actions helpers: `mcnux-fetch.th` fetch list, `download.sh`, `build.sh` — consumed by `.github/workflows/ci.yml`, never run locally). Builds and the Cactus regression harness run only inside the sandbox (or CI).
- `TODO.md` — the durable build plan (written by the plan orchestrator only).
- `.research/` — transient plan-loop findings; regenerated every plan run.
- Reference checkouts are siblings: `../amrex`, `../warpx`, `../CarpetX`, `../WeakLibInterp` (flesh only inside the sandbox).

## Build & run

- Build: `./agent_scripts/build.sh` from the repo root, inside the sandbox only (needs `$CACTUSX`/`$ETKCFG`, set by `agent_scripts/sandbox/setup.sh`; arrangement symlinks must exist). Incremental rebuilds take seconds.
- Tests: `./agent_scripts/test.sh` (Cactus regression harness, sandbox only). Compile-time selftests (static_asserts in shared headers) run on every build via inclusion from `MCNuX/src/stub.cxx`.
- Non-harness gates: `./agent_scripts/gate_paramcheck.sh` (sandbox only, after `build.sh`) — exact-criterion check for deliberately-aborting runs the harness cannot express; parfiles in `agent_scripts/gates/*.par`, transient output under `$CACTUSX/GATE/mcnux/`, never in the repo. Precedent for later gate tasks (T25b/T25c). Any MCNuX parfile needs `Cactus::presync_mode = "mixed-error"` to pass CarpetX's own paramcheck.
- Spec-corpus linter: `bash specs/tools/validate_specs.sh` from the repo root, sandbox only (no build needed) — the `flesh/` reference root resolves solely through `MCNX_FLESH_ROOT`, persisted by `agent_scripts/sandbox/setup.sh`; outside the sandbox every `flesh/...` citation check fails spuriously (environment defect, never a corpus defect).
