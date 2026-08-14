# MCNuX

## Mission

MCNuX is a GPU Monte Carlo neutrino transport code implemented as a Cactus thorn on the CarpetX driver, with transport packets carried by the AMReX particle module and weak-interaction opacities/EOS quantities evaluated through the WeakLibInterp library. The canonical mission and correctness contract is `specs/README.md`; the 12 coded specs under `specs/` are the acceptance source of truth, and `specs/README.md` governs on any conflict. The build is driven by a Ralph loop reading `TODO.md`.

## Layout

- `specs/` — the spec corpus: `README.md` (canonical) + 12 coded leaf specs; `specs/tools/validate_specs.sh` is the corpus linter.
- `MCNuX/` — the production thorn (source root): `interface.ccl`, `param.ccl`, `configuration.ccl` (`REQUIRES CarpetX`, `REQUIRES WeakLibInterp`), `schedule.ccl`, `src/` (`make.code.defn`). Shared headers live flat in `MCNuX/src/`, named `mcnux_*.hxx`, contents in `namespace MCNuX` (first instance: `mcnux_units.hxx`, constants/units/species). `MCNuX/test/` (Cactus test home: `test.ccl`, parfiles, golden data) does not exist yet — created by the first benchmark task.
- `TestMCNuX/` — the sibling test thorn (initial-data writers for benchmarks): `interface.ccl`, `param.ccl`, `configuration.ccl`, `schedule.ccl`, `src/make.code.defn` (no sources yet).
- `agent_scripts/` — `build.sh`, `test.sh`, and `sandbox/` (`setup.sh`, `mcnux.th` pinned ThornList, `README.md`). Builds and the Cactus regression harness run only inside the sandbox.
- `TODO.md` — the durable build plan (written by the plan orchestrator only).
- `.research/` — transient plan-loop findings; regenerated every plan run.
- Reference checkouts are siblings: `../amrex`, `../warpx`, `../CarpetX`, `../WeakLibInterp` (flesh only inside the sandbox).

## Build & run

- Build: `./agent_scripts/build.sh` from the repo root, inside the sandbox only (needs `$CACTUSX`/`$ETKCFG`, set by `agent_scripts/sandbox/setup.sh`; arrangement symlinks must exist). Incremental rebuilds take seconds.
- Tests: `./agent_scripts/test.sh` (Cactus regression harness, sandbox only). Compile-time selftests (static_asserts in shared headers) run on every build via inclusion from `MCNuX/src/stub.cxx`.
