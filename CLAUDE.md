# MCNuX

## Mission

MCNuX is a GPU Monte Carlo neutrino transport code implemented as a Cactus thorn on the CarpetX driver, with transport packets carried by the AMReX particle module and weak-interaction opacities and EOS quantities evaluated through the WeakLibInterp library. The canonical mission and correctness contract live in `specs/README.md` — it is the source of truth if anything here conflicts with it. Work is driven by `TODO.md` (the durable build plan) against the acceptance specs in `specs/`.

## Layout

- `specs/` — the spec corpus: `README.md` (canonical) + exactly 12 coded acceptance specs; `specs/tools/validate_specs.sh` is the spec-text linter (it checks spec files only, never code — a green run is not implementation progress).
- `MCNuX/` — the production transport thorn: `interface.ccl`, `param.ccl`, `schedule.ccl`, `configuration.ccl` (CCL contracts); source in `MCNuX/src/` (every new file must be added to `MCNuX/src/make.code.defn` `SRCS` or it silently does not compile); acceptance tests belong under `MCNuX/test/` (`test.ccl`, `.par` parfiles, golden data) — this directory does not exist yet and is created by the harness-bootstrap task.
- `TestMCNuX/` — the sibling test thorn (analytic initial-data writers only, no transport physics): same CCL + `TestMCNuX/src/` + `TestMCNuX/src/make.code.defn` structure.
- `agent_scripts/` — build/test entry points: `build.sh`, `test.sh`, and `sandbox/` (environment `setup.sh`, pinned ThornList `mcnux.th`).
- Shared C++ utilities (conventions header, RNG, tetrad transform, kernel idiom) live in `MCNuX/src/` as single-owner headers — write once, reuse; per-consumer copies are a test failure.
- Reference repos (read-only, sibling checkouts): `../amrex`, `../warpx`, `../CarpetX`, `../WeakLibInterp` (see `specs/README.md` for env-var overrides; the Cactus `flesh` repo is NOT checked out here).
