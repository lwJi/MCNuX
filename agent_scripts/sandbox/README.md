# MCNuX sandbox environment for Docker Sandboxes (sbx)

MCNuX reuses the **CarpetX sandbox template** (`lwji/sandbox-templates:claude-carpetx`, recipe committed at `CarpetX/agent_scripts/sandbox/`) — the ETK dependency stack, the baked pristine Cactus tree at `$CACTUSX`, `Compile-ETK`, and the `ubuntu-arm64.cfg` option list all carry over unchanged. This directory adds only the MCNuX-specific wiring, applied at runtime by `setup.sh`; no separate image is needed.

## Contents

| File | Purpose |
| --- | --- |
| `mcnux.th` | The pinned test ThornList of `specs/build-and-integration.md` [MCNX-BLD-04] — a compile list over the baked Cactus checkout |
| `setup.sh` | One-time in-sandbox setup: arrangement symlinks (MCNuX, WeakLibInterp), AMReX build from the mounted source, weaklib-table symlink, validator env |

## Running

From the MCNuX repo root on the host:

```bash
cd ~/docker-workspace/repos/MCNuX
sbx run -t lwji/sandbox-templates:claude-carpetx --name claude-MCNuX claude . \
  ../amrex:ro ../warpx:ro ../CarpetX:ro ../WeakLibInterp:ro ../thornado:ro \
  /Users/liwei/Datas/wl_tables/use_for_production:ro
```

Mount roles: `amrex` (AMReX built from source by setup.sh, shared with the carpetx configuration), `CarpetX` (driver + base thorns; the baked tree's `arrangements/CarpetX` symlink resolves to it), `WeakLibInterp` (the opacity/EOS provider thorn builds the library from this checkout into the configuration's scratch space — the read-only mount is fine by design), `warpx`/`thornado` (reference reading + spec-validator citation roots only; not built), and the weaklib tables (runtime data for table-mode tests).

Inside the sandbox, once per sandbox (idempotent):

```bash
agent_scripts/sandbox/setup.sh
```

then the usual

```bash
./agent_scripts/build.sh   # Cactus configuration `mcnux`; first build configures from scratch, later builds are incremental
./agent_scripts/test.sh    # `make mcnux-testsuite` over MCNuX/test/
bash specs/tools/validate_specs.sh   # spec-corpus linter (setup.sh persists MCNX_FLESH_ROOT for it)
```

Named sandboxes persist across `sbx stop`/restart, so `configs/mcnux`, `exe/`, `TEST/`, and the AMReX install are cached and only the first build is slow. The `mcnux` configuration coexists with the CarpetX repo's `carpetx` configuration in the same Cactus tree (separate `configs/<name>` directories and log files).

## Table paths

The canonical table location in parameter files is `$HOME/Datas/wl_tables/use_for_production/wl-*.h5`:

- **on the host** (`$HOME = /Users/liwei`) that is the real directory;
- **in the sandbox** (`$HOME = /home/agent`) `setup.sh` symlinks `$HOME/Datas/wl_tables` to the mounted `/Users/liwei/Datas/wl_tables`.

So the same committed parameter file works in both places.
