#!/bin/bash

# One-time setup inside an sbx sandbox created from the CarpetX sandbox
# template (lwji/sandbox-templates:claude-carpetx) with MCNuX as the primary
# workspace (idempotent, safe to re-run):
#   1. point the baked Cactus tree's arrangements/MCNuX at this repo
#   2. point arrangements/WeakLibInterp at the mounted WeakLibInterp repo's
#      cactus/thorns directory (the WeakLibInterp arrangement)
#   3. build AMReX from the mounted (read-only) amrex workspace into
#      $HOME/cactus/amrex-lib (shared with the carpetx configuration; skipped
#      if already installed); without an amrex mount, fall back to the AMReX
#      release baked into the image at /usr/local
#   4. link $HOME/Datas/wl_tables at the mounted weaklib table directory, so
#      table paths written as $HOME/Datas/wl_tables/use_for_production/... work
#      identically on the host (where $HOME/Datas/wl_tables is the real
#      directory) and in the sandbox
#   5. persist MCNX_FLESH_ROOT=$CACTUSX so specs/tools/validate_specs.sh
#      resolves flesh/ citations against the baked Cactus tree
#
# Usage: agent_scripts/sandbox/setup.sh [--force-amrex]

set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [ -z "${CACTUSX:-}" ] || [ ! -d "$CACTUSX" ]; then
  echo "✗ CACTUSX is unset or does not exist (is this the sandbox image?)" >&2
  exit 1
fi

# 1. Link this repo into the Cactus tree as the MCNuX arrangement
arr="$CACTUSX/arrangements/MCNuX"
if [ "$(readlink "$arr" 2>/dev/null)" != "$repo" ]; then
  rm -rf "$arr"
  ln -s "$repo" "$arr"
fi
echo "✓ $arr -> $repo"

# 2. Link the WeakLibInterp arrangement (thorn dirs live under cactus/thorns)
wli_src="${WLI_SRC:-$(dirname "$repo")/WeakLibInterp}"
wli_arr_src="$wli_src/cactus/thorns"
if [ ! -f "$wli_arr_src/WeakLibInterp/configuration.ccl" ]; then
  echo "✗ WeakLibInterp arrangement not found at $wli_arr_src (mount ../WeakLibInterp:ro, or set WLI_SRC)" >&2
  exit 1
fi
arr="$CACTUSX/arrangements/WeakLibInterp"
if [ "$(readlink "$arr" 2>/dev/null)" != "$wli_arr_src" ]; then
  rm -rf "$arr"
  ln -s "$wli_arr_src" "$arr"
fi
echo "✓ $arr -> $wli_arr_src"

# 3. AMReX (same recipe and install location as CarpetX's setup.sh, so the
#    carpetx and mcnux configurations share one build)
amrex_src="${AMREX_SRC:-$(dirname "$repo")/amrex}"
amrex_lib="$HOME/cactus/amrex-lib"
amrex_build="$HOME/cactus/amrex-build"

if [ ! -d "$amrex_src" ]; then
  if [ ! -e "$amrex_lib" ]; then
    ln -s /usr/local "$amrex_lib"
  fi
  echo "✓ no amrex workspace mounted; using baked AMReX ($amrex_lib -> $(readlink "$amrex_lib" 2>/dev/null || echo "$amrex_lib"))"
elif [ -d "$amrex_lib/lib" ] && [ ! -L "$amrex_lib" ] && [ "${1:-}" != "--force-amrex" ]; then
  echo "✓ AMReX already installed in $amrex_lib (rerun with --force-amrex to rebuild)"
else
  echo "Building AMReX from $amrex_src ..."
  [ -L "$amrex_lib" ] && rm "$amrex_lib"
  rm -rf "$amrex_build" "$amrex_lib"
  log="$HOME/cactus/amrex-build.log"

  # Same options as the host's Install-AMReX; out-of-source build because
  # the amrex workspace is mounted read-only
  if (
    export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 &&
    cmake -B "$amrex_build" -S "$amrex_src" \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX="$amrex_lib" \
      -DAMReX_FORTRAN=OFF \
      -DAMReX_FORTRAN_INTERFACES=OFF \
      -DAMReX_OMP=ON \
      -DAMReX_PARTICLES=ON \
      -DAMReX_PRECISION=DOUBLE &&
    cmake --build "$amrex_build" -j"$(nproc)" --target install
  ) > "$log" 2>&1; then
    echo "✓ AMReX installed in $amrex_lib"
  else
    echo "--- last 40 lines ---"
    tail -40 "$log"
    echo "✗ AMReX build failed — full log: $log" >&2
    exit 1
  fi
fi

# 4. WeakLib opacity/EOS tables at the canonical $HOME/Datas/wl_tables path
tables_src="${WL_TABLES_SRC:-/Users/liwei/Datas/wl_tables}"
tables="$HOME/Datas/wl_tables"
probe="use_for_production/wl-EOS-SFHo-15-25-50.h5"
if [ ! -r "$tables/$probe" ]; then
  if [ ! -r "$tables_src/$probe" ]; then
    echo "✗ weaklib tables not found at $tables_src/$probe (mount /Users/liwei/Datas/wl_tables/use_for_production:ro, or set WL_TABLES_SRC)" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$tables")"
  ln -sfn "$tables_src" "$tables"
fi
echo "✓ tables at $tables/use_for_production"

# 5. Spec-validator flesh root (the baked Cactus tree carries lib/sbin/*.pl);
#    /etc/sandbox-persistent.sh is sourced before every shell in the sandbox
persist=/etc/sandbox-persistent.sh
if [ -w "$persist" ]; then
  if ! grep -q '^export MCNX_FLESH_ROOT=' "$persist"; then
    echo "export MCNX_FLESH_ROOT=$CACTUSX" >> "$persist"
  fi
  echo "✓ MCNX_FLESH_ROOT=$CACTUSX persisted in $persist"
else
  echo "! $persist not writable — export MCNX_FLESH_ROOT=$CACTUSX yourself for specs/tools/validate_specs.sh" >&2
fi
