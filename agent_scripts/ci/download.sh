#!/bin/bash

# CI: assemble the Cactus tree at $CACTUSX ($PWD/Cactus by default).
#
# GetComponents fetches the flesh, CactusBase/IOUtil, and the ExternalLibraries
# providers (agent_scripts/ci/mcnux-fetch.th); the CarpetX, WeakLibInterp, and
# MCNuX arrangements are symlinked to the checkouts the workflow put next to
# this repo — the same wiring agent_scripts/sandbox/setup.sh applies in the
# sandbox. Modeled on CarpetX/scripts/download.sh.
#
# Expected layout (all under the workflow workspace, which is $PWD):
#   MCNuX/          this repo             (MCNUX_SRC overrides)
#   CarpetX/        pinned CarpetX        (CARPETX_SRC overrides)
#   WeakLibInterp/  pinned WeakLibInterp  (WLI_SRC overrides)
#   Cactus/         created here

set -ex

mcnux_src="${MCNUX_SRC:-$PWD/MCNuX}"
carpetx_src="${CARPETX_SRC:-$PWD/CarpetX}"
wli_src="${WLI_SRC:-$PWD/WeakLibInterp}"
cactusx="${CACTUSX:-$PWD/Cactus}"

# GetComponents writes to ./Cactus (the fetch list's !DEFINE ROOT)
test "$cactusx" = "$PWD/Cactus"
test -f "$carpetx_src/CarpetX/configuration.ccl"
test -f "$wli_src/cactus/thorns/WeakLibInterp/configuration.ccl"

wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$mcnux_src/agent_scripts/ci/mcnux-fetch.th"

# Arrangement symlinks for the three repos CI checks out directly
mkdir -p "$cactusx/arrangements"
ln -sfn "$carpetx_src" "$cactusx/arrangements/CarpetX"
ln -sfn "$wli_src/cactus/thorns" "$cactusx/arrangements/WeakLibInterp"
ln -sfn "$mcnux_src" "$cactusx/arrangements/MCNuX"
