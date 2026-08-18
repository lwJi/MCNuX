#!/bin/bash

# CI: configure and build the `mcnux` Cactus configuration inside the tree
# assembled by agent_scripts/ci/download.sh.
#
# The ThornList is the pinned compile list of specs/build-and-integration.md
# [MCNX-BLD-04] (agent_scripts/sandbox/mcnux.th) — the same list the sandbox
# builds. The option list is CarpetX's committed CI optionlist for the
# lwji/carpetx:cpu-real64 image, taken from the pinned CarpetX checkout so the
# two stay in lockstep. Plain `make` (the Compile-ETK vocabulary without the
# wrapper); no SimFactory.

set -ex

mcnux_src="${MCNUX_SRC:-$PWD/MCNuX}"
carpetx_src="${CARPETX_SRC:-$PWD/CarpetX}"
cactusx="${CACTUSX:-$PWD/Cactus}"
cfg="${MCNX_CI_CFG:-$carpetx_src/scripts/actions-cpu-real64.cfg}"
thornlist="$mcnux_src/agent_scripts/sandbox/mcnux.th"

test -f "$cfg"
test -f "$thornlist"

cd "$cactusx"
make mcnux-config PROMPT=no options="$cfg" THORNLIST="$thornlist"
time make -j"$(nproc)" mcnux

# The harness and the paramcheck gate both drive this executable
test -x exe/cactus_mcnux
