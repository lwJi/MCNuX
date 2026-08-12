#!/bin/bash

# Build the `mcnux` Cactus configuration (pinned test ThornList of
# specs/build-and-integration.md [MCNX-BLD-04]) inside the Cactus tree at
# $CACTUSX. Lives alongside the CarpetX repo's `carpetx` configuration in the
# same tree; logs and configs are kept separate by name.

if [ -z "${CACTUSX:-}" ]; then
  echo "✗ missing required env var: CACTUSX" >&2
  exit 1
fi

# Option list: explicit ETKCFG (the sandbox interface) takes precedence over
# the host interface (ETKGUIDE). The thornlist is this repo's committed pinned
# list — the same everywhere by design; MCNXTHORNLIST overrides for experiments.
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cfg="${ETKCFG:-${ETKGUIDE:+$ETKGUIDE/macos.cfg}}"
thornlist="${MCNXTHORNLIST:-$repo/agent_scripts/sandbox/mcnux.th}"
if [ -z "$cfg" ]; then
  echo "✗ missing required env vars: set ETKGUIDE or ETKCFG" >&2
  exit 1
fi

log="$CACTUSX/last-build-mcnux.log"
if (
  cd "$CACTUSX" &&
  Compile-ETK -e mcnux -j8 \
    -c "$cfg" \
    -t "$thornlist"
) > "$log" 2>&1; then
  echo "✓ build"
else
  echo "--- first error lines ---"
  grep -n -m 12 -E 'error:|Error|\*\*\*' "$log" || echo "(no obvious error lines)"
  echo "--- last 40 lines ---"
  tail -40 "$log"
  echo "✗ build failed — full log: $log" >&2
  exit 1
fi
