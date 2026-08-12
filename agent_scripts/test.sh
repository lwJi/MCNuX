#!/bin/bash

# Run the Cactus regression testsuite for the `mcnux` configuration. Tests are
# discovered per ThornList thorn from MCNuX/test/ (specs/build-and-integration.md
# [MCNX-BLD-05]: the whole suite is keyed to the MCNuX thorn).

if [ -z "${CACTUSX:-}" ]; then
  echo "✗ missing required env var: CACTUSX" >&2
  exit 1
fi

log="$CACTUSX/last-test-mcnux.log"
summary="$CACTUSX/TEST/mcnux/summary.log"

# PROMPT=no makes RunTest.pl take the default answer for every prompt,
# so the testsuite runs non-interactively regardless of TTY
if ! (
  cd "$CACTUSX" &&
  make mcnux-testsuite PROMPT=no
) > "$log" 2>&1; then
  echo "--- last 40 lines ---"
  tail -40 "$log"
  echo "✗ testsuite run failed — full log: $log" >&2
  exit 1
fi

# A zero-test pass is reported distinctly: "Number failed -> 0" with zero
# tests run would otherwise read as green when discovery is broken (or, as
# during the stub phase, when MCNuX/test/ does not exist yet)
if [ -f "$summary" ] && grep -Eq 'Number failed *-> *0' "$summary"; then
  npassed=$(grep -E 'Number of tests passed' "$summary" | grep -oE '[0-9]+' || echo 0)
  if [ "${npassed:-0}" -eq 0 ]; then
    echo "✓ test (but 0 tests ran — no test/*.par discovered for any ThornList thorn)"
  else
    echo "✓ test ($npassed passed)"
  fi
  exit 0
fi

# List failed tests as Thorn/test from the "Tests failed:" section
failed=$(awk '/Tests failed:/{f=1;next} f && /^=+/{exit} f && NF {gsub(/[()]/,""); print $3 "/" $1}' "$summary" 2>/dev/null)
if [ -z "$failed" ]; then
  echo "--- last 40 lines ---"
  tail -40 "$log"
  echo "✗ test failed (could not parse $summary) — full log: $log" >&2
  exit 1
fi

n=$(printf '%s\n' "$failed" | grep -c .)
echo "✗ test: $n failed"
printf '%s\n' "$failed" | head -20
[ "$n" -gt 20 ] && echo "... and $((n - 20)) more (see $summary)"

# Bounded excerpts for the first few failures; read the full per-test
# logs at $CACTUSX/TEST/mcnux/<Thorn>/<test>.log for more
printf '%s\n' "$failed" | head -3 | while IFS=/ read -r thorn test; do
  tlog="$CACTUSX/TEST/mcnux/$thorn/$test.log"
  diffs="$CACTUSX/TEST/mcnux/$thorn/$test.diffs"
  echo "--- $tlog (last 20 lines) ---"
  tail -20 "$tlog" 2>/dev/null || echo "(missing)"
  if [ -s "$diffs" ]; then
    echo "--- $diffs (first 15 lines) ---"
    head -15 "$diffs"
  fi
done

echo "✗ test failed — per-test logs: \$CACTUSX/TEST/mcnux/<Thorn>/<test>.log, run log: $log" >&2
exit 1
