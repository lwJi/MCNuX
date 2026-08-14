#!/bin/bash

# Non-harness gate for the paramcheck guard of
# specs/carpetx-thorn-integration.md [MCNX-CTX-06] ("Paramcheck guard
# (exact)"): a parameter file with CarpetX::use_subcycling = yes must abort at
# CCTK_PARAMCHECK with MCNuX's error; the identical file with `no` must run.
#
# This cannot be a Cactus regression test (specs/verification-suite-design.md:
# a deliberately-aborting run would be scored as a failure by the harness), so
# it drives the built executable directly. Run it after ./agent_scripts/build.sh,
# inside the sandbox. All transient output goes to $CACTUSX, never the repo.

set -u

if [ -z "${CACTUSX:-}" ]; then
  echo "✗ missing required env var: CACTUSX" >&2
  exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
gates="$repo/agent_scripts/gates"
exe="$CACTUSX/exe/cactus_mcnux"
outdir="$CACTUSX/GATE/mcnux"

if [ ! -x "$exe" ]; then
  echo "✗ missing executable: $exe (run ./agent_scripts/build.sh first)" >&2
  exit 1
fi

# The two parfiles are the same run except for the switch under test; drift
# between them would invalidate the "identical file" clause of the criterion.
if [ "$(grep -cv '^#\|^$' "$gates/subcycling_yes.par")" \
     != "$(grep -cv '^#\|^$' "$gates/subcycling_no.par")" ] ||
   [ "$(grep -v '^#\|^$\|use_subcycling' "$gates/subcycling_yes.par")" \
     != "$(grep -v '^#\|^$\|use_subcycling' "$gates/subcycling_no.par")" ]; then
  echo "✗ gate parfiles differ in more than CarpetX::use_subcycling" >&2
  exit 1
fi

mkdir -p "$outdir"
rc=0

run_case() { # <tag>; leaves output in $outdir/<tag>.log, returns exit status
  local tag="$1"
  ( cd "$CACTUSX" && "$exe" "$gates/subcycling_$tag.par" ) \
    > "$outdir/$tag.log" 2>&1
}

# Case 1: use_subcycling = yes must abort, and the abort must be MCNuX's own
# (a nonzero status from some other paramcheck routine would not prove the
# guard fired).
run_case yes
status=$?
if [ "$status" -eq 0 ]; then
  echo "✗ use_subcycling=yes: run completed (expected a paramcheck abort)"
  rc=1
elif ! grep -q 'in thorn MCNuX' "$outdir/yes.log" ||
     ! grep -q 'MCNuX requires CarpetX::use_subcycling = no' "$outdir/yes.log"; then
  echo "✗ use_subcycling=yes: aborted (status $status) but without MCNuX's error"
  grep -m 6 -A 3 'WARNING level 0' "$outdir/yes.log"
  rc=1
else
  echo "✓ use_subcycling=yes: aborted at paramcheck with MCNuX's error (status $status)"
fi

# Case 2: the identical file with `no` must run to completion.
run_case no
status=$?
if [ "$status" -ne 0 ]; then
  echo "✗ use_subcycling=no: run failed (status $status)"
  tail -20 "$outdir/no.log"
  rc=1
elif ! grep -q 'MCNuX_ParamCheck' "$outdir/no.log"; then
  echo "✗ use_subcycling=no: run completed but MCNuX_ParamCheck was not scheduled"
  rc=1
else
  echo "✓ use_subcycling=no: ran to completion with the guard scheduled"
fi

if [ "$rc" -eq 0 ]; then
  echo "✓ gate paramcheck (MCNX-CTX-06)"
else
  echo "✗ gate paramcheck failed — logs: $outdir/{yes,no}.log" >&2
fi
exit "$rc"
