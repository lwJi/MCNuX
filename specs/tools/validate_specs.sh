#!/usr/bin/env bash
#
# validate_specs.sh — spec-set linter / harness for the MCNuX `specs/` set.
#
# PURPOSE
#   Mechanically enforce that every registered spec is a complete, self-contained
#   Ralph-loop contract and that the README index is internally consistent. The check is
#   fully CI-reproducible: it relies only on files committed under specs/ plus read-only
#   probes into the five local source-of-truth repositories (amrex, warpx, CarpetX, the
#   Cactus flesh, WeakLibInterp).
#
# CHECKS
#   For each registered spec:
#     (a) the 7 mandated section headers are present, in order;
#     (b) it names >= 1 concrete numeric tolerance (a scientific-notation literal);
#     (c) every cited amrex/ warpx/ CarpetX/ flesh/ WeakLibInterp/ source path resolves
#         to a real file under the corresponding repository root;
#     (d) its registered fixed-string claims are present.
#   For the README:
#     (e) the index table links every registered spec, every linked file exists on disk,
#         the on-disk spec count equals the registry size (no orphans, no missing links).
#   Corpus closure:
#     (f) every correctness-requirement id (MCNX-<TAG>-<NN>) declared by a registered
#         spec appears in verification-suite-design.md's coverage matrix, and every id
#         appearing there is declared by a registered spec (bidirectional closure). The
#         scan runs over whichever specs exist on disk, so it holds at every phase of the
#         corpus's growth.
#
# REPOSITORY ROOTS
#   Overridable via MCNUX_AMREX_ROOT, MCNUX_WARPX_ROOT, MCNUX_CARPETX_ROOT,
#   MCNUX_FLESH_ROOT, MCNUX_WEAKLIBINTERP_ROOT. When unset, sibling checkouts and the
#   standard local locations are probed; an unresolvable root fails loudly.
#
# EXTENDING (later phases)
#   Add one entry per new spec to REGISTERED_SPECS and its fixed-string claims to
#   SPEC_REQUIRE_IN_SPEC. New requirement ids and matrix rows are picked up by the
#   closure check automatically. No other part of this script needs editing.
#
# EXIT: 0 if every check passes; 1 otherwise (all failures in a run are reported).

set -u

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPECS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SPECS_DIR/.." && pwd)"
README="$SPECS_DIR/README.md"
SUITE_SPEC="verification-suite-design.md"

# Source-of-truth repository roots. Overridable via environment; if unset, probe a few
# candidate locations and pick the first containing the probe file (so a symlinked or
# non-adjacent checkout still resolves).
resolve_root() {
  # $1 = override value (may be empty); $2 = probe file relative to root; $3.. = candidates.
  local override="$1" probe="$2"; shift 2
  if [ -n "$override" ] && [ -f "$override/$probe" ]; then echo "$override"; return; fi
  local cand
  for cand in "$@"; do
    [ -n "$cand" ] || continue
    if [ -f "$cand/$probe" ]; then echo "$cand"; return; fi
  done
  # Nothing resolved: echo the override (or first candidate) for a meaningful error.
  echo "${override:-$1}"
}

PARENT="$(dirname "$REPO_ROOT")"
AMREX_ROOT="$(resolve_root "${MCNUX_AMREX_ROOT:-}" "Src/Base/AMReX_GpuContainers.H" \
  "$PARENT/amrex" "$HOME/docker-workspace/repos/amrex")"
WARPX_ROOT="$(resolve_root "${MCNUX_WARPX_ROOT:-}" "Source/Evolve/WarpXEvolve.cpp" \
  "$PARENT/warpx" "$HOME/docker-workspace/repos/warpx")"
CARPETX_ROOT="$(resolve_root "${MCNUX_CARPETX_ROOT:-}" "CarpetX/src/driver.cxx" \
  "$PARENT/CarpetX" "$HOME/docker-workspace/repos/CarpetX")"
FLESH_ROOT="$(resolve_root "${MCNUX_FLESH_ROOT:-}" "lib/sbin/RunTestUtils.pl" \
  "$PARENT/flesh" "$HOME/docker-workspace/EinsteinToolkit/Cactus/repos/flesh")"
WEAKLIBINTERP_ROOT="$(resolve_root "${MCNUX_WEAKLIBINTERP_ROOT:-}" "src/eos/wli_eos.H" \
  "$PARENT/WeakLibInterp" "$HOME/docker-workspace/repos/WeakLibInterp")"

FAILURES=0
fail()  { echo "FAIL: $*" >&2; FAILURES=$((FAILURES + 1)); }
pass()  { echo "ok:   $*"; }
info()  { echo "      $*"; }

# --------------------------------------------------------------------------------------
# The 7 mandated section headers, in order. Every registered spec must contain these
# `## ` headers in exactly this sequence (extra subsections may appear, but these seven
# anchors must appear and be ordered).
# --------------------------------------------------------------------------------------
SECTION_HEADERS=(
  "Purpose & scope"
  "Source of truth"
  "Inputs & outputs"
  "Correctness requirements"
  "Verification"
  "Implementation freedom"
  "Open questions / assumptions"
)

# --------------------------------------------------------------------------------------
# REGISTRY
# --------------------------------------------------------------------------------------
# REGISTERED_SPECS : spec basenames that must follow the 7-section template AND be linked
# from the README index. Later phases append to this array (and to the claims below).
REGISTERED_SPECS=(
  "conventions-and-units.md"
  "rng-and-statistical-acceptance.md"
  "packet-representation-and-sampling.md"
  "geodesic-propagation.md"
  "opacity-eos-evaluation.md"
  "neutrino-matter-interactions.md"
  "verification-suite-design.md"
)

# Per-spec fixed-string claims, expressed as "specfile|||needle". The needle (a fixed
# string, grep -F) must appear somewhere in the spec file.
SPEC_REQUIRE_IN_SPEC=(
  "conventions-and-units.md|||(−,+,+,+)"
  "conventions-and-units.md|||G = c = M_sun = 1"
  "conventions-and-units.md|||g = 4"
  "conventions-and-units.md|||1.160451812e10"
  "conventions-and-units.md|||6.175828e17"
  "conventions-and-units.md|||log10"
  "conventions-and-units.md|||ADMBaseX"
  "conventions-and-units.md|||HydroBaseX"
  "conventions-and-units.md|||lepton number"
  "conventions-and-units.md|||1e-14"
  "conventions-and-units.md|||Global correctness contract"
  "conventions-and-units.md|||X ∈ {E, P, S}"
  "conventions-and-units.md|||wli_eos_inversion.H"
  "rng-and-statistical-acceptance.md|||counter-based"
  "rng-and-statistical-acceptance.md|||(S, q, e, k)"
  "rng-and-statistical-acceptance.md|||(seed, packet id, event counter)"
  "rng-and-statistical-acceptance.md|||20260721"
  "rng-and-statistical-acceptance.md|||4σ"
  "rng-and-statistical-acceptance.md|||6.33e-5"
  "rng-and-statistical-acceptance.md|||1e-12"
  "packet-representation-and-sampling.md|||E_tot(cell, s) = α √γ ΔV Δt η_s"
  "packet-representation-and-sampling.md|||N_k = E_p / ε_k"
  "packet-representation-and-sampling.md|||equal-energy"
  "packet-representation-and-sampling.md|||isotropic in the fluid frame"
  "packet-representation-and-sampling.md|||six uniform draws"
  "packet-representation-and-sampling.md|||g = 4"
  "packet-representation-and-sampling.md|||ε = −p_μ u^μ"
  "packet-representation-and-sampling.md|||p_t = β^i p_i − α √(γ^{ij} p_i p_j)"
  "packet-representation-and-sampling.md|||1e-10"
  "packet-representation-and-sampling.md|||1e-14"
  "packet-representation-and-sampling.md|||4σ"
  "geodesic-propagation.md|||p^t = √(γ^{ij} p_i p_j) / α"
  "geodesic-propagation.md|||dx^i/dt = γ^{ij} p_j / p^t − β^i"
  "geodesic-propagation.md|||dp_i/dt = −α p^t ∂_i α + p_j ∂_i β^j − (p_j p_k / (2 p^t)) ∂_i γ^{jk}"
  "geodesic-propagation.md|||β^i p_i − α √(γ^{ij} p_i p_j)"
  "geodesic-propagation.md|||Kerr–Schild"
  "geodesic-propagation.md|||one cell width"
  "geodesic-propagation.md|||1e-14"
  "geodesic-propagation.md|||2e-2"
  "geodesic-propagation.md|||1e-2"
  "geodesic-propagation.md|||1e-12"
  "opacity-eos-evaluation.md|||already \`log10\`'d by the caller"
  "opacity-eos-evaluation.md|||1.160451812e10"
  "opacity-eos-evaluation.md|||| 0 |"
  "opacity-eos-evaluation.md|||| 01 |"
  "opacity-eos-evaluation.md|||| 02 |"
  "opacity-eos-evaluation.md|||| 03 |"
  "opacity-eos-evaluation.md|||| 10 |"
  "opacity-eos-evaluation.md|||| 11 |"
  "opacity-eos-evaluation.md|||| 13 |"
  "opacity-eos-evaluation.md|||T = 0"
  "opacity-eos-evaluation.md|||1.66054e3"
  "opacity-eos-evaluation.md|||3.16409e15"
  "opacity-eos-evaluation.md|||permissive"
  "opacity-eos-evaluation.md|||ResidentTable"
  "opacity-eos-evaluation.md|||TableView"
  "opacity-eos-evaluation.md|||nOpacities = 2"
  "opacity-eos-evaluation.md|||REQUIRES WeakLibInterp"
  "opacity-eos-evaluation.md|||Electron Antineutrino"
  "opacity-eos-evaluation.md|||28/3"
  "opacity-eos-evaluation.md|||4π E³ / (hc)³"
  "opacity-eos-evaluation.md|||μ_ν,1 = −μ_ν,0"
  "opacity-eos-evaluation.md|||1e-12"
  "opacity-eos-evaluation.md|||1e-14"
  "neutrino-matter-interactions.md|||Δt_a = −ln(r_a) · p^t / (κ_a ε)"
  "neutrino-matter-interactions.md|||Δt_s = −ln(r_s) · p^t / (κ_s ε)"
  "neutrino-matter-interactions.md|||memoryless"
  "neutrino-matter-interactions.md|||T^6"
  "neutrino-matter-interactions.md|||T_low = 0.5 MeV"
  "neutrino-matter-interactions.md|||re-isotropized in the fluid frame"
  "neutrino-matter-interactions.md|||1.476625e5"
  "neutrino-matter-interactions.md|||κ_a/(κ_a + κ_s)"
  "neutrino-matter-interactions.md|||u_eq = η_tot/κ_a"
  "neutrino-matter-interactions.md|||4σ"
  "neutrino-matter-interactions.md|||1e-14"
  "neutrino-matter-interactions.md|||1e-10"
  "neutrino-matter-interactions.md|||3e-2"
  "neutrino-matter-interactions.md|||1e-2"
  "verification-suite-design.md|||ABSTOL = RELTOL = 1e-12"
  "verification-suite-design.md|||out_tsv_vars"
  "verification-suite-design.md|||out_norm_vars"
  "verification-suite-design.md|||RunTestUtils.pl"
  "verification-suite-design.md|||Coverage matrix"
  "verification-suite-design.md|||m < 0"
  "verification-suite-design.md|||SKIP"
)

# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------

# Extract the ordered list of `## ` header titles (trimmed) from a markdown file.
spec_headers() {
  grep -E '^## ' "$1" | sed -E 's/^##[[:space:]]+//; s/[[:space:]]+$//'
}

# Check the 7 mandated headers appear in order: the spec's `## ` headers, filtered to the
# mandated set, must equal the mandated sequence.
check_section_order() {
  local file="$1" base; base="$(basename "$file")"
  local -a found=()
  local h want
  while IFS= read -r h; do
    for want in "${SECTION_HEADERS[@]}"; do
      if [ "$h" = "$want" ]; then found+=("$h"); break; fi
    done
  done < <(spec_headers "$file")

  if [ "${#found[@]}" -ne "${#SECTION_HEADERS[@]}" ]; then
    fail "$base: expected ${#SECTION_HEADERS[@]} mandated sections, found ${#found[@]} (${found[*]:-none})"
    return
  fi
  local i
  for i in "${!SECTION_HEADERS[@]}"; do
    if [ "${found[$i]}" != "${SECTION_HEADERS[$i]}" ]; then
      fail "$base: section #$((i+1)) is '${found[$i]}', expected '${SECTION_HEADERS[$i]}' (out of order)"
      return
    fi
  done
  pass "$base: 7 mandated sections present and in order"
}

# Check the spec names at least one concrete numeric tolerance (a scientific-notation
# literal such as 1e-12 or 4.9e-6).
check_has_tolerance() {
  local file="$1" base; base="$(basename "$file")"
  if grep -Eq '[0-9](\.[0-9]+)?[eE][-+]?[0-9]+' "$file"; then
    pass "$base: names a concrete numeric tolerance"
  else
    fail "$base: no concrete numeric tolerance found"
  fi
}

# Resolve every source-of-truth path the spec cites. We scan the whole file for tokens
# that look like a local-repo source path (a known repo prefix, then a relative path with
# a recognized extension) and check each resolves under its repository root.
PATH_TOKEN_RE='(amrex|warpx|CarpetX|flesh|WeakLibInterp)/[A-Za-z0-9_./-]+\.(H|hxx|hpp|h|cxx|cpp|c|ccl|md|F90|f90|par|pl|sh|tex|h5ls|txt)'

check_source_of_truth_paths() {
  local file="$1" base; base="$(basename "$file")"
  local any=0 ok=1
  local tok root rel resolved
  while IFS= read -r tok; do
    [ -z "$tok" ] && continue
    any=1
    case "$tok" in
      amrex/*)         root="$AMREX_ROOT";         rel="${tok#amrex/}" ;;
      warpx/*)         root="$WARPX_ROOT";         rel="${tok#warpx/}" ;;
      CarpetX/*)       root="$CARPETX_ROOT";       rel="${tok#CarpetX/}" ;;
      flesh/*)         root="$FLESH_ROOT";         rel="${tok#flesh/}" ;;
      WeakLibInterp/*) root="$WEAKLIBINTERP_ROOT"; rel="${tok#WeakLibInterp/}" ;;
      *) continue ;;
    esac
    resolved="$root/$rel"
    if [ -f "$resolved" ]; then
      info "$base: source-of-truth resolves: $tok"
    else
      fail "$base: source-of-truth path does not resolve: $tok (looked at $resolved)"
      ok=0
    fi
  done < <(grep -oE "$PATH_TOKEN_RE" "$file" | sort -u)

  if [ "$any" -eq 0 ]; then
    fail "$base: no local-repo source-of-truth path cited"
  elif [ "$ok" -eq 1 ]; then
    pass "$base: all cited source-of-truth paths resolve"
  fi
}

# Per-spec fixed-string assertions.
check_require_in_spec() {
  local entry="$1" file needle base
  file="${entry%%|||*}"; needle="${entry##*|||}"
  base="$file"; file="$SPECS_DIR/$file"
  if [ ! -f "$file" ]; then fail "$base: missing (required claim: $needle)"; return; fi
  if grep -Fq "$needle" "$file"; then
    pass "$base: contains required claim: $needle"
  else
    fail "$base: missing required claim: $needle"
  fi
}

# --------------------------------------------------------------------------------------
# README index integrity: every registered spec is linked, every link resolves, and the
# README links EXACTLY the set of spec .md files on disk (no orphans, no missing links).
# --------------------------------------------------------------------------------------
readme_linked_specs() {
  grep -oE '\]\(\.?/?[A-Za-z0-9_-]+\.md\)' "$README" | sed -E 's/^\]\(\.?\/?//; s/\)$//' | sort -u
}

check_readme_links() {
  if [ ! -f "$README" ]; then fail "README.md missing at specs/README.md"; return; fi
  local before="$FAILURES"

  local -a linked=()
  local target
  while IFS= read -r target; do
    [ -n "$target" ] && linked+=("$target")
  done < <(readme_linked_specs)

  # (e1) every registered spec is linked and present on disk.
  local s l hit
  for s in "${REGISTERED_SPECS[@]}"; do
    hit=0
    for l in "${linked[@]}"; do [ "$l" = "$s" ] && hit=1; done
    if [ "$hit" -eq 1 ]; then
      info "README links registered spec: $s"
    else
      fail "README index does not link registered spec: $s"
    fi
    if [ ! -f "$SPECS_DIR/$s" ]; then
      fail "registered spec file is missing on disk: $s"
    fi
  done

  # (e2) every link target exists on disk (no broken intra-repo links).
  for l in "${linked[@]}"; do
    if [ ! -f "$SPECS_DIR/$l" ]; then
      fail "README links a nonexistent file: $l"
    fi
  done

  if [ "$FAILURES" -eq "$before" ]; then
    pass "README index: all registered specs linked, all links resolve"
  fi
}

check_readme_exact_registered() {
  if [ ! -f "$README" ]; then fail "README.md missing at specs/README.md"; return; fi
  local before="$FAILURES"

  # Spec .md files present on disk (README.md itself is the index, not a spec).
  local -a on_disk=()
  local f
  while IFS= read -r f; do
    f="$(basename "$f")"
    [ "$f" = "README.md" ] && continue
    on_disk+=("$f")
  done < <(find "$SPECS_DIR" -maxdepth 1 -name '*.md' | sort)

  local -a linked=()
  local target
  while IFS= read -r target; do
    [ -n "$target" ] && linked+=("$target")
  done < <(readme_linked_specs)

  # exactly as many spec files on disk as there are registered specs.
  local want="${#REGISTERED_SPECS[@]}"
  if [ "${#on_disk[@]}" -eq "$want" ]; then
    pass "exactly $want spec files present on disk (excluding README)"
  else
    fail "expected exactly $want spec files on disk, found ${#on_disk[@]}: ${on_disk[*]}"
  fi

  # every on-disk spec is linked (no orphans).
  local l hit
  for f in "${on_disk[@]}"; do
    hit=0
    for l in "${linked[@]}"; do [ "$l" = "$f" ] && hit=1; done
    [ "$hit" -eq 1 ] || fail "orphan spec file on disk not linked from README: $f"
  done

  # every README spec link points at a file on disk (no dangling spec links).
  for l in "${linked[@]}"; do
    hit=0
    for f in "${on_disk[@]}"; do [ "$f" = "$l" ] && hit=1; done
    [ "$hit" -eq 1 ] || fail "README links a spec not present on disk: $l"
  done

  if [ "$FAILURES" -eq "$before" ]; then
    pass "README links exactly the $want on-disk spec files (no orphans, no missing)"
  fi
}

# --------------------------------------------------------------------------------------
# Requirement-id closure: every MCNX-<TAG>-<NN> id declared by a registered spec (other
# than the suite spec itself) appears in verification-suite-design.md's coverage matrix,
# and every id appearing there is declared by a registered spec. Runs over whichever
# specs exist on disk, so it holds at every phase boundary.
# --------------------------------------------------------------------------------------
REQID_RE='MCNX-[A-Z]+-[0-9]+'

check_reqid_closure() {
  local suite="$SPECS_DIR/$SUITE_SPEC"
  if [ ! -f "$suite" ]; then fail "$SUITE_SPEC missing"; return; fi

  # Requirement ids declared across the registered specs (excluding the suite spec).
  local -a declared=()
  local s id
  for s in "${REGISTERED_SPECS[@]}"; do
    [ "$s" = "$SUITE_SPEC" ] && continue
    [ -f "$SPECS_DIR/$s" ] || continue
    while IFS= read -r id; do
      [ -n "$id" ] && declared+=("$id")
    done < <(grep -ohE "$REQID_RE" "$SPECS_DIR/$s" | sort -u)
  done

  if [ "${#declared[@]}" -eq 0 ]; then
    fail "no requirement ids ($REQID_RE) declared by any registered spec"
    return
  fi

  # Forward closure: every declared id has a coverage-matrix appearance.
  local -a declared_u=()
  while IFS= read -r id; do
    [ -n "$id" ] && declared_u+=("$id")
  done < <(printf '%s\n' "${declared[@]}" | sort -u)

  for id in "${declared_u[@]}"; do
    if grep -Fq "$id" "$suite"; then
      pass "coverage matrix covers requirement id: $id"
    else
      fail "$SUITE_SPEC coverage matrix omits declared requirement id: $id"
    fi
  done

  # Reverse closure: every id the suite spec mentions is declared by a registered spec.
  local hit d
  while IFS= read -r id; do
    [ -n "$id" ] || continue
    hit=0
    for d in "${declared_u[@]}"; do [ "$d" = "$id" ] && hit=1; done
    if [ "$hit" -eq 1 ]; then
      info "matrix id is declared by a spec: $id"
    else
      fail "$SUITE_SPEC references requirement id declared by no registered spec: $id (stale row?)"
    fi
  done < <(grep -ohE "$REQID_RE" "$suite" | sort -u)
}

# --------------------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------------------
echo "=== MCNuX spec-set validator ==="
echo "specs dir:          $SPECS_DIR"
echo "amrex root:         $AMREX_ROOT"
echo "warpx root:         $WARPX_ROOT"
echo "CarpetX root:       $CARPETX_ROOT"
echo "flesh root:         $FLESH_ROOT"
echo "WeakLibInterp root: $WEAKLIBINTERP_ROOT"
echo

# Repository roots must resolve before path checks mean anything.
[ -d "$AMREX_ROOT" ]         || fail "amrex root does not resolve (set MCNUX_AMREX_ROOT)"
[ -d "$WARPX_ROOT" ]         || fail "warpx root does not resolve (set MCNUX_WARPX_ROOT)"
[ -d "$CARPETX_ROOT" ]       || fail "CarpetX root does not resolve (set MCNUX_CARPETX_ROOT)"
[ -d "$FLESH_ROOT" ]         || fail "flesh root does not resolve (set MCNUX_FLESH_ROOT)"
[ -d "$WEAKLIBINTERP_ROOT" ] || fail "WeakLibInterp root does not resolve (set MCNUX_WEAKLIBINTERP_ROOT)"

echo "--- per-spec: section order, tolerance, source-of-truth ---"
for s in "${REGISTERED_SPECS[@]}"; do
  f="$SPECS_DIR/$s"
  if [ ! -f "$f" ]; then fail "registered spec not found: $s"; continue; fi
  check_section_order "$f"
  check_has_tolerance "$f"
  check_source_of_truth_paths "$f"
done

echo
echo "--- per-spec: registered fixed-string claims ---"
for e in "${SPEC_REQUIRE_IN_SPEC[@]}"; do check_require_in_spec "$e"; done

echo
echo "--- README index integrity ---"
check_readme_links

echo
echo "--- closure: README links exactly the registered spec files ---"
check_readme_exact_registered

echo
echo "--- closure: requirement ids vs the coverage matrix (both directions) ---"
check_reqid_closure

echo
if [ "$FAILURES" -eq 0 ]; then
  echo "ALL CHECKS PASSED ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
  exit 0
else
  echo "VALIDATION FAILED: $FAILURES check(s) failed"
  exit 1
fi
