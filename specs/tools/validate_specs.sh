#!/usr/bin/env bash
#
# validate_specs.sh — spec-set linter / harness for the MCNuX `specs/` set.
#
# PURPOSE
#   Mechanically enforce that every spec is a complete, self-contained Ralph-loop
#   contract and that the README index is internally consistent. Structural port of
#   WeakLibInterp/specs/tools/validate_specs.sh with MCNuX-specific registries and checks.
#
# CHECKS (all default-mode, CI-reproducible; no tables, no build needed)
#   Per coded spec:
#     (a) the 7 mandated section headers are present, in order;
#     (b) it names >= 1 concrete numeric tolerance;
#     (c) every repo-prefixed reference path it cites resolves under one of the five
#         reference roots (amrex, warpx, CarpetX, flesh, WeakLibInterp);
#     (d) registered fixed-string claims (SPEC_REQUIRE_IN_SPEC) are present;
#     (e) registered structural claims are anchored in WeakLibInterp's committed
#         specs/fixtures/*.h5ls snapshots under the resolved WeakLibInterp root
#         (cited in place — MCNuX vendors no fixture copies);
#     (f) citation-registry closure: every inline-cited bare arXiv id (and DOI) has a
#         version-pinned row in the spec's "## Source of truth" registry, that row
#         appears in the README master registry, and the pinned versions agree;
#     (g) requirement-id closure: every [MCNX-XXX-NN] id occurring in a coded spec has
#         a known owning spec (its code), is declared in the owner's
#         "## Correctness requirements" section (code-ownership rule), and is covered
#         by verification-suite-design.md's coverage matrix.
#   Corpus-wide:
#     (h) each of the five reference roots is cited at least once (five-root closure);
#     (i) README link closure: every registered spec linked, every link resolves, and
#         the on-disk spec set is EXACTLY the 13-file registry (no orphans, no extras).
#
# ROOT RESOLUTION
#   Five env overrides — MCNX_AMREX_ROOT, MCNX_WARPX_ROOT, MCNX_CARPETX_ROOT,
#   MCNX_FLESH_ROOT, MCNX_WEAKLIB_INTERP_ROOT — each honored only if it contains the
#   per-repo probe file; otherwise candidates (sibling checkout, then a fixed $HOME
#   layout) are searched, so a stale override self-heals on this checkout.
#
# EXTENDING
#   The registries below (REGISTERED_SPECS, CODED_SPECS, SPEC_REQUIRE_IN_SPEC,
#   SPEC_REQUIRE_IN_SNAPSHOT) are the sole edit points for growing coverage.
#   Every appended assertion must be proven live by an (uncommitted) red-run:
#   copy specs/ to scratch, inject the matching defect, require a nonzero exit.
#   This script locates specs/ relative to its own path, so a scratch copy is
#   self-contained and runnable in place.
#
# EXIT: 0 iff every check passes; all failures within a run are reported before exiting.

set -u

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPECS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SPECS_DIR/.." && pwd)"
README="$SPECS_DIR/README.md"

# Reference source-of-truth repos. Overridable via environment for portability; an
# override is honored only if it actually contains the probe file, otherwise the
# candidate list (sibling checkout, then fixed $HOME layout) is searched.
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
AMREX_ROOT="$(resolve_root "${MCNX_AMREX_ROOT:-}" "Src/Particle/AMReX_ParticleTile.H" \
  "$PARENT/amrex" "$HOME/docker-workspace/repos/amrex")"
WARPX_ROOT="$(resolve_root "${MCNX_WARPX_ROOT:-}" "Source/Particles/WarpXParticleContainer.H" \
  "$PARENT/warpx" "$HOME/docker-workspace/repos/warpx")"
CARPETX_ROOT="$(resolve_root "${MCNX_CARPETX_ROOT:-}" "CarpetX/src/driver.hxx" \
  "$PARENT/CarpetX" "$HOME/docker-workspace/repos/CarpetX")"
FLESH_ROOT="$(resolve_root "${MCNX_FLESH_ROOT:-}" "lib/sbin/RunTestUtils.pl" \
  "$PARENT/flesh" "$HOME/docker-workspace/EinsteinToolkit/Cactus/repos/flesh")"
WEAKLIB_INTERP_ROOT="$(resolve_root "${MCNX_WEAKLIB_INTERP_ROOT:-}" "specs/fixtures/tables.provenance" \
  "$PARENT/WeakLibInterp" "$HOME/docker-workspace/repos/WeakLibInterp")"

# WeakLibInterp's committed table snapshots, cited in place (MCNuX vendors no copies).
FIXTURES_DIR="$WEAKLIB_INTERP_ROOT/specs/fixtures"

FAILURES=0
fail()  { echo "FAIL: $*" >&2; FAILURES=$((FAILURES + 1)); }
pass()  { echo "ok:   $*"; }
info()  { echo "      $*"; }

# --------------------------------------------------------------------------------------
# The 7 mandated section headers, in order. Every coded spec must contain these `## `
# headers in exactly this sequence (extra `## ` sections may appear; these seven anchors
# must appear and be ordered). README.md is the one file with a different anatomy.
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
# REGISTERED_SPECS : the complete corpus — exactly these 13 files, pinned from day one.
#                    The on-disk spec set must equal this registry exactly.
REGISTERED_SPECS=(
  "README.md"
  "conventions-and-units.md"
  "rng-and-statistical-acceptance.md"
  "packet-representation-and-sampling.md"
  "geodesic-propagation.md"
  "neutrino-matter-interactions.md"
  "trapped-regime-treatment.md"
  "opacity-eos-evaluation.md"
  "hydro-coupling-source-terms.md"
  "particle-container-and-gpu.md"
  "carpetx-thorn-integration.md"
  "verification-suite-design.md"
  "build-and-integration.md"
)

# CODED_SPECS : the 12 coded specs, "CODE|||basename". The code owns the [MCNX-CODE-NN]
#               requirement-id namespace (code-ownership rule).
CODED_SPECS=(
  "CNV|||conventions-and-units.md"
  "RNG|||rng-and-statistical-acceptance.md"
  "PKT|||packet-representation-and-sampling.md"
  "GEO|||geodesic-propagation.md"
  "INT|||neutrino-matter-interactions.md"
  "TRP|||trapped-regime-treatment.md"
  "OPA|||opacity-eos-evaluation.md"
  "HYD|||hydro-coupling-source-terms.md"
  "GPU|||particle-container-and-gpu.md"
  "CTX|||carpetx-thorn-integration.md"
  "VER|||verification-suite-design.md"
  "BLD|||build-and-integration.md"
)

spec_by_code() {
  case "$1" in
    CNV) echo "conventions-and-units.md" ;;
    RNG) echo "rng-and-statistical-acceptance.md" ;;
    PKT) echo "packet-representation-and-sampling.md" ;;
    GEO) echo "geodesic-propagation.md" ;;
    INT) echo "neutrino-matter-interactions.md" ;;
    TRP) echo "trapped-regime-treatment.md" ;;
    OPA) echo "opacity-eos-evaluation.md" ;;
    HYD) echo "hydro-coupling-source-terms.md" ;;
    CTX) echo "carpetx-thorn-integration.md" ;;
    GPU) echo "particle-container-and-gpu.md" ;;
    VER) echo "verification-suite-design.md" ;;
    BLD) echo "build-and-integration.md" ;;
    *)   echo "" ;;
  esac
}

# Per-spec fixed-string assertions, "specfile|||needle". The needle (a fixed string,
# grep -F) must appear somewhere in the spec file. Grows by append as each spec lands;
# every appended needle is proven live by a red-run (see header).
SPEC_REQUIRE_IN_SPEC=(
  # README: global contract anchors.
  "README.md|||1e-14"
  "README.md|||1e-12"
  "README.md|||1e-10"
  "README.md|||4σ"
  "README.md|||ABSTOL=RELTOL=1e-12"
  # conventions-and-units: pinned defining constants and conversions.
  "conventions-and-units.md|||2.99792458e10"
  "conventions-and-units.md|||1.380649e-16"
  "conventions-and-units.md|||1.602176634e-19"
  "conventions-and-units.md|||6.67430e-8"
  "conventions-and-units.md|||1.3271244e26"
  "conventions-and-units.md|||1.160451812e10"
  "conventions-and-units.md|||8.617333262e-11"
  "conventions-and-units.md|||g = 4"
  # rng-and-statistical-acceptance: generator pin, KAT anchors, statistical bar.
  "rng-and-statistical-acceptance.md|||Philox4x32-10"
  "rng-and-statistical-acceptance.md|||0xD2511F53"
  "rng-and-statistical-acceptance.md|||0xCD9E8D57"
  "rng-and-statistical-acceptance.md|||0x9E3779B9"
  "rng-and-statistical-acceptance.md|||0xBB67AE85"
  "rng-and-statistical-acceptance.md|||6627e8d5"
  "rng-and-statistical-acceptance.md|||408f276d"
  "rng-and-statistical-acceptance.md|||d16cfe09"
  "rng-and-statistical-acceptance.md|||2^-53"
  "rng-and-statistical-acceptance.md|||4σ"
  "rng-and-statistical-acceptance.md|||1048576"
  "rng-and-statistical-acceptance.md|||1296518744"
  "rng-and-statistical-acceptance.md|||20260730"
  "rng-and-statistical-acceptance.md|||10.1145/2063384.2063405"
  # opacity-eos-evaluation: Kelvin at every temperature argument (incl. NES/Pair),
  # the νx mapping literals, Brem weights, and the g-placement prohibition.
  "opacity-eos-evaluation.md|||× 1.160451812e10"
  "opacity-eos-evaluation.md|||including NES/Pair"
  "opacity-eos-evaluation.md|||κ_a(νx, E) = 0"
  "opacity-eos-evaluation.md|||η(νx, E) = 0"
  "opacity-eos-evaluation.md|||κ_s(νx, E) = ½ [κ_s(νe, E) + κ_s(ν̄e, E)]"
  "opacity-eos-evaluation.md|||28/3"
  "opacity-eos-evaluation.md|||never in per-species transport coefficients"
  # geodesic-propagation: restated-equation anchors (Eqs. 25-26 + null closure).
  "geodesic-propagation.md|||dx^i/dt = γ^{ij} p_j/p^t − β^i"
  "geodesic-propagation.md|||dp_i/dt = −α (∂_i α) p^t + (∂_i β^k) p_k − ½ (∂_i γ^{jk}) p_j p_k/p^t"
  "geodesic-propagation.md|||p^t = √(γ^{ij} p_i p_j)/α"
  # packet-representation-and-sampling: restated-equation anchors (Eqs. 18-19) and
  # the g = 4 exactly-once-at-creation rule.
  "packet-representation-and-sampling.md|||N_p = g_s √−g ΔV Δt η_b(s)/E_p"
  "packet-representation-and-sampling.md|||p^{μ′} = ν (1, sinθ cosφ, sinθ sinφ, cosθ)"
  "packet-representation-and-sampling.md|||g = 4"
  "packet-representation-and-sampling.md|||exactly once, at packet creation"
  # particle-container-and-gpu: the pure-SoA layout pin and the anti-pattern warning.
  "particle-container-and-gpu.md|||ParticleContainerPureSoA"
  "particle-container-and-gpu.md|||are an anti-pattern for GPU particle operators and are not adopted"
  # carpetx-thorn-integration: the hard driver dependency, the observable cadence
  # contract, the single-rate pin, and the same anti-pattern warning (required
  # verbatim in BOTH technical specs).
  "carpetx-thorn-integration.md|||REQUIRES CarpetX"
  "carpetx-thorn-integration.md|||advance exactly once per coarsest-level Δt"
  "carpetx-thorn-integration.md|||CarpetX::use_subcycling = no"
  "carpetx-thorn-integration.md|||are an anti-pattern for GPU particle operators and are not adopted"
  # neutrino-matter-interactions: the restated Eq. 24 anchor (in the corpus's [0,1)
  # draw mapping), the discrete-absorption phrase, and the elastic-outcome phrase.
  "neutrino-matter-interactions.md|||Δt_{s,a} = −ln(1 − u_{s,a}) p^t/(κ_{s,a} ν)"
  "neutrino-matter-interactions.md|||removed entirely (no continuous weight decay)"
  "neutrino-matter-interactions.md|||at fixed fluid-frame energy ν"
  # trapped-regime-treatment: the D8-mandated needles in the binding 2021 α
  # convention, EXACT codepoints (′ U+2032 PRIME, − U+2212 MINUS SIGN — an ASCII
  # lookalike must not satisfy these), plus the flip-warning and limit phrases.
  "trapped-regime-treatment.md|||η′ = αη"
  "trapped-regime-treatment.md|||κ_a′ = ακ_a"
  "trapped-regime-treatment.md|||κ_s′ = κ_s + (1−α)κ_a"
  "trapped-regime-treatment.md|||κ_a′Δt"
  "trapped-regime-treatment.md|||the letter's α is this corpus's 1−α"
  "trapped-regime-treatment.md|||α → 1 reproduces the explicit scheme"
  # hydro-coupling-source-terms: the zero-then-add protocol phrase, the corpus-derived
  # 8-point cell-to-vertex average, and the MeV interface requirement on the partner.
  "hydro-coupling-source-terms.md|||zeroed at one schedule point, accumulated into by contributors, and only then read"
  "hydro-coupling-source-terms.md|||eT(v) += Σ_{d∈{0,1}³} q(v − d)/8"
  "hydro-coupling-source-terms.md|||populates HydroBaseX::temperature in MeV"
)

# Per-spec snapshot-structure assertions, "specfile|||snapshot.h5ls|||needle": the needle
# must appear in the named committed snapshot under the resolved WeakLibInterp root AND
# the spec must cite that snapshot's table by name. Grows by append as specs land.
SPEC_REQUIRE_IN_SNAPSHOT=(
  # opacity-eos-evaluation: the structural claims of the five-family contract and the
  # νx mapping — EOS chemical potentials (Kirchhoff closure), exactly-two electron-type
  # species datasets (EmAb + Iso), the NES/Pair EtaGrid + species-less Kernels, and
  # Brem's single S_sigma dataset with (rho, T) trailing axes.
  "opacity-eos-evaluation.md|||wl-EOS-SFHo-15-25-50.h5ls|||/DependentVariables/Electron\\ Chemical\\ Potential Dataset {30, 81, 185}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-EmAb.h5ls|||/EmAb/Electron\\ Neutrino Dataset {30, 81, 185, 40}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-EmAb.h5ls|||/EmAb/Electron\\ Antineutrino Dataset {30, 81, 185, 40}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-Iso.h5ls|||/Scat_Iso_Kernels/Electron\\ Neutrino Dataset {30, 81, 185, 2, 40}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-NES.h5ls|||/EtaGrid/Values          Dataset {120}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-NES.h5ls|||/Scat_NES_Kernels/Kernels Dataset {120, 81, 4, 40, 40}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-Pair.h5ls|||/Scat_Pair_Kernels/Kernels Dataset {120, 81, 4, 40, 40}"
  "opacity-eos-evaluation.md|||wl-Op-SFHo-15-25-50-E40-Brem.h5ls|||/Scat_Brem_Kernels/S_sigma Dataset {81, 185, 1, 40, 40}"
)

# Repo-prefixed reference-path token pattern (extension-anchored, plus bare README files).
PATH_RE='(amrex|warpx|CarpetX|flesh|WeakLibInterp)/[A-Za-z0-9_.+/-]+\.(H|hxx|cxx|cpp|h|c|ccl|F90|f90|pl|par|md|th|h5ls|provenance)|(amrex|warpx|CarpetX|flesh|WeakLibInterp)/[A-Za-z0-9_+/-]+/README'

# Requirement-id pattern.
ID_RE='\[MCNX-[A-Z]{3}-[0-9]{2}\]'

# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------

# Extract the ordered list of `## ` header titles (trimmed) from a markdown file.
spec_headers() {
  grep -E '^## ' "$1" | sed -E 's/^##[[:space:]]+//; s/[[:space:]]+$//'
}

# Body of one "## " section of one spec (by basename): lines after the exact header up to
# (not including) the next "## " header.
section_text() {
  awk -v h="## $2" '$0==h{f=1;next} /^## /{f=0} f' "$SPECS_DIR/$1"
}

# Check the 7 mandated headers appear in order.
check_section_order() {
  local file="$1" base; base="$(basename "$file")"
  local -a found=()
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

# Check the spec names at least one concrete numeric tolerance.
check_has_tolerance() {
  local file="$1" base; base="$(basename "$file")"
  if grep -Eiq '[0-9](\.[0-9]+)?[eE][-+]?[0-9]+' "$file"; then
    pass "$base: names a concrete numeric tolerance"
  else
    fail "$base: no concrete numeric tolerance found"
  fi
}

# Resolve every repo-prefixed reference path the spec cites against the five roots.
check_cited_paths() {
  local file="$1" base; base="$(basename "$file")"
  local any=0 ok=1
  while IFS= read -r tok; do
    [ -z "$tok" ] && continue
    any=1
    local root rel resolved
    case "$tok" in
      amrex/*)         root="$AMREX_ROOT";          rel="${tok#amrex/}" ;;
      warpx/*)         root="$WARPX_ROOT";          rel="${tok#warpx/}" ;;
      CarpetX/*)       root="$CARPETX_ROOT";        rel="${tok#CarpetX/}" ;;
      flesh/*)         root="$FLESH_ROOT";          rel="${tok#flesh/}" ;;
      WeakLibInterp/*) root="$WEAKLIB_INTERP_ROOT"; rel="${tok#WeakLibInterp/}" ;;
      *) continue ;;
    esac
    resolved="$root/$rel"
    if [ -f "$resolved" ]; then
      info "$base: cited path resolves: $tok"
    else
      fail "$base: cited path does not resolve: $tok (looked at $resolved)"
      ok=0
    fi
  done < <(grep -oE "$PATH_RE" "$file" | sort -u)

  if [ "$any" -eq 0 ]; then
    info "$base: cites no repo-prefixed reference paths"
  elif [ "$ok" -eq 1 ]; then
    pass "$base: all cited reference paths resolve"
  fi
}

# Corpus-wide: each of the five reference roots is cited at least once across the
# on-disk registered specs (five-root closure — no repo goes silently uncited).
check_all_roots_cited() {
  local prefix s hits
  for prefix in amrex warpx CarpetX flesh WeakLibInterp; do
    hits=0
    for s in "${REGISTERED_SPECS[@]}"; do
      [ -f "$SPECS_DIR/$s" ] || continue
      if grep -oE "$PATH_RE" "$SPECS_DIR/$s" | grep -q "^$prefix/"; then hits=1; break; fi
    done
    if [ "$hits" -eq 1 ]; then
      pass "corpus cites at least one $prefix/ reference path"
    else
      fail "corpus cites no $prefix/ reference path in any on-disk spec (five-root closure)"
    fi
  done
}

# Per-spec fixed-string assertions.
check_require_in_spec() {
  local entry="$1" file needle base
  file="${entry%%|||*}"; needle="${entry##*|||}"
  base="$file"; file="$SPECS_DIR/$file"
  if [ -f "$file" ] && grep -Fq -- "$needle" "$file"; then
    pass "$base: contains required claim: $needle"
  else
    fail "$base: missing required claim: $needle"
  fi
}

# Per-spec snapshot-structure assertions against WeakLibInterp's committed snapshots.
check_require_in_snapshot() {
  local entry="$1" file snap needle base table
  file="${entry%%|||*}"
  local rest="${entry#*|||}"
  snap="${rest%%|||*}"
  needle="${rest##*|||}"
  base="$file"
  table="${snap%.h5ls}.h5"

  local snapfile="$FIXTURES_DIR/$snap"
  if [ ! -f "$snapfile" ]; then
    fail "$base: WeakLibInterp snapshot not found under resolved root: specs/fixtures/$snap (root: $WEAKLIB_INTERP_ROOT)"
    return
  fi
  if ! grep -Fq -- "$needle" "$snapfile"; then
    fail "$base: snapshot $snap does not contain documented structure: $needle"
    return
  fi
  if ! grep -Fq -- "$table" "$SPECS_DIR/$file" 2>/dev/null; then
    fail "$base: spec does not cite its reference table by name: $table"
    return
  fi
  pass "$base: documented structure '$needle' confirmed in snapshot $snap"
}

# --------------------------------------------------------------------------------------
# Citation-registry closure (per coded spec):
#   every inline-cited bare arXiv id must have exactly one version-pinned row in the
#   spec's "## Source of truth" registry; the README master registry must pin the same
#   id; and the two pinned versions must agree. Non-arXiv sources close on their DOI.
# --------------------------------------------------------------------------------------
check_citation_registry() {
  local base="$1" f="$SPECS_DIR/$1"
  local ids id esc spec_ver readme_ver nspec nreadme ok=1 any=0

  ids="$(grep -oE 'arXiv:[0-9]{4}\.[0-9]{4,5}' "$f" | sort -u || true)"
  if [ -n "$ids" ]; then
    while IFS= read -r id; do
      [ -n "$id" ] || continue
      any=1
      esc="$(printf '%s' "$id" | sed 's/\./\\./g')"
      spec_ver="$(section_text "$base" "Source of truth" | grep -oE "${esc}v[0-9]+" | sort -u || true)"
      nspec="$(printf '%s' "$spec_ver" | grep -c . || true)"
      if [ "$nspec" -eq 0 ]; then
        fail "$base: cites $id but its Source of truth registry has no version-pinned row (${id}vN)"
        ok=0; continue
      elif [ "$nspec" -gt 1 ]; then
        fail "$base: Source of truth registry pins $id at multiple versions: $(echo $spec_ver)"
        ok=0; continue
      fi
      readme_ver="$(grep -oE "${esc}v[0-9]+" "$README" 2>/dev/null | sort -u || true)"
      nreadme="$(printf '%s' "$readme_ver" | grep -c . || true)"
      if [ "$nreadme" -eq 0 ]; then
        fail "$base: $id is not version-pinned in the README master registry"
        ok=0; continue
      elif [ "$nreadme" -gt 1 ]; then
        fail "README master registry pins $id at multiple versions: $(echo $readme_ver)"
        ok=0; continue
      fi
      if [ "$spec_ver" != "$readme_ver" ]; then
        fail "$base: $id pinned as $spec_ver in spec registry but $readme_ver in README master registry (version disagreement)"
        ok=0
      fi
    done <<< "$ids"
  fi

  local dois doi
  dois="$(grep -oE '10\.[0-9]{4,9}/[A-Za-z0-9._/-]+' "$f" | sort -u || true)"
  if [ -n "$dois" ]; then
    while IFS= read -r doi; do
      [ -n "$doi" ] || continue
      any=1
      if ! section_text "$base" "Source of truth" | grep -Fq -- "$doi"; then
        fail "$base: cites DOI $doi but its Source of truth registry has no matching row"
        ok=0
      fi
      if ! grep -Fq -- "$doi" "$README" 2>/dev/null; then
        fail "$base: DOI $doi is not in the README master registry"
        ok=0
      fi
    done <<< "$dois"
  fi

  if [ "$any" -eq 0 ]; then
    info "$base: cites no published sources"
  elif [ "$ok" -eq 1 ]; then
    pass "$base: citation registry closed (inline ids ⊂ spec registry ⊂ README master registry, versions agree)"
  fi
}

# --------------------------------------------------------------------------------------
# Requirement-id closure:
#   every [MCNX-XXX-NN] occurrence in any coded spec must (1) carry a known spec code,
#   (2) be declared in the owning spec's "## Correctness requirements" section
#   (code-ownership rule — an id's prefix always names its owner), and (3) appear in
#   verification-suite-design.md's coverage matrix (coverage closure). A matrix id with
#   no declaration fails (2), so phantom matrix rows are caught by the same loop.
# --------------------------------------------------------------------------------------
check_id_closure() {
  local entry f all_ids id code owner
  all_ids="$( { for entry in "${CODED_SPECS[@]}"; do
                  f="$SPECS_DIR/${entry##*|||}"
                  [ -f "$f" ] && grep -hoE "$ID_RE" "$f"
                done; } | sort -u || true)"

  if [ -z "$all_ids" ]; then
    info "no [MCNX-XXX-NN] requirement ids found in on-disk coded specs"
    return
  fi

  while IFS= read -r id; do
    [ -n "$id" ] || continue
    code="$(printf '%s' "$id" | sed -E 's/^\[MCNX-([A-Z]{3}).*/\1/')"
    owner="$(spec_by_code "$code")"
    if [ -z "$owner" ]; then
      fail "$id: unknown spec code '$code' (no owning spec)"
      continue
    fi
    if [ ! -f "$SPECS_DIR/$owner" ]; then
      fail "$id: owning spec $owner is missing on disk"
    elif section_text "$owner" "Correctness requirements" | grep -Fq -- "$id"; then
      pass "$id: declared in ## Correctness requirements of $owner"
    else
      fail "$id: not declared in ## Correctness requirements of $owner (code-ownership rule)"
    fi
    if grep -Fq -- "$id" "$SPECS_DIR/verification-suite-design.md" 2>/dev/null; then
      info "$id: covered by verification-suite-design.md"
    else
      fail "coverage: $id not referenced in verification-suite-design.md coverage matrix"
    fi
  done <<< "$all_ids"
}

# --------------------------------------------------------------------------------------
# README index integrity
# --------------------------------------------------------------------------------------
check_readme_links() {
  if [ ! -f "$README" ]; then fail "README.md missing at specs/README.md"; return; fi

  local -a linked=()
  while IFS= read -r target; do
    target="${target#./}"
    linked+=("$target")
  done < <(grep -oE '\]\(\.?/?[A-Za-z0-9_-]+\.md\)' "$README" | sed -E 's/^\]\(\.?\/?//; s/\)$//' | sort -u)

  # (i1) every registered spec is linked (README itself is the index, not a link target).
  local s
  for s in "${REGISTERED_SPECS[@]}"; do
    [ "$s" = "README.md" ] && continue
    local hit=0 l
    for l in "${linked[@]}"; do [ "$l" = "$s" ] && hit=1; done
    if [ "$hit" -eq 1 ]; then
      info "README links registered spec: $s"
    else
      fail "README index does not link registered spec: $s"
    fi
  done

  # (i2) every link target exists on disk (no broken intra-repo links).
  local l
  for l in "${linked[@]}"; do
    if [ ! -f "$SPECS_DIR/$l" ]; then
      fail "README links a nonexistent file: $l"
    fi
  done

  # (i3) every registered spec must be on disk.
  for s in "${REGISTERED_SPECS[@]}"; do
    if [ ! -f "$SPECS_DIR/$s" ]; then
      fail "registered spec file is missing on disk: $s"
    fi
  done

  if [ "$FAILURES" -eq 0 ]; then pass "README index: all registered specs linked, all links resolve, no orphans"; fi
}

# README links EXACTLY the registered spec set: on-disk .md count equals the registry
# size (13, including README.md itself), no orphans, no dangling links.
check_readme_exact_registered() {
  if [ ! -f "$README" ]; then fail "README.md missing at specs/README.md"; return; fi

  local -a on_disk=()
  while IFS= read -r f; do
    on_disk+=("$(basename "$f")")
  done < <(find "$SPECS_DIR" -maxdepth 1 -name '*.md' | sort)

  local -a linked=()
  while IFS= read -r target; do
    target="${target#./}"
    linked+=("$target")
  done < <(grep -oE '\]\(\.?/?[A-Za-z0-9_-]+\.md\)' "$README" | sed -E 's/^\]\(\.?\/?//; s/\)$//' | sort -u)

  local want="${#REGISTERED_SPECS[@]}"
  if [ "${#on_disk[@]}" -eq "$want" ]; then
    pass "exactly $want spec files present on disk (including README.md)"
  else
    local missing="" extra="" s f hit
    for s in "${REGISTERED_SPECS[@]}"; do
      hit=0
      for f in "${on_disk[@]}"; do [ "$f" = "$s" ] && hit=1; done
      [ "$hit" -eq 1 ] || missing="$missing $s"
    done
    for f in "${on_disk[@]}"; do
      hit=0
      for s in "${REGISTERED_SPECS[@]}"; do [ "$s" = "$f" ] && hit=1; done
      [ "$hit" -eq 1 ] || extra="$extra $f"
    done
    fail "expected exactly $want spec files on disk (including README.md), found ${#on_disk[@]}; missing:${missing:- none}; unregistered:${extra:- none}"
  fi

  # every on-disk spec (other than README) is linked from the README (no orphans).
  local f l hit
  for f in "${on_disk[@]}"; do
    [ "$f" = "README.md" ] && continue
    hit=0
    for l in "${linked[@]}"; do [ "$l" = "$f" ] && hit=1; done
    [ "$hit" -eq 1 ] || fail "orphan spec file on disk not linked from README: $f"
  done

  # every README spec link points at an on-disk file (no dangling links).
  for l in "${linked[@]}"; do
    hit=0
    for f in "${on_disk[@]}"; do [ "$f" = "$l" ] && hit=1; done
    [ "$hit" -eq 1 ] || fail "README links a spec not present on disk: $l"
  done

  if [ "$FAILURES" -eq 0 ]; then pass "README links exactly the $want registered on-disk spec files (no orphans, no missing)"; fi
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
echo "WeakLibInterp root: $WEAKLIB_INTERP_ROOT"
echo

echo "--- per-spec: section order, tolerance, cited paths ---"
for s in "${REGISTERED_SPECS[@]}"; do
  f="$SPECS_DIR/$s"
  if [ ! -f "$f" ]; then fail "registered spec not found: $s"; continue; fi
  if [ "$s" != "README.md" ]; then check_section_order "$f"; fi
  check_has_tolerance "$f"
  check_cited_paths "$f"
done

echo
echo "--- corpus: five-root citation closure ---"
check_all_roots_cited

echo
echo "--- per-spec: registered fixed-string claims ---"
for e in "${SPEC_REQUIRE_IN_SPEC[@]}"; do check_require_in_spec "$e"; done

echo
echo "--- per-spec: documented structure vs WeakLibInterp committed snapshots ---"
if [ "${#SPEC_REQUIRE_IN_SNAPSHOT[@]}" -gt 0 ]; then
  for e in "${SPEC_REQUIRE_IN_SNAPSHOT[@]}"; do check_require_in_snapshot "$e"; done
else
  info "no snapshot-structure assertions registered"
fi

echo
echo "--- per-spec: citation-registry closure ---"
for entry in "${CODED_SPECS[@]}"; do
  s="${entry##*|||}"
  [ -f "$SPECS_DIR/$s" ] || continue   # absence already reported above
  check_citation_registry "$s"
done

echo
echo "--- closure: requirement-id ownership, declaration, coverage ---"
check_id_closure

echo
echo "--- README index integrity ---"
check_readme_links

echo
echo "--- closure: README links exactly the registered spec files ---"
check_readme_exact_registered

echo
if [ "$FAILURES" -eq 0 ]; then
  echo "ALL CHECKS PASSED ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
  exit 0
else
  echo "VALIDATION FAILED: $FAILURES check(s) failed"
  exit 1
fi
