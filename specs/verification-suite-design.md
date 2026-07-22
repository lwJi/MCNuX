# Verification-suite design (technical, the runnable coverage matrix)

> Technical spec. It binds every MCNuX spec's "Verification" section into one coherent, runnable suite and fixes the pass/fail discipline, but defers the correctness arithmetic to the owning specs and the shared contracts. It references `conventions-and-units.md` (notation/units) and `rng-and-statistical-acceptance.md` (seeds, counts, the 4σ recipe) and restates only what the suite as a whole must guarantee. `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope

This spec defines the verification suite as a *contract*: the harness formats (Cactus regression tests with golden TSV/norm data, standalone unit checks), the design by which statistical N-sigma acceptance reduces to numeric golden data, the pass/fail discipline, and the corpus-wide **coverage matrix** in which every correctness-requirement id declared by any MCNuX spec must appear. It is the single place a fresh agent learns **what the Ralph loop gates on**, so no spec's verification section is an isolated island.

In scope:

- The two check families and their formats: (A) Cactus regression tests, (B) standalone unit checks.
- The statistical-margin design: N-sigma checks computed by an analysis routine, gated at runtime, and pinned as golden norm data.
- The coverage matrix: columns, the requirement-id scheme, and the closure rule (grown as the corpus grows; complete over whatever specs exist on disk).
- The pass/fail discipline: every check asserts and fails the suite on violation; no print-only checks.

Out of scope:

- The physics and numerics each requirement asserts (owned by the spec that declares the requirement id; this file only binds each id to a benchmark, tier, and pass/fail number).
- The tolerance-tier definitions (inherited from `README.md`) and the RNG/seed/count contract (see `rng-and-statistical-acceptance.md`).
- The test framework, directory layout, and assertion library — explicitly implementation freedom (see below).
- Implementing the suite (this corpus specifies; the Ralph loop implements).

## Source of truth

- `flesh/lib/sbin/RunTestUtils.pl` — the Cactus test harness: a test is `<thorn>/test/<name>.par` with a same-named sibling directory of golden output files; discovery requires the directory to exist; comparison is a numeric line-by-line diff with default `ABSTOL = RELTOL = 1e-12`, overridable per-thorn/per-test via `test.ccl`.
- `CarpetX/CarpetX/src/io_tsv.cxx` — the TSV line-cut output (`CarpetX::out_tsv_vars`) used as golden reference data; one file per variable group per axis with header `# 1:iteration 2:time 3:patch 4:level 5:i 6:j 7:k 8:x 9:y 10:z 11:<var>`.
- `CarpetX/CarpetX/src/io_norm.cxx` — the norm-table output (`CarpetX::out_norm_vars`): per-group `norms/<name>.tsv` files carrying min/max/sum/avg/stddev/L1norm/L2norm/maxabs per iteration.
- `CarpetX/WaveToyX/test/standing.par` (and its sibling golden directory, including `norms/`) — the working in-repo example of exactly this test shape.
- `WeakLibInterp/specs/regression-suite-design.md` — pattern provenance for the coverage-matrix/closure design (a proven Ralph-loop suite contract in this ecosystem); nothing in it is normative for MCNuX.

## Inputs & outputs

This spec does not define a callable surface. Its "inputs" are the correctness-requirement ids declared by the MCNuX specs; its "output" is a binary pass/fail per check and an aggregate pass/fail for the suite.

### Check family A — Cactus regression tests (golden data)

- A test is a parameter file `test/<name>.par` in the MCNuX thorn plus a same-named directory `test/<name>/` of committed golden output, diffed numerically by the Cactus harness at `ABSTOL = RELTOL = 1e-12` (per `flesh/lib/sbin/RunTestUtils.pl`); a test may override tolerances via `test.ccl` only with a rationale recorded in this spec's matrix row.
- Golden reference outputs are CarpetX TSV line-cuts (`CarpetX::out_tsv_vars`) and norm tables (`CarpetX::out_norm_vars`) of MCNuX grid variables and per-packet diagnostic output.
- Deterministic regression tests run in **single-rank CPU mode at pinned seeds**, where per-packet event sequences are bitwise reproducible (requirement MCNX-RNG-03 in `rng-and-statistical-acceptance.md`); grid-tally goldens from parallel configurations are compared as norms at the relaxed tier with the rationale stated in the matrix row.

### Check family B — standalone unit checks

- Host-runnable checks with no Cactus executable and no grid: constants, species metadata, RNG key-tuple sweeps, sampling-kernel moments, single-packet pushes on analytic metrics. One executable (or target) per concern, runnable in CI; each asserts its tier and exits nonzero on violation.
- Unit checks are the natural home for machine-tier and exact-tier cells; they must not silently skip — a check that cannot run in an environment reports a distinct SKIP state, never a pass.

### The statistical-margin design (how 4σ reduces to golden data)

Statistical acceptance cells are computed by an **analysis routine** (a scheduled MCNuX analysis-bin routine or a unit check, per the matrix row) that reduces the assertion to numbers:

1. The routine computes the estimator Q̂, the target Q, the standard error s (per the recipe in `rng-and-statistical-acceptance.md`), and the margin `m = 4·s − |Q̂ − Q|`, at the pinned seed and packet count.
2. The routine **hard-fails at runtime** (aborts the run with an error) if `m < 0` — the acceptance criterion gates even if golden data were ever regenerated from a failing run.
3. The routine writes (Q̂, Q, s, m) as grid-scalar/norm output, which is committed as golden data; the harness diff then also catches silent drift of the statistical outcome at the pinned seed (a regression bound far tighter than the 4σ acceptance bound).

So every statistical cell is simultaneously an *acceptance* check (sign of m, runtime-gated) and a *regression* check (golden values at pinned seed, harness-gated).

### The requirement-id scheme and the coverage matrix

Every MCNuX spec declares its hard correctness requirements with ids of the form `MCNX-<TAG>-<NN>` (tag unique per spec, two-digit sequence). The coverage matrix below binds each id to a benchmark, a check family, a tolerance tier, and a concrete pass/fail criterion. **Closure rule:** every id declared in any spec on disk appears in this matrix, and every id in this matrix is declared by a spec on disk; `specs/tools/validate_specs.sh` enforces both directions mechanically, so the matrix is complete at every phase of the corpus's growth.

#### Coverage matrix

| Requirement id | Owning spec | Benchmark / check | Family | Tier | Pass/fail criterion |
|---|---|---|---|---|---|
| MCNX-CNV-01 | [conventions-and-units](./conventions-and-units.md) | recompute every conversion factor from the five defining constants; compare against implementation values and against the spec's quoted decimals | B (unit) | machine | derived factors agree at `rtol = 1e-14`, `atol = 1e-30`; quoted decimals at `rtol = 1e-6` |
| MCNX-CNV-02 | [conventions-and-units](./conventions-and-units.md) | species metadata table; synthetic νx emission/tally round trip for the single-g-factor rule | B (unit) | exact + relaxed | indices/g/ℓ integer-equal; νx tally carries exactly one factor g = 4, round trip at `rtol = 1e-10` |
| MCNX-RNG-01 | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | ≥ `1e3` key tuples evaluated serial/reverse/shuffled/threaded/GPU/2-rank | B (unit) | exact | bitwise-identical doubles across all execution orders and devices |
| MCNX-RNG-02 | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | uniform-moment battery (mean, mean of u², lag-1 correlation, cross-stream correlation) at `N = 1e6` draws, seeds `S0/S1/S2` | B (unit) | 4σ | margin `m ≥ 0` for all four statistics at all three pinned seeds |
| MCNX-RNG-03 | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | fixed-seed transport problem run twice, single-rank CPU; per-packet TSV vs. committed golden | A (regression) | golden | zero violations at `ABSTOL = RELTOL = 1e-12` |
| MCNX-RNG-04 | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | same problem ≥ 3 threaded/GPU runs; pairwise grid-tally norm comparison | A (regression, norms) | relaxed | all pairwise norm differences ≤ `1e-10` relative |
| MCNX-RNG-05 | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | enumerate all packet ids in a 2-rank run with mid-run checkpoint/restart | B (unit) / A | exact | zero duplicate ids |

(The matrix grows by appended rows as domain and technical specs join the corpus; the closure rule above keeps it complete at every phase boundary.)

### Tolerance tiers (referenced, not redefined)

The tiers — machine `~1e-14`, golden/parity `1e-12` (the Cactus harness default `ABSTOL = RELTOL = 1e-12`), relaxed `1e-10`, exact/no-tolerance, and the 4σ statistical bar — and the mixed comparison form `|got − expected| <= rtol·|expected| + atol` (`atol = 1e-30` default) are defined in `README.md`; seeds and packet counts in `rng-and-statistical-acceptance.md`. Each matrix cell asserts at the tier its owning spec dictates.

## Correctness requirements

- **Closure (both directions).** Every requirement id declared by any spec on disk has exactly one matrix row here; every matrix row's id is declared by a spec on disk. Enforced mechanically by `specs/tools/validate_specs.sh`; a spec adding a requirement id without a matrix row fails the gate, as does a stale row.
- **Every row is concrete.** A matrix row names its benchmark, its check family, its tier, and a numeric (or exact/integer) pass/fail criterion. "Verified by inspection" is not a valid row.
- **Statistical rows are pinned and double-gated.** Every 4σ row runs at the pinned seed(s) and count(s) of `rng-and-statistical-acceptance.md`, hard-fails at runtime on `m < 0`, and commits (Q̂, Q, s, m) as golden data.
- **Assertions gate; nothing print-only.** Every check in both families asserts its criterion and propagates failure to the suite's aggregate exit status (a failing Cactus test fails the harness run; a failing unit check exits nonzero). A check that cannot run reports SKIP distinctly from pass.
- **Perturbation meta-check.** The suite must include at least one deliberate-perturbation check per family: inject a known error (e.g. scale one expected value by 10, corrupt one golden number) and confirm the corresponding check **fails** — proving the assertions are thresholded and wired, not decorative.
- **Golden data provenance.** Every committed golden directory records, in a header or sibling provenance note, the seed, packet count, and configuration that produced it; regenerating golden data for a statistical row is valid only from a run whose runtime margin gate passed.

## Verification

A fresh agent confirms the suite design itself is honored by these checks:

1. **Mechanical closure.** `bash specs/tools/validate_specs.sh` passes: this file carries the 7 mandated sections in order; every requirement id declared across the on-disk specs appears in the matrix above and vice versa; cited `flesh/`, `CarpetX/`, and `WeakLibInterp/` paths resolve.
2. **Matrix realization (once implementation exists).** For every matrix row there is a runnable check whose name or registration traceably references the requirement id, and running the suite executes every row's check (no permanently-skipped rows in the reference CI environment).
3. **Perturbation.** Executing the perturbation meta-checks demonstrates a failing exit status in each family.
4. **Harness conformance.** A regression test added per family-A rules is discovered and diffed by the Cactus harness exactly as `CarpetX/WaveToyX/test/standing.par` is (same file shape, same tolerance mechanics).

## Implementation freedom

- **The test framework, assertion library, harness glue, and directory layout** for family-B unit checks (CTest, a hand-rolled driver, GoogleTest, …), and how checks map to executables or targets.
- The transport problems used as fixed-seed regression fixtures (grid sizes, packet counts within the pinned options, which variables are output) — provided each matrix row's benchmark description is realized and seeds/counts are the pinned ones.
- How SKIP is reported, how golden provenance is recorded, and how the perturbation meta-checks are implemented — provided the observable guarantees above hold.
- Whether family-B checks live in the thorn source tree, a `test/` subtree, or a separate utility target.

## Open questions / assumptions

- **The matrix grows with the corpus (by design, durable).** At any moment the matrix covers exactly the requirement ids of the specs on disk; rows for future specs (packet sampling, geodesics, interactions, trapped regime, coupling, container, thorn integration, build) are appended when those specs land. The closure rule makes "complete" a checkable property at every phase, so this growth does not weaken the contract.
- **GPU rows depend on a device build existing (assumption, non-blocking).** Rows that mention GPU execution (MCNX-RNG-01's device sweep, MCNX-RNG-04's GPU tallies) run in CPU-only form where no device build is available and report the device leg as SKIP — loudly, never silently.
- **Norm-based comparison granularity (assumption).** Norm tables aggregate per group per iteration; where a future row needs spatial localization beyond norms, it uses TSV line-cuts instead. No new output format is introduced without a spec change here.
