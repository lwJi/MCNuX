# RNG and statistical acceptance (cross-cutting stochastic-correctness contract)

> Cross-cutting spec. Every MCNuX spec that samples anything (emission, event distances, scattering angles) or asserts a Monte Carlo outcome references this file for the random-number contract, the reproducibility policy, and the 4σ acceptance-band recipe, and restates only the draws it consumes. If a leaf spec and this file ever disagree, `README.md` is the canonical arbiter; absent that, this file governs stochastic correctness.

## Purpose & scope

This spec defines the stochastic-correctness contract: the stateless counter-based random-number generator keyed on (seed, packet id, event counter), the reproducibility guarantees (bitwise per-packet event sequences in single-rank CPU test mode; floating-point-reduction tolerance for grid tallies), the 4σ statistical acceptance-band recipe, and the pinned packet counts and seeds every statistical check in the corpus uses.

In scope:

- The generator contract: key structure, output range, order-independence, device-independence.
- Packet-id uniqueness and event-counter semantics (the observable RNG identity of a packet).
- The reproducibility policy and which comparisons it licenses (golden diffs vs. norm bounds).
- The 4σ acceptance-band recipe: standard-error estimation, pass condition, false-failure budget, pinned-seed discipline.
- The pinned packet counts and pinned seeds (normative for the whole corpus).

Out of scope:

- *Which* random draws each physics process consumes and what distributions they feed (owned by the packet-representation and interaction specs, which cite this file for the draw primitive).
- The verification harness mechanics that turn these bounds into runnable tests (see `verification-suite-design.md`).
- Units, species, and notation (see `conventions-and-units.md`).

## Source of truth

The RNG contract is **derived and pinned in this spec** (no single external reference suffices); provenance and reference artifacts:

- Counter-based ("stateless") parallel RNG design: Salmon, Moraes, Dror & Shaw 2011, *Parallel random numbers: as easy as 1, 2, 3* (SC'11, doi:10.1145/2063384.2063405) — the Random123 family (Philox, Threefry) demonstrating crush-resistant pure-function generators keyed on (counter, key). This spec pins the *contract* such generators satisfy, not a specific family.
- `amrex/Src/Base/AMReX_Random.H` — AMReX's default per-thread **stateful** generators. Cited as the counter-example: a per-thread stateful stream's output depends on execution order and thread assignment and therefore does **not** satisfy this contract; an implementation may not use it for any draw that affects physics outcomes.
- `flesh/lib/sbin/RunTestUtils.pl` — the Cactus test harness numeric comparison (default `ABSTOL = RELTOL = 1e-12`) that the single-rank CPU golden-output guarantee is designed to satisfy.
- `CarpetX/CarpetX/src/io_norm.cxx` — the norm-table output used for reduction-tolerance and statistical-margin comparisons.

## Inputs & outputs

### The generator contract

The generator is a pure function

```text
R : (S, q, e, k) → u ∈ [0, 1)
```

- `S` — 64-bit unsigned global seed (a run parameter).
- `q` — 64-bit unsigned packet id, unique per packet within a run (see below).
- `e` — 64-bit unsigned event counter, scoped to the packet: the index of the RNG-consuming event in that packet's life (creation is event 0; each subsequent event that consumes randomness increments it).
- `k` — 32-bit unsigned draw index within the event (0, 1, 2, … for the first, second, third draw the event consumes).
- `u` — an IEEE-754 `double` in [0, 1), uniformly distributed, with at least 53 significant bits drawn from the generator's output word(s) and an exact (error-free) integer-to-double scaling, so the same key tuple yields the same `double` on every platform.

Draws for other distributions (isotropic angles, exponential optical depths, …) are constructed from these uniform draws; the construction is implementation freedom, but each physics spec states *how many* uniform draws its events consume and in what order (`k` assignment), because that ordering is part of the reproducible event sequence.

### Packet id and event counter (the RNG identity)

- Every packet carries its RNG identity `(q, e)` as part of its state (nomenclature reserved in `conventions-and-units.md`).
- `q` is assigned at packet creation and is **unique for the lifetime of the run**: no two packets ever created in a run share a `q`, including across ranks, levels, restarts from checkpoint, and packets created by splitting/russian-roulette. Uniqueness is exact (integer), not statistical.
- `e` starts at 0 at creation and increments by exactly 1 per RNG-consuming event; it never repeats or rewinds for a given packet (checkpoint/restart must preserve `(q, e)`).

### Pinned seeds and packet counts (normative for every statistical check)

| Name | Value | Use |
|---|---|---|
| Primary seed `S0` | `20260721` | every statistical check's default seed |
| Alternate seed `S1` | `42` | second-seed confirmation runs |
| Alternate seed `S2` | `271828` | third-seed confirmation runs |
| Smoke count `N_smoke` | `1e4` packets | quick statistical smoke checks |
| Standard count `N_standard` | `1e6` packets (or `1e6` draws for generator-level checks) | the gating statistical acceptance runs |

A statistical check in any spec runs at `(S0, N_standard)` unless that spec explicitly states otherwise; a check may *add* runs at `S1`/`S2` or `N_smoke`, never substitute a different unpinned seed or count.

### The 4σ acceptance-band recipe

For an estimator Q̂ of a target value Q with standard error s:

1. Obtain s as (in order of preference): (a) the analytically known standard error of the estimator; else (b) `s = σ̂/√N`, with σ̂ the sample standard deviation of the per-packet (or per-draw) contributions; else (c) a batch-means estimate with ≥ 16 batches.
2. The check **passes** iff `|Q̂ − Q| ≤ 4·s` (two-sided).
3. If `s = 0` (the outcome is deterministic at the pinned seed), the check degrades to the deterministic tier its owning spec states (golden `1e-12` or exact).
4. Every statistical check reports (Q̂, Q, s) and the margin `m = 4·s − |Q̂ − Q|`; the pass condition is `m ≥ 0`.

False-failure budget: a two-sided Gaussian 4σ band has per-check false-failure probability `6.33e-5`; a suite of M ≤ 200 statistical cells therefore expects ≤ 0.013 spurious failures, and because seeds are pinned, a pass or failure is *reproducible*, not merely probable.

**Pinned-seed discipline.** Seeds are pinned before results are computed. If a pinned-seed check fails its 4σ band, the response is investigation (of the implementation, the target value, or the standard-error estimate) — never replacing the seed, widening the band, or rerunning until it passes. Changing a pinned seed or count is a spec change to this file.

## Correctness requirements

- **[MCNX-RNG-01] Order- and device-independence (exact).** `R(S, q, e, k)` returns the identical `double` for a given key tuple regardless of the order in which draws are evaluated, the number of MPI ranks and how packets are distributed over them, tiling/iteration order, and CPU vs. GPU execution. Consequently a packet's entire random stream is a pure function of `(S, q)` and its event history. Per-thread stateful generation (the `amrex/Src/Base/AMReX_Random.H` default) violates this requirement and is excluded from all physics-affecting draws. Pass criterion: bitwise equality of draws for identical key tuples across all execution configurations tested.
- **[MCNX-RNG-02] Statistical quality at the pinned seed (4σ).** With seed `S0 = 20260721` and `N_standard = 1e6` consecutive draws (fixed key-tuple enumeration: q = 0, e = 0, k = 0…N−1), each of the following sample statistics lies within its 4σ band:
  - sample mean of u vs. 1/2, standard error `s = 1/√(12N)`;
  - sample mean of u² vs. 1/3, standard error `s = √(4/45)/√N` (Var(u²) = 1/5 − 1/9 = 4/45);
  - lag-1 serial correlation vs. 0, standard error `s = 1/√N`;
  - cross-stream correlation between streams `q = 0` and `q = 1` (paired by k) vs. 0, standard error `s = 1/√N`.
  The same four statistics must also pass at seeds `S1` and `S2`. This is a necessary smoke bar, not a substitute for choosing a generator family with published TestU01/BigCrush results (see Implementation freedom).
- **[MCNX-RNG-03] Single-rank CPU bitwise reproducibility (golden 1e-12).** In single-rank CPU test mode, two runs of an identical configuration (same parameter file, same seed) produce bitwise-identical per-packet event sequences — for every packet: the ordered list of (event counter, event type, draws consumed, resulting state) — and therefore identical per-packet TSV output. Pass criterion: the Cactus harness numeric diff of per-packet output at its default `ABSTOL = RELTOL = 1e-12` shows zero violations (the outputs are in fact bitwise identical; `1e-12` is the harness bar the golden scheme rides on).
- **[MCNX-RNG-04] Grid-tally reproducibility bound (relaxed 1e-10).** Repeated identical runs in any parallel configuration (multi-thread, GPU atomics, multi-rank) produce grid-tally fields whose norms agree within relative `1e-10`. Rationale for the relaxed tier: deposition order is unordered (atomic adds, MPI reduction order), so tallies differ by floating-point reassociation only; with ≤ `1e4` contributions per cell in the pinned test configurations the reassociation error is bounded by a few times `n·ulp ≈ 1e4 · 2.2e-16 ≈ 2e-12` relative, and `1e-10` carries two orders of margin. No bitwise cross-run requirement exists for tallies; any check needing bitwise stability must use per-packet outputs under MCNX-RNG-03.
- **[MCNX-RNG-05] Uniqueness of q (exact).** Packet-id uniqueness as defined above is an invariant: a run that ever assigns a duplicate `q` is incorrect regardless of downstream statistics.
- **Statistical checks corpus-wide.** Every statistical acceptance in every MCNuX spec uses the recipe, seeds, and counts of this file; a leaf spec may tighten (more packets, more seeds) but not loosen, and may not introduce unpinned randomness into a gating check.

## Verification

- **MCNX-RNG-01**: a unit check evaluates a fixed set of ≥ `1e3` key tuples (spanning small/large S, q, e, k, including k = 0 and the maximum draw index used by any event) serially, in reverse, shuffled, multi-threaded, and (where a device build exists) on GPU, and asserts bitwise-identical `double` outputs across all orders and devices; plus a two-rank MPI run asserting the same tuples give the same values on both ranks. Exact equality, no tolerance.
- **MCNX-RNG-02**: a unit check draws the pinned sequences and asserts the four moment statistics inside their 4σ bands at `S0`, `S1`, `S2`; it reports each margin `m`.
- **MCNX-RNG-03**: a Cactus regression test (see `verification-suite-design.md`) runs a small fixed-seed transport problem twice in single-rank CPU mode and diffs per-packet TSV output against committed golden data at `ABSTOL = RELTOL = 1e-12`.
- **MCNX-RNG-04**: the same fixed-seed problem run ≥ 3 times in a threaded (and, where available, GPU) configuration; grid-tally norms compared pairwise at relative `1e-10`.
- **MCNX-RNG-05**: a check enumerating all packet ids created in a multi-rank run (including a mid-run checkpoint/restart) and asserting zero duplicates.
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `amrex/`, `flesh/`, and `CarpetX/` paths, and contains the required claim strings (counter-based keying, the pinned seed, the 4σ recipe).

## Implementation freedom

- **The generator family** — Philox, Threefry, a keyed PCG variant, or any other pure-function generator — provided it satisfies MCNX-RNG-01/02 and has published results from a recognized statistical battery (TestU01 Crush/BigCrush or equivalent) for the chosen parameterization.
- How the four key components are packed into the generator's counter/key words, and how output words are combined into the [0, 1) `double` — provided the scaling is exact and platform-independent.
- How packet ids are generated (rank-prefixed counters, reserved-range allocation, AMReX particle id machinery) — provided uniqueness-for-the-run holds exactly, including across restarts.
- Transformations from uniform draws to other distributions (inversion, rejection with a stated draw budget, Box–Muller, …) — provided the consuming spec's distributional requirements and draw-order statement are met.
- Caching, batching, or vectorizing draw evaluation — provided observable values equal the pure-function definition.

## Open questions / assumptions

- **Draw budget per event (assumption, non-blocking).** The 32-bit draw index k assumes no single event consumes more than `2^32` uniform draws; all currently specified event types consume ≤ 8. If a future rejection-sampling scheme could exceed a bounded draw count, its spec must state the bound and the fallback (e.g. advancing e).
- **Counter-space exhaustion (assumption, non-blocking).** 64-bit q and e cannot plausibly wrap at target packet counts (≤ `1e9` packets, ≤ `1e6` events per packet); no wrap-handling behavior is specified.
- **GPU bitwise identity across vendors (assumption).** MCNX-RNG-01's device-independence rests on integer arithmetic and exact int-to-double scaling, which IEEE-conforming CPU and GPU backends both guarantee; exotic backends without conforming 64-bit integer or double support are out of scope.
