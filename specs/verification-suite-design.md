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
| MCNX-PKT-01 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | fixed-seed emission fixture: every created packet validated against the state contract (types, ranges, species, p_t consistency, six-draw creation budget via the logged event sequence) | B (unit) | exact + machine | integer fields exact; zero range/finiteness violations; p_t identity at `rtol = 1e-14`; draws k = 0…5 of event e = 0 in order |
| MCNX-PKT-02 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | synthetic-η emission fixture: per (cell, species, step), Σ N_k ε_k vs E_tot = α √γ ΔV Δt η_s (νx includes g = 4 once) | B (unit) | relaxed | agreement at `rtol = 1e-10` for every (cell, species, step) |
| MCNX-PKT-03 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | same fixture: equal fluid-frame energy within (cell, species, step); N_k = E_p / ε_k | B (unit) | machine | pairwise N_k ε_k equality at `rtol = 1e-14` |
| MCNX-PKT-04 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | ≥ 8-group synthetic spectrum spanning ≥ 3 decades; group occupation f_b vs η_{s,b}/η_s at (S0, N_standard) | B (unit) | 4σ | margin m ≥ 0 for every group, s_b = √(P(b)(1−P(b))/N) |
| MCNX-PKT-05 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | moving-fluid cell (v = (0.3, 0.2, −0.1)): fluid-frame direction moment battery; per-packet boost identities −p_μu^μ = ε and g^{μν}p_μp_ν = 0 | B (unit) | 4σ + machine | all moment margins m ≥ 0; identities at `rtol = 1e-14` for every packet |
| MCNX-PKT-06 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | single-cell position sampling at (S0, N_standard): per-coordinate mean vs cell center, variance vs w²/12 | B (unit) | 4σ | margin m ≥ 0 per coordinate for mean and variance |
| MCNX-PKT-07 | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | homogeneous 16³ fixture, N_tgt = 64, 4 steps: created-packet counts per (cell, species, step) | B (unit) | exact | every N_p ∈ [32, 128]; N_p ≥ 1 wherever E_tot > 0 |
| MCNX-GEO-01 | [geodesic-propagation](./geodesic-propagation.md) | null-condition identity along every scored Schwarzschild-benchmark trajectory, every substep | B (unit) | machine | \|g^{μν}p_μp_ν\|/(α p^t)² ≤ `1e-14` |
| MCNX-GEO-02 | [geodesic-propagation](./geodesic-propagation.md) | substep-halving p_t-drift convergence on analytically evaluated Kerr–Schild (Δt_sub ∈ {M/5, M/10, M/20}) | B (unit) | numeric | e(Δt)/e(Δt/2) ≥ 3.5 for both halvings |
| MCNX-GEO-03 | [geodesic-propagation](./geodesic-propagation.md) | fixed-seed single-rank CPU regression with time-varying lapse α(t) = 1 − `1e-2`·sin(t); per-packet TSV vs committed golden | A (regression) | golden | zero violations at `ABSTOL = RELTOL = 1e-12` |
| MCNX-GEO-04 | [geodesic-propagation](./geodesic-propagation.md) | instrumented fixed-seed runs: step/event/cell-crossing caps and escape accounting | B (unit) / A | exact | zero cap violations; escape records match removals exactly (totals at `rtol = 1e-10`); event-cap leg SKIP until interactions land |
| MCNX-GEO-05 | [geodesic-propagation](./geodesic-propagation.md) | flat-spacetime straight-line propagation: ≥ 100 isotropic packets, 256 substeps to t = 1.5 | B (unit) | machine | position error ≤ `1e-12` code units per component; \|Δp_i\| ≤ `1e-14`·\|p\|; \|Δp_t\| ≤ `1e-14`·\|p_t\| |
| MCNX-GEO-06 | [geodesic-propagation](./geodesic-propagation.md) | Kerr–Schild Schwarzschild p_t drift: N_smoke packets from (10M, 0, 0), t = 30M, Δx ∈ {M/4, M/8, M/16}, capture excision r < 3M | B (unit) | numeric + convergence | at Δx = M/8: all δ ≤ `2e-2`, ≥ 99% δ ≤ `1e-2`; median δ monotone decreasing with median(M/16) ≤ median(M/4)/1.5 |
| MCNX-OPA-01 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | argument-convention parity: ≥ `1e3` probes per channel through the MCNuX path vs direct WeakLibInterp `_Point` reference calls on synthetic tables (production-table leg env-gated, loud SKIP) | B (unit) | golden/parity | every probe agrees at `rtol = 1e-12`, `atol = 1e-30` |
| MCNX-OPA-02 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | range-enforcement sweep: out-of-range probes on every axis/side plus NaN probes | B (unit) | exact | below-range ρ → interaction-free with zero table calls; clamped calls bitwise-equal edge evaluations; non-finite input aborts; diagnostic counters integer-exact |
| MCNX-OPA-03 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | inversion error protocol: in-bounds, each out-of-bounds axis, NaN, no-root cases | B (unit) | exact | exact code set {0, 01, 02, 03, 10, 11, 13}; `T = 0` on every failure; failed T never consumed; pinned fallback selected |
| MCNX-OPA-04 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | species-axis mapping on synthetic two-dataset tables with distinguishable content | B (unit) | exact + machine | νe/ν̄e dataset selection exact; κ_a(νx) ≡ 0 and η_νx ≡ 0; νx Iso mean at `rtol = 1e-14` |
| MCNX-OPA-05 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | detailed-balance emissivity identity η_s(E)/(c κ_a) vs (4π E³/(hc)³) f_eq over a (s, E, T, μ) probe lattice | B (unit) | machine | identity at `rtol = 1e-14` wherever κ_a > 0; μ_ν,1 = −μ_ν,0 exact |
| MCNX-OPA-06 | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | residency instrumentation: read/upload counts; fixed probe set at step 0 vs step k and on 1 vs 2 ranks | B (unit) | exact | one read + one upload per table per run; probe evaluations bitwise identical |
| MCNX-INT-01 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | homogeneous-cell exponential event-time law: pure-absorption, pure-scattering, and mixed legs at (S0, N_standard) | B (unit) | 4σ | margins m ≥ 0 for interval mean and variance per leg and for the absorbed fraction κ_a/(κ_a + κ_s) |
| MCNX-INT-02 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | logged per-packet event sequences of fixed-seed fixtures, single-rank CPU | B (unit) | exact | draw budgets (2/2/0) and k-order exact; e increments by 1 per RNG-consuming event only |
| MCNX-INT-03 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | collimated beam scattered once in a moving-fluid cell (v = (0.3, 0.2, −0.1)), N_standard events | B (unit) | machine + 4σ | per-event \|−p′_μu^μ − ε\| ≤ `1e-14` ε and null identity ≤ `1e-14` ε²; outgoing-direction isotropy margins m ≥ 0 |
| MCNX-INT-04 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | instrumented fixed-seed absorption accounting on the INT-01/INT-06 fixtures | B (unit) | exact + relaxed | removal semantics and counts/species integer-exact; tally sums at `rtol = 1e-10`; ℓ = {+1, −1, 0} exact |
| MCNX-INT-05 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | low-T rule sweep T ∈ {0, 0.1, 0.3, 0.49, 0.5, 0.51, 5} MeV on a synthetic spectrum | B (unit) | machine + exact | η̃ = η(T_low)·(T/T_low)^6 at `rtol = 1e-14` below T_low; bitwise η̃ = η at/above; η̃ = 0 at T = 0 |
| MCNX-INT-06 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | 0-D cooling (ledger recurrence + continuum decay) and thermalization (census vs u_eq(1 − e^{−κ_a t})) benchmarks | B (unit) | relaxed + numeric | ledger recurrence at `rtol = 1e-10`; continuum decay within `3e-2` at t = 3 t_cool; thermalization within `1e-2` of u_eq at t = {1, 2, 4}/κ_a with ≥ `1e5` packets |
| MCNX-INT-07 | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | two-slab beam attenuation (κ_1 d_1 = 1.0, κ_2 d_2 = 0.5) at (S0, N_standard) | B (unit) | 4σ | transmitted fraction vs e^{−1.5}, binomial margin m ≥ 0 |
| MCNX-TRP-01 | [trapped-regime-treatment](./trapped-regime-treatment.md) | instrumented effective-coefficient stage over κ_a Δt_c ∈ [`1e-3`, `1e3`]: explicit selector everywhere, relabeling selector where κ_a Δt_c ≤ ξ | B (unit) | exact | bitwise pass-through of (η′, κ_a′, κ_s′); zero diffusion dispositions under `explicit` |
| MCNX-TRP-02 | [trapped-regime-treatment](./trapped-regime-treatment.md) | (κ_a, κ_s, η̃, Δt_c, ξ) probe lattice spanning ≥ 6 decades of κ_a Δt_c: relabeling equations (2021 α convention) and invariants recomputed independently | B (unit) | machine | η′ = αη̃, κ_a′ = ακ_a, κ_s′ = κ_s + (1−α)κ_a, κ_a′ + κ_s′ = κ_a + κ_s, η′/κ_a′ = η̃/κ_a each at `rtol = 1e-14` |
| MCNX-TRP-03 | [trapped-regime-treatment](./trapped-regime-treatment.md) | same lattice + instrumented runs: α = min(1, ξ/(κ_a Δt_c)) and the enforced budget | B (unit) | machine + exact | α at `rtol = 1e-14` (clamp branch bitwise 1); zero κ_a′Δt_c > ξ violations |
| MCNX-TRP-04 | [trapped-regime-treatment](./trapped-regime-treatment.md) | α → 1 equivalence: one fixed-seed single-rank CPU problem run under both selectors (κ_a Δt_c ≤ ξ everywhere) against the same committed golden per-packet TSV | A (regression) | golden | both diffs pass with zero violations at `ABSTOL = RELTOL = 1e-12` |
| MCNX-TRP-05 | [trapped-regime-treatment](./trapped-regime-treatment.md) | optically-thick thermalization with relabeling on: single gray cell, κ_a Δt_c = 10, ξ = 1 (α = 0.1), seed S0, census ≥ `1e5` packets | B (unit) | numeric | \|u_rad − u_eq(1 − e^{−κ_a′ t})\| ≤ `1e-2` u_eq at t = {1, 2, 4}/κ_a′, u_eq the unmodified η̃_tot/κ_a |
| MCNX-TRP-06 | [trapped-regime-treatment](./trapped-regime-treatment.md) | diffusion-regime fidelity: uniform pure-scattering medium (κ_s = 32, Δt = 0.5, t_end = 8, S0, ≥ `1e5` packets), diffusion path vs τ_diff = ∞ explicit path | B (unit) | numeric + 4σ + exact | position variance vs 2 D t_end (D = 1/(3 κ_s′)) within `5e-2`; means within 4σ, variances within `5e-2` between paths; per-leg energy at `rtol = 1e-14`; 8-draw budget and trigger discipline exact |
| MCNX-SRC-01 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | structural check of `rad_force`/`lep_source` (names, `CCC` centering, time levels) + instrumented write-phase discipline (zero → adds only → read-only) | B (unit) | exact | declared shape matches; zero non-accumulating writes outside the transport phase |
| MCNX-SRC-02 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | per-(cell, step) assembly identity from logged per-event tallies on fixed-seed fixtures | B (unit) | relaxed | S_μ and S_L match the assembly formulas at `rtol = 1e-10`; event-to-cell attribution integer-exact |
| MCNX-SRC-03 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | cadence/lag regression: fixed seed, single-rank CPU, ≥ 4 steps with mid-run checkpoint/restart; source variables at the read point vs committed golden TSV | A (regression) | golden | step-n read equals step-(n−1) tallies, step-0 read exactly zero; zero violations at `ABSTOL = RELTOL = 1e-12` |
| MCNX-SRC-04 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | global conservation ledger (deposits + census change + escapes) for four-momentum and lepton number on instrumented fixed-seed runs | B (unit) | relaxed + exact | both closures at `rtol = 1e-10` of gross exchange; νx lepton contributions integer-exact zero |
| MCNX-SRC-05 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | optional Tmunu contribution: synthetic packet populations with analytic T_rad^{μν} on an analytic-metric grid; off/on behavior, seeded-prior accumulation, 8-point cell-to-vertex average | B (unit) | exact + relaxed | off ⇒ eT* untouched; on ⇒ `+=` preserves prior exactly, zero packets ⇒ zero, vertex values at `rtol = 1e-10` |
| MCNX-SRC-06 | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | two-region lepton diffusion (nubhlight-style): closed optically-thick box, Ye = 0.35/0.15 halves, scripted consumer applying dYe/dt = m_b S_L/ρ, S0, ≥ `1e5` census | B (unit) | numeric + relaxed | total lepton number conserved at `rtol = 1e-10` every step; region-mean ΔYe strictly decreasing over 8 pinned outputs; every cell within `5e-2`·ΔYe(0) of Ye_eq at t_end |
| MCNX-PAR-01 | [particle-container-and-gpu](./particle-container-and-gpu.md) | synthetic ≥ `1e4`-packet population with distinguishable per-component values pushed through store/read-back, 1- and 2-rank redistribution, and checkpoint write/read | B (unit) | exact | every schema component of every packet bitwise-identical after each operation; s and (q, e) integer-exact |
| MCNX-PAR-02 | [particle-container-and-gpu](./particle-container-and-gpu.md) | two-level AMR fixture (refined central box), ≥ 4 fixed-seed transport steps at 1 and 2 ranks; exhaustive owner sweep after each step | B (unit) | exact | zero packets away from their owning (rank, level, grid, tile); owning level = finest level covering the position; step-created packets included |
| MCNX-PAR-03 | [particle-container-and-gpu](./particle-container-and-gpu.md) | instrumented per-packet displacement between redistributions + local-vs-global equivalence rerun of one step from an identical checkpointed state | B (unit) | exact | zero displacement-bound violations (bound = substeps × one cell width, and ≤ any local distance parameter); local and global paths yield identical assignments and bitwise-identical populations |
| MCNX-PAR-04 | [particle-container-and-gpu](./particle-container-and-gpu.md) | per-event-logged deposition on (a) a single-level fixture with emitters within one cell of box boundaries, (b) the two-level fixture, (c) ≥ 3 repeated threaded/GPU runs | B (unit) | relaxed | grid totals = per-event ledger at `rtol = 1e-10` (incl. guard-cell folding); fine-mask hierarchy sum and restricted coarse values at `rtol = 1e-10`; repeat-run tally norms within relative `1e-10` (MCNX-RNG-04); GPU leg SKIP where no device build |
| MCNX-PAR-05 | [particle-container-and-gpu](./particle-container-and-gpu.md) | fixed-seed single-rank CPU run with mid-run checkpoint/restart vs. uninterrupted golden; dual-grid restart leg (different rank count and max grid size) | A (regression) + B | golden + exact | per-packet TSV zero violations at `ABSTOL = RELTOL = 1e-12`; dual-grid leg: census multiset equal, residency invariant holds post-restart, zero duplicate q |
| MCNX-CTX-01 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | structural audit of the thorn tree: four ccl files, `REQUIRES CarpetX`, declared centerings, no Loop-layer accessor on transport/deposition paths | B (unit) | exact | all fixed-string/shape assertions pass |
| MCNX-CTX-02 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | full fixed-seed run under driver enforcement: presync mode, `poison_undefined_values = yes`, checksums on, Tmunu leg enabled | A (regression) | exact | clean completion with zero validity errors, zero poison reads, zero checksum violations; output diff passes |
| MCNX-CTX-03 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | schedule-declaration audit + per-routine invocation counters for every `MFIter`/`ParIter`-walking routine | B (unit) | exact | every such routine is level/global mode; invoked once per (level or run) per traversal, never per tile |
| MCNX-CTX-04 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | cadence counters on a ≥ 4-step instrumented run + the MCNX-SRC-03 regression run against the pinned schedule realization | A + B | exact + golden | transport phase executes exactly once per `CCTK_EVOL` traversal; packet/source state changes only in transport phases; SRC-03 diff passes at `ABSTOL = RELTOL = 1e-12` |
| MCNX-CTX-05 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | guard launch with `CarpetX::use_subcycling = yes` plus a control launch with `no` | B (unit) | exact | guarded launch aborts at `CCTK_PARAMCHECK` with the diagnostic and zero transport side effects; control launch proceeds |
| MCNX-CTX-06 | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | linear probe fields f = a + b_i x^i: vertex-to-cell-center bridge at every cell center; gather-at-position at ≥ `1e3` interior probe positions | B (unit) | machine | both bridges reproduce f at `rtol = 1e-14`; repeated evaluation bitwise deterministic |

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

- **The matrix grows with the corpus (by design, durable).** At any moment the matrix covers exactly the requirement ids of the specs on disk; rows for future specs (build-and-integration) are appended when those specs land. The closure rule makes "complete" a checkable property at every phase, so this growth does not weaken the contract.
- **GPU rows depend on a device build existing (assumption, non-blocking).** Rows that mention GPU execution (MCNX-RNG-01's device sweep, MCNX-RNG-04's GPU tallies) run in CPU-only form where no device build is available and report the device leg as SKIP — loudly, never silently.
- **Norm-based comparison granularity (assumption).** Norm tables aggregate per group per iteration; where a future row needs spatial localization beyond norms, it uses TSV line-cuts instead. No new output format is introduced without a spec change here.
