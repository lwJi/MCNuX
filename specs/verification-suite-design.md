# Verification suite design

> Technical leaf spec. Self-contained: an agent can implement the MCNuX verification
> suite — the `TestMCNuX` sibling test thorn, the analytic opacity mode, the benchmark
> set in the Cactus regression-harness format, and the requirement coverage matrix —
> from this file alone. It restates only the conventions it uses and references
> [conventions-and-units](./conventions-and-units.md) (units, constants, species) and
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (seeds, packet
> counts, the 4σ bar) for the rest. `README.md` is canonical if any restated convention
> here conflicts with it.

## Purpose & scope

This spec designs how MCNuX is verified — the concrete test infrastructure and the
closure of that infrastructure over every correctness requirement in the corpus:

- The **`TestMCNuX` sibling test thorn**: analytic backgrounds provided through the
  base thorns' native parameter-extension idiom (`SHARES:` + `EXTENDS KEYWORD`), with
  the exact background formulas — Schwarzschild spacetime, uniform-sphere and
  equilibration-box hydro profiles — restated in the corpus's own 3+1 notation.
- The **analytic (gray) opacity mode**, a first-class MCNuX physics capability
  (parameter-selected opacity source), with its binding per-species formulas. It is the
  νx coverage mechanism (tabulated νx emission is identically zero,
  [opacity-eos-evaluation](./opacity-eos-evaluation.md)) and the trapped-regime
  equilibration harness.
- The **deterministic benchmark set** in the Cactus regression-harness format —
  `<thorn>/test/<name>.par` plus a same-named directory of golden TSV/norm data, diffed
  at the harness defaults `ABSTOL=RELTOL=1e-12` — including Schwarzschild `p_t`
  conservation, fixed-seed single-rank reproducibility runs, and the α → 1 limit check
  of [trapped-regime-treatment](./trapped-regime-treatment.md).
- The **statistical benchmark reduction**: every N-sigma Monte Carlo check is computed
  by an analysis routine into numeric norm/TSV outputs, so statistical acceptance also
  reduces to deterministic golden data at the pinned seed.
- The **coverage matrix** — the closure list over every `[MCNX-XXX-NN]` requirement id
  declared anywhere in the corpus, mapping each id to the check(s) that verify it.

Out of scope:

- The correctness requirements themselves — each is owned by the spec whose code its
  id carries (see the coverage matrix; the id prefix names the owner).
- The build and dependency contract, the pinned test ThornList, and the
  production/test activation split — [build-and-integration](./build-and-integration.md)
  (this spec pins *what* the tests are; that spec pins how the configuration that runs
  them is assembled, and that golden data lives under `MCNuX/test/`).
- Implementing MCNuX or `TestMCNuX` (this corpus specifies; the Ralph loop implements),
  and modifying the flesh harness (consumed as-is).
- The RNG algorithm, seeds, packet counts, and the 4σ bar
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)); the
  tabulated opacity contract ([opacity-eos-evaluation](./opacity-eos-evaluation.md)).

## Source of truth

Harness mechanics, output formats, and the test-thorn idiom are pinned from the flesh
and the CarpetX stack; the background formulas are corpus-derived closed forms
(standard textbook solutions restated in CNV notation, normative as written per
`README.md`). Reference artifacts:

- `flesh/lib/sbin/RunTestUtils.pl` — the regression harness: test discovery
  (`<thorn>/test/<name>.par` + same-named golden directory), whitespace-tokenized
  per-column diffing with strong failure at
  `|new − old| ≥ max(ABSTOL, RELTOL·max(|old|,|new|))`, defaults
  `ABSTOL = RELTOL = 1e-12`, tolerance resolution global → thorn `test.ccl` → per-test
  `TEST` block; `flesh/lib/sbin/RunTest.pl` — the runner.
- `CarpetX/CarpetX/test/test.ccl` — the single global `EXTENSIONS tsv` declaration
  that makes `.tsv` golden files recognized for every thorn in the ThornList.
- `CarpetX/CarpetX/src/io_tsv.cxx`, `CarpetX/CarpetX/src/io_norm.cxx` — the
  deterministic TSV line-cut and norm-table writers (`out_tsv_vars`, `out_norm_vars`)
  whose outputs are the golden formats.
- `CarpetX/WaveToyX/test/standing.par` — the worked harness example this suite's
  parfiles mirror; `CarpetX/TestODESolvers/test/test.ccl` — the per-test tolerance
  block precedent (any non-default tolerance carries an in-file justification).
- `CarpetX/ADMBaseX/param.ccl` — ADMBaseX owns the RESTRICTED keywords `initial_data`
  (default `"Cartesian Minkowski"`), `initial_lapse` (default `"one"`), `initial_shift`
  (default `"zero"`); `CarpetX/ADMBaseX/schedule.ccl` — the `ADMBaseX_InitialData` and
  `ADMBaseX_InitialGauge` hook groups and the keyword-conditional scheduling of the
  base thorn's own initializers. Minkowski needs no new code.
- `CarpetX/HydroBaseX/param.ccl` — HydroBaseX owns the RESTRICTED keyword
  `initial_hydro` (sole base choice `"vacuum"`); `CarpetX/HydroBaseX/schedule.ccl` —
  the `HydroBaseX_InitialData` hook group `TestMCNuX` schedules into;
  `CarpetX/HydroBaseX/interface.ccl` — the cell-centered fluid variables the profiles
  populate.
- `flesh/lib/sbin/parameter_parser.pl` and `flesh/lib/sbin/ImpParamConsistency.pl` —
  the flesh's `EXTENDS KEYWORD` machinery (a sharing thorn adds ranges to another
  implementation's RESTRICTED keyword); external-extension evidence in the stack:
  `CarpetX/CarpetX/par/brill-lindquist-write.par` sets
  `ADMBaseX::initial_data = "Brill-Lindquist"`, a value absent from ADMBaseX's own
  choice list (see Open questions).

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; precedent for analytic-opacity validation problems and the source of the restated transport equations the benchmarks exercise |
| arXiv:1903.09273v1 | Miller, Ryan, Dolence 2019 — νbhlight; analytic-opacity cooling/equilibration test precedent |
| arXiv:1505.05119v1 | Ryan, Dolence, Gammie 2015, ApJ 807, 31 — bhlight; optically-thin analytic cooling test precedent |
| arXiv:1507.03606v2 | Richers et al. 2015, ApJ 813, 38 — Sedonu; blackbody emission/equilibration test precedent |

## Inputs & outputs

### Suite artifact shape

| Artifact | Location | Content |
|---|---|---|
| benchmark parameter files | `MCNuX/test/<name>.par` (the MCNuX thorn's `test/` directory; [build-and-integration](./build-and-integration.md) owns the location pin) | one per benchmark below; ActiveThorns line activates `TestMCNuX` where a non-trivial background is needed |
| golden data | `MCNuX/test/<name>/` | committed TSV line-cuts and `norms/*.tsv` tables produced by the CarpetX writers at the primary seed |
| tolerance config | `MCNuX/test/test.ccl` | `NPROCS 1` for every test; per-test `TEST` blocks only where a non-default tolerance is justified in-file |
| `TestMCNuX` thorn | sibling thorn in the MCNuX repository | initial-data extensions only ([MCNX-VER-02]); no transport physics, no output machinery of its own |
| non-harness gates | CI scripts / compile-time asserts (uncommitted to `specs/`) | build gates, inspection gates, the paramcheck-abort check (a deliberately-aborting run cannot be a harness test — the harness treats nonzero exit as failure) |

Outputs: harness verdicts (every benchmark passes against its golden data), plus the
gate results. Statistical checks output their standardized deviations as numeric data
(see [MCNX-VER-07]); nothing in the suite requires human judgment to score.

### Analytic opacity mode parameters (per species s ∈ {νe, ν̄e, νx})

| Parameter | Meaning | Units / range |
|---|---|---|
| `kappa_a0(s)` | constant absorption opacity | cm⁻¹, ≥ 0 |
| `kappa_s0(s)` | constant elastic scattering opacity | cm⁻¹, ≥ 0 |
| `eta_scale(s)` | emissivity scale on the Kirchhoff form | dimensionless, ≥ 0 (1 = equilibrium Kirchhoff; 0 = absorption-only medium) |

## Correctness requirements

- **[MCNX-VER-01] Harness-format contract.** Every runnable MCNuX verification
  benchmark is a Cactus regression test: `MCNuX/test/<name>.par` plus the same-named
  golden directory of CarpetX TSV/norm outputs, discovered and diffed by the flesh
  harness (`flesh/lib/sbin/RunTestUtils.pl`) at the defaults `ABSTOL=RELTOL=1e-12`;
  `.tsv` recognition comes from the CarpetX `EXTENSIONS tsv` declaration. Every test
  runs `NPROCS 1` (single rank — the exact-reproducibility mode of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)). Golden data
  is generated at the primary seed `1296518744`; a tolerance different from the
  1e-12 defaults is admissible only via a per-test `TEST` block carrying an in-file
  justification (the `CarpetX/TestODESolvers/test/test.ccl` precedent).

- **[MCNX-VER-02] `TestMCNuX` thorn design (the native extension idiom).**
  `TestMCNuX` is a sibling thorn in the MCNuX repository whose sole job is analytic
  initial data; it is specified here (no 14th spec file) and never activated outside
  test parameter files ([build-and-integration](./build-and-integration.md) owns that
  pin). Binding design:

  - `param.ccl` uses the base thorns' native extension idiom — `SHARES: ADMBaseX` with
    `EXTENDS KEYWORD initial_data` adding the value `"Schwarzschild"` and
    `EXTENDS KEYWORD initial_lapse` adding `"Schwarzschild"`; `SHARES: HydroBaseX`
    with `EXTENDS KEYWORD initial_hydro` adding `"uniform sphere"` and
    `"equilibration box"`. The quoted keyword values are binding (parameter files and
    golden data reference them); `TestMCNuX`'s own parameter names (masses, radii,
    profile values) are free but must cover every symbol in [MCNX-VER-03] and
    [MCNX-VER-04]. The Schwarzschild shift needs no extension: ADMBaseX's default
    `initial_shift = "zero"` is the exact solution's shift.
  - `schedule.ccl` schedules the writers into the base thorns' empty hook groups,
    conditional on the keyword values: the metric/curvature writer
    `IN ADMBaseX_InitialData`, the lapse writer `IN ADMBaseX_InitialGauge`, the hydro
    writers `IN HydroBaseX_InitialData` — with complete WRITES declarations
    (`(everywhere)`), per [carpetx-thorn-integration](./carpetx-thorn-integration.md)
    [MCNX-CTX-04].
  - `TestMCNuX` contains no transport physics, no output machinery, and writes only
    ADMBaseX/HydroBaseX variables through the extension groups above. Minkowski
    backgrounds use no `TestMCNuX` code at all: `ADMBaseX::initial_data` defaults to
    `"Cartesian Minkowski"` with lapse-one/shift-zero initializers scheduled by the
    base thorn itself.

- **[MCNX-VER-03] Schwarzschild background (exact formulas, normative).** The
  `"Schwarzschild"` initial data is the Schwarzschild solution of mass `M` in
  **isotropic coordinates** (corpus-derived closed form in CNV 3+1 notation; static,
  `∂_t g_μν = 0`):

  ```text
  r = √(x² + y² + z²) ,   ψ = 1 + M/(2r)
  γ_ij = ψ⁴ δ_ij ,   K_ij = 0 ,   β^i = 0
  α = (1 − M/(2r)) / (1 + M/(2r))
  ```

  with `M` a `TestMCNuX` parameter (default 1, geometrized units). The horizon sits at
  isotropic `r = M/2` (areal radius 2M); the formulas are valid and the writers are
  exercised only for `r > M/2` — every pinned benchmark places its computational
  domain entirely in `r > M/2` (an offset box; behavior at `r ≤ M/2` is unspecified
  and unreachable in the pinned parfiles). Because the metric is static and `β^i = 0`,
  the covariant energy `p_t = −α² p^t` of every packet is exactly conserved along
  continuum geodesics ([geodesic-propagation](./geodesic-propagation.md)
  [MCNX-GEO-05]) — the property the `schwarzschild-pt` benchmark measures. Writers
  must reproduce these closed forms at every grid point to the machine tier (`~1e-14`
  relative).

- **[MCNX-VER-04] Hydro test profiles (exact formulas, normative).** Both profiles
  set the fluid at rest (`velˣ = velʸ = velᶻ = 0`, i.e. `u^μ = (1/α, 0, 0, 0)`
  normalized against the active metric) and populate `HydroBaseX::temperature` **in
  MeV** — the test-side exercise of the interface requirement of
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-04].
  With `(ρ₀, T₀, Yₑ₀)` and `(ρ_amb, T_amb, Yₑ_amb)` `TestMCNuX` parameters
  (ρ in g/cm³, T in MeV, Yₑ dimensionless) and sphere center `x₀^i` and radius `R_s`
  in geometrized coordinates:

  ```text
  uniform sphere:      (ρ, T, Yₑ)(x) = (ρ₀, T₀, Yₑ₀)      if |x − x₀| ≤ R_s
                                       (ρ_amb, T_amb, Yₑ_amb)  otherwise
  equilibration box:   (ρ, T, Yₑ)(x) = (ρ₀, T₀, Yₑ₀)      everywhere
  ```

  All other HydroBaseX variables the writer touches (`eps`, `press`, `entropy`,
  `Bvec`, vector potentials) are set to zero — MCNuX transport reads only
  `(ρ, T, Yₑ, vel)`. Discontinuity convention at the sphere surface: the cell-centered
  value is the profile evaluated at the cell center (no volume-fraction smoothing).

- **[MCNX-VER-05] Analytic (gray) opacity mode (first-class capability; formulas
  normative).** MCNuX carries a parameter-selected opacity source: the production
  weaklib tables ([opacity-eos-evaluation](./opacity-eos-evaluation.md)) or the
  analytic mode below, feeding the identical per-species coefficient interface
  (κ_a, κ_s, η per single species, g = 1). This is a physics capability, not test
  scaffolding — constant-opacity benchmarks are how the reference codes validate
  (arXiv:1708.08452 §3; arXiv:1903.09273; arXiv:1505.05119; arXiv:1507.03606) — and
  it is internal by necessity: opacity evaluation happens inside MCNuX's device
  kernels, and the weaklib table schema cannot express νx emissivity. The binding
  formulas, per species s, with E and T in MeV:

  ```text
  κ_a(s, E) = κ_a0(s)
  κ_s(s, E) = κ_s0(s)
  η(s, E)   = η_scale(s) · c κ_a0(s) · 4πE³/(hc)³ · 1/(exp(E/T) + 1)
  ```

  (units as the coefficient interface: κ in cm⁻¹, η in MeV cm⁻³ s⁻¹ MeV⁻¹; constants
  `c = 2.99792458e10 cm/s`, `hc = 1.239841984e-10 MeV cm` — the same pins as
  [opacity-eos-evaluation](./opacity-eos-evaluation.md)). `η_scale = 1` gives exact
  Kirchhoff equilibrium at μ_ν = 0 (equilibrium spectral energy density
  `J_eq(E) = η/κ_a · 1/c = 4πE³/(hc)³ · f_eq(E)`, the analytic target of the
  equilibration benchmarks); `η_scale = 0` gives an absorption-only medium (beam
  benchmarks). The analytic mode is the **νx coverage mechanism**: `κ_a0(νx) > 0`
  and `η_scale(νx) = 1` exercise νx emission, absorption, and the g = 4 creation
  factor, none of which the tables can reach; it is likewise the **trapped-regime
  harness** — the relabeling of
  [trapped-regime-treatment](./trapped-regime-treatment.md) applies to the analytic
  coefficients identically (a pure post-map on the coefficient interface). Everything
  downstream of the coefficient interface (sampling, events, tallies, RNG) is
  identical in both modes.

- **[MCNX-VER-06] Deterministic benchmark set (binding).** The suite contains at
  least the following harness benchmarks, each fixed-seed (`1296518744`) and
  single-rank, each archiving golden TSV/norm data at `ABSTOL=RELTOL=1e-12`:

  | Benchmark | Setup | Primary checks |
  |---|---|---|
  | `unit-selftest` | Minkowski, no transport; `MCNuX::run_selftests = yes` runs the self-test battery at initial time and writes one row per named check (measured error, pass flag) as TSV | the battery of [MCNX-VER-08]'s matrix: RNG KAT vectors and [0,1) closed forms, conversion-factor recomputation and round-trips, νx mapping and relabeling identities, monotonicity tripwire, tetrad/null identities, Tmunu closed forms, sign fixtures, range-policy fixtures, inversion-protocol fixtures, `sizeof(amrex::Real) == 8` |
  | `minkowski-freestream` | Minkowski, pinned initial packets, zero opacities | straight lines at `\|dx/dt\| = 1`, `p_i` constant to machine tier per step, flat-spacetime exactness of [geodesic-propagation](./geodesic-propagation.md) [MCNX-GEO-04] |
  | `schwarzschild-pt` | `TestMCNuX` `"Schwarzschild"` (M = 1), domain in `r > M/2`, no interactions; pinned trajectory family at isotropic `r ∈ [6M, 10M]`, integrated for `T = 100M` at resolution `h = M/16`; per-packet `p_t` written as TSV | `p_t` drift `\|p_t(T) − p_t(0)\|/\|p_t(0)\| ≤ 1e-5` for every trajectory at the pinned resolution (design bound — see Open questions); companion legs at `Δt, Δt/2, Δt/4` show convergence order ≥ 2 of the drift ([MCNX-GEO-04], [MCNX-GEO-05]) |
  | `emission-fixedseed` | equilibration box, analytic mode, one emission-only step | exact packet set (ids, states) golden; bin-center energies; count law; creation draw-map audit ([packet-representation-and-sampling](./packet-representation-and-sampling.md) [MCNX-PKT-05]) |
  | `interactions-fixedseed` | uniform sphere, analytic mode, pinned small population, several steps; per-packet event log written as TSV | episode structure and draw-map audits ([MCNX-INT-06]; the constant five-draw diffusion audit of [MCNX-TRP-07]), the diffusion leg (prelude handoff and regime decision on Δt′_rem per [MCNX-TRP-06], causal cap, two-leg position update), discreteness, elasticity identities, id retirement, ledger closure per step |
  | `ownership-multibox` | Minkowski, ≥ 2 boxes per rank, packets crossing box boundaries | after-step ownership invariant, bounded-motion audit, one-invocation-per-(patch, level) scheduling counter ([MCNX-GPU-05], [MCNX-GPU-06], [MCNX-CTX-03]) |
  | `cadence-lag` | uniform sphere, analytic mode, with ODESolvers evolving a trivial state; transport-step counter and source-tally checksums as norms | transport advances exactly once per `CCTK_EVOL` iteration; injected fluid step-change at iteration n first affects source reads at n+1; checkpoint/recovery leg reproduces partner-visible values ([MCNX-CTX-05], [MCNX-HYD-01], [MCNX-HYD-03]) |
  | `trp-alpha-limit-explicit` / `trp-alpha-limit-relabeled` | equilibration box, analytic mode; the second parfile enables the relabeling with the fixed-α override pinned to α = 1 | the α → 1 limit check of [trapped-regime-treatment](./trapped-regime-treatment.md) [MCNX-TRP-05]: both parfiles are diffed against the **same golden numbers** (generated from the explicit run), so any divergence at α = 1 fails the harness at 1e-12 |

  Two guards are non-harness by necessity: the paramcheck-abort check
  (`CarpetX::use_subcycling = yes` must abort at `CCTK_PARAMCHECK`,
  [carpetx-thorn-integration](./carpetx-thorn-integration.md) [MCNX-CTX-06]) and the
  configure-without-CarpetX failure ([MCNX-CTX-02]) — both run as CI/script gates,
  because the harness scores any nonzero exit as failure.

- **[MCNX-VER-07] Statistical checks reduce to norm-based numeric golden data.**
  Every statistical (4σ) check in the corpus is computed by an MCNuX analysis routine
  that writes, per check: the estimate, the expected value, the standard error σ, and
  the standardized deviation `z = (estimate − expected)/σ`, as grid scalars/norm
  outputs in the same TSV/norm formats. Acceptance is `|z| ≤ 4` at the pinned packet
  count and seed — verified once when the golden data is created — and thereafter the
  archived z values are themselves golden numbers diffed at 1e-12 (a statistical check
  at pinned seed is deterministic; a code change that moves any estimator shows up as
  a golden diff even when `|z| ≤ 4` still holds). The binding statistical benchmark
  set: `stats-emission` (count law, emitted energy, isotropy, and the νx g = 4 leg —
  νx emitted energy = 4× the single-species analytic law), `stats-beam` (transmitted
  fraction `exp(−κ_a L)` through an `η_scale = 0` slab), `stats-scatterbox` (Poisson
  event counts and post-scatter isotropy), `stats-equilibration-explicit` /
  `stats-equilibration-relabeled` (spectrum and energy density against the analytic
  `J_eq` of [MCNX-VER-05], the relabeled leg at α < 1 with a measured event-count
  reduction), `stats-diffusion` (atom frequency vs `exp(−κ_s′ Δt′_rem)`;
  continuous-branch r̃ moments vs the normalized truncated analytic moments; sampled
  cos θ₂ distribution vs the restated θ₂ law of [MCNX-TRP-10]), and
  `trp-overlap-consistency` (the trapped-regime harness pinned just above and just
  below `τ_diff`; statistical comparison of the two populations' end-of-step
  position and direction distributions — both share the discrete prelude of
  [MCNX-TRP-06], so they differ only in how the step remainder is finished).
  Defaults `N_p = 1048576`, seed `1296518744` per
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md); a benchmark
  pinning different values documents them in its parfile.

- **[MCNX-VER-08] Coverage matrix is the closure list.** The coverage matrix in the
  Verification section below lists **every** `[MCNX-XXX-NN]` requirement id declared
  in any coded spec's Correctness requirements section, and maps each to the named
  check(s) — benchmarks of [MCNX-VER-06]/[MCNX-VER-07], `unit-selftest` battery
  entries, or build/inspection/mechanical gates — that verify it. The matrix may not
  contain an id no spec declares (a phantom row), and no declared id may be absent
  (uncovered). `tools/validate_specs.sh` enforces both directions mechanically on
  every run; coverage grows by extending the matrix when a spec adds an id.

### Conventions restated (the subset this leaf uses)

- Signature `(−,+,+,+)`; lapse `α`, shift `β^i`, spatial metric `γ_ij`;
  `p_t = −α² p^t + β^i p_i`; geometrized units `G = c = M_sun = 1` for coordinates,
  M, Δt, and packet momenta. In [MCNX-VER-03] `α` is the lapse; the trapped-regime
  relabeling fraction α appears only in the `trp-alpha-limit` / relabeled-benchmark
  rows.
- Microphysics: E and T in MeV, ρ in g/cm³, Yₑ dimensionless; `HydroBaseX::temperature`
  is populated in MeV by every test background; constants `c` and `hc` exactly as
  pinned by [conventions-and-units](./conventions-and-units.md) and
  [opacity-eos-evaluation](./opacity-eos-evaluation.md).
- Tolerance tiers per `README.md`: machine `~1e-14`, golden parity `1e-12`
  (`ABSTOL=RELTOL=1e-12` harness defaults), relaxed `1e-10` with stated rationale,
  exact for KATs/enumerations/event logs; statistical acceptance 4σ at `N_p = 1048576`,
  primary seed `1296518744`, secondary seed `20260730`.

## Verification

This spec's own correctness is verified by the suite existing in exactly the pinned
shape and by mechanical closure:

- **Harness discovery (exact).** `flesh/lib/sbin/RunTestUtils.pl` discovers exactly
  the benchmarks of [MCNX-VER-06]/[MCNX-VER-07] under `MCNuX/test/`; every parfile has
  its same-named golden directory; every test runs `NPROCS 1` and passes at the
  1e-12 defaults (or a justified per-test block).
- **Background closed forms (machine tier, `~1e-14`).** A `unit-selftest` leg
  evaluates the `TestMCNuX` writers' outputs at pinned grid points against the closed
  forms of [MCNX-VER-03]/[MCNX-VER-04]; on the Schwarzschild data it additionally
  checks the algebraic identity `α (1 + M/(2r)) = 1 − M/(2r)` and `γ_ij = ψ⁴ δ_ij`
  component-wise, each to `1e-14` relative at every sampled point.
- **Analytic-mode identities (machine tier).** `η(s, E)/κ_a(s, E)` equals
  `η_scale · c · 4πE³/(hc)³ · f_eq(E)` to `1e-14` relative wherever `κ_a0 > 0`; the
  mode switch changes no downstream code path (verified by the α = 1-style identical
  golden data of a table-mode vs analytic-mode run pinned to identical coefficient
  values on a matched state).
- **Statistical-reduction audit (exact).** Every 4σ check named anywhere in the
  corpus appears as a z-score row in exactly one `stats-*` benchmark's golden data;
  `|z| ≤ 4` holds in every committed golden file.
- **α → 1 golden identity (exact / 1e-12).** The two `trp-alpha-limit` parfiles diff
  against byte-identical golden numbers, per [MCNX-VER-06].

### Coverage matrix (the closure list)

Every requirement id declared in the corpus, with the check(s) that verify it. Check
names: benchmarks per [MCNX-VER-06]/[MCNX-VER-07]; `selftest:` legs of
`unit-selftest`; `gate:` non-harness build/inspection/CI gates; `mech:` checks run by
`tools/validate_specs.sh` on every validator invocation.

| Id | Verified by |
|---|---|
| [MCNX-CNV-01] | selftest: tetrad/null/`p_t` identity legs (notation consistency exercised wherever equations are evaluated); mech: restated-equation needles |
| [MCNX-CNV-02] | selftest: unit round-trips (cgs ↔ code, MeV ↔ K); gate: no-runtime-unit-parameter source sweep |
| [MCNX-CNV-03] | selftest: factor recomputation from the five defining constants; gate: single-constant-set source sweep |
| [MCNX-CNV-04] | selftest: every derived factor recomputed to `1e-14` relative |
| [MCNX-CNV-05] | selftest: MeV→K applied at every WeakLibInterp call site (call-boundary fixtures, incl. NES/Pair); stats-equilibration-explicit |
| [MCNX-CNV-06] | selftest: species enumeration/lepton-number fixtures; stats-emission (νx g = 4 leg) |
| [MCNX-RNG-01] | selftest: Philox4x32-10 KAT vectors, host and device |
| [MCNX-RNG-02] | selftest: the three KAT vectors, exact word equality |
| [MCNX-RNG-03] | selftest: counter/key packing fixtures (pinned tuples → pinned words) |
| [MCNX-RNG-04] | selftest: [0, 1) closed forms (`u64 = 0, 2^11 − 1, 2^11, 2^64 − 1`) |
| [MCNX-RNG-05] | gate: inspection — no stateful-RNG call in any physics path; emission-fixedseed and interactions-fixedseed draw-map audits |
| [MCNX-RNG-06] | emission-fixedseed, interactions-fixedseed (run-twice bitwise identity of event logs); every deterministic golden benchmark |
| [MCNX-RNG-07] | all `stats-*` benchmarks (z-score reduction per [MCNX-VER-07]) |
| [MCNX-RNG-08] | gate: no bitwise cross-configuration assertion exists in the suite (inspection); norm-tolerance comparisons in `stats-*` cross-checks |
| [MCNX-PKT-01] | emission-fixedseed (exact packet-set golden); interactions-fixedseed (id/weight immutability) |
| [MCNX-PKT-02] | stats-emission (count law and g = 4 leg); emission-fixedseed (floor/Bernoulli determinism) |
| [MCNX-PKT-03] | emission-fixedseed (bin-center exactness, weight `E_p/ν`, position/time draws) |
| [MCNX-PKT-04] | selftest: static-flat tetrad identity, null and ν-recovery identities; stats-emission isotropy |
| [MCNX-PKT-05] | emission-fixedseed draw-map audit (exact `(q, e, k)` trace) |
| [MCNX-PKT-06] | stats-emission |
| [MCNX-GEO-01] | schwarzschild-pt (drift + convergence legs); minkowski-freestream |
| [MCNX-GEO-02] | selftest: null-closure identity `α²(p^t)² = γ^{ij} p_i p_j`; schwarzschild-pt |
| [MCNX-GEO-03] | selftest: linear-field gather exactness (values and derivatives) |
| [MCNX-GEO-04] | minkowski-freestream (flat exactness); schwarzschild-pt Δt-convergence legs (order ≥ 2) |
| [MCNX-GEO-05] | schwarzschild-pt (`p_t` drift ≤ 1e-5 at pinned resolution; convergent at order ≥ 2) |
| [MCNX-INT-01] | stats-beam (`exp(−κ_a L)`); interactions-fixedseed (episode kinematics audit) |
| [MCNX-INT-02] | interactions-fixedseed (episode/competition audit); stats-scatterbox (Poisson counts) |
| [MCNX-INT-03] | interactions-fixedseed (discreteness: weights never change, whole-packet removal); stats-beam |
| [MCNX-INT-04] | selftest: elasticity/null identities per scattering event; stats-scatterbox isotropy |
| [MCNX-INT-05] | selftest: coefficient provenance parity (event-side κ equals coefficient-interface value bitwise) |
| [MCNX-INT-06] | interactions-fixedseed draw-map audit |
| [MCNX-TRP-01] | trp-alpha-limit-explicit; stats-equilibration-explicit |
| [MCNX-TRP-02] | selftest: relabeling formulas and the two invariants (`η′/κ_a′ = η/κ_a`, `κ_a′ + κ_s′ = κ_a + κ_s`) |
| [MCNX-TRP-03] | selftest: convention tripwire (monotonicity of κ_a′ and κ_s′ in α) |
| [MCNX-TRP-04] | selftest: α-selection rule reproduction; per-cell `κ_a′Δt_c ≤ ξ` diagnostic in stats-equilibration-relabeled |
| [MCNX-TRP-05] | trp-alpha-limit-explicit vs trp-alpha-limit-relabeled (identical golden numbers) |
| [MCNX-TRP-06] | interactions-fixedseed diffusion leg (prelude handoff, regime decision on Δt′_rem); trp-overlap-consistency |
| [MCNX-TRP-07] | interactions-fixedseed diffusion draw-map audit (constant five draws) |
| [MCNX-TRP-08] | stats-diffusion (atom frequency; r̃ moments); interactions-fixedseed (causal cap) |
| [MCNX-TRP-09] | interactions-fixedseed diffusion leg (two-leg position update); machine-tier f_free-limit identities |
| [MCNX-TRP-10] | stats-diffusion (cos θ₂ distribution); machine-tier atom identity (cos θ₂ = 1) |
| [MCNX-OPA-01] | selftest: kernel parity — MCNuX values vs directly-invoked WeakLibInterp `_Point` kernels, bitwise, all five families |
| [MCNX-OPA-02] | selftest: Kelvin call-boundary fixtures at every family (incl. NES/Pair) |
| [MCNX-OPA-03] | selftest: log10-ownership fixtures (raw vs pre-logged argument parity) |
| [MCNX-OPA-04] | selftest: Kirchhoff-closure identity; stats-equilibration-explicit (spectrum → f_eq) |
| [MCNX-OPA-05] | selftest: νx mapping identities (zeros exact; ½-mean to `1e-14`) |
| [MCNX-OPA-06] | selftest: range-policy fixtures (transparency floor, clamping, no-NaN, clamp counters) |
| [MCNX-OPA-07] | selftest: inversion-protocol fixtures (`Error ≠ 0 ⇒ T = 0` surfaced; test-mode hard failure) |
| [MCNX-OPA-08] | mech: snapshot anchoring against the committed `.h5ls` fixtures, re-run every validator invocation |
| [MCNX-HYD-01] | cadence-lag (checkpoint/recovery leg); gate: validity-machinery run (complete READS/WRITES) |
| [MCNX-HYD-02] | selftest: sign fixtures (single-event ledger entries); interactions-fixedseed ledger closure |
| [MCNX-HYD-03] | cadence-lag (zero/add/read ordering); gate: contributor-split meta-test + poison/checksum run |
| [MCNX-HYD-04] | cadence-lag (one-step lag); MeV interface exercised by every `TestMCNuX` background ([MCNX-VER-04]); gate: coupled-configuration checklist |
| [MCNX-HYD-05] | per-step five-channel ledger-closure diagnostic in every deterministic benchmark (relative 1e-13) |
| [MCNX-HYD-06] | selftest: Tmunu closed forms (one-packet cell value; uniform and linear `q` fields; vertex-total identity) |
| [MCNX-HYD-07] | gate: output-surface audit (MCNuX writes only its own variables + the optional TmunuBaseX contribution); mech: D15 no-GRMHD-citation sweep |
| [MCNX-GPU-01] | gate: compile-time layout pin (`ParticleContainerPureSoA` static_assert; trivially-copyable tile view) |
| [MCNX-GPU-02] | gate: GPU-backend build of every per-packet operator as device kernels + by-value POD capture inspection |
| [MCNX-GPU-03] | gate: compile-time schema check (one SoA component per physical component, binary64); emission-fixedseed |
| [MCNX-GPU-04] | interactions-fixedseed id-contract leg (uniqueness, immutability, retirement, bit 63 clear) |
| [MCNX-GPU-05] | ownership-multibox (`OK()`-style placement check after every step) |
| [MCNX-GPU-06] | ownership-multibox (bounded-motion audit: measured ≤ derived; regrid leg redistributes before next transport op) |
| [MCNX-GPU-07] | gate: inspection — no host per-particle loops, no AoS container in MCNuX; mech: anti-pattern needle in both technical specs |
| [MCNX-CTX-01] | gate: build + inspection (grid access via `ghext`/AmrCore/MultiFab/`Array4`/`indextype` only); validity-machinery run |
| [MCNX-CTX-02] | gate: configure without CarpetX fails on `REQUIRES CarpetX`; pinned ThornList builds clean |
| [MCNX-CTX-03] | ownership-multibox scheduling-mode counter (one invocation per (patch, level), never per tile) |
| [MCNX-CTX-04] | gate: poison/presync validity run + deliberately-under-declared meta-test fires |
| [MCNX-CTX-05] | cadence-lag (once-per-EVOL counter; one-step lag golden TSV) |
| [MCNX-CTX-06] | gate: paramcheck-abort check (`use_subcycling = yes` aborts; `no` runs) |
| [MCNX-CTX-07] | gate: inspection (interpolator execution model absent from MCNuX); mech: anti-pattern needle in both technical specs |
| [MCNX-VER-01] | harness discovery gate: every benchmark is a discovered `<name>.par` + golden dir, `NPROCS 1`, 1e-12 defaults |
| [MCNX-VER-02] | gate: `TestMCNuX` ccl inspection (extension idiom, keyword values, hook-group scheduling); backgrounds exercised by schwarzschild-pt and the hydro benchmarks |
| [MCNX-VER-03] | selftest: Schwarzschild writer closed-form leg; schwarzschild-pt |
| [MCNX-VER-04] | selftest: hydro-profile closed-form leg; stats-equilibration-explicit, stats-beam (profiles in use) |
| [MCNX-VER-05] | selftest: analytic-mode identity legs; stats-emission νx leg; stats-equilibration pair |
| [MCNX-VER-06] | harness run: the pinned benchmark set exists and passes; α-limit golden identity |
| [MCNX-VER-07] | statistical-reduction audit: every 4σ check appears as a z-score row with `\|z\| ≤ 4` in committed golden data |
| [MCNX-VER-08] | mech: requirement-id closure (ownership, declaration, coverage; phantom rows rejected) on every validator run |
| [MCNX-BLD-01] | gate: repository/arrangement layout audit; pinned ThornList compiles both thorns |
| [MCNX-BLD-02] | gate: configure-time dependency checks (`REQUIRES CarpetX`, `REQUIRES WeakLibInterp`); build without either fails |
| [MCNX-BLD-03] | selftest: `sizeof(amrex::Real) == 8` and `sizeof(amrex::ParticleReal) == 8`; gate: CPU and GPU configurations both build |
| [MCNX-BLD-04] | gate: the pinned test ThornList configures, builds, and runs the full harness suite clean |
| [MCNX-BLD-05] | harness discovery gate: all benchmarks (and only they) discovered under `MCNuX/test/`; `TestMCNuX` ships no `test/` directory |
| [MCNX-BLD-06] | gate: production-parfile audit — no non-test parfile activates `TestMCNuX` |

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the `TestMCNuX` idiom literals `SHARES: ADMBaseX`,
`EXTENDS KEYWORD initial_data`, `SHARES: HydroBaseX`, and
`EXTENDS KEYWORD initial_hydro`, and the matrix identity phrase `the closure list`;
that its cited `flesh/` and `CarpetX/` paths resolve; and — the closure this file
exists to carry — that every `[MCNX-XXX-NN]` id declared in any coded spec appears in
the coverage matrix above, that every id in the matrix is declared by its owning spec
(no phantom rows), and that this file's own `[MCNX-VER-NN]` ids are declared in its
Correctness requirements section.

## Implementation freedom

- All `TestMCNuX` and `MCNuX` routine names, source layout, and languages; the
  self-test battery's internal organization and reporting format beyond "one named
  row per check with measured error and pass flag".
- `TestMCNuX` parameter names, defaults, and any additional analytic backgrounds or
  profile parameters — provided the pinned keyword values and formulas of
  [MCNX-VER-02]–[MCNX-VER-04] are among them.
- Benchmark grid extents, box placements, packet counts for the deterministic (non
  `stats-*`) runs, output cadences, and which variables beyond the pinned diagnostics
  are archived — once golden data is committed, these become pinned per benchmark by
  that data.
- Adding benchmarks and matrix rows beyond the binding set (the matrix grows by
  append; removing or weakening a pinned benchmark is a spec change).
- How the non-harness gates are scripted (CI harness, shell, ctest) — their pass/fail
  criteria above are binding, their packaging is not.
- How the analysis routines compute σ (analytic propagation vs documented sample
  estimate) — the choice is documented per check and frozen by its golden data.

## Open questions / assumptions

- **Assumption: `EXTENDS KEYWORD` across repositories (recorded).** No thorn in the
  CarpetX checkout itself uses `EXTENDS KEYWORD`, but the flesh implements it
  (`flesh/lib/sbin/parameter_parser.pl`, range merging and the RESTRICTED-only rule in
  `flesh/lib/sbin/ImpParamConsistency.pl`), the extended keywords are RESTRICTED in
  the base thorns' param.ccl, and the stack's own parameter files reference externally
  extended values (`CarpetX/CarpetX/par/brill-lindquist-write.par` sets
  `ADMBaseX::initial_data = "Brill-Lindquist"`). This spec relies on that mechanism;
  if it regressed in the flesh, `TestMCNuX` would need a fallback (own keyword +
  `"none"` base values), which would be a spec change here.
- **Assumption: the `schwarzschild-pt` drift bound is a design bound.** The
  `1e-5`-relative bound at the pinned resolution (`h = M/16`, `T = 100M`,
  `r ∈ [6M, 10M]`) is set from the expected error budget of an order-≥ 2 integrator on
  a trilinear-gathered analytic metric; the binding, resolution-independent
  requirements are the convergence order and the golden parity. If golden-data
  creation measures a convergent drift above the bound, the pinned resolution (never
  the order requirement) is re-pinned with the measurement documented here.
- **Assumption: statistical golden data is seed-frozen.** The z-score reduction of
  [MCNX-VER-07] deliberately freezes each statistical check at the primary seed; the
  secondary seed `20260730` exists for non-golden robustness cross-checks. A false
  4σ failure at the pinned seed is a fixed property discovered at golden-data
  creation and handled by the documented re-pin protocol of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) — never by
  seed shopping.
- **Assumption: cell-centered sharp sphere surface.** The uniform-sphere profile's
  cell-center evaluation (no volume-fraction smoothing) makes the realized sphere
  resolution-dependent at the one-cell level; benchmarks that integrate over the
  sphere (beam attenuation path lengths) pin their grids so the surface cells are a
  ≤ 1% volume effect, absorbed into the expected-value computation of the analysis
  routine.
- **Open question: GPU-backend golden comparisons.** Golden data is generated
  single-rank CPU (the exact-reproducibility mode). GPU runs of the same benchmarks
  are compared through norm-based tolerances only (unordered atomics,
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)
  [MCNX-RNG-08]); pinning concrete GPU-vs-golden norm tolerances requires measured
  reduction spreads and is deferred until the suite first runs on a GPU backend —
  recorded here so the gap is a decision, not an oversight.
