# Trapped-regime treatment

> Leaf spec. Self-contained: an agent can implement the optically-thick transport
> scheme — the explicit baseline, the opacity-relabeling extension with its stability
> control, and the diffusion-approximation advection — from this file alone. It restates
> only the conventions it uses and references
> [conventions-and-units](./conventions-and-units.md) (units, species),
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (draws, 4σ bar),
> [opacity-eos-evaluation](./opacity-eos-evaluation.md) (the unprimed coefficients), and
> [neutrino-matter-interactions](./neutrino-matter-interactions.md) (the event machinery
> the relabeling feeds) for the rest. `README.md` is canonical if any restated
> convention here conflicts with it.

## Purpose & scope

In optically thick regions the explicit event loop becomes both stiff (absorption
timescale ≪ Δt couples radiation and fluid explosively at finite step size) and
expensive (packet event counts scale with optical depth). This spec pins MCNuX's
treatment:

- The **explicit scheme as the baseline correctness contract**: unmodified
  (η, κ_a, κ_s) driving the emission law of
  [packet-representation-and-sampling](./packet-representation-and-sampling.md) and the
  event machinery of [neutrino-matter-interactions](./neutrino-matter-interactions.md).
- The **opacity-relabeling extension** (arXiv:2103.16588 Eq. 43, restated below in the
  **binding 2021 α convention**), selected by parameter, with the boxed warning against
  the flipped convention of arXiv:2008.08089.
- The **stability bound** (arXiv:2103.16588 Eqs. 48–49, restated) and the practical
  run-time control `κ_a′Δt_c ≤ ξ` with its binding per-cell α selection rule.
- The named **limit-check requirement** that **α → 1 reproduces the explicit scheme**.
- The **diffusion-approximation advection** replacing discrete scattering in
  scattering-dominated cells, with its sampling law, causal truncation, and RNG map.

Out of scope:

- The unprimed coefficients η, κ_a, κ_s themselves (tabulated assembly, νx mapping,
  analytic mode) — [opacity-eos-evaluation](./opacity-eos-evaluation.md); the
  relabeling is a pure post-map on that interface's outputs.
- The event-sampling law, episode structure, and event outcomes —
  [neutrino-matter-interactions](./neutrino-matter-interactions.md); the relabeling
  changes only which coefficient values that machinery consumes.
- The emission count law — [packet-representation-and-sampling](./packet-representation-and-sampling.md);
  the relabeling substitutes η′ for η in it, nothing else.
- Source-term tallies of the resulting events —
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md).
- Geodesic propagation ([geodesic-propagation](./geodesic-propagation.md)); device
  execution ([particle-container-and-gpu](./particle-container-and-gpu.md)); concrete
  equilibration benchmarks and golden data
  ([verification-suite-design](./verification-suite-design.md)).

## Source of truth

The relabeling, its stability analysis, and the practical control are restated in full
below (normative per `README.md`) from the pinned 2021 implementation paper. The
diffusion-approximation sampling law is restated from the methods paper with the 2021
paper's equivalent forms as co-provenance; this corpus restates it in dimensionally
explicit geometrized form (spec text normative, per `README.md`).

| Pinned id | Informational gloss |
|---|---|
| arXiv:2103.16588v2 | Foucart et al. 2021, ApJ 920, 82 — SpEC MC implementation paper; Eq. 43 relabeling, Eqs. 44–50 stability analysis, Eqs. 40–42 diffusion sampling. **The binding α convention.** |
| arXiv:2008.08089v3 | Foucart et al. 2020, ApJL 902, L27 — the letter that introduced the relabeling, with the **flipped** α convention (comparison/warning only; never a source of formulas) |
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; Eqs. 27–30 are the diffusion-approximation displacement sampling restated below |
| arXiv:1806.02349v1 | Foucart et al. 2018, PRD 98, 063007 — cited here only for the α symbol-collision note (its α is a spectral-fit index, a third unrelated meaning) |

## Inputs & outputs

| Quantity | Role | Type / units |
|---|---|---|
| `η(s, E), κ_a(s, E), κ_s(s, E)` | unprimed coefficients from [opacity-eos-evaluation](./opacity-eos-evaluation.md) | per that spec; κ converted to geometrized code units downstream |
| `ξ` | run-time stability-control parameter of the relabeling | dimensionless > 0, typically `ξ = 1` |
| scheme selector | run-time parameter: explicit or relabeled; plus a fixed-α override for verification runs | — |
| `τ_diff` | run-time diffusion-regime threshold (per-step scattering optical depth) | dimensionless > 0 |
| `Δt`, `Δt_c` | transport step; cell light-crossing time (minimum cell width, c = 1) | geometrized code units |
| `u^μ`, tetrad | cell fluid four-velocity and the pinned tetrad of [packet-representation-and-sampling](./packet-representation-and-sampling.md) | geometrized |

Outputs: the **active coefficient triple** consumed by emission and event sampling —
`(η, κ_a, κ_s)` under the explicit scheme, `(η′, κ_a′, κ_s′)` under the relabeling —
and, in diffusion-regime cells, the packet's step outcome (absorbed, or advected and
redirected) replacing the discrete-scattering event loop.

## Correctness requirements

- **[MCNX-TRP-01] Explicit scheme is the baseline.** With the relabeling deselected
  (or α = 1 everywhere), MCNuX transport is exactly the explicit scheme: emission per
  [packet-representation-and-sampling](./packet-representation-and-sampling.md) with
  the unprimed η, events per
  [neutrino-matter-interactions](./neutrino-matter-interactions.md) with the unprimed
  κ_a, κ_s. This baseline is a complete, correct (if expensive and conditionally
  stable) configuration; every extension in this spec is measured against it.

- **[MCNX-TRP-02] The opacity relabeling (restated, normative; provenance
  arXiv:2103.16588 Eq. 43; binding 2021 α convention).** When the relabeling is
  selected, the active coefficients are, per (cell, species, energy bin),

  ```text
  η′ = αη ,   κ_a′ = ακ_a ,   κ_s′ = κ_s + (1−α)κ_a ,        α ∈ (0, 1]
  ```

  substituted for (η, κ_a, κ_s) in the emission count law and the event-sampling law —
  and nowhere else; tallies remain the actual packet events. Physical basis (recorded):
  absorption followed by re-emission at equilibrium is practically indistinguishable
  from elastic scattering, so a fraction (1−α) of absorption is relabeled as
  scattering. Two identities are exact consequences and must hold to machine tier
  (`~1e-14` relative) at every evaluated state:

  ```text
  η′/κ_a′ = η/κ_a            (equilibrium energy density unmodified, κ_a > 0)
  κ_a′ + κ_s′ = κ_a + κ_s    (total interaction rate / diffusion timescale unmodified)
  ```

- **[MCNX-TRP-03] α-convention guard.**

  > **WARNING — flipped α convention in the literature (the corpus's
  > highest-consequence trap).** The 2020 letter arXiv:2008.08089 (its Eq. 4) uses the
  > same symbol α with the **opposite** meaning: **the letter's α is this corpus's
  > 1−α**. The 2021 paper does not remark on the change. This corpus pins the 2021
  > convention (arXiv:2103.16588 Eq. 43) everywhere; never take an α formula from the
  > letter without translating α ↦ 1−α.
  >
  > | Quantity | arXiv:2008.08089 Eq. 4 (do **not** use) | arXiv:2103.16588 Eq. 43 (**binding**) |
  > |---|---|---|
  > | Modified emissivity | η′ = (1−α)η | η′ = αη |
  > | Modified absorption | κ_a′ = (1−α)κ_a | κ_a′ = ακ_a |
  > | Modified scattering | κ_s′ = κ_s + ακ_a | κ_s′ = κ_s + (1−α)κ_a |
  > | Explicit-scheme limit | α → 0 | α → 1 |
  > | Stability inequality direction | α **>** bound (its Eq. 5) | α **<** bound (Eq. 49) |

  Binding requirement: every α-dependent quantity in MCNuX (parameters, diagnostics,
  formulas, documentation) uses the 2021 convention — α is the **retained** fraction:
  increasing α toward 1 means *less* relabeling. Monotonicity signature (checkable):
  for fixed unprimed coefficients, `κ_a′` is increasing and `κ_s′` decreasing in α.

- **[MCNX-TRP-04] Stability bound and the practical control (restated, normative;
  provenance arXiv:2103.16588 Eqs. 44, 48–50).** The forward-Euler stability analysis
  of the coupled radiation–fluid exchange bounds the relabeled absorption per step:

  ```text
  κ_a′Δt ≤ 1/(1 + β̃)                       (Eq. 48)
  α < [1/(1 + β̃)] · [1/(κ_a Δt)]           (Eq. 49)
  β = max( |du_ν/du_fl| , m_p c² |dn_ν/d(ρYₑ)_fl| )    (Eq. 44)
  ```

  with `u_ν = η/κ_a` the equilibrium neutrino energy density, `β̃` the per-channel
  stiffness (`du_ν/du_fl` for energy at constant (ρ, Yₑ); `m_p dn_ν/d(ρYₑ)` for lepton
  number at constant (ρ, T)), and the paper's imposed run-time form
  `α < [1/(1 + β)] · [1/(2C κ_a Δt_c)]` (Eq. 50, safety factor C ≥ 1). The **binding
  practical control**, per the paper's stated practice: α is selected per (cell,
  species, energy bin) as

  ```text
  α = min(1, ξ/(κ_a Δt_c))        so that κ_a′Δt_c ≤ ξ by construction
  ```

  with `Δt_c` the cell light-crossing time and `ξ` a runtime parameter, typically
  `ξ = 1`. In optically thin cells this yields α = 1 (no relabeling) automatically. A
  fixed-α override parameter (verification use, [MCNX-TRP-05]) may replace the
  selection rule; `κ_a′Δt_c ≤ ξ` must hold as a per-cell diagnostic whenever the
  selection rule is active. Enforcement of the β-dependent Eq. 50 bound is an open
  question below.

- **[MCNX-TRP-05] α → 1 limit check (named requirement; the built-in limit check of
  the design).** **α → 1 reproduces the explicit scheme.** At α = 1 the relabeling is
  the identity (η′ = η, κ_a′ = κ_a, κ_s′ = κ_s, exactly), and a run of the relabeled
  scheme with α pinned to 1 (fixed-α override) must reproduce the explicit scheme's
  run at the same fixed seed **bitwise** on a single rank: identical per-packet event
  sequences, identical packet populations after every step, identical golden outputs
  at the harness tier (`ABSTOL=RELTOL=1e-12`). This is exact-tier equivalence, not
  statistical agreement: at α = 1 both configurations must consume identical draws.

- **[MCNX-TRP-06] Diffusion-approximation advection (restated, normative; provenance
  arXiv:1708.08452 Eqs. 27–30; arXiv:2103.16588 Eqs. 40–42).** In cells whose per-step
  scattering optical depth exceeds the threshold — `κ_s Δt′ > τ_diff`, with `κ_s` here
  the **active** scattering opacity (κ_s′ under the relabeling) and `Δt′ = Δt/u^t` the
  fluid proper time elapsed over the step — the discrete-scattering event loop and the
  geodesic push are replaced, once per (packet, step), by:

  1. **Absorption test**: the packet is absorbed iff `u_a < 1 − exp(−κ_a Δt′)` (active
     κ_a at the packet's ν; a diffusing packet's fluid-frame path length over Δt′ is
     Δt′, c = 1). Absorbed packets are removed and tallied at the cell per
     [neutrino-matter-interactions](./neutrino-matter-interactions.md) /
     [hydro-coupling-source-terms](./hydro-coupling-source-terms.md).
  2. **Displacement sampling**: the surviving packet's fluid-frame displacement
     magnitude `r_d` over Δt′ is drawn from the diffusion distribution

     ```text
     f(r̃) = (4/√π) r̃² exp(−r̃²) ,      r̃ = r_d √(3κ_s/(4Δt′))
     F(x) = ∫₀^x f dr̃ = erf(x) − (2/√π) x exp(−x²)
     ```

     by inverting the **causally truncated** CDF: draw P uniform and solve

     ```text
     F(r̃(P)) = P · F(r̃_max) ,      r̃_max = √(3κ_s Δt′/4)    (⇔ r_d ≤ Δt′)
     ```

     so no sampled displacement exceeds the light-travel bound. The displacement
     direction `n̂` is drawn isotropically in the fluid frame.
  3. **Position update**: `x^i ← x^i + (u^i/u^t) Δt + ê^i_(j′) n̂^{j′} r_d` —
     advection with the fluid plus the fluid-frame displacement mapped through the
     spatial legs of the pinned tetrad of
     [packet-representation-and-sampling](./packet-representation-and-sampling.md).
  4. **Momentum redraw**: the outgoing momentum is redrawn isotropically in the fluid
     frame at fixed fluid-frame energy ν and tetrad-transformed — the elastic-scattering
     outcome of [neutrino-matter-interactions](./neutrino-matter-interactions.md),
     applied once. The net momentum exchange of the step is tallied as one effective
     scattering event.

  The regime decision uses the packet's cell at the start of the step and holds for the
  whole step. Numerical inversion of F is admissible with absolute error ≤ 1e-12 in r̃.

- **[MCNX-TRP-07] RNG consumption map (trapped regime).** Relabeled event sampling
  consumes draws exactly per
  [neutrino-matter-interactions](./neutrino-matter-interactions.md) [MCNX-INT-06]
  (primed coefficients change values, never the draw structure). A diffusion-regime
  step consumes one episode: increment the packet's event counter `e` by 1, then

  ```text
  k = 0 → u_a (absorption test)
  k = 1 → P  (radial CDF)            k = 2 → cosθ_d ,  k = 3 → φ_d   (displacement direction)
  k = 4 → cosθ_p ,  k = 5 → φ_p     (outgoing momentum direction)
  ```

  with `k = 1…5` consumed only if the packet survives step 1, in exactly this order,
  and no other draws.

### Conventions restated (the subset this leaf uses)

- **α in this spec is exclusively the trapped-regime relabeling fraction** (the 2021
  convention); the lapse never appears in this spec's formulas (the fluid time
  dilation is written via `u^t`). The spectral-fit α of arXiv:1806.02349 is a third,
  unrelated symbol ([conventions-and-units](./conventions-and-units.md) carries the
  three-way collision warning).
- Primes on coefficients (η′, κ_a′, κ_s′) denote the relabeled values — Unicode
  U+2032 PRIME, as throughout this corpus; minus signs in formulas are U+2212.
- Geometrized units: κ in code units (`κ[code] = κ[cm⁻¹] × 1.476625038e5`), Δt, Δt_c,
  Δt′, r_d in code units, c = 1; coefficient lookups take E = ν in MeV.
- Species `{νe, ν̄e, νx}` indices (0, 1, 2); fluid-frame energy `ν = −p_μ u^μ`.
- All draws are `u(S, q, e, k)` evaluations in `[0, 1)`; statistical checks use the 4σ
  bar at pinned counts and seeds of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md).

## Verification

- **Relabeling identities (machine tier, `~1e-14`).** Over a pinned sweep of states
  and α values: the computed (η′, κ_a′, κ_s′) equal the restated formulas, and the two
  invariants of [MCNX-TRP-02] hold to `1e-14` relative.
- **Convention tripwire (exact).** For fixed unprimed coefficients and α₁ < α₂ in
  (0, 1]: `κ_a′(α₁) < κ_a′(α₂)` and `κ_s′(α₁) > κ_s′(α₂)` — the monotonicity
  signature of the 2021 convention; the flipped convention reverses both and fails
  this check.
- **Control enforcement (exact).** In every step of every relabeled benchmark, the
  per-cell diagnostic verifies `κ_a′Δt_c ≤ ξ` and (selection rule active) reproduces
  `α = min(1, ξ/(κ_a Δt_c))` bitwise from the cell's κ_a.
- **α → 1 limit (exact tier / golden parity).** The [MCNX-TRP-05] check: relabeled
  run at fixed α = 1 vs explicit run, same seed, single rank — bitwise-identical
  per-packet event logs and packet populations; golden outputs identical at
  `ABSTOL=RELTOL=1e-12`.
- **Equilibration box (statistical, 4σ).** In the uniform equilibration-box benchmark
  ([verification-suite-design](./verification-suite-design.md)): the relabeled scheme
  at α < 1 reaches the same equilibrium spectrum and energy density as the explicit
  scheme within the 4σ bar, with a measured reduction in event count; the equilibrium
  energy density matches `η/κ_a` (the invariant the relabeling preserves).
- **Diffusion displacement distribution (statistical, 4σ).** Sampled r̃ over many
  draws: first and second moments within 4σ of the truncated-`f` analytic moments at
  the benchmark's `κ_s Δt′`.
- **Overlap consistency (statistical, 4σ).** In a moderate-depth cell (`κ_s Δt′` near
  `τ_diff`), runs with the threshold set just above and just below produce
  end-of-step packet position and direction distributions agreeing within 4σ — the
  discrete event loop and the diffusion approximation must join smoothly where both
  are valid.
- **Draw-map audit (exact).** A `(q, e, k)` trace of a pinned diffusion-regime step
  matches [MCNX-TRP-07] exactly (one draw for an absorbed packet, six for a survivor).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the 2021-convention relabeling literals `η′ = αη`, `κ_a′ = ακ_a`, and
`κ_s′ = κ_s + (1−α)κ_a`, the stability-control string `κ_a′Δt`, the flip-warning
phrase `the letter's α is this corpus's 1−α`, and the limit phrase
`α → 1 reproduces the explicit scheme` — all with the exact codepoints (`′` U+2032,
`−` U+2212); and that its `[MCNX-TRP-NN]` ids are declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- How and where the primed coefficients are computed (fused into the coefficient
  evaluation, cached per cell, recomputed per packet) — observable values must equal
  the pure post-map of [MCNX-TRP-02].
- The F-inversion method for the diffusion CDF (Newton, bisection, precomputed
  monotone interpolant) within the stated 1e-12 accuracy; caching of `r̃_max` and
  `F(r̃_max)` per cell.
- Parameter names and defaults for ξ, τ_diff, the scheme selector, and the fixed-α
  override (the *semantics* above are binding; `ξ = 1` typical is guidance, not a
  frozen default).
- Whether diffusion-regime packets are processed in the same kernel as free-streaming
  packets or a separate pass, and all batching/fusion choices — subject to the draw
  maps.

## Open questions / assumptions

- **Assumption: α range (recorded finding — the source states no range).** Two
  independent full-text passes of arXiv:2103.16588 found **no explicit α range**
  (such as 0 < α ≤ 1) anywhere; α is constrained only implicitly by Eq. 49/50 and
  three qualitative validity conditions, and the paper reports α = 0.9 in production.
  This corpus **assumes and enforces α ∈ (0, 1]** (paramcheck-rejected otherwise):
  α > 1 would drive κ_s′ below κ_s (negative relabeled scattering contribution) and
  α ≤ 0 would extinguish emission entirely — both outside the scheme's physical basis.
- **Open question: enforcing the β-dependent bound (Eqs. 49–50).** The baseline
  enforces only the practical control `κ_a′Δt_c ≤ ξ` (the paper's stated practice).
  Computing β per Eq. 44 requires derivatives of tabulated η and κ_a with respect to
  the fluid state; consuming those (via the WeakLibInterp differentiate kernels) and
  enforcing Eq. 50 with a safety factor C would be a purely additive extension of
  [MCNX-TRP-04]. Until then, stiff configurations rely on choosing ξ conservatively.
- **Assumption: dimensional restatement of the diffusion scaling.** The sources write
  the displacement scaling in forms tied to their internal unit conventions; this
  corpus restates it dimensionally explicitly (`r̃ = r_d √(3κ_s/(4Δt′))`, r_d and Δt′
  both in code units, c = 1), which reproduces the standard random-walk variance
  `⟨r_d²⟩ = 6DΔt′ = 2Δt′/κ_s` (D = 1/(3κ_s), c = 1) and the causal cap r_d ≤ Δt′. Per
  `README.md`, the restated form is normative; a PDF spot-check discrepancy updates
  this spec.
- **Assumption: one diffusion episode per step.** The regime decision and the
  coefficients are frozen at step start; a packet whose sampled displacement carries it
  into a cell of the other regime completes the step under the starting cell's rules.
  Accepted at baseline; a sub-step hybrid scheme would be a spec change.
- **Assumption: relabeling and lepton number.** Relabeled "scattering" (the (1−α)κ_a
  component) exchanges no lepton number, while the absorption + re-emission it stands
  in for would exchange lepton number with zero net expectation at equilibrium — the
  scheme's validity conditions (well inside the neutrinosphere, near-equilibrium)
  are what make this admissible; the equilibration benchmark is the arbiter.
