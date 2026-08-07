# Trapped-regime treatment

> Leaf spec. Self-contained: an agent can implement the optically-thick transport
> scheme — the explicit baseline, the opacity-relabeling extension with its stability
> control, and the diffusion-approximation treatment — from this file alone. It restates
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
- The **diffusion-approximation treatment** of scattering-dominated cells — entered
  after a performed scattering of a step, finishing the step with a two-branch
  displacement draw (causally truncated continuous kernel plus a discrete no-scatter
  atom), a two-leg advection/free-streaming position update, the correlated
  θ₂/B(f_free) momentum model, and a fixed five-draw RNG map.

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
diffusion-approximation scheme is restated with the 2021 paper's §2.6 as the
procedural baseline and the 2018 methods paper supplying the pieces the 2021 paper
cites without restating — the θ₂/B(f_free) momentum model (its Eqs. 34–35) and the
fluid-advection provenance (its Eqs. 31–33); this corpus restates the scheme in
dimensionally explicit geometrized form (spec text normative, per `README.md`). The
θ₂ law and its B(f_free) fit were certified from the typeset arXiv:1708.08452v2 PDF
per the README's escalation path; the verbatim transcription record is
`.tasks/fix-trapped-regime-treatment/06-source-certification-b-ffree.md`.

| Pinned id | Informational gloss |
|---|---|
| arXiv:2103.16588v2 | Foucart et al. 2021, ApJ 920, 82 — SpEC MC implementation paper; Eq. 43 relabeling, Eqs. 44–50 stability analysis, §2.6 (Eqs. 36–42) diffusion procedure. **The binding α convention.** |
| arXiv:2008.08089v3 | Foucart et al. 2020, ApJL 902, L27 — the letter that introduced the relabeling, with the **flipped** α convention (comparison/warning only; never a source of formulas) |
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; its §2.3.6 is the diffusion-approximation scheme restated below (Eqs. 27–30 displacement sampling, Eqs. 31–33 fluid advection, Eqs. 34–35 the θ₂/B(f_free) momentum model, certified from the typeset PDF) |
| arXiv:1806.02349v1 | Foucart et al. 2018, PRD 98, 063007 — cited here only for the α symbol-collision note (its α is a spectral-fit index, a third unrelated meaning) |

## Inputs & outputs

| Quantity | Role | Type / units |
|---|---|---|
| `η(s, E), κ_a(s, E), κ_s(s, E)` | unprimed coefficients from [opacity-eos-evaluation](./opacity-eos-evaluation.md) | per that spec; κ converted to geometrized code units downstream |
| `ξ` | run-time stability-control parameter of the relabeling | dimensionless > 0, typically `ξ = 1` |
| scheme selector | run-time parameter: explicit or relabeled; plus a fixed-α override for verification runs | — |
| `τ_diff` | run-time diffusion-regime threshold (scattering optical depth of the step remainder, [MCNX-TRP-06]) | dimensionless > 0, default 3 |
| `Δt`, `Δt_c` | transport step; cell light-crossing time (minimum cell width, c = 1) | geometrized code units |
| `Δt_rem`, `Δt′_rem` | coordinate-time remainder of the step after a performed scattering ([MCNX-TRP-06]); its fluid-frame counterpart `Δt′_rem = Δt_rem/u^t` | geometrized code units |
| `u^μ`, tetrad | cell fluid four-velocity and the pinned tetrad of [packet-representation-and-sampling](./packet-representation-and-sampling.md) | geometrized |

Outputs: the **active coefficient triple** consumed by emission and event sampling —
`(η, κ_a, κ_s)` under the explicit scheme, `(η′, κ_a′, κ_s′)` under the relabeling —
and, when the diffusion criterion of [MCNX-TRP-06] holds after a performed
scattering, the packet's remainder-of-step outcome: the two-leg position update and
the correlated momentum redraw that finish the step.

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

- **[MCNX-TRP-06] Diffusion-regime entry and the first-scattering prelude (restated,
  normative; provenance arXiv:1708.08452 §2.3.6; arXiv:2103.16588 §2.6).** In every
  cell, whatever its optical depth, a packet's step begins in the ordinary discrete
  event loop of [neutrino-matter-interactions](./neutrino-matter-interactions.md)
  ([MCNX-INT-01..06]) on the active coefficients. After **each performed scattering
  event** ([MCNX-INT-04]), the diffusion criterion is evaluated in the cell where
  that scattering occurred, over the time remaining in the step:

  ```text
  κ_s′ Δt′_rem > τ_diff ,        Δt′_rem = Δt_rem/u^t
  ```

  with `κ_s′` the **active** scattering opacity of that cell at the packet's
  fluid-frame energy ν (the relabeled value when the relabeling is selected; under
  the explicit scheme the active value is the unprimed κ_s — pinning the active
  coefficient in the criterion is a recorded deviation below), `Δt_rem` the
  coordinate time remaining in the step, and `Δt′_rem` its fluid-frame counterpart.
  If the criterion holds, the packet's step is finished by **exactly one diffusion
  episode** over `Δt_rem`: displacement per [MCNX-TRP-08], position update per
  [MCNX-TRP-09], outgoing momentum per [MCNX-TRP-10], draws per [MCNX-TRP-07] — the
  post-scatter sampling episode of [MCNX-INT-02] does not begin; the diffusion
  episode is the packet's next episode. If it does not hold, the discrete loop
  continues normally (new sampling episode per [MCNX-INT-02]). `τ_diff` is a runtime
  parameter with default 3 (the 2021 paper's fixed value; the 2018 paper leaves it a
  runtime parameter, this corpus's form). The active coefficients, `u^μ`, and the
  tetrad of the scattering cell are pinned for the whole remainder (recorded
  deviation below). Consequences: at most one diffusion episode per (packet, step) —
  the episode consumes all remaining step time — and a packet that performs no
  scattering during the step never enters the diffusion regime; packets created
  mid-step participate through the discrete loop from their creation time per
  [MCNX-INT-02].

- **[MCNX-TRP-08] Two-branch displacement law (restated, normative; provenance
  arXiv:1708.08452 Eqs. 27–30; arXiv:2103.16588 Eqs. 36–42).** A single uniform draw
  `P` samples the packet's fluid-frame displacement magnitude `r_d` over `Δt′_rem`
  from the papers' full displacement distribution — a normalized continuous
  diffusion kernel truncated at the causal bound, plus a discrete no-scatter outcome
  whose mass is exactly the zero-scattering probability:

  ```text
  p_ns = exp(−κ_s′ Δt′_rem)           (probability of zero scatterings over Δt′_rem)

  f(r̃) = (4/√π) r̃² exp(−r̃²) ,      r̃ = r_d √(3κ_s′/(4Δt′_rem))
  F(x) = ∫₀^x f dr̃ = erf(x) − (2/√π) x exp(−x²)

  P < 1 − p_ns  →  continuous branch: solve  F(r̃(P)) = P · F(r̃_max) / (1 − p_ns)
  P ≥ 1 − p_ns  →  no-scatter atom:   r_d = Δt′_rem    (r̃ = r̃_max , f_free = 1)

  r̃_max = √(3κ_s′ Δt′_rem/4)         (⇔ r_d ≤ Δt′_rem)
  ```

  The branch is decided by comparison on the one draw `P` — no second draw — and no
  sampled displacement exceeds the light-travel bound `r_d ≤ Δt′_rem`. The atom is
  the papers' no-scatter/free-streaming branch: the probability mass beyond
  `1 − p_ns` is exactly the packets that experience no scattering over the
  remainder, consistent with the discrete loop the episode hands off from.
  Presentation note: the sources leave F as an implicit integral to be tabulated
  (2018 Eqs. 29–30; 2021 Eq. 41); the closed form above is its exact antiderivative
  — the same function, no probability changed, pinned here in closed form. Numerical
  inversion of F on the continuous branch is admissible with absolute error ≤ 1e-12
  in r̃.

- **[MCNX-TRP-09] Two-leg position update (restated, normative; provenance
  arXiv:1708.08452 Eqs. 31–33, restated in 4-velocity form — recorded deviation
  below; the coordinate-time split is likewise a recorded deviation).** The sampled
  displacement fixes the free-streaming fraction of the remainder,

  ```text
  f_free = r_d / Δt′_rem              (∈ [0, 1]; the no-scatter atom gives f_free = 1)
  ```

  and the remainder splits into two legs, executed in this order:

  ```text
  advection leg        x^i ← x^i + (u^i/u^t)(1 − f_free) Δt_rem
                       (the packet is stationary in the fluid frame; this is the
                       coordinate distance the fluid covers in that time)
  free-streaming leg   geodesic push over the coordinate time f_free Δt_rem along
                       the direction n̂_f, drawn isotropically in the fluid frame at
                       fixed ν (cosθ_f = 2u₁ − 1, φ_f = 2πu₂, with u₁, u₂ the
                       k = 1, 2 draws of [MCNX-TRP-07]), tetrad-mapped per
                       [packet-representation-and-sampling](./packet-representation-and-sampling.md)
                       [MCNX-PKT-04] and advanced by the pusher of
                       [geodesic-propagation](./geodesic-propagation.md) with
                       Δt_push = f_free Δt_rem
  ```

  Both legs consume coordinate time and their durations sum to exactly `Δt_rem`, so
  the packet lands exactly on the step boundary (recorded deviation below). The
  limits are exact: `f_free = 1` (the atom) is a pure geodesic push over `Δt_rem`;
  `f_free → 0` is pure fluid advection.

- **[MCNX-TRP-10] Correlated momentum model (restated, normative; provenance
  arXiv:1708.08452 Eqs. 34–35, certified from the typeset PDF — see Source of
  truth).** The packet's momentum at the end of the diffusion episode is drawn
  **correlated with the free-streaming direction**, at fixed fluid-frame energy ν:
  in a spherical polar system whose polar axis is `n̂_f`, the outgoing fluid-frame
  direction sits at polar angle θ₂ and azimuth φ₂, with `φ₂` uniform in [0, 2π) (by
  symmetry), `r` uniform in [0, 1), and

  ```text
  cos θ₂ = B − (1 + B) · ((B−1)/(B+1))^r ,        B = B(f_free)
  ```

  (the paper's `exp[r ln((B − 1)/(B + 1))]` written in the algebraically identical
  power form — the same function for every tabulated B > 1; a presentation note in
  the same class as [MCNX-TRP-08]'s closed-form F). `B(f_free)` is the 2018 paper's
  calibrated fit, restated in full: **linear interpolation in f_free** through the
  21 tabulated values `B_i = B(0.05i)`, `i = 0, 1, ..., 20`,

  ```text
  B_i = ( 1000, 18.74, 7.52, 4.75, 3.51, 2.78, 2.32, 2.00,
          1.77, 1.60, 1.47, 1.36, 1.28, 1.21, 1.15, 1.10,
          1.07, 1.04, 1.019, 1.0027, 1.0000001 )
  ```

  with the fit's boundary anchors `B(0) → ∞` (isotropy at f_free = 0; the paper's
  B_0 = 1000 is its finite stand-in) and `B(1) = 1` (forward at f_free = 1; the
  paper's B_20 = 1.0000001 sits just above 1, keeping the law finite over the whole
  table). **Atom branch pinned by rule**: at `f_free = 1` (the no-scatter atom of
  [MCNX-TRP-08]) the outgoing direction is `n̂_f` **exactly** — `cos θ₂ = 1`, the
  azimuth irrelevant — with the draws `r` and `φ₂` consumed but inert
  ([MCNX-TRP-07]). The paper instead samples its regularized endpoint B_20; the pin
  is a corpus decision, noted here inline: it makes "never scattered" mean exactly
  free streaming and avoids legislating the floating-point corner of
  `((B−1)/(B+1))^r` at B = 1 (IEEE `pow(0, 0) = 1` would flip cos θ₂ to −1 at the
  reachable draw r = 0). The outgoing fluid-frame momentum `ν(1, n̂_out)` is
  tetrad-transformed per
  [packet-representation-and-sampling](./packet-representation-and-sampling.md)
  [MCNX-PKT-04]; weight and species are unchanged, and the null and ν-recovery
  identities of [MCNX-INT-04] hold to machine tier (`~1e-14` relative). The net
  momentum exchange of the diffusion episode is tallied as **one effective
  scattering event** per
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) (corpus-defined;
  the papers' diffusion sections are silent on the tally).

- **[MCNX-TRP-07] RNG consumption map (trapped regime).** Relabeled event sampling —
  the discrete prelude of [MCNX-TRP-06] — consumes draws exactly per
  [neutrino-matter-interactions](./neutrino-matter-interactions.md) [MCNX-INT-06]
  (primed coefficients change values, never the draw structure). A diffusion episode
  consumes exactly **five draws, unconditionally**: increment the packet's event
  counter `e` by 1, then

  ```text
  k = 0 → P                         (displacement draw / branch selector, [MCNX-TRP-08])
  k = 1 → cosθ_f ,  k = 2 → φ_f     (free-streaming direction n̂_f, [MCNX-TRP-09])
  k = 3 → r ,       k = 4 → φ₂      (θ₂ law and azimuth about n̂_f, [MCNX-TRP-10])
  ```

  All five draws are consumed in **both** branches of [MCNX-TRP-08] — the no-scatter
  atom still needs `n̂_f`, and `r` and `φ₂` are consumed but inert per
  [MCNX-TRP-10] — in exactly this order, and no other draws are consumed. RNG
  consumption therefore never diverges across packets within a diffusion episode.

### Conventions restated (the subset this leaf uses)

- **α in this spec is exclusively the trapped-regime relabeling fraction** (the 2021
  convention); the lapse never appears in this spec's formulas (the fluid time
  dilation is written via `u^t`). The spectral-fit α of arXiv:1806.02349 is a third,
  unrelated symbol ([conventions-and-units](./conventions-and-units.md) carries the
  three-way collision warning).
- Primes on coefficients (η′, κ_a′, κ_s′) denote the relabeled values — Unicode
  U+2032 PRIME, as throughout this corpus; minus signs in formulas are U+2212.
- **Step-remainder notation (owned by this spec)**: `Δt_rem` is the coordinate time
  remaining in the current transport step when the first-scattering prelude hands
  off ([MCNX-TRP-06]); `Δt′_rem = Δt_rem/u^t` is its fluid-frame (proper-time)
  counterpart — the same `u^t` time-dilation convention as `Δt′ = Δt/u^t`.
- Geometrized units: κ in code units (`κ[code] = κ[cm⁻¹] × 1.476625038e5`), Δt, Δt_c,
  Δt_rem, Δt′_rem, r_d in code units, c = 1; coefficient lookups take E = ν in MeV.
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
- **Diffusion sampling law (`stats-diffusion`, statistical, 4σ — three clauses from
  one run).** At the benchmark's pinned `κ_s′ Δt′_rem`: (1) the no-scatter atom's
  frequency is within 4σ of `exp(−κ_s′ Δt′_rem)`; (2) the continuous branch's
  sampled r̃ first and second moments are within 4σ of the analytic moments of the
  normalized truncated distribution `f(r̃)/F(r̃_max)` on `[0, r̃_max)`; (3) the
  sampled cos θ₂ distribution, conditioned on f_free bins, is within 4σ of the
  restated θ₂ law of [MCNX-TRP-10] with the interpolated B(f_free).
- **Scheme identities (exact / machine tier `~1e-14`).** In a pinned deterministic
  diffusion run: every atom-branch episode ends with the outgoing direction equal to
  `n̂_f` and `cos θ₂ = 1` exactly; at `f_free = 1` the position update is a pure
  geodesic push over `Δt_rem` (advection leg exactly zero); as `f_free → 0` the
  update reduces to pure fluid advection `x^i ← x^i + (u^i/u^t) Δt_rem` (machine
  tier); and the two legs' coordinate times sum to `Δt_rem` exactly, so every packet
  ends the step on the step boundary.
- **Overlap consistency (`trp-overlap-consistency`, statistical, 4σ).** In a
  moderate-depth cell (`κ_s′ Δt′_rem` near `τ_diff`), runs with the threshold set
  just above and just below produce end-of-step packet position and direction
  distributions agreeing within 4σ. Both populations share the identical discrete
  prelude up to and including the first performed scattering, so they differ only in
  how the step remainder is finished — the discrete event loop and the diffusion
  episode must join smoothly where both are valid.
- **Draw-map audit (exact, constant).** A `(q, e, k)` trace of a pinned
  diffusion-regime run matches [MCNX-TRP-07] exactly: exactly five draws per
  diffusion episode, both branches, in the pinned order, with the discrete prelude's
  episodes audited per [MCNX-INT-06].

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the 2021-convention relabeling literals `η′ = αη`, `κ_a′ = ακ_a`, and
`κ_s′ = κ_s + (1−α)κ_a`, the stability-control string `κ_a′Δt`, the flip-warning
phrase `the letter's α is this corpus's 1−α`, and the limit phrase
`α → 1 reproduces the explicit scheme`; contains the diffusion-scheme needles
`κ_s′ Δt′_rem > τ_diff`, `k = 0 → P`, `F(r̃(P)) = P · F(r̃_max) / (1 − p_ns)`,
`f_free = r_d / Δt′_rem`, and `((B−1)/(B+1))^r` — all with the exact codepoints
(`′` U+2032, `−` U+2212, `→` U+2192, `·` U+00B7); and that its `[MCNX-TRP-NN]` ids
are declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- How and where the primed coefficients are computed (fused into the coefficient
  evaluation, cached per cell, recomputed per packet) — observable values must equal
  the pure post-map of [MCNX-TRP-02].
- The F-inversion method for the diffusion CDF (Newton, bisection, precomputed
  monotone interpolant) within the stated 1e-12 accuracy; caching of `r̃_max`,
  `F(r̃_max)`, and `p_ns` per cell.
- Storage and lookup of the `B_i` table (constant array, precomputed interpolation
  slopes) — the 21 values and the linear interpolation of [MCNX-TRP-10] are binding.
  Likewise the azimuth reference direction for `φ₂` about `n̂_f` (φ₂ is uniform, so
  any fixed convention is statistically identical; the pinned choice is frozen per
  benchmark by its golden data).
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
  corpus restates it dimensionally explicitly (`r̃ = r_d √(3κ_s′/(4Δt′_rem))`, r_d
  and Δt′_rem both in code units, c = 1, κ_s′ the active scattering opacity), which
  reproduces the standard random-walk variance `⟨r_d²⟩ = 6DΔt′_rem = 2Δt′_rem/κ_s′`
  (D = 1/(3κ_s′), c = 1) and the causal cap r_d ≤ Δt′_rem. Per `README.md`, the
  restated form is normative; a PDF spot-check discrepancy updates this spec.
- **Resolution logged: coordinate-time split and step-boundary landing (deviation
  from the source paper's letter, spec text normative per `README.md`).** The 2018
  paper converts the free-streaming leg's fluid-frame duration to simulation time
  via `Δt_free = f_free Δt [p^t/(ν u^t)]` and states that, because of this
  weighting, a packet "may be propagated to a time larger or smaller than initially
  requested." That per-packet time drift is incompatible with this corpus's
  architecture: packet state carries no time coordinate, the container contract
  places every packet at its owning cell after every step
  ([particle-container-and-gpu](./particle-container-and-gpu.md) [MCNX-GPU-05]), and
  the ledger closure of
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) [MCNX-HYD-05]
  balances synchronized end-of-step snapshots over one global Δt. [MCNX-TRP-09]
  instead applies the split to the coordinate-time remainder — the free-streaming
  leg consumes exactly `f_free Δt_rem` — so every packet lands exactly on the step
  boundary. The two schemes agree in the comoving limit and differ only by the
  relativistic aberration weighting on the free-streaming leg.
- **Resolution logged: coefficients and tetrad pinned to the scattering cell
  (deviation from the source paper's letter, spec text normative per `README.md`).**
  The papers do not address cell crossings during the diffusion displacement. This
  corpus pins the active coefficients, `u^μ`, and the tetrad of the cell where the
  triggering scattering occurred for the entire remainder ([MCNX-TRP-06]): a
  displacement that carries the packet into a cell of the other regime (or of
  different coefficients) completes under the scattering cell's rules. The error is
  O(cell size), the same tier as the diffusion approximation's own locality
  assumption; a sub-step re-evaluation would be a spec change.
- **Resolution logged: advection leg restated in 4-velocity form (deviation from the
  source paper's letter, spec text normative per `README.md`).** The 2018 paper
  implements fluid advection by constructing a null momentum `p^μ = A t^μ + B u^μ`
  (its Eqs. 31–33, with a singular-fluid guard at small `δ_ij u^i u^j`) and pushing
  along it for a derived time so the packet covers exactly the coordinate distance
  the fluid covers. [MCNX-TRP-09] restates the leg directly as that transport —
  `x^i ← x^i + (u^i/u^t)(1 − f_free) Δt_rem` — the same displacement with no
  null-momentum construction and no singular limit (`u^i/u^t` is regular at
  `u^i = 0`). Same style as the dimensional-restatement entry above: the sources'
  coordinate algebra is replaced by the dimensionally explicit form it encodes.
- **Resolution logged: active κ_s′ in the regime criterion (deviation from the
  source paper's letter, spec text normative per `README.md`).** The papers write
  the diffusion criterion with plain κ_s; neither links the Eq. 43 relabeling to the
  criterion. This corpus pins the **active** scattering opacity in [MCNX-TRP-06]:
  the relabeled component (1−α)κ_a stands in for absorption + re-emission that traps
  packets exactly as scattering does, the relabeling preserves the total interaction
  rate (`κ_a′ + κ_s′ = κ_a + κ_s`, [MCNX-TRP-02]), and the criterion must see the
  same coefficient that drives the displacement sampling. Under the explicit scheme
  (or α = 1) the two readings coincide.
- **Assumption: relabeling and lepton number.** Relabeled "scattering" (the (1−α)κ_a
  component) exchanges no lepton number, while the absorption + re-emission it stands
  in for would exchange lepton number with zero net expectation at equilibrium — the
  scheme's validity conditions (well inside the neutrinosphere, near-equilibrium)
  are what make this admissible; the equilibration benchmark is the arbiter.
