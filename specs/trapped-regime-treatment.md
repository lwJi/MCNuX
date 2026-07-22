# Trapped-regime treatment (domain leaf: transport where the fluid is optically thick)

> Domain leaf spec. Self-contained: an agent can implement the optically-thick transport treatment — the explicit baseline, the opacity-relabeling extension with its α control and stability bound, and the diffusion-regime advection for scattering-dominated cells — from this file alone, referencing `conventions-and-units.md` for units/species/packet nomenclature, `rng-and-statistical-acceptance.md` for the uniform-draw primitive, pinned seeds/counts, and the 4σ band recipe, `neutrino-matter-interactions.md` for the event physics the effective coefficients feed, and `opacity-eos-evaluation.md` for how the unmodified coefficients are assembled. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact, 4σ) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

**⚠ α-convention warning (read before touching any equation).** Two Foucart papers describe the same opacity-relabeling scheme with the **same symbol α carrying the opposite meaning**. The 2020 ApJL letter ([arXiv:2008.08089](https://arxiv.org/abs/2008.08089), its Eq. 4) writes κ_a′ = (1−α)κ_a, κ_s′ = κ_s + ακ_a, η′ = (1−α)η — its α is the *removed* absorption fraction. The 2021 methods paper ([arXiv:2103.16588](https://arxiv.org/abs/2103.16588), its Eq. 43) writes η′ = αη, κ_a′ = ακ_a, κ_s′ = κ_s + (1−α)κ_a — its α is the *kept* absorption fraction. The two papers' α are related by α ↔ 1−α. **This corpus pins the 2021 convention exclusively**: α ∈ (0, 1] is the kept fraction, α = 1 means "no relabeling" (the explicit scheme), smaller α means more absorption relabeled to scattering. Any equation, parameter, diagnostic, or test in MCNuX using the symbol α with the opposite meaning is a correctness bug. Do not "fix" the equations below against the 2020 letter.

## Purpose & scope

This spec defines how MCNuX transports packets through optically thick fluid, where the explicit absorption–emission timescale 1/κ_a becomes small against the transport step Δt and naive explicit coupling is both unstable and prohibitively expensive: the **explicit scheme** (unmodified coefficients) as the baseline correctness contract; the **opacity-relabeling scheme** (Foucart 2021 convention) as a required, separately-verified extension selected by run parameter, with its per-cell α control rule, the invariants relabeling must preserve, and the restated stability bound that motivates it; the **diffusion-regime advection** replacing individual scattering events in scattering-dominated cells; and the built-in consistency requirement that the relabeling scheme with α = 1 reproduces the explicit scheme exactly.

In scope:

- The scheme selector and the explicit baseline (identity transformation, normative).
- The restated relabeling equations (2021 convention), the per-(cell, species, energy) α control rule κ_a′Δt_c ≤ ξ, and the restated stability bound.
- The algebraic invariants of relabeling: total interaction opacity preserved; equilibrium ratio η′/κ_a′ = η/κ_a preserved.
- The composition order with the low-temperature emissivity rule of `neutrino-matter-interactions.md`.
- The diffusion-regime trigger, behavior contract, and draw budget for scattering-dominated cells.
- The α → 1 equivalence requirement and the optically-thick benchmarks with concrete pass/fail numbers.

Out of scope:

- Assembly of the unmodified coefficients κ_a(s, ε), κ_s(s, ε), η_s(E) from the tables (`opacity-eos-evaluation.md`); this spec transforms coefficients it is handed.
- The event physics consuming the effective coefficients — event-time sampling, absorption/scattering outcomes, tallies (`neutrino-matter-interactions.md`; every equation there takes the effective coefficients produced here).
- Emission sampling mechanics (`packet-representation-and-sampling.md`; the sampler receives the effective η′).
- Geodesic propagation between events (`geodesic-propagation.md`); the diffusion-regime advection defined here *replaces* the geodesic push for the affected packet-steps but remains subject to that spec's substep caps.
- The source-term assembly (`hydro-coupling-source-terms.md`); relabeling changes *which* events occur, never the tally definitions.
- Implicit Monte Carlo (Fleck–Cummings), Fleck–Canfield random walks on frozen backgrounds, and moment-closure hybrids (discarded alternatives; see Open questions).

## Source of truth

The equations restated under "Inputs & outputs" are **normative**; citations are provenance only (restate-and-pin, per `README.md`):

- Opacity relabeling, stability bound, and practical control: Foucart et al. 2021 [arXiv:2103.16588](https://arxiv.org/abs/2103.16588) Eq. 43 (relabeling), Eqs. 48–49 (forward-Euler coupling stability), §2.7 (practical choice κ_a′Δt ≲ 1); the scheme is IMC-inspired (Fleck & Cummings 1971 provenance via the 2021 paper).
- The flipped-convention counter-reference: Foucart et al. 2020 [arXiv:2008.08089](https://arxiv.org/abs/2008.08089) Eqs. 4–5 — cited **only** as the warning above; nothing from it is normative here.
- Diffusion-approximation advection for scattering-dominated cells: Foucart et al. 2021 §2.8 (packets advected with the fluid plus a sampled random-walk displacement; the paper's fitted angular/displacement functions are provenance for the freedom, not restated — the observable contract below governs).
- The explicit-scheme baseline (resolve the physical timestep, no modification): Miller, Ryan & Dolence 2019 (nubhlight) [arXiv:1903.09273](https://arxiv.org/abs/1903.09273) §2 (fully explicit treatment).
- `WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` — the entry points producing the κ_a and κ_s this spec relabels (consumption contract owned by `opacity-eos-evaluation.md`).
- `CarpetX/HydroBaseX/interface.ccl` — the cell-centered fluid state whose per-cell values the per-cell α is computed from.
- `conventions-and-units.md` — packet nomenclature (ε = −p_μ u^μ in MeV, p^t), species enumeration, the length conversion `1.476625e5` cm per code unit under which κ [cm⁻¹] becomes inverse code length.
- `rng-and-statistical-acceptance.md` — the uniform-draw primitive R(S, q, e, k), event-counter semantics, pinned seeds/counts, the 4σ recipe.

## Inputs & outputs

### Interface

The treatment is a per-cell, per-species, per-energy transformation applied once per transport step, after coefficient assembly (and after the low-T rule) and before any sampling:

- **Scheme selector** (run parameter, keyword): `explicit` (default) or `relabeling`. The active scheme is fixed for a run.
- **ξ** (run parameter, real > 0, default `ξ = 1`): the relabeling scheme's per-step effective-absorption budget (used only when the relabeling scheme is active).
- **τ_diff** (run parameter, real > 0, default `τ_diff = 10`): the per-step scattering-optical-depth threshold above which a packet enters the diffusion regime (relabeling scheme only; the explicit scheme never uses the diffusion path).
- Inputs per cell: the transport step Δt_c (coordinate time, code units; the same Δt as the transport cadence), and the assembled coefficients — κ_a(s, ε), κ_s(s, ε) in inverse code length and the emissivity spectrum η̃_{s,b} (already low-T-modified per `neutrino-matter-interactions.md`).
- Outputs per cell: the **effective coefficients** η′_{s,b}, κ_a′(s, ε), κ_s′(s, ε) handed to the emission sampler and the event sampler, the per-(cell, species, group) diagnostic α, and per-packet diffusion-regime dispositions.

Downstream, `neutrino-matter-interactions.md` and `packet-representation-and-sampling.md` operate on the effective coefficients with **no other change**: event times use κ_a′ and κ_s′ in their −ln(r) equations, emission uses η′, and all tally definitions are untouched.

### Restated equations (normative)

**1. Explicit scheme (baseline; identity transformation).** With the `explicit` selector, the effective coefficients equal the inputs bitwise:

```text
η′ = η̃ ,   κ_a′ = κ_a ,   κ_s′ = κ_s      (no relabeling, no diffusion regime)
```

This is the nubhlight-style fully explicit treatment: correctness is unconditional, cost and coupling stability are the run's responsibility (the step must resolve 1/κ_a′ where coupling matters). Every physics requirement of the interaction and sampling specs is stated against exactly this scheme.

**2. Opacity relabeling (2021 convention, normative).** With the `relabeling` selector, per (cell, species s, energy):

```text
η′ = αη̃ ,   κ_a′ = ακ_a ,   κ_s′ = κ_s + (1−α)κ_a ,      α ∈ (0, 1]
```

Provenance: Foucart 2021 Eq. 43. A fraction (1−α) of absorption is relabeled as elastic scattering, and emission is reduced by the same factor: packets in equilibrium then scatter (cheaply, conservatively) instead of being destroyed and re-created at unresolvable rates. The relabeling is applied pointwise at the (s, ε) of each evaluation, with α from rule 3 evaluated at that (cell, s, ε); for emission, η′_{s,b} = α(cell, s, E_b) η̃_{s,b} with E_b the group's representative energy (the same one the sampler uses).

**3. The α control rule (normative).** α is a deterministic per-(cell, species, energy) function of the local unmodified absorption opacity and the step:

```text
α = min(1, ξ/(κ_a Δt_c))      with κ_a in inverse code time along the ray ≡ κ_a(s, ε) ε/p^t
```

evaluated with the same κ_a(s, ε), ε, p^t combination that enters the event-time equations of `neutrino-matter-interactions.md` (for cell-level quantities such as emission, the fluid-frame form κ_a Δt_c with ε/p^t → 1 is used — the static-fluid limit; the choice is pinned: cell-level α uses the fluid-frame κ_a Δt_c, packet-level relabeled opacities use the same cell-level α so that a single α governs a (cell, s, group)). Consequences, both normative:

```text
κ_a′Δt_c ≤ ξ         (everywhere, every step — the enforced invariant)
α = 1  exactly  wherever κ_a Δt_c ≤ ξ      (well-resolved cells are untouched)
```

With ξ = 1 (typical, the default), the effective absorption optical depth per step never exceeds one. Provenance: Foucart 2021 §2.7.

**4. Stability bound (restated provenance for rule 3).** The forward-Euler operator-split coupling of radiation to fluid internal energy is linearly stable only if

```text
κ_a′Δt ≤ 1/(1 + β̃) ,   equivalently   α < 1/((1 + β̃) κ_a Δt) ,      β̃ ≡ dU_ν/dU_fl
```

with β̃ the equilibrium stiffness ratio (change of neutrino equilibrium energy density per change of fluid energy density). Provenance: Foucart 2021 Eqs. 48–49. Rule 3 with ξ = 1 satisfies this whenever β̃ ≪ 1; in regimes where β̃ is not small, ξ must be reduced to ξ ≤ 1/(1 + β̃) (run-configuration responsibility; the spec pins the invariant κ_a′Δt_c ≤ ξ, and the bound here is the recipe for choosing ξ).

**5. Algebraic invariants of relabeling (normative).** Wherever the relabeling is applied:

```text
κ_a′ + κ_s′ = κ_a + κ_s               (total interaction opacity preserved exactly)
η′/κ_a′ = η̃/κ_a   wherever κ_a > 0   (the equilibrium ratio — hence u_eq — is preserved)
```

Both follow algebraically from equation 2; they are stated as independent requirements because they are what makes relabeling physically admissible (same total event rate, same equilibrium state, only the absorbed-vs-scattered split changes).

**6. Composition with the low-T rule (restated from `neutrino-matter-interactions.md`, binding).** The relabeling applies to the already-low-T-modified emissivity: η′ = αη̃ where η̃ is the output of the T^6 rule. The composition order low-T-first, relabeling-second is normative in both specs.

**7. Diffusion-regime advection (relabeling scheme only; behavior contract).** A packet-step leg is in the **diffusion regime** iff, in its current cell,

```text
τ_s ≡ κ_s′ ε Δt_leg / p^t ≥ τ_diff    and    κ_s′ ≥ 10 κ_a′
```

(Δt_leg the leg's coordinate-time extent; the factor-10 scattering-dominance guard is pinned, not a parameter). For such a leg, individual scattering events are **not** sampled; instead, in one RNG-consuming event (e increments by 1) consuming exactly **eight draws** of R(S, q, e, k):

- k = 0 — absorption test: the packet is absorbed during the leg with probability `1 − exp(−κ_a′ ε Δt_leg / p^t)`;
- k = 1 — absorption-time draw (always consumed; used, when absorbed, to place the removal time in the leg by the truncated-exponential inverse);
- k = 2…5 — fluid-frame random-walk displacement: three Gaussian components with per-dimension standard deviation `√(2 D t_fl)`, `D = 1/(3 κ_s′)` (fluid-frame diffusion coefficient, c = 1), t_fl the leg's fluid-frame duration, constructed from these four uniforms (the uniform→Gaussian map is implementation freedom within this budget);
- k = 6, 7 — exit direction: isotropic in the fluid frame (μ = 2u − 1, φ = 2πu, as elastic scattering).

A surviving packet ends the leg co-moving-displaced (advected with the fluid plus the sampled displacement, mapped to coordinates), with its fluid-frame energy preserved exactly (elastic chain) and its four-momentum rebuilt from the exit direction; the substep caps of `geodesic-propagation.md` (step end, |Δx^i| ≤ one cell width) still bind, and the scattering momentum-transfer tally ΔP_scat,μ += N (p_μ^before − p_μ^after) is deposited exactly as for discrete scattering. An absorbed packet is removed and tallied per the absorption rules of `neutrino-matter-interactions.md`. The sampled displacement/direction *distributions* beyond the Gaussian/isotropic forms above (e.g. the fitted finite-τ corrections of Foucart 2021 §2.8) are implementation freedom bounded by MCNX-TRP-06.

## Correctness requirements

- **[MCNX-TRP-01] Scheme selection and pass-through (exact/bitwise).** With the `explicit` selector, the effective coefficients are bitwise equal to the inputs for every (cell, species, energy), no α is applied anywhere, and no packet ever takes the diffusion path. With the `relabeling` selector in cells where κ_a Δt_c ≤ ξ, α = 1 exactly and the effective coefficients are again bitwise equal to the inputs. Reference: equations 1–3.
- **[MCNX-TRP-02] Relabeling identities (machine).** Wherever α < 1 is applied: η′ = αη̃, κ_a′ = ακ_a, and κ_s′ = κ_s + (1−α)κ_a each hold at `rtol = 1e-14` against independently recomputed values; the invariants κ_a′ + κ_s′ = κ_a + κ_s and η′/κ_a′ = η̃/κ_a (wherever κ_a > 0) hold at `rtol = 1e-14`. Reference: equations 2 and 5; provenance Foucart 2021 Eq. 43 — **2021 convention only** (see the warning box).
- **[MCNX-TRP-03] α control and the enforced budget (machine + exact).** α equals min(1, ξ/(κ_a Δt_c)) at `rtol = 1e-14` over the probe sweep of "Verification" (and exactly 1, bitwise, on the clamped branch); consequently κ_a′Δt_c ≤ ξ holds for every (cell, species, group) at every step — zero violations in instrumented runs (exact; a `1e-14`-relative grace above ξ is permitted for the equality case only). Reference: equations 3–4; provenance Foucart 2021 Eqs. 48–49 and §2.7.
- **[MCNX-TRP-04] α → 1 reproduces the explicit scheme (golden `1e-12`).** A fixed-seed single-rank CPU transport run with the `relabeling` selector on a fixture where κ_a Δt_c ≤ ξ everywhere (so α = 1 in every cell at every step) produces per-packet output identical to the same run with the `explicit` selector: the Cactus-harness numeric diff at `ABSTOL = RELTOL = 1e-12` shows zero violations (the runs are in fact bitwise identical — same coefficients, same draws, same event sequences). This is the built-in consistency limit of the extension. Reference: equations 1–3.
- **[MCNX-TRP-05] Optically-thick thermalization with relabeling on (numeric).** In the single-cell gray fixture of "Verification" with κ_a Δt_c = 10, κ_s = 0, ξ = 1 (so α = 0.1, κ_a′Δt_c = 1): the packet-census radiation energy density follows u_rad(t) = u_eq (1 − e^{−κ_a′ t}) with the **unmodified** equilibrium u_eq = η̃_tot/κ_a = η′_tot/κ_a′, within `1e-2` u_eq at each of t = {1, 2, 4}/κ_a′, at seed S0 with a census of ≥ `1e5` packets. This exercises both invariants of equation 5 end-to-end: the equilibrium is unchanged, the relaxation rate is the relabeled κ_a′. Reference: equations 2–5.
- **[MCNX-TRP-06] Diffusion-regime fidelity (numeric + 4σ).** In the uniform pure-scattering fixture of "Verification": (a) the per-dimension variance of packet positions at t_end matches the analytic diffusion value 2 D t_end, D = 1/(3 κ_s′), within `5e-2` relative; (b) the diffusion-path run and the explicit-scattering run (τ_diff = ∞) agree — per-dimension census position means within 4σ of each other and variances within `5e-2` relative; (c) every diffusion-regime leg preserves the fluid-frame energy at `rtol = 1e-14` and consumes exactly the eight draws of equation 7 in k-order (exact, via the logged event sequence); (d) zero diffusion-regime legs occur where the trigger condition is not met (exact). Reference: equation 7; provenance Foucart 2021 §2.8.

## Verification

- **MCNX-TRP-01**: family-B unit check: instrument the effective-coefficient stage on a fixed-seed fixture spanning κ_a Δt_c from `1e-3` to `1e3`; assert bitwise pass-through for the `explicit` selector everywhere and for the `relabeling` selector wherever κ_a Δt_c ≤ ξ; assert zero diffusion dispositions under `explicit`.
- **MCNX-TRP-02 / TRP-03**: family-B unit check sweeping a probe lattice of (κ_a, κ_s, η̃, Δt_c, ξ) spanning ≥ 6 decades of κ_a Δt_c: recompute α, η′, κ_a′, κ_s′ independently and assert the `1e-14` identities, the bitwise α = 1 clamp branch, and κ_a′Δt_c ≤ ξ at every probe.
- **MCNX-TRP-04**: family-A Cactus regression: one fixed-seed transport problem (single-rank CPU, seed S0) run under both selectors against the **same** committed golden per-packet TSV; both diffs pass at `ABSTOL = RELTOL = 1e-12` (per the harness contract in `verification-suite-design.md`).
- **MCNX-TRP-05**: family-B unit check: single cell, flat static metric, constant gray η̃_tot and κ_a with κ_a Δt_c = 10, κ_s = 0, ξ = 1, seed S0, census ≥ `1e5` packets; assert |u_rad − u_eq(1 − e^{−κ_a′ t})| ≤ `1e-2` u_eq at t = {1, 2, 4}/κ_a′ and report the margins (statistical legs per the margin design of `verification-suite-design.md`).
- **MCNX-TRP-06**: family-B unit check: uniform static-fluid flat-metric medium, κ_a = 0, κ_s = 32 (inverse code length), monoenergetic packets launched from the origin at S0 (≥ `1e5`), transport steps Δt = 0.5 (per-step τ_s = 16 ≥ τ_diff = 10) to t_end = 8 (analytic per-dimension variance 2 t_end/(3 κ_s′) = 1/6); a second run with τ_diff = ∞ forces explicit scatterings. Assert (a), (b), (c), (d) of MCNX-TRP-06.
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `WeakLibInterp/` and `CarpetX/` paths, and contains the required claim strings (the three 2021-convention relabeling equations, the α ↔ 1−α warning, the κ_a′Δt_c ≤ ξ budget, the stability bound, the equilibrium invariant, each benchmark tolerance).

## Implementation freedom

- **Where the relabeling is applied** — precomputed per-cell effective-coefficient tables per step, or on-the-fly at each evaluation — provided the observable values satisfy MCNX-TRP-02/03 and identical inputs give identical outputs.
- The representation of α (stored per (cell, species, group), recomputed on demand) and of the diffusion disposition — provided the diagnostics of MCNX-TRP-03/06 are observable.
- The uniform→Gaussian construction inside draws k = 2…5 of the diffusion event (Box–Muller, inverse-CDF, …) and any finite-τ corrections to the displacement/direction distributions (e.g. the fitted functions of Foucart 2021 §2.8) — bounded by MCNX-TRP-06's variance, agreement, and draw-budget requirements.
- How the fluid-frame displacement is mapped to coordinates in curved/moving backgrounds — bounded by the caps of `geodesic-propagation.md` and the static-flat benchmark here (see Open questions).
- Kernel structure, batching, and fusion with coefficient assembly — provided draw budgets and event-sequence reproducibility (MCNX-RNG-03) hold.

## Open questions / assumptions

- **β̃ is not computed by the baseline (assumption, run-configuration burden).** Rule 3 enforces κ_a′Δt_c ≤ ξ but does not estimate β̃; the default ξ = 1 assumes β̃ ≪ 1 (radiation energy density a small perturbation to the fluid's). A future automatic ξ control from an on-the-fly β̃ estimate (via dU_ν/dT and dU_fl/dT from the EOS) would be a spec change here, not a silent addition.
- **Diffusion advection is verified in static-flat fixtures only (assumption).** MCNX-TRP-06 pins the diffusion contract in a static, flat, uniform medium. Curved/moving-background accuracy of the displacement mapping is unverified by this corpus's benchmarks; production use in such cells relies on the displacement being ≪ cell size (mean free path 1/κ_s′ tiny by construction of the trigger). A curved-background diffusion benchmark is a recorded future spec change.
- **The scattering-dominance guard (κ_s′ ≥ 10 κ_a′) is pinned, not derived (assumption).** It ensures the diffusion path is taken only where scattering genuinely dominates within a leg; cells failing it fall back to discrete events, which is always correct (just slower). Retuning the constant is a spec change here.
- **Per-group α granularity (assumption).** α is per (cell, species, group), matching the granularity of the tabulated coefficients; a single per-cell α (cheaper, coarser) would satisfy the invariants only approximately across the spectrum and is not permitted without a spec change.
- **Discarded alternatives (recorded).** Fully-implicit or IMC (Fleck–Cummings) coupling and the Sedonu-style Fleck–Canfield random walk on frozen backgrounds were considered and rejected in design (poor fit for a time-dependent operator-split thorn); the explicit baseline + 2021 relabeling pair was chosen for its built-in α → 1 consistency limit (MCNX-TRP-04).
- **Published-PDF equation numbering (shared corpus assumption).** The 2021/2020 equation numbers cited as provenance follow the arXiv renderings (two independent renderings cross-checked); the restated equations are normative regardless.
