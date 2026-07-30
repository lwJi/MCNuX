# Packet representation and sampling

> Leaf spec. Self-contained: an agent can implement the Monte Carlo packet state and
> the emission sampling of energies, angles, and species from this file alone. It
> restates only the conventions it uses and references
> [conventions-and-units](./conventions-and-units.md) (units, constants, species) and
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (the random-draw
> primitive and the statistical bar) for the rest. `README.md` is canonical if any
> restated convention here conflicts with it.

## Purpose & scope

This spec pins:

- The **packet state**: the physical components every Monte Carlo packet carries
  (position, coordinate-frame four-momentum, weight, species, RNG identity), their
  units, meanings, and valid ranges.
- **Emission sampling**: how many packets are created per cell, step, species, and
  energy bin; how each packet's energy, position, creation time, and direction are
  drawn; and the tetrad transform from the fluid frame to grid coordinates
  (arXiv:1708.08452 Eqs. 18–21, restated below).
- The rule that the νx degeneracy factor **g = 4 enters exactly once, at packet
  creation** — in the emission count — and nowhere else.
- The RNG consumption map for creation: every stochastic draw in this spec is an
  evaluation of the pure function `u(S, q, e, k)` of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md), with the
  `(q, e, k)` assignment pinned here.

Out of scope:

- The in-memory layout of these components (pure-SoA container, runtime components,
  id storage/uniqueness) — [particle-container-and-gpu](./particle-container-and-gpu.md).
- Propagation of created packets ([geodesic-propagation](./geodesic-propagation.md))
  and post-creation interaction events
  ([neutrino-matter-interactions](./neutrino-matter-interactions.md)).
- The emissivity η_b itself — tabulated assembly and the νx mapping
  ([opacity-eos-evaluation](./opacity-eos-evaluation.md)); the analytic-mode formulas
  ([verification-suite-design](./verification-suite-design.md)).
- The trapped-regime modification of emission (η′ = αη relabeling) —
  [trapped-regime-treatment](./trapped-regime-treatment.md); this spec's count law
  applies to whatever η the active scheme supplies.
- The fluid-interface contract that provides `(ρ, T, Yₑ, u^μ)`
  ([hydro-coupling-source-terms](./hydro-coupling-source-terms.md)).

## Source of truth

The sampling equations are restated in full below (normative per `README.md`) from the
pinned methods paper; the bin-center energy convention carries its own provenance.

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; Eqs. 18–21 are the emission sampling restated below |
| arXiv:1706.06187v2 | Richers et al. 2017, ApJ 847, 133 — provenance of the bin-center packet-energy choice adopted by the methods paper |

## Inputs & outputs

### Packet state (binding components; layout is not prescribed here)

| Component | Symbol | Type | Units | Valid range / meaning |
|---|---|---|---|---|
| position | `x^i` | 3 × double | geometrized coordinates | inside the computational domain |
| momentum | `p_i` | 3 × double | geometrized | lower spatial components of the coordinate-frame four-momentum; not all zero |
| weight | `N` | double | dimensionless count | > 0; number of physical neutrinos represented (for νx: summed over the four lumped heavy-lepton states) |
| species | `s` | integer | — | ∈ {0 (νe), 1 (ν̄e), 2 (νx)} |
| packet id | `q` | uint64 | — | unique, immutable for the packet's life (the AMReX `idcpu` identity; uniqueness contract in [particle-container-and-gpu](./particle-container-and-gpu.md)) |
| event counter | `e` | uint32 | — | 0 at creation; incremented per stochastic event by the consuming specs |

Derived (not independent state): `p^t = √(γ^{ij} p_i p_j)/α` (null closure, owned by
[geodesic-propagation](./geodesic-propagation.md)); the fluid-frame energy
`ν = −p_μ u^μ > 0`, converted to MeV at the microphysics boundary via the pinned
energy factor of [conventions-and-units](./conventions-and-units.md).

### Emission inputs

Per cell, transport step, species `s`, and tabulated energy bin `b`: the bin-integrated
emissivity `η_b(s)` (MeV cm⁻³ s⁻¹ from the active opacity source), the cell fluid state
`(ρ, T, Yₑ, u^μ)` (cell-centered), the metric at creation positions (gathered per
[geodesic-propagation](./geodesic-propagation.md)'s gather contract), the cell
coordinate volume ΔV and step Δt, and the run-time target packet energy `E_p` [MeV].
Output: the set of newly created packets with the state above.

## Correctness requirements

- **[MCNX-PKT-01] Packet state schema.** Every packet carries exactly the physical
  content of the state table above (however stored). `q` never changes; `N`, `s` are
  set at creation and change only as sibling specs prescribe (baseline: never — an
  absorbed packet is removed, not down-weighted); `e` starts at 0 and is
  non-decreasing.

- **[MCNX-PKT-02] Emission count — g enters exactly once (restated, normative;
  provenance arXiv:1708.08452 Eq. 18).** The expected number of packets created in a
  cell during one transport step, per species `s` and energy bin `b`, is

  ```text
  N_p = g_s √−g ΔV Δt η_b(s)/E_p ,        g_s = (1, 1, 4) for (νe, ν̄e, νx)
  ```

  with `√−g = α√γ` and the product evaluated in one documented unit system (the
  dimensionless result is unit-system independent to machine tier). Exactly `⌊N_p⌋`
  packets are created, plus one more with probability `N_p − ⌊N_p⌋` (one uniform draw
  per (cell, step, s, b), pinned in [MCNX-PKT-05]). The νx degeneracy factor g = 4
  enters the corpus **exactly once, at packet creation**, in this count law; it appears
  never in per-species transport coefficients
  ([opacity-eos-evaluation](./opacity-eos-evaluation.md)) and never again downstream
  (`N` already counts all lumped states).

- **[MCNX-PKT-03] Per-packet initialization (restated; provenance arXiv:1708.08452
  Eq. 18 prose; bin-center choice arXiv:1706.06187).** Every packet created in bin
  `b = [E_{b−1}, E_b]` is initialized with:

  - fluid-frame energy `ν = 0.5 (E_{b−1} + E_b)` (the bin center, exact — energies are
    never drawn continuously in the baseline);
  - weight `N = E_p/ν` (so each packet carries fluid-frame energy `E_p`);
  - position drawn uniformly in the cell's coordinate volume (three independent uniform
    draws, one per coordinate);
  - creation time drawn uniformly in the step interval `[t_n, t_n + Δt)` (one draw);
    the packet participates in transport from its creation time to the end of the step.

- **[MCNX-PKT-04] Direction sampling and tetrad transform (restated, normative;
  provenance arXiv:1708.08452 Eqs. 19–21).** The packet's four-momentum is drawn
  isotropically in the fluid rest frame:

  ```text
  p^{μ′} = ν (1, sinθ cosφ, sinθ sinφ, cosθ)
  cosθ = 2u₁ − 1 ∈ [−1, 1) ,   φ = 2π u₂ ∈ [0, 2π)
  ```

  with `u₁, u₂` uniform draws, and transformed to grid coordinates with the
  orthonormal tetrad `ê^μ_(λ′)`:

  ```text
  p^t = ê^t_(λ′) p^{λ′} ,      p_i = g_{iμ} ê^μ_(λ′) p^{λ′}
  ```

  The tetrad is pinned: timelike leg `ê^μ_(0′) = u^μ`; spatial legs obtained by
  Gram–Schmidt orthonormalization of the coordinate basis vectors against `g_μν`,
  **in fixed coordinate order (x, y, z)** (the paper leaves the order open; it is
  pinned here for golden-data determinism). Consistency (machine tier, `~1e-14`
  relative): the transformed momentum is null (`g^{μν} p_μ p_ν = 0` normalized by
  `(αp^t)²`), and `−p_μ u^μ` recovers the drawn ν.

- **[MCNX-PKT-05] RNG consumption map (creation).** All draws above are evaluations of
  the pure function `u(S, q, e, k)` of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md); AMReX's
  stateful device RNG is forbidden for them. The pinned assignment:

  - **Cell-level count draw** (the Bernoulli of [MCNX-PKT-02]): `q = K_cell`,
    `e = n` (the transport-step index, counted from 0), `k = 3b + s` (bin b, species
    s). The cell key is

    ```text
    K_cell = 2⁶³ + patch·2⁵⁶ + level·2⁴⁸ + (k_c−k₀)·2³² + (j_c−j₀)·2¹⁶ + (i_c−i₀)
    ```

    with `(i_c, j_c, k_c)` the cell's integer coordinates on its level,
    `(i₀, j₀, k₀)` the level's domain lower corner, and the 2⁶³ flag bit
    disambiguating cell keys from packet ids (see Open questions for the capacity
    assumption).
  - **Per-packet creation draws** (event `e = 0` of the new packet's own id `q`):
    `k = 0, 1, 2` → position x, y, z; `k = 3` → creation time; `k = 4` → `u₁`
    (cosθ); `k = 5` → `u₂` (φ). Draws are consumed in exactly this order; the
    packet's event counter is 0 during creation and creation consumes no other draws.

- **[MCNX-PKT-06] Emission statistics.** At the pinned packet count and seeds of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md), over a pinned
  emission-only benchmark: the realized packet count per (cell, s, b) has mean `N_p`
  (within 4σ of the binomial/floor construction); the total emitted fluid-frame energy
  per (cell, step, s) is within 4σ of `g_s √−g ΔV Δt Σ_b η_b(s)`; and the
  fluid-frame angular distribution is isotropic (first moments of `cosθ` and `φ`
  within 4σ of their nulls).

### Conventions restated (the subset this leaf uses)

- Species `{νe, ν̄e, νx}` with indices (0, 1, 2), degeneracies g = (1, 1, 4), lepton
  numbers (+1, −1, 0); g = 4 for νx appears in this corpus only in [MCNX-PKT-02].
- Signature `(−,+,+,+)`; lapse `α`, spatial metric `γ_ij`, `√−g = α√γ`; fluid-frame
  energy `ν = −p_μ u^μ`; `u^μ` normalized `g_μν u^μ u^ν = −1` (checked to machine
  tier at consumption).
- Geometrized units for `x^i`, `p_i`, ΔV, Δt; MeV for ν, `E_p`, and η_b's energy
  content; conversions use exactly the pinned factor set of
  [conventions-and-units](./conventions-and-units.md).
- All stochastic draws are `u(S, q, e, k)` evaluations; uniform variates lie in
  `[0, 1)`.

## Verification

- **Static-flat closed form (machine tier, `~1e-14`).** On Minkowski with a static
  fluid (`u^μ = (1, 0, 0, 0)`), the tetrad is the identity: `p^t = ν` and
  `p_i = ν n̂_i` with `n̂_i` the drawn unit vector, to machine tier; the null and
  ν-recovery identities of [MCNX-PKT-04] hold at every created packet.
- **Bin-center exactness (exact tier).** Every created packet's fluid-frame ν equals
  its bin center bitwise; every weight equals `E_p/ν` to machine tier.
- **Count-law determinism (exact tier).** At fixed seed on a single rank, two runs of
  an emission-only step create identical packet sets (ids, states) bitwise — the
  reproducibility contract of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) applied to
  creation.
- **Emission statistics (4σ at pinned count and seed).** The checks of [MCNX-PKT-06]
  on the emission benchmark of
  [verification-suite-design](./verification-suite-design.md), including the analytic
  opacity mode's νx channel verifying the g = 4 factor (νx emitted energy = 4× the
  single-species analytic law, within 4σ).
- **Draw-map audit (exact tier).** A trace of `(q, e, k)` tuples consumed during a
  pinned small emission step matches the map of [MCNX-PKT-05] exactly (no reordering,
  no extra draws).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the restated-equation anchors
`N_p = g_s √−g ΔV Δt η_b(s)/E_p` and `p^{μ′} = ν (1, sinθ cosφ, sinθ sinφ, cosθ)`, the
literal `g = 4`, and the phrase `exactly once, at packet creation`; and that its
`[MCNX-PKT-NN]` ids are declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- The storage layout, component naming, and any internal rescaling of stored momenta by
  a documented global constant — provided observables (ν in MeV, geometrized momenta at
  the interfaces, golden outputs) are as specified.
- Whether the creation loop is organized per cell, per bin, per species, or fused;
  packet-buffer growth strategy; where ids are assigned within creation (subject to the
  uniqueness contract of [particle-container-and-gpu](./particle-container-and-gpu.md)
  and single-rank determinism).
- Caching of `√−g`, tetrad legs, and bin constants; batching of draws (values must
  equal the pure function).
- Whether `ν` is stored or recomputed from `p_μ` and `u^μ` at the point of use
  (observable values must agree to machine tier with the derived definition).

## Open questions / assumptions

- **Assumption: cell-key capacity.** The `K_cell` packing supports patch < 128,
  level < 256, and per-dimension level extents < 65536 cells, and relies on packet
  `idcpu` values never having bit 63 set (AMReX packs 40-bit id + 24-bit cpu; bit 63
  requires cpu ≥ 2²³ ranks). Exceeding any of these would require a deliberate,
  golden-data-invalidating re-pin of the packing.
- **Assumption: uniform in-cell sampling.** Positions are drawn uniformly in
  *coordinate* volume (per the source paper), which under-resolves proper-volume
  variation across a cell on strongly curved/coarse grids; accepted at baseline, and a
  higher-order sampling would be a spec change, not an implementation choice.
- **Assumption: emission uses cell-centered fluid state.** η_b and the tetrad's `u^μ`
  are the creation cell's cell-centered values (no sub-cell fluid interpolation at
  emission); the metric at the drawn position uses the gather contract of
  [geodesic-propagation](./geodesic-propagation.md).
- **Open question: packet-count control.** `E_p` is the only packet-count control
  pinned at baseline (per the source paper). Weight windowing, splitting/rouletting,
  or per-region `E_p` would modify [MCNX-PKT-02]'s realized counts and must be
  specified here before adoption.
