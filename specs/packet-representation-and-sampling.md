# Packet representation and sampling (domain leaf: what a packet is and how emission creates it)

> Domain leaf spec. Self-contained: an agent can implement the packet state contract and emission sampling from this file alone, referencing `conventions-and-units.md` for the shared notation, unit systems, and species enumeration, and `rng-and-statistical-acceptance.md` for the uniform-draw primitive, the pinned seeds and packet counts, and the 4σ band recipe. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact, 4σ) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines the Monte Carlo packet — the fundamental sample of the neutrino distribution function — and the emission-sampling process that creates packets from the fluid state: the packet state contract (components, types, units, valid ranges, invariants), the cell-integrated emission energy, the packet-count and equal-energy weighting rule, fluid-frame energy-spectrum sampling, fluid-frame isotropic angle sampling with the boost to the coordinate frame, in-cell position sampling, and the RNG draw budget of the creation event.

In scope:

- The packet state contract: position x^i, covariant coordinate-frame four-momentum components p_t and p_i, particle-count weight N, species s, RNG identity (q, e) — with types, units, valid ranges, and consistency invariants.
- The restated emission-sampling equations: cell-integrated emission energy, packet count, equal-energy weighting, spectrum sampling, angle sampling and boost, position sampling.
- The creation event's uniform-draw budget and draw ordering (part of the reproducible event sequence).

Out of scope:

- Propagation of packets between events (see `geodesic-propagation.md`; the p_t consistency identity is restated here because p_t is packet state).
- Absorption and scattering events, optical-depth sampling, and packet destruction (the neutrino–matter interactions spec).
- The optically-thick modification of emissivities and opacities (the trapped-regime spec); the equations here take the effective η actually handed to the sampler.
- How η is evaluated from the weak-interaction tables — units at the call, `log10` conventions, range enforcement (the opacity/EOS evaluation spec). This spec treats the emissivity spectrum as a given input with the units stated below.
- In-memory storage layout, container choice, and AMR residency (the particle-container spec). This spec pins observable components and invariants, not layout.

## Source of truth

The equations restated under "Inputs & outputs" are **normative**; citations are provenance only (restate-and-pin, per `README.md`):

- Cell-integrated emission energy: Foucart 2018 [arXiv:1708.08452](https://arxiv.org/abs/1708.08452) Eq. 18; Foucart et al. 2021 [arXiv:2103.16588](https://arxiv.org/abs/2103.16588) Eq. 20.
- Packet-as-sample representation (a packet carries a particle count N and a four-momentum p_μ; the distribution function is a sum of such weighted samples): Foucart et al. 2021 Eq. 3.
- Equal-energy weighting (weight w = C/ν so every packet carries equal energy, with the normalization set by a target packet count): Miller, Ryan & Dolence 2019 (nubhlight) [arXiv:1903.09273](https://arxiv.org/abs/1903.09273) Eqs. 28–30.
- `conventions-and-units.md` — packet component nomenclature (names, index placement, units), the species enumeration {νe, ν̄e, νx} with degeneracy g and lepton number ℓ, the geometrized↔cgs/MeV conversion factors used below.
- `rng-and-statistical-acceptance.md` — the uniform-draw primitive R(S, q, e, k), packet-id uniqueness (MCNX-RNG-05), the pinned seeds/counts, the 4σ recipe.
- `CarpetX/ADMBaseX/interface.ccl` and `CarpetX/HydroBaseX/interface.ccl` — the grid variables emission consumes: vertex-centered α, β^i, γ_ij; cell-centered rho, temperature, Ye, vel. HydroBaseX's `vel` is the Valencia 3-velocity v^i measured by the Eulerian observer (HydroBaseX README: Valencia formulation).
- Emissivity provenance: η values originate in the weak-interaction tables consumed through WeakLibInterp (`WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` and companions); the evaluation contract is owned by the opacity/EOS evaluation spec.

## Inputs & outputs

### Packet state contract

| Component | Symbol | Type | Units | Valid range / invariant |
|---|---|---|---|---|
| position | x^i | 3 × `double` | code (geometrized) length | finite; inside the computational domain (boundary handling owned by `geodesic-propagation.md`) |
| four-momentum | p_t, p_i | 4 × `double` | MeV | finite; γ^{ij} p_i p_j > 0; p_t < 0; p_t consistency identity below |
| weight | N | `double` | dimensionless (physical neutrinos represented) | N > 0, finite; for νx the degeneracy g = 4 is already folded in (exactly once, at creation) |
| species | s | integer | — | s ∈ {0, 1, 2} = {νe, ν̄e, νx} per the binding enumeration in `conventions-and-units.md` |
| RNG identity | (q, e) | 2 × `uint64` | — | q unique for the run (MCNX-RNG-05, exact); e starts at 0 at creation and increments by 1 per RNG-consuming event |

Derived (not independent state): the fluid-frame energy ε = −p_μ u^μ with u^μ the fluid four-velocity, in MeV; the contravariant time component p^t from the null condition (owned by `geodesic-propagation.md`).

**p_t consistency identity (machine tier).** The stored p_t is redundant given p_i and the local metric; wherever p_t is written or reported, it must satisfy, with the same metric values (α, β^i, γ^{ij}) used at that point,

```text
p_t = β^i p_i − α √(γ^{ij} p_i p_j)
```

to relative `1e-14`. (This is the null condition with future-directed p^t > 0; the algebra is restated and owned by `geodesic-propagation.md`.)

### Emission-sampling interface

Emission runs once per transport step per cell per species. Inputs, per (cell, species s):

- Cell geometry: coordinate widths (Δx, Δy, Δz), coordinate volume ΔV = Δx Δy Δz, in code units; cell-center metric values α and √γ (vertex-centered ADMBaseX data bridged to the cell center; the centering-bridge obligation is owned by the CarpetX-integration spec).
- Fluid state at the cell: the Valencia 3-velocity v^i (HydroBaseX `vel`), from which the fluid four-velocity is built as `W = (1 − γ_ij v^i v^j)^(−1/2)`, `u^t = W/α`, `u^i = W (v^i − β^i/α)`.
- The fluid-frame emissivity spectrum η_{s,b} ≥ 0 for energy groups b = 1…B with group edges [E_b^lo, E_b^hi) in MeV: the energy emission rate of species s per unit proper volume in group b, in `MeV cm^-3 s^-1`, evaluated at the cell's fluid state on the current time slice. For s = νx, η_{s,b} entering the equations below is the per-heavy-lepton-species table value multiplied by g = 4, exactly once (per `conventions-and-units.md`; no downstream factor).
- The transport step Δt (coordinate time, code units), the global seed S, and the emission-granularity run parameter N_tgt (target created packets per emitting (cell, species) per step).

Output: N_p ≥ 0 new packets appended to the packet population, each satisfying the state contract above, with fresh unique ids q and e = 0. Given identical inputs and seed, the output is identical (creation is deterministic at fixed key tuples).

### Restated emission-sampling equations (normative)

**1. Cell-integrated emission energy.** The total fluid-frame energy emitted by a cell in species s during one transport step is

```text
E_tot(cell, s) = α √γ ΔV Δt η_s ,    η_s = Σ_b η_{s,b}
```

where α √γ ΔV Δt is the invariant four-volume of the (cell × step) slab (√−g = α √γ). Units: with η_s in `MeV cm^-3 s^-1` and ΔV, Δt in code units, E_tot in MeV is obtained by multiplying by the cube of the length factor `1.476625e5` cm and the time factor `4.925491e-6` s from `conventions-and-units.md`. Provenance: Foucart 2018 Eq. 18 / Foucart 2021 Eq. 20.

**2. Packet count and equal-energy weighting.** Whenever E_tot > 0, the sampler creates N_p ≥ 1 packets (no silently dropped emission). Every packet created in the same (cell, species, step) carries the same fluid-frame energy

```text
E_p = E_tot / N_p
```

and, given its sampled fluid-frame neutrino energy ε_k, the particle-count weight

```text
N_k = E_p / ε_k    (so N_k ε_k = E_p for every k)
```

Provenance for the equal-energy convention: nubhlight Eqs. 28–30. The count N_p is a deterministic function of the cell inputs and run parameters, tracking N_tgt (requirement MCNX-PKT-07); the adaptivity policy itself is implementation freedom.

**3. Energy-spectrum sampling.** Packet k's fluid-frame energy ε_k is sampled so that the probability of landing in group b is the energy-weighted fraction

```text
P(b) = η_{s,b} / η_s
```

and ε_k lies inside the selected group's edges [E_b^lo, E_b^hi). Placement *within* the group is implementation freedom (bin center, uniform, or table-informed), provided it is deterministic at fixed key tuples and consumes only the draw budgeted below.

**4. Angle sampling and boost.** Directions are isotropic in the fluid frame: with uniform draws u_μ, u_φ,

```text
μ = 2 u_μ − 1,   φ = 2π u_φ,   Ω̂ = (√(1−μ²) cos φ, √(1−μ²) sin φ, μ)
```

The coordinate-frame four-momentum is built from the fluid-frame null momentum via an orthonormal tetrad {u^μ, e_a^μ} (g_μν e_a^μ e_b^ν = δ_ab, g_μν u^μ e_a^ν = 0):

```text
p^μ = ε_k (u^μ + Ω̂^a e_a^μ),   then   p_μ = g_μν p^ν
```

The tetrad construction is implementation freedom; the observable requirements are the exact identities ε = −p_μ u^μ = ε_k and g^{μν} p_μ p_ν = 0 (machine tier) and fluid-frame isotropy (4σ, MCNX-PKT-05).

**5. Position sampling.** The creation position is uniform in the cell's coordinate volume: per direction, x = x_lo + u · Δx with an independent uniform draw per coordinate.

**6. Draw budget of the creation event.** Creation is the packet's event e = 0 and consumes exactly **six uniform draws** of R(S, q, 0, k), in this k-order: k = 0, 1, 2 → position x, y, z; k = 3 → energy (group selection and, if used, in-group placement from the single draw); k = 4 → μ; k = 5 → φ. This ordering is part of the reproducible per-packet event sequence (MCNX-RNG-03).

## Correctness requirements

- **[MCNX-PKT-01] State validity and draw budget (exact + machine).** Every packet, at creation and at every subsequent observation point, satisfies the state-contract table: s ∈ {0, 1, 2} and (q, e) semantics exact (integer); N > 0, ε > 0, γ^{ij} p_i p_j > 0, p_t < 0, all components finite (no NaN/Inf); the p_t consistency identity holds at relative `1e-14`; and the creation event consumes exactly the six draws in the stated k-order (observable through the reproducible event sequence of MCNX-RNG-03). Reference: the state contract and draw budget above.
- **[MCNX-PKT-02] Cell-integrated emission energy identity (relaxed `1e-10`).** For every (cell, species, step) with E_tot > 0: Σ_k N_k ε_k over the packets created there equals E_tot(cell, s) = α √γ ΔV Δt η_s at relative `1e-10` (deterministic bookkeeping identity; relaxed tier for the summation). For s = νx this includes g = 4 exactly once, inside η. Reference: equations 1–2 above; provenance Foucart 2018 Eq. 18 / 2021 Eq. 20.
- **[MCNX-PKT-03] Equal-energy weighting (machine).** Within one (cell, species, step), the products N_k ε_k are pairwise equal at relative `1e-14`, and each equals E_p = E_tot/N_p; equivalently N_k = E_p / ε_k. Reference: equation 2; provenance nubhlight Eqs. 28–30.
- **[MCNX-PKT-04] Spectrum fidelity (4σ).** In a fixture with a known synthetic spectrum η_{s,b}, at seed S0 and N_standard = `1e6` created packets, the observed fraction f_b of packets in each group matches P(b) = η_{s,b}/η_s within 4σ, with standard error s_b = √(P(b)(1 − P(b))/N_p) per group. Seeds, counts, and the band recipe per `rng-and-statistical-acceptance.md`.
- **[MCNX-PKT-05] Fluid-frame isotropy and boost identities (4σ + machine).** In a fixture with a moving fluid (W > 1): (a) the fluid-frame direction moments over N_standard packets satisfy ⟨Ω̂^a⟩ = 0 for each component (standard error `1/√(3N)`), ⟨(Ω̂^a)²⟩ = 1/3 for each component (standard error `√(4/45)/√N`), and ⟨Ω̂^a Ω̂^b⟩ = 0 for a ≠ b (sample-based standard error, recipe (b) of the RNG spec), each within 4σ; (b) for **every** created packet, |−p_μ u^μ − ε_k| ≤ `1e-14` ε_k and |g^{μν} p_μ p_ν| ≤ `1e-14` ε_k² (machine tier), i.e. the boost preserves the fluid-frame energy and nullness exactly. Reference: equation 4 above.
- **[MCNX-PKT-06] In-cell position uniformity (4σ).** Over N_standard packets created in one cell: per coordinate, the sample mean matches the cell center within 4σ (standard error `w/√(12N)`, w the cell width) and the sample variance matches w²/12 within 4σ (sample-based standard error). Reference: equation 5 above.
- **[MCNX-PKT-07] Emission granularity tracks the target (exact integer band).** In a homogeneous fixture (uniform η, uniform fluid, flat metric) with N_tgt = 64 and default parameters, every (cell, species, step) creates N_p ∈ [32, 128] packets (i.e. within a factor 2 of N_tgt), and N_p ≥ 1 wherever E_tot > 0. Integer check, no tolerance. Reference: equation 2 and the interface contract above.

## Verification

- **MCNX-PKT-01**: a unit check runs a fixed-seed emission fixture and validates every created packet against the state table (integer fields exact; finiteness; sign constraints; p_t identity at `rtol = 1e-14` using the fixture's metric values), and asserts from the logged event sequence that creation consumed draws k = 0…5 of event e = 0 in order.
- **MCNX-PKT-02 / PKT-03**: the same fixture with synthetic η: per (cell, species) assert Σ_k N_k ε_k against the analytically computed E_tot at `rtol = 1e-10`, and pairwise equality of N_k ε_k at `rtol = 1e-14`; a νx cell asserts the single g = 4 factor (cross-check with MCNX-CNV-02).
- **MCNX-PKT-04**: unit check at (S0, N_standard) on a ≥ 8-group synthetic spectrum spanning ≥ 3 orders of magnitude in η_{s,b}; per-group margin m = 4·s_b − |f_b − P(b)| reported and asserted m ≥ 0.
- **MCNX-PKT-05**: unit check on a cell with v = (0.3, 0.2, −0.1) (W ≈ 1.08): moment battery at 4σ; per-packet boost identities at `rtol = 1e-14`.
- **MCNX-PKT-06**: unit check at (S0, N_standard) in a single cell; mean and variance margins asserted per coordinate.
- **MCNX-PKT-07**: unit check on a 16³-cell homogeneous fixture; assert the integer band for every (cell, species, step) over 4 steps.
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` and `WeakLibInterp/` paths, and contains the required claim strings (the E_tot equation, the equal-energy weight rule, the isotropy statement, the six-draw budget).

## Implementation freedom

- The packet-count adaptivity policy (how N_p is chosen from E_tot, N_tgt, and history) — provided MCNX-PKT-07's band, determinism at fixed seed, N_p ≥ 1 for E_tot > 0, and the equal-energy identities hold.
- In-group energy placement (bin center, uniform-in-group, table-informed) — provided groups are honored, the placement is deterministic, and only draw k = 3 is consumed.
- The orthonormal tetrad construction used for the boost (Gram–Schmidt seeded from any coordinate triad, or closed-form) — provided the machine-tier identities and 4σ isotropy of MCNX-PKT-05 hold.
- Storage layout, precision of intermediates, batching of creation over cells, and whether ε or p_μ is the stored primary (the observable state contract governs).
- Emission-loop order over cells/species and any parallelization — provided per-packet streams are pure functions of (S, q) per the RNG contract.

## Open questions / assumptions

- **In-group spectral placement is unconstrained beyond determinism (assumption).** The 4σ spectrum requirement is stated at group granularity, matching the group structure of the tabulated emissivities; no sub-group spectral shape is asserted. If a future consumer (e.g. pair annihilation) needs sub-group fidelity, MCNX-PKT-04 gains a refined form in a spec change here.
- **Packet splitting and Russian roulette are deferred (assumption, non-blocking).** No packet is ever split or stochastically terminated by this spec; if the interactions or trapped-regime specs later introduce weight control, they must preserve MCNX-PKT-01/02/03 semantics and assign fresh unique q to any packet they create.
- **Emission floor (assumption: none by default).** N_p ≥ 1 whenever E_tot > 0 means arbitrarily small E_tot still emits one (low-energy-weight) packet per species. A documented floor parameter that skips cells below a threshold — with the skipped energy accounted in the conservation ledger — is an allowed future extension recorded here, default off.
- **HydroBaseX `vel` is the Valencia 3-velocity (assumption, shared with `conventions-and-units.md`).** The W and u^μ formulas above assume it; if the producing GRMHD code supplies a different velocity variable, only the u^μ construction at the read boundary changes.
