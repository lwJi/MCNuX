# Neutrino–matter interactions (domain leaf: sampling and outcomes of absorption and elastic scattering)

> Domain leaf spec. Self-contained: an agent can implement the collision physics — event-time sampling, absorption, elastic scattering, the low-temperature emissivity rule, and the per-event tallies — from this file alone, referencing `conventions-and-units.md` for units/species/packet nomenclature, `rng-and-statistical-acceptance.md` for the uniform-draw primitive, pinned seeds/counts, and the 4σ band recipe, `packet-representation-and-sampling.md` for the packet state contract, and `geodesic-propagation.md` for how the next-event time caps propagation. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact, 4σ) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines what happens to a packet between creation and removal, other than free streaming: the exponential (−ln r) optical-depth sampling of absorption and elastic-scattering event times, the event-time resampling policy and its statistical exactness, the absorption outcome (packet removal plus exact tallying), the elastic-scattering outcome (fluid-frame re-isotropization at exactly preserved fluid-frame energy), the low-temperature emissivity extrapolation (∝ T^6) applied to the emissivity spectrum before emission sampling, the RNG draw budgets of each event type, and the per-event tally definitions the hydro-coupling source-term spec consumes.

In scope:

- The restated event-time sampling equations and the resampling (memorylessness) policy.
- Absorption and elastic-scattering outcomes as observable state transitions.
- The low-T emissivity rule η ∝ T^6 below the pinned threshold.
- Event draw budgets and ordering (part of the reproducible per-packet event sequence).
- Per-event tally definitions (emission, absorption, scattering) at cell granularity, with their bookkeeping identities.
- The 0-D cooling/thermalization and two-region attenuation benchmarks with concrete pass/fail numbers.

Out of scope:

- How κ_a(s, ε), κ_s(s, ε), and η_s(E) are evaluated from the weak-interaction tables — units at the call, `log10` conventions, range enforcement, species→dataset mapping (see `opacity-eos-evaluation.md`). This spec consumes the assembled coefficients.
- Emission sampling mechanics — packet counts, weights, spectrum/angle/position draws (see `packet-representation-and-sampling.md`); this spec only modifies the η handed to that sampler at low T.
- Geodesic motion between events and the substep caps (see `geodesic-propagation.md`; this spec supplies the next-event time t_ev that spec honors).
- The optically-thick relabeling η′ = αη, κ_a′ = ακ_a, κ_s′ = κ_s + (1−α)κ_a (the trapped-regime spec). Every equation below takes the *effective* coefficients actually handed to it; the baseline is the unmodified coefficients.
- Inelastic scattering (NES), pair annihilation/production, and bremsstrahlung as transport processes (deferred; see `opacity-eos-evaluation.md` Open questions).
- The assembly of tallies into the radiation four-force G^μ and lepton source (the hydro-coupling source-term spec); this spec defines the raw per-event tallies only.

## Source of truth

The equations restated under "Inputs & outputs" are **normative**; citations are provenance only (restate-and-pin, per `README.md`):

- Exponential optical-depth event sampling in coordinate time: Foucart et al. 2021 [arXiv:2103.16588](https://arxiv.org/abs/2103.16588) Eq. 24; Miller, Ryan & Dolence 2019 (nubhlight) [arXiv:1903.09273](https://arxiv.org/abs/1903.09273) Eqs. 31–32.
- Low-temperature emissivity extrapolation ∝ T^6: Foucart et al. 2021 Eq. 19.
- Elastic re-isotropization in the fluid frame at preserved fluid-frame energy: Foucart et al. 2021 §2 (isoenergetic scattering); Richers et al. 2015 (Sedonu) [arXiv:1507.03606](https://arxiv.org/abs/1507.03606) (elastic fluid-frame re-isotropization).
- 0-D cooling/thermalization benchmark design: nubhlight Eqs. 58–59 and §4 (artificial cooling against analytic decay; late-time deviation ≲ 1%).
- `conventions-and-units.md` — species enumeration and lepton numbers ℓ = {+1, −1, 0}, packet nomenclature (ε = −p_μ u^μ in MeV), the length conversion factor `1.476625e5` cm per code unit used to convert opacities.
- `rng-and-statistical-acceptance.md` — the uniform-draw primitive R(S, q, e, k), event-counter semantics, pinned seeds/counts, the 4σ recipe.
- `packet-representation-and-sampling.md` — the packet state contract, the fluid four-velocity construction u^μ from the HydroBaseX Valencia 3-velocity, and the tetrad/boost identities the scattering outcome reuses.
- `CarpetX/HydroBaseX/interface.ccl` — the cell-centered fluid variables (rho, temperature, Ye, vel) whose cell values define the interaction state a packet sees.

## Inputs & outputs

### Interface

Interactions operate per packet within one transport step, interleaved with geodesic propagation. Inputs:

- Packet state (x^i, p_i, p_t, N, s, (q, e)) per `packet-representation-and-sampling.md`; the fluid-frame energy ε = −p_μ u^μ (MeV) evaluated with the current cell's fluid state.
- The effective coefficients of the packet's current cell at its current fluid-frame energy: κ_a(s, ε) and κ_s(s, ε) in cm⁻¹ from `opacity-eos-evaluation.md`, converted to inverse code length by multiplying by `1.476625e5` (cm per code length unit) before use in the equations below.
- The uniform-draw primitive R(S, q, e, k); each draw below uses r = 1 − u ∈ (0, 1] so that ln r is finite (u ∈ [0, 1) per the RNG contract).
- The transport-step end time t_n + Δt and the cell geometry (for crossing detection, owned by `geodesic-propagation.md`).

Outputs: the packet's next-event coordinate time t_ev and event type; on absorption, packet removal plus tallies; on scattering, updated (p_i, p_t) plus tallies; the modified emissivity spectrum η̃ handed to the emission sampler at low T.

### Restated equations (normative)

**1. Event-time sampling (competing exponentials).** Whenever event times are (re)sampled, draw two independent uniforms r_a, r_s and form the candidate coordinate-time intervals

```text
Δt_a = −ln(r_a) · p^t / (κ_a ε) ,    Δt_s = −ln(r_s) · p^t / (κ_s ε)
```

with κ in inverse code length, ε and p^t in MeV (their ratio is dimensionless, so Δt is in code time). A vanishing opacity gives Δt = ∞ (no candidate). The next event is the earlier candidate:

```text
t_ev = t_cur + min(Δt_a, Δt_s) ,   event type = absorption if Δt_a ≤ Δt_s else scattering
```

If t_ev exceeds the earlier of the step end and the next cell crossing, no event occurs in that leg. Provenance: Foucart 2021 Eq. 24 (the invariant optical depth along a null trajectory is dτ = κ ε dt / p^t, so τ = −ln r gives the intervals above); nubhlight Eqs. 31–32.

**2. Resampling policy (memorylessness, normative).** Event times are resampled — a fresh RNG-consuming event, both draws — (a) at the start of each transport step, (b) whenever the packet enters a new cell, and (c) immediately after a scattering event. Between resamplings, κ_a, κ_s, and ε are held frozen at the values used for the draw (consistent with the frozen-metric-per-cell policy of `geodesic-propagation.md`). Because the exponential distribution is memoryless, discarding an unexpired event time at a cell boundary and resampling with the new cell's coefficients is *statistically exact*, not an approximation; the frozen-ε/κ within a cell leg is a bounded discretization approximation that vanishes with resolution (provenance: Foucart 2021 §2).

**3. Absorption outcome.** The packet is removed from transport at t_ev: it is never propagated past t_ev, contributes the absorption tallies below, and consumes **zero** uniform draws.

**4. Elastic-scattering outcome.** At t_ev the packet's fluid-frame energy is preserved exactly and its direction is re-isotropized in the fluid frame: with two uniform draws u_μ, u_φ,

```text
μ′ = 2 u_μ − 1 ,   φ′ = 2π u_φ ,   Ω̂′ = (√(1−μ′²) cos φ′, √(1−μ′²) sin φ′, μ′)
p′^μ = ε (u^μ + Ω̂′^a e_a^μ) ,   p′_μ = g_μν p′^ν
```

using the same orthonormal-tetrad construction {u^μ, e_a^μ} as emission (`packet-representation-and-sampling.md`; the tetrad is implementation freedom, the identities are not): −p′_μ u^μ = ε and g^{μν} p′_μ p′_ν = 0 at machine tier. N, s, and q are unchanged; transport then resumes with freshly sampled event times (rule 2c).

**5. Low-temperature emissivity rule (∝ T^6).** For a cell whose fluid temperature satisfies T < T_low, the emissivity spectrum handed to the emission sampler is

```text
η̃_s(E; T) = η_s(E; T_low) · (T / T_low)^6 ,    T_low = 0.5 MeV (run parameter, default 0.5; T_low ≥ the table floor 0.1 MeV)
```

i.e. the spectrum is evaluated once at T_low and scaled by (T/T_low)^6; for T ≥ T_low the rule is inactive and η̃ = η exactly. Absorption and scattering opacities are *not* modified by this rule (they use the clamped table evaluation of `opacity-eos-evaluation.md`). Provenance: Foucart 2021 Eq. 19 — charged-current emissivities fall steeply (∝ T^6) below ≈ 0.5 MeV where table entries are unreliable; the rule suppresses spurious low-T emission smoothly to zero as T → 0.

**6. Per-event tallies (cell granularity; the raw inputs to the hydro-coupling spec).** Per (cell, species, step), accumulated over events occurring in that cell during that step, with ℓ_s the lepton number and every sum weighted by the packet weight N:

```text
Emission   (per created packet k):   E_em  += N_k ε_k ,   L_em  += ℓ_s N_k ,   P_em,μ  += N_k p_μ,k
Absorption (per absorbed packet):    E_abs += N ε ,       L_abs += ℓ_s N ,     P_abs,μ += N p_μ
Scattering (per scattering event):   ΔP_scat,μ += N ( p_μ^before − p_μ^after )     (μ = t, x, y, z)
```

E-tallies are fluid-frame energies (MeV); P-tallies are coordinate-frame covariant components (MeV). The assembly of these tallies into the radiation four-force density and lepton source — including unit conversion to geometrized grid quantities and the deposition protocol — is owned by the hydro-coupling source-term spec; the definitions above are binding for it.

**7. Draw budgets (part of the reproducible event sequence).** Each RNG-consuming event increments the packet's event counter e by exactly 1 and consumes exactly these draws of R(S, q, e, k), in k-order:

| Event type | Draws | k-order |
|---|---|---|
| event-time sampling (rules 2a/2b/2c) | 2 | k = 0 → r_a, k = 1 → r_s |
| elastic scattering outcome | 2 | k = 0 → u_μ, k = 1 → u_φ |
| absorption | 0 | — (not an RNG-consuming event; e unchanged) |

(Packet creation consumes six draws per `packet-representation-and-sampling.md`.)

## Correctness requirements

- **[MCNX-INT-01] Exponential event-time law (4σ).** In a homogeneous static-fluid flat-metric fixture with constant gray κ_a and κ_s = 0, the sampled absorption intervals Δt_a over N_standard = `1e6` packets at seed S0 are exponential with mean m = p^t/(κ_a ε): the sample mean matches m within 4σ (standard error m/√N) and the sample variance matches m² within 4σ (sample-based standard error, recipe (b) of the RNG spec). The same holds for Δt_s with κ_a = 0, κ_s > 0, and for the event-type split with both nonzero: the absorbed fraction matches κ_a/(κ_a + κ_s) within 4σ (binomial standard error). Reference: equation 1; provenance Foucart 2021 Eq. 24.
- **[MCNX-INT-02] Draw budgets and event-sequence discipline (exact).** Every event consumes exactly the draws of the budget table in the stated k-order; e increments by exactly 1 per RNG-consuming event and never for absorption; the per-packet event sequence is bitwise reproducible in single-rank CPU test mode (rides on MCNX-RNG-03). Integer/bitwise check, no tolerance. Reference: equation 7.
- **[MCNX-INT-03] Elastic scattering: energy preservation and re-isotropization (machine + 4σ).** Per scattering event, the fluid-frame energy is preserved: |−p′_μ u^μ − ε| ≤ `1e-14` ε, and |g^{μν} p′_μ p′_ν| ≤ `1e-14` ε² (machine tier, every event). Over N_standard scattering events of an initially collimated beam in a moving-fluid cell (v = (0.3, 0.2, −0.1), as MCNX-PKT-05), the outgoing fluid-frame direction moments satisfy the isotropy battery — ⟨Ω̂′^a⟩ = 0 (standard error `1/√(3N)`), ⟨(Ω̂′^a)²⟩ = 1/3 (standard error `√(4/45)/√N`), ⟨Ω̂′^a Ω̂′^b⟩ = 0 for a ≠ b (sample-based standard error) — each within 4σ, demonstrating no memory of the incoming direction. Reference: equation 4.
- **[MCNX-INT-04] Absorption accounting (exact + relaxed).** In instrumented fixed-seed runs: every packet whose event resolves to absorption is removed at t_ev and never observed later; the absorption tallies equal the summed content of the removed packets — counts and species exactly (integer), E_abs/L_abs/P_abs sums at `rtol = 1e-10` (relaxed tier for summation order); L-tallies carry exactly ℓ = {+1, −1, 0} per species. Reference: equations 3 and 6.
- **[MCNX-INT-05] Low-T rule (machine + exact).** For fixture cells with T < T_low: the sampler's η̃_{s,b} equals η_{s,b}(T_low) · (T/T_low)^6 at `rtol = 1e-14` for every group; for T ≥ T_low: η̃ = η bitwise (the rule is inactive, exact). At T = 0, η̃ = 0 exactly. Reference: equation 5; provenance Foucart 2021 Eq. 19.
- **[MCNX-INT-06] 0-D cooling and thermalization benchmarks (numeric + relaxed).** On the single-cell fixtures of "Verification":
  - *Cooling*: with synthetic gray emissivity η_tot = u_fl/t_cool tied to a fluid internal-energy ledger and κ_a = κ_s = 0, the ledger follows the discrete operator-split recurrence u_{n+1} = u_n (1 − Δt/t_cool) at `rtol = 1e-10` (deterministic bookkeeping — emission energy per step is deterministic), and matches the continuum decay u(t) = u(0) e^{−t/t_cool} within `3e-2` relative at t = 3 t_cool with Δt = t_cool/64 (the O(Δt) splitting bias is the dominant term).
  - *Thermalization*: with constant synthetic η_tot and gray κ_a (κ_s = 0), the packet-census radiation energy density follows u_rad(t) = u_eq (1 − e^{−κ_a t}) with u_eq = η_tot/κ_a: at each of t = {1, 2, 4}/κ_a, |u_rad − exact| ≤ `1e-2` u_eq at seed S0 with a census of ≥ `1e5` packets (statistical noise ≈ 0.3%; provenance for the ≲ 1% late-time bar: nubhlight Eqs. 58–59, §4).
- **[MCNX-INT-07] Two-region attenuation (4σ).** A collimated monoenergetic beam launched through two adjacent homogeneous slabs with distinct gray absorption opacities κ_1 ≠ κ_2 (κ_s = 0, static fluid, flat metric) is transmitted with probability exp(−κ_1 d_1 − κ_2 d_2) (d_i the path lengths in each slab, opacities in inverse code length): over N_standard packets at S0, the transmitted fraction matches within 4σ (binomial standard error). This exercises the resampling-at-crossing policy: any double-counting or dropped optical depth at the boundary biases the result far beyond 4σ. Reference: equations 1–2.

## Verification

- **MCNX-INT-01**: family-B unit check on a single homogeneous cell (periodic/large domain so no packet escapes early); three legs (pure absorption, pure scattering, mixed) at (S0, N_standard); margins m = 4s − |Q̂ − Q| reported and asserted m ≥ 0 for each statistic.
- **MCNX-INT-02**: asserted from the logged per-packet event sequences of the fixed-seed fixtures (creation, event-time, scattering events with their k-draws), single-rank CPU mode; integer equality.
- **MCNX-INT-03**: family-B unit check: a beam of N_standard identical packets in the moving-fluid cell, each forced through exactly one scattering event; per-event machine-tier identities and the 4σ moment battery on outgoing directions.
- **MCNX-INT-04**: instrumented fixed-seed runs of the INT-01 and INT-06 fixtures; assert removal semantics, integer counts, and tally sums at `rtol = 1e-10`.
- **MCNX-INT-05**: family-B unit check sweeping T ∈ {0, 0.1, 0.3, 0.49, 0.5, 0.51, 5} MeV against a synthetic η spectrum; assert the scaling identity at `rtol = 1e-14` below T_low, bitwise equality at/above, and zero at T = 0.
- **MCNX-INT-06**: family-B checks (statistical legs per the margin design of `verification-suite-design.md`): single cell, flat static metric; cooling leg with the ledger recurrence (relaxed `1e-10`) and continuum bound (`3e-2` at t = 3 t_cool); thermalization leg with census assertions (`1e-2` of u_eq at the three pinned times, ≥ `1e5` packets, seed S0).
- **MCNX-INT-07**: family-B unit check: two-slab fixture, e.g. κ_1 d_1 = 1.0 and κ_2 d_2 = 0.5 (transmission e^{−1.5} ≈ 0.223); binomial 4σ margin asserted at (S0, N_standard).
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` path, and contains the required claim strings (both event-time equations, the memorylessness policy, the T^6 rule, the isotropy statement, the draw budgets, each benchmark tolerance).

## Implementation freedom

- **Event-loop structure** — whether event times are sampled per cell leg on entry (as stated) via per-packet time stepping, kernel fusion with the geodesic push, batching by cell or by species, and host/device split — provided the observable draws, budgets, and outcomes match equations 1–7.
- The representation of "no candidate" (∞), and how min/tie-breaking is computed — provided the tie rule Δt_a ≤ Δt_s → absorption is honored exactly (it is part of the reproducible sequence).
- The tetrad construction for the scattering boost (shared freedom with emission) — bounded by the machine-tier identities and 4σ isotropy of MCNX-INT-03.
- How tallies are accumulated (per-thread partials, atomics + `SumBoundary`-style reduction) — bounded by MCNX-INT-04's `1e-10` bookkeeping identity and the grid-tally reproducibility contract MCNX-RNG-04.
- Where η̃ is applied (at spectrum assembly or inside the sampler) — provided MCNX-INT-05's identities hold on the sampler's actual input.
- Diagnostics, counters, and event logging format — provided the event sequence needed by MCNX-INT-02 is observable in single-rank CPU test mode.

## Open questions / assumptions

- **Frozen coefficients within a cell leg (bounded approximation, by design).** κ and ε are frozen between resamplings; in strong-gradient or high-velocity cells this is an O(Δx) discretization error that vanishes with resolution (same character as the frozen-metric policy of `geodesic-propagation.md`). No additional in-cell subdivision is required; if profiling ever shows this dominating error budgets, a sub-cell resampling rule would be a spec change here.
- **No biasing/oversampling in the baseline (assumption, non-blocking).** nubhlight's per-channel biasing factors (its Eqs. 33–35) oversample rare processes; the baseline corpus samples unbiased. Introducing biasing would change draw budgets and weights and is a spec change here plus `packet-representation-and-sampling.md`.
- **Absorption consumes no draws (pinned; consequence for sequences).** Because absorption is draw-free, packets absorbed at their first event have exactly one RNG-consuming event after creation (the event-time sampling). Any future absorption-side sampling (e.g. partial absorption via weight reduction) changes MCNX-INT-02 and is a spec change.
- **T_low interaction with the trapped regime (assumption).** The T^6 rule and the trapped-regime relabeling compose as: relabeling applies to the already-modified η̃ (η̃′ = αη̃). This ordering is asserted here as the binding composition; the trapped-regime spec must restate it when it lands.
- **Fluid state is per-cell constant (shared assumption).** The interaction state a packet sees is its current cell's cell-centered HydroBaseX values (no sub-cell interpolation of ρ, T, Yₑ, v^i); consistent with the emission sampler and the provenance codes. Sub-cell fluid interpolation would be a coordinated spec change here and in `packet-representation-and-sampling.md`.
