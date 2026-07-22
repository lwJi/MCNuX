# Hydro-coupling source terms (domain leaf: closing the loop back to GRMHD)

> Domain leaf spec. Self-contained: an agent can implement the radiation→fluid coupling surface — the MCNuX-owned cell-centered source-term grid variables, their assembly from the per-event tallies, the zero-then-add schedule protocol and one-step-lag cadence, the optional TmunuBaseX contribution, and the conservation budgets — from this file alone, referencing `conventions-and-units.md` for units/species/conversion factors, `neutrino-matter-interactions.md` for the binding per-event tally definitions consumed here, and `rng-and-statistical-acceptance.md` for the grid-tally reproducibility bound and pinned seeds/counts. The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact, 4σ) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict. The GRMHD evolution code itself is out of scope corpus-wide: this spec pins only the interface MCNuX exposes.

## Purpose & scope

This spec defines how the energy, momentum, and lepton number exchanged between neutrinos and matter reach the fluid solver: the MCNuX-owned, cell-centered source-term grid variables (the densitized radiation four-force and lepton source) with their names, centering, units, and sign convention; the normative discrete assembly of those variables from the per-event tallies of `neutrino-matter-interactions.md`; the zero-then-add accumulation protocol and the observable cadence contract (one-step lag); the optional accumulation of the Monte Carlo radiation stress-energy into TmunuBaseX for spacetime sourcing, with the required cell-to-vertex bridge; and the corpus's conservation budgets with explicit numeric tolerances.

In scope:

- The exposed grid variables: group membership, names, centering (cell, `CCC`), units (geometrized code units), sign convention, and the consumer's READS obligation.
- The restated continuum coupling equations and their normative discrete realization from the tallies.
- The zero-then-add protocol (modeled on `CarpetX/TmunuBaseX/schedule.ccl`) and the cadence contract, stated observably.
- The optional TmunuBaseX contribution: the restated MC stress-energy formula and the 8-point cell-to-vertex average.
- Conservation budgets — deterministic tally identities and the two-region lepton benchmark — with concrete pass/fail numbers.

Out of scope:

- The per-event tally *definitions* (E_em, L_em, P_em,μ; E_abs, L_abs, P_abs,μ; ΔP_scat,μ) — owned by `neutrino-matter-interactions.md` and restated below only as consumed.
- The GRMHD code's use of the sources (its formulation, densitization, and RHS mechanics); MCNuX stays ignorant of the consumer beyond the interface pinned here.
- The trapped-regime relabeling (`trapped-regime-treatment.md`); relabeling changes which events occur, never the tally or assembly definitions here.
- Deposition mechanics — atomics, guard cells, `SumBoundary`-style reduction, AMR level reconciliation (the particle-container spec); this spec pins what the deposited values must equal, to the reproducibility bound of MCNX-RNG-04.
- Literal schedule-file text and schedule-bin placement (the CarpetX-integration spec); this spec states the protocol and cadence as observable behavior.
- Neutrino pair annihilation as a distinct momentum-deposition channel (deferred; see Open questions).

## Source of truth

The equations restated under "Inputs & outputs" are **normative**; citations are provenance only (restate-and-pin, per `README.md`):

- Radiation four-force from packet-event tallies, operator-split once per fluid step: Foucart et al. 2021 [arXiv:2103.16588](https://arxiv.org/abs/2103.16588) Eqs. 5–10 (fluid-coupling moments tallied from packet events).
- Lepton-number source and the discrete per-event Ye increment: Miller, Ryan & Dolence 2019 (nubhlight) [arXiv:1903.09273](https://arxiv.org/abs/1903.09273) Eqs. 4, 23 (G_Ye), 36 (discrete per-event increment); its Eq. 2 for the four-force placement in the fluid equations.
- MC radiation stress-energy reconstructed from packets: Foucart 2018 [arXiv:1708.08452](https://arxiv.org/abs/1708.08452) Eq. 16.
- `CarpetX/TmunuBaseX/schedule.ccl` — the zero-then-add accumulation protocol this spec's source-term protocol mirrors (an explicit zeroing routine, then an empty accumulation group consumers schedule `IN`, with `TmunuBaseX_ZeroTmunu` guaranteed first); `CarpetX/TmunuBaseX/interface.ccl` — the vertex-centered eTtt/eTti/eTij the optional contribution accumulates into.
- `CarpetX/HydroBaseX/interface.ccl` — the cell centering (`CCC`) the source-term variables match, and the fluid variables (rho, Ye, temperature, vel) the coupling ultimately feeds.
- `flesh/lib/sbin/RunTestUtils.pl` — the Cactus regression harness the cadence and lag goldens ride on.
- `conventions-and-units.md` — the species lepton numbers ℓ = {+1, −1, 0}, the MeV→code-energy conversion (1 code energy unit = `1.115416e60` MeV), packet nomenclature, and the centering facts (metric vertex, hydro cell).
- `neutrino-matter-interactions.md` — the binding per-event tally definitions (its equation 6) this spec consumes unchanged.

## Inputs & outputs

### The exposed grid variables (the interface, binding)

MCNuX declares and owns these grid-function groups, all **cell-centered** (`CCC`, matching `CarpetX/HydroBaseX/interface.ccl`), single time level, checkpointed (see the cadence contract):

| Group | Members | Content | Units |
|---|---|---|---|
| `MCNuX::rad_force` | `rad_srct`, `rad_srcx`, `rad_srcy`, `rad_srcz` | S_μ — covariant components of the **densitized radiation four-force on the fluid**, S_μ = √−g F_μ | code energy per code volume per code time |
| `MCNuX::lep_source` | `lep_src` | S_L — the **densitized lepton-number source on the fluid**, S_L = √−g F_L | number per code volume per code time |

Sign convention (binding): S_μ and S_L are the **fluid's gain**. Positive S_t-component energy flow heats the fluid; a νe-absorbing cell has S_L > 0. The continuum statement being discretized (restated; provenance nubhlight Eq. 2, Foucart 2021 §2.1):

```text
∇_ν T_fl^{μν} = F^μ ,    ∇_ν T_rad^{μν} = −F^μ ,    ∇_μ (n_e u^μ) = F_L
```

with F^μ the radiation four-force density on the fluid and n_e the fluid's net electron(−lepton) number density. The stored variables are the √−g-densitized covariant components, so a consumer evolving per-coordinate-volume densitized conserved variables adds S_μ (and m_b-scaled S_L, m_b its baryon-mass convention) to its RHS directly, with no further metric factors.

**Consumer obligation:** the GRMHD thorn declares `READS: MCNuX::rad_force(interior)` and `READS: MCNuX::lep_source(interior)` in the RHS routines that consume them, and never writes them. MCNuX declares the matching WRITES; the driver's validity machinery enforces both sides.

### Normative discrete assembly (from the binding tallies)

Consumed inputs, per (cell, step n), restated from `neutrino-matter-interactions.md` equation 6 (binding there): the N-weighted per-event tallies

```text
P_em,μ  = Σ N_k p_μ,k   (created packets) ,    P_abs,μ = Σ N p_μ   (absorbed packets) ,
ΔP_scat,μ = Σ N (p_μ^before − p_μ^after)   (scattering events) ,
L_em = Σ ℓ_s N_k ,    L_abs = Σ ℓ_s N
```

with p_μ in MeV, L dimensionless, and events attributed to the cell containing the event position. The grid variables for step n are then, per cell (normative):

```text
S_μ(cell, n) = c_E · ( P_abs,μ + ΔP_scat,μ − P_em,μ ) / (ΔV Δt)
S_L(cell, n) =        ( L_abs − L_em ) / (ΔV Δt)
```

with ΔV the cell's coordinate volume and Δt the transport step (code units), and c_E = 1/`1.115416e60` (MeV per code energy unit, `conventions-and-units.md`; the sole unit conversion in the assembly). No √−g division occurs: the tallies count events in the coordinate cell-step slab, so the quotient *is* the densitized (√−g-weighted) source. Summed over cells and steps, Σ S_μ ΔV Δt is exactly the net coordinate-frame four-momentum (code units) transferred to the fluid — the ledger identity of MCNX-SRC-04. Provenance: Foucart 2021 Eqs. 5–10; nubhlight Eqs. 23/36 for the per-event lepton increments (ℓ_s per the species table of `conventions-and-units.md`; νx events contribute exactly zero to S_L).

### Zero-then-add protocol and cadence (observable contract)

Modeled on the `CarpetX/TmunuBaseX/schedule.ccl` pattern (an explicit zeroing routine followed by accumulation, in a fixed order within one group):

1. **Zero.** At a fixed point inside MCNuX's per-step transport phase — after any consumer read of the previous values, before any tallying — MCNuX zeroes `rad_force` and `lep_source` everywhere.
2. **Add.** Transport accumulates (only `+=`, never `=`) event tallies into the variables as events occur during step n's transport; by the end of the transport phase the variables hold exactly S_μ(·, n), S_L(·, n).
3. **Read.** The consumer reads the variables outside MCNuX's transport phase. **Cadence contract (normative, stated observably): the source terms read by the fluid at step n are those tallied during step n−1's transport.** Equivalently: at every point outside MCNuX's transport phase, the variables hold the most recent *completed* transport step's sources, and they change exactly once per coarsest-level Δt.
4. **Restart.** The variables are checkpointed and restored, so the one-step lag survives checkpoint/restart bit-for-bit in single-rank CPU mode (they are *state*, not scratch — unlike TmunuBaseX's `checkpoint="no"` groups).

The first transport step of a run (and of a restart-less initial slice) begins from zeroed sources: the fluid's step-0 read sees zero, matching the "no completed transport step yet" reading of the cadence contract.

### Optional TmunuBaseX contribution (parameter-selected; default off)

When the boolean run parameter `rad_tmunu` is enabled, MCNuX additionally accumulates the Monte Carlo radiation stress-energy into TmunuBaseX so the spacetime solver sees the neutrinos, via a routine scheduled `IN TmunuBaseX_AddToTmunu` (the empty accumulation group of `CarpetX/TmunuBaseX/schedule.ccl`, guaranteed to run after `TmunuBaseX_ZeroTmunu`). The cell-centered MC stress-energy, restated (normative; provenance Foucart 2018 Eq. 16):

```text
T_rad^{μν}(cell) = c_E / (√−g ΔV) · Σ_{k ∈ cell} N_k p_k^μ p_k^ν / p_k^t
```

summed over the packets resident in the cell on the current slice (p^μ raised with the local metric; units: code energy density). Because eTtt/eTti/eTij are **vertex**-centered (`CarpetX/TmunuBaseX/interface.ccl`) while the packet census is tallied per cell, the contribution is bridged by the **8-point cell-to-vertex average**: each vertex receives the unweighted mean of the values of the (up to) 8 surrounding cells, and the routine *accumulates* (`+=`) the lowered components T_μν into eTtt/eTti/eTij — never overwrites (the group contract guarantees prior zeroing and other contributors). With `rad_tmunu` off (the default), MCNuX schedules nothing in the group and eT* is untouched by MCNuX.

## Correctness requirements

- **[MCNX-SRC-01] Interface shape and accumulation discipline (exact).** The two groups exist with the stated names, members, cell centering (`CCC`), and single time level; between the per-step zeroing and the end of the transport phase the variables change only by accumulation of event contributions; outside the transport phase they are never written. Verified structurally and by instrumented runs (any non-accumulating write or centering mismatch is a hard failure). Reference: the interface tables and protocol above.
- **[MCNX-SRC-02] Assembly identity (relaxed `1e-10`).** For every (cell, step) of an instrumented fixed-seed run, the deposited S_μ and S_L equal the assembly formulas evaluated from the independently logged per-event tallies — `S_μ = c_E (P_abs,μ + ΔP_scat,μ − P_em,μ)/(ΔV Δt)`, `S_L = (L_abs − L_em)/(ΔV Δt)` — at `rtol = 1e-10` (deterministic tally identity; relaxed tier for summation order), with event-to-cell attribution exact (each event deposits into the cell containing its position, integer check). Reference: the assembly equations above; provenance Foucart 2021 Eqs. 5–10, nubhlight Eqs. 23/36.
- **[MCNX-SRC-03] Cadence and lag (golden `1e-12`).** In a fixed-seed single-rank CPU regression, the source variables sampled at the consumer's read point at step n are bitwise those assembled from step n−1's logged tallies, for every step including across a mid-run checkpoint/restart; the step-0 read is exactly zero. Harness diff at `ABSTOL = RELTOL = 1e-12` against committed golden TSV (the values are bitwise-reproducible under MCNX-RNG-03). Reference: the cadence contract above.
- **[MCNX-SRC-04] Global conservation ledger (relaxed `1e-10` + exact).** Over any fixed-seed run, for each μ and for lepton number, the closed ledger holds:

  ```text
  Σ_cells,steps S_μ ΔV Δt + c_E ( Δ census Σ N p_μ + Σ_escaped N p_μ ) = 0
  Σ_cells,steps S_L ΔV Δt + Δ census Σ ℓ_s N + Σ_escaped ℓ_s N = 0
  ```

  each at `rtol = 1e-10` relative to the gross exchange (Σ of the absolute tallied contributions); νx packets contribute exactly zero to every lepton term (ℓ = 0, exact). What the fluid gains, the radiation field (census + escapes) loses — a bookkeeping identity independent of sampling noise. Sampled *physical* budgets ride the 4σ benchmarks (MCNX-SRC-06, MCNX-TRP-05). Reference: the assembly equations; conservation contract of `README.md`.
- **[MCNX-SRC-05] Tmunu contribution (exact + relaxed).** With `rad_tmunu` off, MCNuX never writes eTtt/eTti/eTij (exact). With it on: the routine runs inside `TmunuBaseX_AddToTmunu` and accumulates — pre-existing contributions are preserved exactly (verified by seeding a known prior value); with zero packets the contribution is exactly zero; with a synthetic packet population of known analytic T_rad^{μν}, each vertex value equals the 8-point cell-to-vertex average of the restated formula at `rtol = 1e-10`, and the metric-lowering uses the local vertex metric (machine-tier identity on analytic metrics). Reference: the Tmunu section above; provenance Foucart 2018 Eq. 16.
- **[MCNX-SRC-06] Two-region lepton diffusion (numeric + relaxed).** In the closed two-region fixture of "Verification" (optically thick νe/ν̄e exchange between two half-domains initialized at Ye_1 = 0.35 and Ye_2 = 0.15, uniform ρ and T, test driver applying dYe/dt = m_b S_L/ρ per cell): (a) total lepton number (fluid + packet census, no escapes) is conserved at `rtol = 1e-10` at every step; (b) the region-mean contrast ΔYe(t) = Ye̅_1 − Ye̅_2 decreases strictly monotonically across the 8 pinned output times; (c) at t_end every cell satisfies |Ye − Ye_eq| ≤ `5e-2` · ΔYe(0) with Ye_eq = (Ye_1 + Ye_2)/2, the fixture designed so the census holds ≤ `1e-2` of the total lepton number at t_end. Provenance for the benchmark design: nubhlight's two-region lepton-diffusion test (its §4). Reference: the S_L assembly above.

## Verification

- **MCNX-SRC-01**: structural check of the declared groups (names, centering, time levels) plus an instrumented fixed-seed run asserting write-phase discipline (zero → adds only → read-only) each step.
- **MCNX-SRC-02**: family-B unit check on the fixed-seed fixtures of `neutrino-matter-interactions.md` (INT-01/INT-06 style) with per-event logging: recompute S_μ, S_L per (cell, step) from the log and assert `rtol = 1e-10` plus exact event-to-cell attribution.
- **MCNX-SRC-03**: family-A Cactus regression: fixed seed, single-rank CPU, ≥ 4 steps with a mid-run checkpoint/restart leg; committed golden TSV of the source variables at the read point diffed at `ABSTOL = RELTOL = 1e-12`; step-0 zero asserted exactly.
- **MCNX-SRC-04**: the same instrumented runs accumulate the full ledger (deposits, census change, escapes) and assert both closures at `rtol = 1e-10` of gross exchange, and integer-exact zero νx lepton contributions.
- **MCNX-SRC-05**: family-B unit check on a small analytic-metric grid: synthetic packet populations (single-cell delta, uniform field) with hand-computable T_rad^{μν}; assert the off/on behaviors, the seeded-prior accumulation, the zero-packet zero, and the 8-point average at `rtol = 1e-10`.
- **MCNX-SRC-06**: family-B check: closed (reflecting/periodic) box split into two half-domains, gray νe/ν̄e opacities sized so each half is optically thick (κ_a L ≥ 10) with κ_a Δt_c ≤ ξ (explicit regime), ≥ `1e5` packets in census at S0, the scripted test driver playing the consumer role (forward-Euler Ye update from `lep_src`, m_b = `1.66054e-24` g); assert (a), (b) at the 8 pinned output times, and (c) at t_end.
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` and `flesh/` paths, and contains the required claim strings (the zero-then-add protocol, the cadence sentence, the assembly combinations, `TmunuBaseX_AddToTmunu`, the 8-point cell-to-vertex average, each conservation tolerance).

## Implementation freedom

- **Deposition mechanics** — per-thread partials, atomics into guard cells with `SumBoundary`-style folding, or serial accumulation — bounded by MCNX-SRC-02's `1e-10` identity and the grid-tally reproducibility bound MCNX-RNG-04 (`rng-and-statistical-acceptance.md`).
- Whether S_μ/S_L are accumulated in tally units (MeV) and converted once at the end of the step, or converted per event — provided the observable end-of-step values satisfy MCNX-SRC-02.
- The internal realization of the zero/add/read phases (one routine or several, kernel fusion with transport) — provided the observable protocol and cadence of MCNX-SRC-01/03 hold; literal schedule text is owned by the CarpetX-integration spec.
- How the Tmunu cell-to-vertex average handles domain-boundary vertices with fewer than 8 neighbor cells (one-sided average over existing cells, or zero-padded) — provided the choice is deterministic, documented in the implementation, and interior vertices satisfy MCNX-SRC-05 exactly as stated.
- The ledger instrumentation (per-step reductions vs. end-of-run) — provided MCNX-SRC-04's closures are observable at the stated tolerance.
- Additional diagnostic outputs (per-species source splits, luminosities) in any units, provided they are not part of the pinned interface.

## Open questions / assumptions

- **The consumer's densitization convention (assumption, to confirm at first coupling).** S_μ = √−g F_μ per coordinate volume and time is the natural source for Valencia-type densitized conserved variables; if the first real GRMHD consumer expects a different weighting (e.g. undensitized F_μ), the conversion happens on the consumer side or this interface gains a documented alternative group — a spec change here, never a silent re-interpretation of the stored values.
- **Baryon-mass convention for Ye (assumption).** S_L is a pure number source; converting it to a Ye source divides by the consumer's baryon number density ρ/m_b, so the m_b convention (amu `1.66054e-24` g in the fixtures) is the consumer's. MCNuX never bakes an m_b into `lep_src`.
- **Fine-coarse source consistency under AMR (deferred to the particle-container spec).** Events tallied on refined levels must reach the level(s) the consumer reads consistently (restriction/reflux of source terms); the single-rate, single-level configurations verified here do not exercise it. The requirement lands with the particle-container and CarpetX-integration specs.
- **Pair annihilation momentum deposition (deferred).** Foucart 2021's effective-κ_p pair-annihilation channel (its Eqs. 27, 31–35) deposits four-momentum without a packet absorption event; it is not in the baseline interaction set (`opacity-eos-evaluation.md` Open questions) and would add a tally term to the S_μ assembly here when specified.
- **Radiation pressure in the fluid's primitive recovery (out of scope, recorded).** The sources exchange four-momentum only; no radiation-pressure contribution to the fluid EOS is exposed or implied. Consumers requiring one need the optional Tmunu surface or a future extension here.
