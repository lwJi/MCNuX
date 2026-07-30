# Hydro coupling and source terms

> Leaf spec. Self-contained: an agent can implement MCNuX's coupling surface — the
> source-term grid variables, the event ledger that fills them, the zero-then-add
> schedule protocol, the interface requirements on the GRMHD partner, and the optional
> stress-energy contribution — from this file alone. It restates only the conventions
> it uses and references [conventions-and-units](./conventions-and-units.md) (units,
> species) and [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)
> (reproducibility limits of atomic tallies) for the rest. `README.md` is canonical if
> any restated convention here conflicts with it.

## Purpose & scope

This spec pins how MCNuX feeds energy, momentum, and lepton-number exchange back to a
fluid evolution — **GRMHD-agnostically**: every contract below is stated against the
CarpetX base thorns plus explicit interface requirements on whatever GRMHD partner
couples to MCNuX. No production GRMHD code is cited anywhere in this corpus.

- The **MCNuX-owned, cell-centered source-term grid variables** (radiation four-force
  density and net lepton-number source) and their binding ledger definition.
- The **zero-then-add schedule protocol** and the **contributor idiom**, as protocol
  requirements modeled on TmunuBaseX.
- The **interface requirements on the GRMHD partner**: `READS:` on the source
  variables, and the **MeV temperature** requirement on `HydroBaseX::temperature`.
- The **conservation ledger** closing packets against grid tallies.
- The **optional Tmunu contract**: the MC stress-energy estimator (arXiv:1708.08452
  Eq. 16, restated) accumulated into TmunuBaseX via the corpus-derived 8-point
  cell-to-vertex average.

Out of scope:

- Which events occur and their outcomes —
  [neutrino-matter-interactions](./neutrino-matter-interactions.md),
  [trapped-regime-treatment](./trapped-regime-treatment.md),
  [packet-representation-and-sampling](./packet-representation-and-sampling.md); this
  spec consumes the event stream those specs define.
- *When* in the step the source variables advance — the observable cadence contract
  (once per coarsest-level Δt, one-step lag) is owned by
  [carpetx-thorn-integration](./carpetx-thorn-integration.md) [MCNX-CTX-05]; this spec
  pins *what* they contain and the intra-step zero/add/read ordering.
- The deposition execution idiom (device atomics, tile buffering) —
  [particle-container-and-gpu](./particle-container-and-gpu.md); this spec pins where
  tallies land and what they must sum to.
- The unit convention the MeV requirement rests on
  ([conventions-and-units](./conventions-and-units.md)); concrete coupled benchmarks
  and golden data ([verification-suite-design](./verification-suite-design.md)).
- The GRMHD partner itself — only the interface MCNuX exposes to it is specified.

## Source of truth

The zero-then-add protocol is derived from the stack's own base-thorn pattern; the
estimator is restated from the pinned methods paper. Reference artifacts:

- `CarpetX/TmunuBaseX/schedule.ccl` — the reference protocol (lines 1–34): an outer
  group `IN ODESolvers_PostStep AFTER ADMBaseX_SetADMVars`, an explicit zeroing routine
  (`WRITES: eTtt(everywhere) eTti(everywhere) eTij(everywhere)`, "Set T_munu to
  zero"), and the empty extension group ordered after it ("Add to T_munu here"):

  ```text
  SCHEDULE GROUP TmunuBaseX_SetTmunuVars IN ODESolvers_PostStep AFTER ADMBaseX_SetADMVars
  SCHEDULE TmunuBaseX_ZeroTmunu IN TmunuBaseX_SetTmunuVars
    { WRITES: eTtt(everywhere) eTti(everywhere) eTij(everywhere) }  "Set T_munu to zero"
  SCHEDULE GROUP TmunuBaseX_AddToTmunu IN TmunuBaseX_SetTmunuVars AFTER TmunuBaseX_ZeroTmunu
    { }  "Add to T_munu here"
  ```

- `CarpetX/TmunuBaseX/interface.ccl` — the `eTtt`, `eTti{eTtx,eTty,eTtz}`,
  `eTij{eTxx…eTzz}` grid functions (covariant components T_μν; no `CENTERING`
  declaration ⇒ driver-default **vertex-centered** — the centering mismatch the
  8-point average of [MCNX-HYD-06] bridges).
- `CarpetX/CarpetX/src/schedule.cxx` — the declare-then-enforce validity machinery
  (READS presync/poison checks, CRC32 checksums on undeclared writes) that makes the
  contributor idiom's `READS:` + `WRITES:` declarations enforceable at runtime.
- `CarpetX/HydroBaseX/interface.ccl` — the cell-centered (`CENTERING={ccc}`) fluid
  variables MCNuX's centering matches, and the `temperature` declaration (line 19)
  carrying **no unit annotation** — the gap the MeV interface requirement closes.

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; its Eq. 16 is the MC stress-energy estimator restated in [MCNX-HYD-06] |
| arXiv:1806.02349v1 | Foucart et al. 2018, PRD 98, 063007 — registered here only for disambiguation: it *also* has an Eq. 16, a spectral fit, **not** the estimator (cite by bare arXiv id, never by author-year) |

The 8-point cell-to-vertex average of [MCNX-HYD-06] is **corpus-derived** (an
unweighted mean needs no external provenance) and is normative as written.

## Inputs & outputs

### MCNuX-owned source-term grid variables (the coupling interface)

All are `TYPE=gf`, **cell-centered** (`CENTERING={ccc}`, matching HydroBaseX), owned
and declared by MCNuX, and pinned by name (they are the published interface a partner's
`READS:` clauses and the golden data reference):

| Group | Variables | Content (geometrized code units) |
|---|---|---|
| `MCNuX::rad_force` | `rf_t`, `rf_x`, `rf_y`, `rf_z` | radiation four-force density on the fluid: `G^t` (coordinate-time component, index **up**) and `G_x, G_y, G_z` (spatial components, index **down**), per unit coordinate volume per unit coordinate time, **undensitized** (no √γ folded in) |
| `MCNuX::lepton_source` | `lep_src` | net electron-lepton **number** source density to the fluid (number per unit coordinate volume per unit coordinate time, same normalization) |

Sign convention: **positive = gained by the fluid**. Inputs: the interaction-event
stream (creations, absorptions, scatterings, diffusion advections) with per-event
packet-content changes; the cell coordinate volume ΔV and transport step Δt; for the
optional Tmunu contract, the metric quantities √γ and the lapse at the cell.

### Event ledger (what accumulates)

Per cell c and transport step, with `ΔP^t`, `ΔP_i`, `ΔL` the net change of total
packet content `Σ_k N_k p_k^t`, `Σ_k N_k p_{i,k}`, `Σ_k ℓ_{s_k} N_k` in c due to
interaction events:

| Event | contribution to ΔP^t, ΔP_i | contribution to ΔL |
|---|---|---|
| packet creation (emission) | `+N p^t`, `+N p_i` | `+ℓ_s N` |
| absorption (discrete or diffusion-test) | `−N p^t`, `−N p_i` | `−ℓ_s N` |
| elastic scattering / diffusion advection | `+N (p^t_out − p^t_in)`, `+N (p_{i,out} − p_{i,in})` | 0 |
| escape through the outer boundary | none (escape tallies, not source terms) | 0 |
| geodesic propagation | none (exchange with spacetime, not fluid) | 0 |

with `ℓ_s = (+1, −1, 0)` for (νe, ν̄e, νx); νx packets never contribute lepton number,
and the degeneracy factor g = 4 appears nowhere here (`N` already counts all lumped
states).

## Correctness requirements

- **[MCNX-HYD-01] Variable ownership and lifetime.** MCNuX declares and owns exactly
  the source-term grid variables of the table above, cell-centered, with complete
  `READS:`/`WRITES:` declarations on every routine touching them (per
  [carpetx-thorn-integration](./carpetx-thorn-integration.md) [MCNX-CTX-04]). Because
  the values tallied at step n are read by the fluid at step n+1, the variables carry
  state across steps: they must be part of checkpoint state, and a
  checkpoint/recovery cycle between transport and the partner's read must not alter
  what the partner reads.

- **[MCNX-HYD-02] Ledger definition (normative).** At the end of each transport step,
  each cell holds

  ```text
  G^t = −ΔP^t/(ΔV Δt) ,   G_i = −ΔP_i/(ΔV Δt) ,   S_ℓ = −ΔL/(ΔV Δt)
  ```

  with the event contributions of the ledger table and every event tallied into the
  cell containing it (the episode cell of
  [neutrino-matter-interactions](./neutrino-matter-interactions.md)). The minus sign
  makes the variables the **fluid-side** sources: what packets gain, the fluid loses.
  Worked signs: emission of a νe packet gives `S_ℓ < 0` (the fluid loses a lepton);
  absorption of a νe packet gives `S_ℓ > 0`. Tallies are accumulated with unordered
  atomics; cross-configuration reproducibility of these grid values is
  floating-point-reduction-tolerance only
  ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)
  [MCNX-RNG-08]).

- **[MCNX-HYD-03] Zero-then-add protocol (protocol requirement).** The source-term
  variables are **zeroed at one schedule point, accumulated into by contributors, and
  only then read** — never overwritten after the zero point. Concretely: at the start
  of MCNuX's transport group each step, one routine zeroes all source variables
  everywhere (`WRITES: …(everywhere)`); every deposition routine thereafter is a
  **contributor** in the TmunuBaseX sense — a self-accumulating `READS:` + `WRITES:`
  pair on the same variables, performing read-modify-write `+=` only, ordered after
  the zeroing (by group ordering, the reference protocol's structure); and no MCNuX
  routine assigns (rather than accumulates) a source variable after the zero point.
  The zeroing must not run between the end of the step's accumulation and the
  partner's read in the next step — the zero point lives at the *start* of the
  transport group, so the tallied values persist across the step boundary for the
  fluid's RHS reads.

- **[MCNX-HYD-04] Interface requirements on the GRMHD partner (stated abstractly, no
  consumer named).** Any GRMHD partner coupling to MCNuX must: (a) declare `READS:`
  on `MCNuX::rad_force` and `MCNuX::lepton_source` in the RHS routines that consume
  them (the driver's validity machinery enforces the declaration); (b) apply its own
  densitization, index raising/lowering, and — for a `ρYₑ`-type variable — its own
  baryon-mass conversion to the documented contents of the table above; (c) populate
  the end-of-step fluid state MCNuX reads (`rho`, `vel*`, `Ye`, and `temperature`),
  and in particular **populates HydroBaseX::temperature in MeV** —
  `CarpetX/HydroBaseX/interface.ccl` declares no unit, and this corpus pins MeV as
  the interface requirement — before MCNuX's transport group runs; and (d) never
  write any MCNuX-owned variable.

- **[MCNX-HYD-05] Conservation ledger closure.** Per transport step and per channel
  `X ∈ {P^t, P_x, P_y, P_z, L}`: the grid total must balance the packet total,

  ```text
  | Σ_cells (source value)·ΔV·Δt + ΔX_packets |  ≤  1e-13 · max(Σ_events |δX|, ε₀)
  ```

  i.e. what all packets gained through interaction events equals minus what the grid
  variables recorded, to relative `1e-13` per step (the `README.md` conservation
  tolerance), with the denominator the gross event activity (so quiet steps with
  near-total cancellation are judged against activity, not the net) and `ε₀` a
  documented tiny floor guarding the all-zero case.

- **[MCNX-HYD-06] Optional Tmunu contract (parameter-selected, default off).** When
  enabled, an MCNuX routine scheduled `IN TmunuBaseX_AddToTmunu` accumulates the MC
  radiation stress-energy into `eTtt`/`eTti`/`eTij`, as a contributor
  ([MCNX-HYD-03]'s idiom, on TmunuBaseX's variables through TmunuBaseX's own
  extension group). The estimator (restated, normative; provenance arXiv:1708.08452
  Eq. 16, index-lowered because TmunuBaseX carries covariant components):

  ```text
  T_{μν}^{(ν),MC} = Σ_{k=1}^{P} N_k p_{μ,k} p_{ν,k} / (√γ α p_k^t) · δ³(x^i − x_k^i)
  ```

  with α here the **lapse** (√γ α = √−g), `p_μ = (p_t, p_i)`,
  `p_t = −α² p^t + β^i p_i`, and the current (end-of-step) packet population.
  Discretized per cell, `δ³` becomes the per-cell sum over contained packets divided
  by the cell coordinate volume:

  ```text
  q_{μν}(c) = (1/ΔV) Σ_{k∈c} N_k p_{μ,k} p_{ν,k} / (√γ α p_k^t)
  ```

  Because `q_{μν}` is cell-centered and the TmunuBaseX variables are vertex-centered,
  the vertex contribution is the corpus-derived **8-point cell-to-vertex average** —
  the unweighted mean over the eight cells sharing the vertex (normative as written;
  an unweighted mean needs no external provenance):

  ```text
  eT(v) += Σ_{d∈{0,1}³} q(v − d)/8
  ```

  where `v − d` indexes, per component, the eight cells adjacent to vertex v. Energy
  units follow [conventions-and-units](./conventions-and-units.md) (geometrized,
  matching TmunuBaseX's convention of the stack). Do **not** confuse the estimator
  with arXiv:1806.02349 Eq. 16 (a spectral fit; see Source of truth).

- **[MCNX-HYD-07] Sole coupling channel (GRMHD-agnosticism, binding).** MCNuX's
  entire grid-variable output surface is: its own source-term variables
  ([MCNX-HYD-01]) plus, when enabled, the TmunuBaseX contribution of [MCNX-HYD-06].
  MCNuX writes **no** variable owned by any other thorn (writing directly into
  consumer-owned variables is forbidden), assumes nothing about the consumer beyond
  the interface requirements of [MCNX-HYD-04], and cites no production GRMHD code.

### Conventions restated (the subset this leaf uses)

- Species lepton numbers `ℓ = (+1, −1, 0)` for (νe, ν̄e, νx); g = 4 enters exactly
  once, at packet creation, never here.
- Signature `(−,+,+,+)`; `√−g = α√γ` with α the lapse (this spec's only α);
  `p_t = −α² p^t + β^i p_i`; packet state carries `(p^t, p_i)` per
  [geodesic-propagation](./geodesic-propagation.md)'s null closure.
- Geometrized code units for all grid variables, ΔV, Δt; conversions use exactly the
  pinned factor set of [conventions-and-units](./conventions-and-units.md).
- Grid tallies deposited by unordered atomics carry no bitwise cross-run/cross-rank/
  cross-device requirement; single-rank fixed-seed runs are exactly reproducible.

## Verification

- **Ledger closure (relative 1e-13 per step).** Every deterministic benchmark
  computes, each step, the five-channel balance of [MCNX-HYD-05] from an event-side
  audit sum and the grid-side total; the bound holds in every step of every benchmark.
- **Zero-then-add protocol (exact + driver-enforced).** A meta-test splits deposition
  into two artificial contributors; the result equals the single-contributor run to
  floating-point-reduction tolerance, and the driver's validity machinery
  (poison/checksum runs, per [carpetx-thorn-integration](./carpetx-thorn-integration.md))
  passes — proving declarations are complete and no assignment follows the zero
  point. Restart-recovery midway between transport and the partner-read point
  reproduces the partner-visible values ([MCNX-HYD-01]).
- **One-step lag (golden parity, 1e-12).** The cadence benchmark of
  [carpetx-thorn-integration](./carpetx-thorn-integration.md) shows a fluid-state
  step-change at iteration n first affecting `MCNuX::rad_force` /
  `MCNuX::lepton_source` reads at iteration n+1, archived as golden TSV at
  `ABSTOL=RELTOL=1e-12`.
- **Sign fixtures (exact).** Single-event unit checks: one νe emission in one cell
  yields `S_ℓ < 0` and `G^t < 0` in exactly that cell; one νe absorption reverses
  both signs; one elastic scattering changes momentum components only.
- **Tmunu closed forms (machine tier, `~1e-14`).** With one packet in one cell on
  Minkowski, the accumulated `q_{μν}` equals `N p_μ p_ν/(√γ α p^t ΔV)` to machine
  tier; a spatially uniform `q` field yields `eT(v) = q` exactly at every interior
  vertex; a linear-in-coordinates `q` field is reproduced exactly at vertices (the
  8-point mean is exact for linear fields); and the vertex-total identity
  `Σ_v (contribution) = Σ_c q(c)` holds to machine tier on a periodic/interior
  region (each cell contributes 1/8 to each of its 8 vertices).
- **MeV interface (test-side).** All corpus test backgrounds
  ([verification-suite-design](./verification-suite-design.md)) document and populate
  `HydroBaseX::temperature` in MeV, exercising the interface as specified; the
  requirement on external partners is contractual ([MCNX-HYD-04]) and enforced by
  documentation plus the coupled-configuration checklist, not by a runtime unit
  probe (none is possible — see Open questions).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the protocol phrase
`zeroed at one schedule point, accumulated into by contributors, and only then read`,
the corpus-derived average `eT(v) += Σ_{d∈{0,1}³} q(v − d)/8`, and the interface
phrase `populates HydroBaseX::temperature in MeV`; that its cited `CarpetX/` paths
resolve; and that its `[MCNX-HYD-NN]` ids are declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Deposition buffering and atomicity strategy (direct atomics vs per-tile buffers
  merged afterward), kernel fusion with the event loop — provided tallies land in the
  owning cells and [MCNX-HYD-05] holds
  ([particle-container-and-gpu](./particle-container-and-gpu.md) owns the idiom).
- Routine and schedule-group names, and how many contributor routines deposition is
  split across — subject to the zero/add/read ordering of [MCNX-HYD-03].
- Whether the event-side audit sums for [MCNX-HYD-05] are computed always or only in
  test configurations (the requirement is on the balance, not the bookkeeping cost).
- The Tmunu-enable parameter's name/default wiring, and whether `q_{μν}` is staged in
  a scratch variable before the vertex pass — observable `eT*` contributions must
  equal [MCNX-HYD-06].
- Internal diagnostics (per-species/per-channel breakdowns, escape-tally outputs)
  beyond the pinned interface variables.

## Open questions / assumptions

- **Assumption: component positioning of the four-force.** The interface stores
  `(G^t, G_x, G_y, G_z)` — time component with index up, spatial components with
  index down — because that matches the packet state the ledger actually sums
  (`p^t` from the null closure, `p_i` evolved) and the conserved-variable form
  typical fluid RHSs consume. The partner reconstructs any other component with its
  own metric ([MCNX-HYD-04]b). A different positioning would be a deliberate,
  golden-data-invalidating interface re-pin.
- **Assumption: undensitized, coordinate-volume normalization.** The stored sources
  are per coordinate volume and coordinate time with no √γ factor; the partner
  applies its own densitization. Recorded because both densitized and proper-volume
  conventions exist in the literature; the ledger equations of [MCNX-HYD-02] are the
  binding definition.
- **Assumption: `HydroBaseX::temperature` unit silence delegates to the coupled
  codes.** The declaration carries no unit; this corpus closes the gap with the MeV
  requirement of [MCNX-HYD-04] (convention pinned in
  [conventions-and-units](./conventions-and-units.md)). No runtime probe can verify a
  partner's unit choice (a number carries no unit tag); the requirement is
  contractual, and a mis-united partner manifests as grossly wrong opacities — the
  coupled-configuration checklist in
  [verification-suite-design](./verification-suite-design.md) carries the sanity
  check.
- **Assumption: one-step-lag splitting error accepted.** The n−1 → n source lag is
  first-order operator splitting, accepted by design per the cadence contract of
  [carpetx-thorn-integration](./carpetx-thorn-integration.md); tightening it (e.g.
  iterated coupling) would be a corpus change, not an implementation choice.
- **Open question: feedback subtraction for the optional Tmunu.** When [MCNX-HYD-06]
  is enabled *and* the partner also applies the four-force sources, the radiation
  four-momentum exchange is represented in both the fluid sources and the spacetime
  stress-energy — double counting at the level of the coupled Einstein–hydro system's
  bookkeeping. The baseline treats [MCNX-HYD-06] as an either/or diagnostic-grade
  capability (validity conditions and any simultaneous-use protocol to be specified
  here before production use with an evolving spacetime).
