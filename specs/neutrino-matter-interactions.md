# Neutrino–matter interactions

> Leaf spec. Self-contained: an agent can implement the absorption and elastic-scattering
> event machinery — interaction-time sampling, event competition, and event outcomes —
> from this file alone. It restates only the conventions it uses and references
> [conventions-and-units](./conventions-and-units.md) (units, constants, species),
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (the random-draw
> primitive and the statistical bar), and
> [opacity-eos-evaluation](./opacity-eos-evaluation.md) (the transport coefficients) for
> the rest. `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope

This spec pins the baseline interaction event set of MCNuX transport — **discrete
absorption** and **elastic (isoenergetic) scattering** — as it acts on propagating
packets:

- The **interaction-time sampling law** (the −ln(r) optical-depth draw of
  arXiv:1708.08452 Eq. 24, restated in full below in this corpus's [0, 1) draw
  convention).
- The **sampling-episode structure**: when interaction times are drawn and re-drawn,
  and how the two channels compete with each other and with the cell/step boundaries.
- The **event outcomes**: absorption removes the packet entirely; elastic scattering
  redraws the fluid-frame direction isotropically at fixed fluid-frame energy ν.
- The **frames and units** of the law, and the pinned zero-opacity behavior.
- The **RNG consumption map** for interaction events, against the pure function
  `u(S, q, e, k)` of [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md).

Out of scope:

- The transport coefficients themselves — tabulated assembly, the νx mapping, the
  analytic mode, range enforcement —
  [opacity-eos-evaluation](./opacity-eos-evaluation.md); this spec consumes
  `κ_a(s, E)`, `κ_s(s, E)` as given.
- **Inelastic channels** (NES, Pair, Brem — energy-changing scattering, pair
  processes): not in the baseline event set. Their consumption contract is pinned and
  their rate assembly is the structured open question of
  [opacity-eos-evaluation](./opacity-eos-evaluation.md); adding them as events would
  extend this spec, not silently appear.
- The **trapped regime**: the opacity relabeling η′/κ_a′/κ_s′ and the
  diffusion-approximation advection that replaces discrete scattering in
  scattering-dominated cells — [trapped-regime-treatment](./trapped-regime-treatment.md).
  When the relabeling is active, this spec's machinery runs unchanged on the primed
  coefficients.
- Free-streaming motion between events ([geodesic-propagation](./geodesic-propagation.md))
  and packet creation ([packet-representation-and-sampling](./packet-representation-and-sampling.md)).
- What events deposit into the source-term grid variables
  ([hydro-coupling-source-terms](./hydro-coupling-source-terms.md)); this spec defines
  the events, that spec defines the ledger.
- Device-kernel execution and container mechanics
  ([particle-container-and-gpu](./particle-container-and-gpu.md)); concrete benchmarks
  and golden data ([verification-suite-design](./verification-suite-design.md)).

## Source of truth

The sampling law and both event outcomes are restated in full below (normative per
`README.md`) from the pinned methods paper. Where this spec pins details the paper
leaves open (the [0, 1) draw mapping, the per-cell episode structure), the spec text is
normative and the deviation is logged under Open questions / assumptions.

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; its Eq. 24 is the interaction-time draw and its §2.3.3 the discrete absorption / elastic scattering outcomes restated below |

## Inputs & outputs

| Quantity | Role | Type / units |
|---|---|---|
| packet state `(x^i, p_i, N, s, q, e)` | input/output | per [packet-representation-and-sampling](./packet-representation-and-sampling.md); `p^t` from the null closure of [geodesic-propagation](./geodesic-propagation.md) |
| `κ_a(s, ν)`, `κ_s(s, ν)` | absorption / elastic-scattering opacity at the packet's cell state and fluid-frame energy | fluid-frame, geometrized code units (converted from cm⁻¹ per [conventions-and-units](./conventions-and-units.md)) |
| `u^μ` | cell-centered fluid four-velocity (tetrad input for the scattering redraw) | geometrized; `g_μν u^μ u^ν = −1` |
| `Δt_rem` | remaining coordinate time in the current transport step for this packet | double ≥ 0 |
| uniform draws | `u(S, q, e, k)` evaluations | `[0, 1)` per [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) |

Outputs: per packet and step, a sequence of zero or more events — each either an
**absorption** (packet removed; exchange recorded per
[hydro-coupling-source-terms](./hydro-coupling-source-terms.md)) or a **scattering**
(packet momentum redrawn; exchange recorded likewise) — interleaved with geodesic
propagation segments owned by [geodesic-propagation](./geodesic-propagation.md).

## Correctness requirements

- **[MCNX-INT-01] Interaction-time draw (restated, normative; provenance
  arXiv:1708.08452 Eq. 24).** At the start of each sampling episode ([MCNX-INT-02]),
  candidate coordinate-time intervals to the next scattering and absorption are drawn
  as

  ```text
  Δt_{s,a} = −ln(1 − u_{s,a}) p^t/(κ_{s,a} ν)
  ```

  with `u_s, u_a` two independent uniform draws in `[0, 1)` ([MCNX-INT-06] pins their
  `(e, k)` assignment). The paper's `−ln(r)` with `r ∈ (0, 1]` is realized as
  `r = 1 − u`: identical in distribution, and the logarithm is finite for every
  representable `u` (no `ln(0)` path exists). Frames and units, binding:

  - `κ_{s,a} = κ_{s,a}(s, ν)` are the **fluid-frame** opacities of the active scheme,
    evaluated at the episode cell's fluid state and at the packet's fluid-frame energy
    `ν = −p_μ u^μ` (converted to MeV for the coefficient lookup), then converted to
    geometrized code units (`κ[code] = κ[cm⁻¹] × 1.476625038e5`).
  - `p^t/ν` is evaluated with both in geometrized units (the ratio is dimensionless);
    it converts the fluid-frame interaction rate `κν` to coordinate time. `Δt_{s,a}`
    is a coordinate-time interval in code units.
  - `p^t`, `ν`, and the opacities are evaluated once at episode start and held fixed
    for the episode (see Open questions).

- **[MCNX-INT-02] Sampling episodes and event competition.** A new sampling episode
  begins, and both channels are re-drawn, exactly at these points of a packet's life:

  1. when the packet begins participating in a transport step (at its creation time,
     or at the step start for pre-existing packets);
  2. immediately after each scattering event;
  3. when the packet enters a new cell during its push.

  Every episode consumes both draws (`u_s` then `u_a`) regardless of the coefficient
  values. Within an episode the candidate event is the **shorter** of `Δt_s` and
  `Δt_a`, measured from the episode start; it fires only if it precedes both the
  packet's exit from the current cell and the end of the step (`Δt_rem`). If the cell
  boundary comes first, the packet crosses and a new episode begins (undrawn remainder
  discarded — statistically exact for piecewise-constant coefficients by
  memorylessness); if the step end comes first, the packet finishes the step with no
  event. A zero-opacity channel (`κ = 0`, e.g. νx absorption in table mode, or the
  transparency floor of [opacity-eos-evaluation](./opacity-eos-evaluation.md)) never
  fires: formally `Δt = +∞`; the implementation must realize this without producing
  NaN, and the channel's draw is still consumed. Ties (`Δt_s = Δt_a` bitwise) resolve
  to scattering.

- **[MCNX-INT-03] Discrete absorption.** When the absorption channel fires, the packet
  is **removed entirely (no continuous weight decay)**: its full content — energy,
  momentum, and lepton number `ℓ_s N` — is transferred to the fluid at the event cell
  per the ledger of [hydro-coupling-source-terms](./hydro-coupling-source-terms.md),
  and its id is retired ([particle-container-and-gpu](./particle-container-and-gpu.md)).
  Packet weights are never reduced by absorption; there is no path-length attenuation
  of `N` anywhere in the baseline.

- **[MCNX-INT-04] Elastic scattering.** When the scattering channel fires, the packet's
  four-momentum is redrawn isotropically in the fluid rest frame
  **at fixed fluid-frame energy ν** (isoenergetic in the fluid frame):

  ```text
  p^{μ′} = ν (1, sinθ cosφ, sinθ sinφ, cosθ) ,
  cosθ = 2u₁ − 1 ,   φ = 2π u₂
  ```

  — the same distribution and the same pinned tetrad transform (timelike leg
  `ê^μ_(0′) = u^μ`, Gram–Schmidt spatial legs in fixed coordinate order (x, y, z)) as
  emission sampling ([packet-representation-and-sampling](./packet-representation-and-sampling.md),
  whose [MCNX-PKT-04] owns the transform). The weight `N` and species `s` are
  unchanged; the recomputed momentum satisfies the null and ν-recovery identities to
  machine tier (`~1e-14` relative); the fluid-frame energy before and after agrees to
  the same tier. The coordinate-frame momentum change `N(p^μ_out − p^μ_in)` is
  exchanged with the fluid per
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md). After the event a
  new sampling episode begins ([MCNX-INT-02]).

- **[MCNX-INT-05] Coefficient provenance.** The opacities consumed by [MCNX-INT-01]
  come exclusively from the coefficient interface of
  [opacity-eos-evaluation](./opacity-eos-evaluation.md) (table or analytic source),
  evaluated at the cell-centered fluid state of the packet's current cell — never
  re-derived, re-interpolated sub-cell, or modified here. When the trapped-regime
  relabeling of [trapped-regime-treatment](./trapped-regime-treatment.md) is active,
  the primed coefficients `κ_a′`, `κ_s′` are substituted for `κ_a`, `κ_s` in
  [MCNX-INT-01] with no other change to this spec's machinery.

- **[MCNX-INT-06] RNG consumption map (interactions).** All draws in this spec are
  evaluations of the pure function `u(S, q, e, k)`; AMReX's stateful device RNG is
  forbidden for them. The pinned assignment: at the start of each sampling episode the
  packet's event counter `e` is incremented by 1 (creation is `e = 0`,
  [packet-representation-and-sampling](./packet-representation-and-sampling.md)); within
  that episode

  ```text
  k = 0 → u_s (scattering time) ,  k = 1 → u_a (absorption time)
  k = 2 → u₁ (cosθ) ,  k = 3 → u₂ (φ)      (consumed only if the episode ends in scattering)
  ```

  in exactly this order, and no other draws are consumed. Episodes that end without an
  event (cell crossing or step end) consume exactly draws `k = 0, 1`; the next episode
  uses the next `e`.

### Conventions restated (the subset this leaf uses)

- Species `{νe, ν̄e, νx}` with indices (0, 1, 2) and lepton numbers (+1, −1, 0); the
  degeneracy factor g = 4 appears nowhere in this spec (it enters exactly once, at
  packet creation).
- Fluid-frame energy `ν = −p_μ u^μ`; `p^t = √(γ^{ij} p_i p_j)/α` (null closure, owned
  by [geodesic-propagation](./geodesic-propagation.md); `α` there is the lapse).
- Geometrized units for `Δt`, `p^t`, `ν` (as used in the ratio), and κ; the microphysics
  lookup takes ν in MeV; conversions use exactly the pinned factor set of
  [conventions-and-units](./conventions-and-units.md).
- All stochastic draws are `u(S, q, e, k)` evaluations with uniform variates in
  `[0, 1)`; statistical checks use the 4σ bar at pinned packet counts and seeds of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md).

## Verification

- **Beam attenuation (statistical, 4σ).** A collimated monoenergetic packet beam
  crossing a uniform pure-absorber slab (analytic opacity mode, `κ_s = 0`): the
  transmitted number fraction equals `exp(−κ_a L)` (L the slab thickness, fluid at
  rest) within 4σ at the pinned packet count and seed.
- **Collision statistics (statistical, 4σ).** In a uniform pure-scattering box
  (`κ_a = 0`), the number of scattering events per packet per fixed path length is
  Poisson with mean `κ_s ℓ_path`: sample mean and variance within 4σ of the Poisson
  values.
- **Scattering isotropy (statistical, 4σ).** Post-scattering fluid-frame directions
  over many events: first moments of `cosθ` and `φ` within 4σ of their isotropic
  nulls (0 and π respectively).
- **Elasticity (machine tier, `~1e-14`).** For every scattering event in a pinned
  deterministic run, the fluid-frame energy recomputed from the outgoing momentum
  equals the incoming ν to `1e-14` relative, and the outgoing momentum satisfies the
  null identity to the same tier.
- **Discreteness (exact).** In every benchmark: no packet's weight `N` ever changes
  after creation; every absorption removes exactly one whole packet; the packet count
  changes only by creations, absorptions, and escapes.
- **Draw-map audit (exact).** A trace of `(q, e, k)` tuples consumed during a pinned
  small run matches [MCNX-INT-06] exactly: two draws per eventless episode, four per
  scattering episode, episode boundaries exactly at creation/step-start, post-scatter,
  and cell crossing.
- **Determinism (exact / golden parity 1e-12).** Fixed-seed single-rank runs of the
  interaction benchmarks reproduce per-packet event sequences bitwise; golden TSV/norm
  outputs archived at the Cactus harness defaults `ABSTOL=RELTOL=1e-12`
  ([verification-suite-design](./verification-suite-design.md) pins the concrete
  benchmarks).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the restated-equation anchor
`Δt_{s,a} = −ln(1 − u_{s,a}) p^t/(κ_{s,a} ν)`, the discrete-absorption phrase
`removed entirely (no continuous weight decay)`, and the elastic-outcome phrase
`at fixed fluid-frame energy ν`; and that its `[MCNX-INT-NN]` ids are declared here
and covered by the [verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Kernel structure: whether episode sampling, the push, and event application are fused
  in one kernel or staged; batching of events; how `Δt_rem` and cell-exit times are
  tracked — provided the episode structure and draw order of [MCNX-INT-02] /
  [MCNX-INT-06] are observed exactly.
- How the `Δt = +∞` zero-opacity case is represented (sentinel, branch, clamped
  inverse) — observable behavior is "channel never fires, draw still consumed, no
  NaN".
- Caching of per-cell coefficient evaluations across the packets in a cell
  (values must equal the per-packet evaluation of [MCNX-INT-05]).
- How removed packets are compacted/retired within a step, subject to the id-retirement
  and container contracts of
  [particle-container-and-gpu](./particle-container-and-gpu.md).

## Open questions / assumptions

- **Resolution logged: per-cell episode structure (deviation from the source paper's
  letter, spec text normative per `README.md`).** arXiv:1708.08452 draws interaction
  times from the opacities at the packet position at the start of the step;
  this corpus pins re-drawing at every cell crossing ([MCNX-INT-02]). For the
  piecewise-constant (cell-centered) coefficient fields MCNuX actually has, the
  exponential's memorylessness makes the per-cell scheme statistically exact, and it
  removes the start-of-step-state bias for packets crossing steep opacity gradients.
  The draw-count consequences are pinned by [MCNX-INT-06]; golden data depends on this
  choice.
- **Assumption: frozen episode kinematics.** `p^t`, ν, and the opacities are evaluated
  at episode start and held fixed for the episode, although `p^t` drifts along the
  geodesic within the cell on curved metrics. The error is O(cell size) and vanishes
  in the convergence limit; re-evaluating mid-episode would be a spec change (it
  alters golden data), not an implementation choice.
- **Assumption: ties resolve to scattering.** `Δt_s = Δt_a` at full double precision
  has probability ~2⁻⁵³ per episode; the tie rule exists so the branch is deterministic,
  not because it is statistically consequential.
- **Assumption: absorption position.** The exchange of an absorption event is tallied
  to the episode cell (the cell whose coefficients drew the event), evaluated at the
  event time; sub-cell position detail does not affect the cell-indexed ledger of
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md).
