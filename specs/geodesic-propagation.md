# Geodesic propagation

> Leaf spec. Self-contained: an agent can implement the packet geodesic push on the
> dynamical 3+1 metric from this file alone. It restates only the conventions it uses
> and references [conventions-and-units](./conventions-and-units.md) for the rest.
> `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope

This spec pins the free-streaming transport of Monte Carlo packets between interaction
events: the null-geodesic equations of motion on the 3+1 metric, the closure of the
evolved state under the null constraint, the contract for gathering the metric (and its
spatial derivatives) from the ADMBaseX grid functions at packet positions, the required
integration accuracy, and the conserved-quantity behavior on stationary backgrounds
that anchors verification.

Out of scope:

- When a push ends early because an interaction event fires — event times and outcomes
  ([neutrino-matter-interactions](./neutrino-matter-interactions.md)); the
  diffusion-approximation advection replacing the push in scattering-dominated cells
  ([trapped-regime-treatment](./trapped-regime-treatment.md)).
- The packet state's components, units, and creation
  ([packet-representation-and-sampling](./packet-representation-and-sampling.md)); this
  spec evolves `(x^i, p_i)` and defines `p^t`.
- How grid data is reached (driver objects, `Array4` views, scheduling) —
  [carpetx-thorn-integration](./carpetx-thorn-integration.md); how the gather/push
  execute as device kernels and the after-step ownership/redistribution invariants —
  [particle-container-and-gpu](./particle-container-and-gpu.md).
- Concrete benchmark parameter files, resolutions, and numeric golden tolerances —
  [verification-suite-design](./verification-suite-design.md) (including the
  Schwarzschild background provided by `TestMCNuX`); this spec pins the invariants and
  convergence orders those benchmarks check.
- The cadence at which the metric is refreshed (once per coarsest-level Δt, end-of-step
  state — the observable cadence contract of
  [carpetx-thorn-integration](./carpetx-thorn-integration.md)). Within one transport
  step the metric is a fixed snapshot; this spec treats it as given.

## Source of truth

The equations of motion are restated in full below (normative per `README.md`) with
provenance arXiv:1708.08452 Eqs. 25–26 (the Hughes et al. 1994 prescription as adopted
by the GR-MC methods paper). Metric-field provenance:

- `CarpetX/ADMBaseX/interface.ccl` — the metric group `gxx…gzz` (`metric`), lapse `alp`
  (`lapse`), and shift `betax…betaz` (`shift`) grid-function declarations MCNuX reads.
  These groups carry **no `CENTERING` declaration**, so they take the CarpetX default:
  **vertex-centered** (`CarpetX/CarpetX/src/driver.cxx`, the "Default: vertex-centred"
  branch of the group-index-type resolution).

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — GR-MC methods paper; Eqs. 25–26 are the geodesic push restated below |

## Inputs & outputs

| Quantity | Role | Type / units |
|---|---|---|
| `x^i` | packet position (input/output) | 3 × double, geometrized coordinates, inside the domain |
| `p_i` | packet lower spatial momentum (input/output) | 3 × double, geometrized units |
| `α, β^i, γ_ij` | lapse, shift, spatial metric at the packet position | gathered from ADMBaseX vertex-centered grid functions |
| `∂_i α, ∂_i β^k, ∂_i γ_jk` | spatial derivatives at the packet position | evaluated per [MCNX-GEO-03] |
| `Δt_push` | coordinate-time push interval | double > 0 (a full transport step or the remainder to an event) |

Outputs: updated `(x^i, p_i)` after propagating for `Δt_push` along the null geodesic.
`p^t` is not independent state: it is recomputed from the null constraint
([MCNX-GEO-02]) wherever it is needed. A packet whose trajectory leaves the outer
domain boundary during a push is removed from transport and accumulated into escape
tallies (per-species energy and number), per the `README.md` boundary policy.

## Correctness requirements

- **[MCNX-GEO-01] Equations of motion (restated, normative).** Between interaction
  events a packet propagates along a null geodesic of the 3+1 metric
  (provenance: arXiv:1708.08452 Eqs. 25–26):

  ```text
  dx^i/dt = γ^{ij} p_j/p^t − β^i
  dp_i/dt = −α (∂_i α) p^t + (∂_i β^k) p_k − ½ (∂_i γ^{jk}) p_j p_k/p^t
  ```

  with `t` the coordinate time, `γ^{ij}` the inverse spatial metric, and all metric
  quantities evaluated at the packet's current position (time-frozen snapshot within
  the step). `∂_i γ^{jk} = −γ^{ja} (∂_i γ_ab) γ^{bk}` relates the inverse-metric
  derivative to the gathered `∂_i γ_ab`.

- **[MCNX-GEO-02] Null closure.** The evolved state is `(x^i, p_i)`; `p^t` is closed by
  the null condition `g^{μν} p_μ p_ν = 0`, i.e.

  ```text
  p^t = √(γ^{ij} p_i p_j)/α
  ```

  recomputed from the current `(p_i, γ^{ij}, α)` at every right-hand-side evaluation
  (and wherever `p^t` is consumed). The null constraint is thus enforced by
  construction and cannot drift; the identity `α² (p^t)² = γ^{ij} p_i p_j` must hold to
  machine tier (`~1e-14` relative) at every evaluation point.

- **[MCNX-GEO-03] Metric gather contract.** `α`, `β^i`, `γ_ij` are interpolated from
  the ADMBaseX vertex-centered grid functions to the packet position by a documented
  interpolation that reproduces constant and linear fields **exactly** (machine tier)
  — trilinear is the minimum admissible order. The spatial derivatives
  `∂_i α, ∂_i β^k, ∂_i γ_jk` are evaluated by a documented scheme (differentiating the
  interpolant, or interpolating finite-difference derivatives of order ≥ 2) that is
  exact on linear fields (constant derivative fields recovered to machine tier). The
  same gathered snapshot (same fields, same scheme) is used for every term of one
  right-hand-side evaluation — no mixing of gather schemes within a step.

- **[MCNX-GEO-04] Integration accuracy.** The time integrator for [MCNX-GEO-01] is
  free but must be at least **second-order convergent** in `Δt` on smooth metrics:
  halving the step size reduces the error in `(x^i, p_i)` against a reference solution
  by a measured factor consistent with order ≥ 2 (within the tolerance the benchmark
  pins). In flat spacetime (`α = 1, β^i = 0, γ_ij = δ_ij`) the right-hand side of
  `p_i` vanishes identically and propagation must be exact to machine tier per step:
  straight lines at coordinate speed `|dx/dt| = 1`.

- **[MCNX-GEO-05] Stationary-background conservation.** On any metric with
  `∂_t g_μν = 0` (Minkowski; the Schwarzschild background of
  [verification-suite-design](./verification-suite-design.md)), the covariant energy

  ```text
  p_t = −α² p^t + β^i p_i
  ```

  is exactly conserved along the continuum geodesic. The implementation must conserve
  each packet's `p_t` to the accuracy implied by its integration order: the drift
  `|p_t(T) − p_t(0)|/|p_t(0)|` over a benchmark trajectory converges to zero at order
  ≥ 2 in `Δt`, and at the benchmark-pinned resolution stays below the numeric
  tolerance pinned per benchmark in
  [verification-suite-design](./verification-suite-design.md). This is the verifiable
  anchor of this spec: it exercises the gather, the derivatives, and the integrator
  together.

### Conventions restated (the subset this leaf uses)

- Signature `(−,+,+,+)`; lapse `α`, shift `β^i`, spatial metric `γ_ij`;
  `γ = det(γ_ij)`; Latin indices spatial, raised/lowered with `γ^{ij}`/`γ_ij` for
  spatial objects; `p^μ` is the coordinate-frame four-momentum.
- Geometrized units `G = c = M_sun = 1` for all quantities in this spec; coordinate
  speed of light is 1.
- The symbol `α` is the lapse throughout this spec (not the trapped-regime relabeling
  fraction, not a spectral index).

## Verification

- **Flat-spacetime exactness (machine tier).** On Minkowski (`ADMBaseX::initial_data =
  "Cartesian Minkowski"`), pushed packets travel straight lines with `p_i` constant
  bitwise-stable per step up to roundoff (`~1e-14` relative over a pinned trajectory)
  and `|dx/dt| = 1` to machine tier.
- **Null-closure identity (machine tier).** At every right-hand-side evaluation of a
  pinned test trajectory, `α² (p^t)² = γ^{ij} p_i p_j` to `1e-14` relative.
- **Linear-field gather exactness (machine tier).** With synthetic metric grid data
  linear in the coordinates, gathered values and their derivatives reproduce the
  closed forms to `1e-14` relative at arbitrary in-cell positions.
- **Convergence study (order ≥ 2).** On the Schwarzschild background, integrate a
  pinned set of trajectories at `Δt, Δt/2, Δt/4`; the observed convergence order of
  the `p_t` drift and of the position error is ≥ 2.
- **Schwarzschild `p_t` conservation (golden parity, 1e-12 harness).** The
  deterministic benchmark of [verification-suite-design](./verification-suite-design.md)
  (fixed initial packets, no interactions, `TestMCNuX` Schwarzschild data) checks the
  `p_t` drift against its per-benchmark numeric tolerance and archives golden TSV
  output at the Cactus harness defaults `ABSTOL=RELTOL=1e-12`.

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the restated-equation anchors `dx^i/dt = γ^{ij} p_j/p^t − β^i`,
`dp_i/dt = −α (∂_i α) p^t + (∂_i β^k) p_k − ½ (∂_i γ^{jk}) p_j p_k/p^t`, and
`p^t = √(γ^{ij} p_i p_j)/α`; that its cited `CarpetX/` paths resolve; and that its
`[MCNX-GEO-NN]` ids are declared here and covered by the
[verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- The integrator (RK2, RK4, symplectic, adaptive substepping) and step-subdivision
  strategy, subject to [MCNX-GEO-04].
- The gather interpolation order beyond the trilinear minimum, and whether derivatives
  come from differentiating the interpolant or from interpolated finite differences,
  subject to [MCNX-GEO-03].
- Caching of gathered metric quantities, inverse-metric computation strategy
  (analytic 3×3 inverse vs decomposition), kernel structure and fusion.
- Whether `p^t` (or `1/p^t`) is cached within a right-hand-side evaluation — observable
  values must equal recomputation from the null closure.

## Open questions / assumptions

- **Assumption: ADMBaseX default centering.** The metric/lapse/shift groups declare no
  `CENTERING`, and the CarpetX driver's documented default is vertex-centered; this
  spec relies on that default (provenance recorded in Source of truth). If a future
  ADMBaseX revision declared an explicit non-vertex centering, the gather contract's
  stencil bookkeeping — not its exactness requirements — would need re-derivation.
- **Assumption: metric smoothness at refinement boundaries.** The gather treats the
  grid data as smooth enough for its interpolation order; prolongation artifacts at
  AMR refinement boundaries degrade accuracy locally but are not a correctness-contract
  violation as long as the convergence benchmarks (run on unigrid backgrounds) pass.
  Redistribution timing across levels is owned by
  [particle-container-and-gpu](./particle-container-and-gpu.md).
- **Assumption: time-frozen metric within a step.** The push uses the end-of-step
  metric snapshot for the whole transport step (operator splitting per the cadence
  contract of [carpetx-thorn-integration](./carpetx-thorn-integration.md)); the
  `O(Δt)` splitting error this introduces on dynamical spacetimes is accepted by
  design and is not counted against [MCNX-GEO-04]'s convergence order, which is
  measured on stationary metrics.
