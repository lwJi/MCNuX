# Geodesic propagation (domain leaf: moving packets through curved spacetime at κ = 0)

> Domain leaf spec. Self-contained: an agent can implement free-streaming packet transport from this file alone, referencing `conventions-and-units.md` for the metric signature, 3+1 variable naming, and packet nomenclature, and `rng-and-statistical-acceptance.md` only for the pinned seeds its benchmarks reuse (propagation itself consumes no random draws). The tolerance tiers named here (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict.

## Purpose & scope

This spec defines how a Monte Carlo packet moves between interaction events: the 3+1 null-geodesic equations advancing (x^i, p_i) in coordinate time, the null-condition enforcement that makes p^t and p_t derived quantities, the metric evaluation and freezing policy, the substep caps (end of transport step, next event time, cell crossing), outer-boundary escape behavior, and the analytic benchmarks (flat spacetime, Schwarzschild p_t conservation) with numeric tolerances.

In scope:

- The restated geodesic equations and the null condition (normative forms below).
- The metric time-slice policy (which slice's metric a transport step uses) and the permitted spatial freezing approximation.
- Substep caps as observable correctness requirements, including the bounded-motion fact the particle-container spec's redistribution contract relies on.
- Outer-boundary escape removal and its exact accounting.
- The flat-spacetime and Schwarzschild free-streaming benchmarks with concrete pass/fail numbers.

Out of scope:

- Emission sampling and the packet state contract's creation-side requirements (see `packet-representation-and-sampling.md`; the p_t identity is shared state algebra, owned here).
- Sampling of absorption/scattering events and their outcomes (the neutrino–matter interactions spec); this spec only honors a supplied next-event time as a substep cap.
- The trapped-regime modification of transport (the trapped-regime spec).
- Container mechanics: how packets are stored, redistributed across ranks/levels, or deposited (the particle-container spec); this spec produces the per-substep bounded-motion guarantee that contract consumes.
- Schedule placement and cadence (the CarpetX-integration spec); this spec takes "one transport step of size Δt on the metric of slice t_n" as given.

## Source of truth

The equations restated under "Inputs & outputs" are **normative**; citations are provenance only (restate-and-pin, per `README.md`):

- 3+1 null-geodesic equations and null-condition enforcement: Foucart et al. 2021 [arXiv:2103.16588](https://arxiv.org/abs/2103.16588) Eqs. 12–14; Foucart 2018 [arXiv:1708.08452](https://arxiv.org/abs/1708.08452) Eqs. 25–26. Note: the 1/(2p^t) factor on the metric-gradient term of dp_i/dt below is required for the equations to be homogeneous of degree one in p_μ (t-parametrized geodesics are invariant under p_μ → λp_μ) and follows from the Hamiltonian H = ½ g^{μν} p_μ p_ν. The arXiv rendering of 2103.16588 Eq. 13 omits this factor (verified against the paper's LaTeX source — a typo); 1708.08452 Eq. 26 includes it. The form written here, with the 1/(2p^t), is normative.
- Metric-frozen-per-cell policy, step caps, and the Schwarzschild free-streaming p_t test design: Foucart et al. 2021 §2 and §3.2 (single emitting region, κ_a = κ_s = 0, per-packet p_t drift "well below 1%" for most packets at merger-typical resolution, decreasing with resolution).
- `CarpetX/ADMBaseX/interface.ccl` — the metric grid variables consumed: vertex-centered `alp` (α), `betax/betay/betaz` (β^i), `gxx…gzz` (γ_ij) on the current time slice.
- `conventions-and-units.md` — signature (−,+,+,+), 3+1 naming, packet nomenclature (x^i in code units, p_μ in MeV), and the fact that the geodesic equations are invariant under constant rescaling of p_μ (so MeV momenta against code-unit coordinates are consistent).
- The Kerr–Schild form of the Schwarzschild metric restated in "Verification" is normative for the benchmark; standard-textbook provenance.

## Inputs & outputs

### Interface

One transport step advances a packet over coordinate-time interval [t_n, t_n + Δt], decomposed into substeps. Inputs:

- Packet state (x^i, p_i) (and stored p_t), per the state contract of `packet-representation-and-sampling.md`; propagation consumes **no** uniform draws and leaves the event counter e unchanged.
- Metric fields α > 0, β^i, γ_ij (positive-definite) of time slice t_n, vertex-centered per `CarpetX/ADMBaseX/interface.ccl`, together with whatever spatial-derivative data the implementation derives from them (finite differences or differentiated interpolants — implementation freedom bounded by MCNX-GEO-02/06).
- The transport step Δt and, when the interactions spec is active, the packet's next-event coordinate time t_ev.

Outputs: updated (x^i, p_i, p_t) at t_n + Δt (or at the event/boundary that ended transport early), with p^t and p_t recomputed from p_i as below; for escaped packets, removal plus an exact escape record (species s, weight N, fluid-frame energy ε, exit time).

Valid ranges: √(γ^{ij} p_i p_j) > 0 (a packet with vanishing momentum is invalid); α > 0 everywhere transport runs; packets are only propagated at positions covered by metric data.

### Restated equations (normative)

The null condition defines p^t from p_i and is **enforced algebraically, never integrated**:

```text
p^t = √(γ^{ij} p_i p_j) / α
```

and the covariant time component follows as

```text
p_t = β^i p_i − α² p^t  =  β^i p_i − α √(γ^{ij} p_i p_j)
```

(p^t > 0 for future-directed packets; p_t < 0 wherever α² > β_i β^i). The coordinate-time evolution equations are

```text
dx^i/dt = γ^{ij} p_j / p^t − β^i
dp_i/dt = −α p^t ∂_i α + p_j ∂_i β^j − (p_j p_k / (2 p^t)) ∂_i γ^{jk}
```

Provenance: Foucart 2021 Eqs. 12–14 / Foucart 2018 Eqs. 25–26, with the homogeneity note under "Source of truth". Only (x^i, p_i) are integrated; p^t (and hence p_t) is recomputed from the null condition at every evaluation point.

### Metric policy

- **Time-slice policy (hard requirement).** All metric data used anywhere in transport step n — including substeps and derivative evaluations — comes from slice t_n; the metric is never updated mid-step and no data from slice t_{n+1} may influence step n's trajectories (MCNX-GEO-03).
- **Spatial evaluation and freezing (bounded approximation).** How α, β^i, γ_ij and their gradients are evaluated at packet positions is implementation freedom: pointwise interpolation, or the frozen-per-cell policy of the provenance (metric and gradients held constant over the cell a packet occupies during a substep). Whatever policy is chosen, the discretization error it induces is bounded observably by the Schwarzschild benchmark and its convergence requirement (MCNX-GEO-06).

### Substep caps and boundary behavior

Each substep of size Δt_sub must satisfy all of:

1. **Step cap** — the substep never advances past t_n + Δt.
2. **Event cap** — the substep never advances past the packet's next-event time t_ev when one is set (contract honored now, exercised once the interactions spec lands).
3. **Cell-crossing cap** — the per-substep coordinate displacement satisfies |Δx^i| ≤ one cell width in every direction i. This keeps the frozen-metric approximation consistent and is the bounded-motion precondition the particle-container spec's local redistribution relies on.

A packet whose substep crosses the outer domain boundary is removed from transport at that substep (escape) and recorded exactly (species, N, ε, exit time); escaped packets are never silently dropped.

## Correctness requirements

- **[MCNX-GEO-01] Null-condition enforcement (machine).** p^t is always computed from p_i via `p^t = √(γ^{ij} p_i p_j) / α` and never carried as an independent integrated variable; consequently, at every substep of every packet, |g^{μν} p_μ p_ν| / (α p^t)² ≤ `1e-14` when evaluated with the metric values used at that substep. Reference: the null condition above; provenance Foucart 2021 Eq. 12.
- **[MCNX-GEO-02] Time-integration order ≥ 2 (numeric convergence).** The substep integrator is at least second-order accurate in Δt_sub: in the analytic-metric halving test of "Verification" (Kerr–Schild evaluated analytically, so grid discretization is excluded), the p_t drift satisfies e(Δt)/e(Δt/2) ≥ 3.5 for both tested halvings (observed order ≥ 1.8, asymptotically 2 for the reference RK2). Provenance: second-order RK per Foucart 2021 §2.
- **[MCNX-GEO-03] Metric time-slice policy (golden `1e-12`).** Transport during step n uses the metric of slice t_n only, never updated mid-step. Verified as a fixed-seed golden regression on a fixture whose lapse varies by O(`1e-2`) per step: any policy deviation (mid-step update, wrong slice) shifts per-packet trajectories orders of magnitude above the harness bar `ABSTOL = RELTOL = 1e-12`.
- **[MCNX-GEO-04] Substep caps, bounded motion, and escape accounting (exact).** In instrumented fixed-seed runs: zero substeps advance past t_n + Δt; zero substeps advance past a set t_ev; zero substeps have |Δx^i| > one cell width in any direction; and the number and total (N, ε) content of escape records equals exactly the packets removed at the outer boundary (integer and bookkeeping identity, escape totals summed at relative `1e-10`). Reference: the caps above.
- **[MCNX-GEO-05] Flat-spacetime propagation (machine).** On a Minkowski grid metric (α = 1, β^i = 0, γ_ij = δ_ij), packets travel straight lines at coordinate speed 1: x^i(t) = x^i(0) + (p_i/|p|) t with p_i and p_t constant. Pass/fail on the fixture of "Verification" (256 substeps, path length ≈ 1.5 code units): final-position error ≤ `1e-12` code units per component; |p_i(end) − p_i(0)| ≤ `1e-14` |p|; |p_t(end) − p_t(0)| ≤ `1e-14` |p_t|. Closed-form exact answer → machine tier (bound is N_steps × a few ULP with margin).
- **[MCNX-GEO-06] Schwarzschild p_t conservation with resolution convergence (numeric).** On the Kerr–Schild Schwarzschild grid-metric benchmark of "Verification" (stationary metric ⇒ p_t exactly conserved along exact geodesics), with δ_k = |p_t(end) − p_t(0)|/|p_t(0)| per scored packet: at the baseline resolution Δx = M/8, **all** scored packets satisfy δ ≤ `2e-2` and ≥ 99% satisfy δ ≤ `1e-2`; and across Δx = M/4 → M/8 → M/16 the median δ decreases monotonically with median(M/16) ≤ median(M/4)/1.5. Provenance for the bound-plus-convergence framing (drift "well below 1%" for most packets, slow but clear convergence): Foucart 2021 §3.2.

## Verification

### Flat-spacetime fixture (MCNX-GEO-05, family B unit check)

Domain [−2, 2]³ with Minkowski metric data on the vertex grid (Δx = 1/16); ≥ 100 packets launched from the origin with unit Eulerian-frame energy and isotropic directions drawn with the corpus generator at seed S0 (test scaffolding, not the emission path); 256 uniform substeps to t = 1.5. Assert the MCNX-GEO-05 bounds per packet. Because the metric is exactly constant, all derivative terms vanish identically and momentum drift beyond a few ULP indicates a wrong term, not rounding.

### Kerr–Schild Schwarzschild benchmark (MCNX-GEO-01/06, family B unit check)

The benchmark metric, restated (normative for the test; M is the mass, code units, M = 1):

```text
H = M/r,   l_i = x_i/r   (r² = x² + y² + z²)
α = 1/√(1 + 2H)
β^i = 2H l^i / (1 + 2H)          (equivalently β_i = 2H l_i)
γ_ij = δ_ij + 2H l_i l_j
γ^{ij} = δ_ij − 2H l_i l_j / (1 + 2H)
√γ = √(1 + 2H)
```

Configuration: domain [−16M, 16M]³; vertex-centered grid metric filled from the expressions above at each of Δx ∈ {M/4, M/8, M/16}; N_smoke = `1e4` packets launched from (10M, 0, 0) with unit Eulerian-frame energy (α p^t = 1 at launch) and isotropic directions with respect to the local γ-orthonormal triad, seed S0; propagate to t = 30M under the substep caps. A packet reaching r < 3M is removed (captured) and excluded from scoring; survivors are scored on δ_k. Assert the MCNX-GEO-06 bounds and, along every scored trajectory at every substep, the MCNX-GEO-01 null identity at `1e-14`.

### Integrator-order fixture (MCNX-GEO-02, family B unit check)

Same spacetime, but α, β^i, γ_ij and their gradients evaluated **analytically** at packet positions (isolating time integration from grid discretization): one packet from (10M, 0, 0), initial direction the local ŷ (tangential) axis of the γ-orthonormal triad, propagated to t = 20M with uniform substeps Δt_sub ∈ {M/5, M/10, M/20}; e(Δt_sub) = |Δp_t|. Assert e(M/5)/e(M/10) ≥ 3.5 and e(M/10)/e(M/20) ≥ 3.5.

### Policy and caps fixtures (MCNX-GEO-03/04)

- MCNX-GEO-03: a family-A Cactus regression: fixed seed, single-rank CPU, a time-dependent analytic lapse α(t) = 1 − `1e-2`·sin(t) fixture; committed golden per-packet TSV diffed at `ABSTOL = RELTOL = 1e-12` (per the harness contract in `verification-suite-design.md`).
- MCNX-GEO-04: instrumented fixed-seed runs (flat and Schwarzschild fixtures above) asserting zero cap violations and exact escape accounting; the event-cap leg activates when the interactions spec lands and is reported SKIP (loudly) until then.

All requirement ids above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).

Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` path, and contains the required claim strings (the null condition, both evolution equations, the p_t expression, the Kerr–Schild benchmark, the drift tolerances).

## Implementation freedom

- **The integrator** — RK2 (the provenance reference), RK3/RK4, or any scheme meeting MCNX-GEO-02 — and the substep-size controller within the caps.
- **Metric spatial evaluation** — interpolation order and stencil, frozen-per-cell vs. pointwise, how gradients are formed (finite differences of grid data or differentiated interpolants) — bounded by MCNX-GEO-06's drift and convergence requirements.
- Whether p_t is stored and refreshed or recomputed on demand (the consistency identity of `packet-representation-and-sampling.md` governs either way), and any caching of per-cell metric quantities.
- Loop structure, batching, kernel fusion with other per-packet work, and precision of intermediates — provided per-packet observables meet the stated tiers.
- How the escape record is represented (per-packet log, per-step reduction), provided the exact accounting of MCNX-GEO-04 is observable.

## Open questions / assumptions

- **Kerr and orbit-count benchmarks are deferred (assumption, non-blocking).** The provenance validates on Schwarzschild only (Foucart 2021 §3.2); this corpus does the same. A Kerr photon-orbit benchmark would strengthen coverage and may be added here as a spec change without affecting existing requirements.
- **Inner excision policy (assumption).** Production spacetimes with black holes need a capture/excision radius inside which packets are removed (the benchmark uses r < 3M). The production excision rule (tied to the apparent horizon or an excision mask) is not yet pinned; until it is, removal-with-exact-accounting mirroring the escape contract is assumed.
- **Convergence bar is deliberately loose (assumption).** median(M/16) ≤ median(M/4)/1.5 encodes the provenance's "slow but clear" convergence; a first-order-accurate frozen-metric policy should beat it (factor ≈ 4 over two halvings). If measured convergence is reliably faster, the bar may be tightened here — never loosened without a documented cause.
- **Time-interpolated metric (recorded improvement).** Using a half-step or time-interpolated metric within a transport step could reduce the O(Δt) freezing error, but would change MCNX-GEO-03's observable policy; it is a spec change here if ever adopted, not an implementation freedom.
- **Subcycled fine-level transport is out of scope here** (single-rate stepping is the verified configuration; the cadence contract and its open question live in the CarpetX-integration spec).
