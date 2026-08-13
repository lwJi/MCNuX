# TODO — MCNuX build plan

Greenfield. The repo is at the stub milestone: implementation is 0% for all 12 coded specs. `MCNuX/src/stub.cxx` is the only translation unit (two `static_assert`s = half of `[MCNX-BLD-03]`); all CCL files are comment-only skeletons except the already-correct `configuration.ccl` REQUIRES lines and the pinned ThornList `agent_scripts/sandbox/mcnux.th` — do not re-touch those. `MCNuX/test/` does not exist. There is no dead code or remediation debt: every task below is new construction, ordered along the dependency spine (list order = priority = build order).

## Standing facts every item inherits

- Tolerance tiers (`specs/README.md`): machine `1e-14`; golden parity `1e-12` (Cactus harness `ABSTOL=RELTOL=1e-12`); relaxed `1e-10` (stated rationale required); exact = bitwise/integer/NaN-equality; statistical = **4σ** at `N_p = 1048576 = 2^20`, primary seed `1296518744`, secondary seed `20260730`; conservation ledger relative `1e-13` per step per channel.
- Any new file under `MCNuX/src/` or `TestMCNuX/src/` must be added to that thorn's `make.code.defn` `SRCS`/`SUBDIRS` or it silently does not compile.
- Every benchmark and `stats-*` parfile must satisfy the activation floor now pinned in `[MCNX-VER-06]`: `CarpetX`, `IOUtil`, `ODESolvers`, `MCNuX` in `ActiveThorns` — `ADMBaseX_InitialData`/`ADMBaseX_InitialGauge` schedule bins only exist when ODESolvers is runtime-active, and the transport group binds `AFTER ODESolvers_Solve`.
- `specs/tools/validate_specs.sh` lints spec text only — it never inspects C++/CCL, so a green run is **not** implementation progress. Conversely, specs pin exact needle strings (e.g. `η′ = αη`, `removed entirely (no continuous weight decay)`, `SHARES: ADMBaseX`) that implementation work must never edit away; run the validator after any spec touch.
- Single-owner shared utilities — write once, reuse everywhere; a per-consumer copy is a spec violation or test failure: conventions header (T02, `[MCNX-CNV-03]` sweep), Philox pure function (T03), coefficient-triple interface (T15, sole provenance per `[MCNX-INT-05]`), tetrad transform (T16, shared verbatim by PKT and INT), device-kernel tile-loop idiom (T08), axis-clamp wrapper (T26, wraps all five WeakLibInterp entry-point families), transport schedule group text (T05 — HYD pins only ordering *within* the group).
- Pinned this planning run by spec repair (all three amendments merged, validator exit 0): the OPA coefficient interface emits **physical units** (κ in cm⁻¹, η in MeV cm⁻³ s⁻¹ MeV⁻¹) and each consumer (INT, PKT) applies `×1.476625038e5` exactly once at its own point of use; foreign-thorn grid access uses `INHERITS: ADMBaseX HydroBaseX` plus fully qualified `READS: ADMBaseX::metric(everywhere)`-style clauses matching the `ghext` `CCTK_GroupIndex` names; `EXTENDS KEYWORD` gets a one-time parse/accept/bind probe with a pre-authorized fallback (own selector keywords + `"none"`/`"vacuum"` base values).

## Tier 0 — test harness and cross-cutting foundations

- [ ] T01 — Test-harness bootstrap + complete the `[MCNX-BLD-03]` precision gate (create `MCNuX/test/` with `test.ccl`, `unit-selftest.par`, golden data; parameter-gated selftest driver routine in thorn MCNuX that later tasks append legs to; AMReX-header `static_assert`s on `sizeof(amrex::Real)`/`sizeof(amrex::ParticleReal)` plus the runtime re-check)
  - spec: specs/build-and-integration.md
  - tests: (a) `agent_scripts/test.sh` discovers ≥1 test and it passes (exits the "0 tests ran" state `test.sh:27-38` special-cases); (b) golden parity `1e-12` on the selftest output file; (c) runtime assertions on `sizeof(amrex::Real)`/`sizeof(amrex::ParticleReal)` — exact; (d) compile-time gate: build fails if any of `double`/`CCTK_REAL`/`amrex::Real`/`amrex::ParticleReal` is not 8 bytes.
  - notes: First task because every other task's verification hangs off this harness — six research slices independently named "no `MCNuX/test/` exists" as their blocking gap. `[MCNX-BLD-05]` pins the parfile/golden location under `MCNuX/test/`. This is also the de-facto first end-to-end run of `agent_scripts/build.sh`/`test.sh` against a real `$CACTUSX` — the toolchain has never been exercised in-repo, so budget for toolchain breakage here. Precedents: `CarpetX/CarpetX/test/test.ccl` (`EXTENSIONS tsv`), `CarpetX/TestODESolvers/test/test.ccl` (per-test tolerance block). Selftest hosted in thorn MCNuX (assumption recorded in Decisions).

- [ ] T02 — Conventions header: pinned constants, nine derived factors, MeV↔Kelvin, species table (`[MCNX-CNV-01..06]`; single shared header, e.g. `MCNuX/src/Conventions.hxx`; factors computed from the five defining constants, not hard-coded; no runtime unit parameter)
  - spec: specs/conventions-and-units.md
  - tests: (a) recompute all nine geometrized↔cgs factors from the defining constants vs the pinned 10-significant-digit values — `1e-14`; (b) cgs↔code and MeV↔K round-trip identities — `1e-14`; (c) species table (index 0/1/2, `g = (1,1,4)`, lepton number `(+1,−1,0)`) — exact; cross-check species/axis naming against the committed WeakLibInterp `.h5ls` fixtures; (d) source-sweep gate (non-harness, exit non-zero on violation): no physical-constant numeric literal outside this header, no runtime unit-selection parameter.
  - notes: Consumed by PKT (g=4 once), OPA (every temperature argument ×1.160451812e10), INT/TRP (κ conversion `1.476625038e5`), HYD. One header only — `[MCNX-CNV-03]`'s no-duplication sweep makes per-consumer copies a test failure.

- [ ] T03 — Philox4x32-10 pure function `u(S,q,e,k)` (`[MCNX-RNG-01..05]`; header-only, host+device-callable, stateless: pinned round constants `M0=0xD2511F53 M1=0xCD9E8D57 W0=0x9E3779B9 W1=0xBB67AE85`, 10 rounds; exact `(S,q,e,k)→(c0..c3,k0,k1)` packing; `[0,1)` via `u64 = x_{2j} + 2^32·x_{2j+1}`, `r = (u64>>11)×2^-53`)
  - spec: specs/rng-and-statistical-acceptance.md
  - tests: (a) 3 known-answer vectors reproduced bitwise; (b) packing check for a pinned `(S,q,e,k)` tuple — exact; (c) double-construction closed forms `u64=0→0.0`, `2^11−1→0.0`, `2^11→2^-53`, `2^64−1→1−2^-53`, plus monotonicity — exact; (d) uniformity smoke: mean/variance/χ² at `N=1048576`, seed `1296518744` — 4σ; (e) gate: source sweep proving no physics path calls AMReX's stateful RNG (`[MCNX-RNG-05]`).
  - notes: API shape/location is implementation freedom — first implementer sets the corpus precedent; no exemplar exists (AMReX's `philox4x32x10` in `AMReX_Random.cpp`/`AMReX_RandomEngine.H` is the *forbidden* stateful device family, not a template). Blocks T17/T18/T19, T24, and every `stats-*` benchmark.

## Tier 1 — thorn wiring

- [ ] T04 — Subcycling paramcheck guard (`[MCNX-CTX-06]`: `MCNuX/param.ccl` `SHARES: CarpetX` + `USES BOOLEAN use_subcycling` (or `CCTK_ParameterGet`, explicitly free); `CCTK_PARAMCHECK` routine aborting when subcycling is on)
  - spec: specs/carpetx-thorn-integration.md
  - tests: non-harness gate: a run with `CarpetX::use_subcycling = yes` aborts at PARAMCHECK with non-zero exit and the pinned message — exact; the default (`no`) run proceeds.
  - notes: Idiom precedent `CarpetX/ODESolvers/param.ccl:23-27`; parameter declared RESTRICTED at `CarpetX/CarpetX/param.ccl:102-104`.

- [ ] T05 — Transport schedule group `AT evol AFTER ODESolvers_Solve` (`[MCNX-CTX-03..05]`: MCNuX-owned group at the pinned placement — `IN ODESolvers_PostStep` and `AT CCTK_POSTSTEP` are explicitly rejected; every container-walking routine `OPTIONS: level`/`global`; complete `READS:`/`WRITES:` per the newly pinned convention (`INHERITS: ADMBaseX HydroBaseX` + fully qualified clauses); `ghext` driver-native grid access)
  - spec: specs/carpetx-thorn-integration.md
  - tests: (a) `cadence-lag`-style benchmark asserting the group runs after `ODESolvers_Solve` and exactly once per step (ordering-counter grid variable or schedule-tree output) — golden parity `1e-12`; (b) inspection gate: no `interp.hxx`/AoS particle container, no host per-particle loop (`[MCNX-CTX-07]`, `[MCNX-GPU-07]`; anti-pattern source `CarpetX/CarpetX/src/interp.hxx:38`).
  - notes: The former open question on foreign-thorn `READS:`/`WRITES:` was resolved this planning run by spec repair — follow the amended `[MCNX-CTX-04]` (qualified `Implementation::group(region)` clauses; read-only on foreign groups except the default-off Tmunu contributor; fallback rungs: per-variable enumeration, then loud failure — never a dropped clause). Access path `ghext->patchdata[p].leveldata[l].groupdata[gi]->mfab[tl].array(mfi)` (`CarpetX/CarpetX/src/driver.hxx:596`); routine-options precedent `CarpetX/ODESolvers/schedule.ccl:39-43`; only in-stack include idiom `CarpetX/ODESolvers/src/solve.hxx:4-7`. The group text is written once, here — HYD (T06/T20) only orders routines within it.

- [ ] T06 — Source-term grid variables + zero-then-add protocol (`[MCNX-HYD-01]`, `[MCNX-HYD-03]`, `[MCNX-HYD-07]`: `PUBLIC:` groups `MCNuX::rad_force {rf_t rf_x rf_y rf_z}` and `MCNuX::lepton_source {lep_src}`, `TYPE=gf CENTERING={ccc}`, checkpointed; zeroing routine `WRITES: …(everywhere)` at the start of T05's group + contributor extension group ordered after it)
  - spec: specs/hydro-coupling-source-terms.md
  - tests: (a) zero-packet benchmark → all five variables identically `0` everywhere after every step — exact; (b) checkpoint/recover parity — `1e-12`; (c) schedule-ordering check that the zero routine precedes every contributor (a stub contributor's `+=` must survive) — exact; (d) gate: MCNuX writes no variable owned by another thorn (`[MCNX-HYD-07]`) — source sweep.
  - notes: Do **not** copy TmunuBaseX's `checkpoint="no"` tag (`CarpetX/TmunuBaseX/interface.ccl:11-15`); extension-group model `CarpetX/TmunuBaseX/schedule.ccl:25-33`. Prerequisite for any `WRITES:` clause in T05 and for the ledger (T20).

## Tier 2 — container and analytic backgrounds

- [ ] T07 — Packet container type + id contract (`[MCNX-GPU-01]`, `[MCNX-GPU-03]`, `[MCNX-GPU-04]`: thin derived `amrex::ParticleContainerPureSoA<...>` realizing the component schema — `x^i`×3, `p_i`×3, weight `N`, species `s`, event counter `e`; packet id `q` = AMReX `idcpu`; keyed to the driver's `AmrCore`/`BoxArray`/`DistributionMapping`)
  - spec: specs/particle-container-and-gpu.md
  - tests: (a) compile-time `static_assert`s: type is-a `ParticleContainerPureSoA`, tile view trivially copyable (`amrex/Src/Particle/AMReX_ParticleTile.H:33,467`) — exact; (b) selftest creating N packets: component count/order matches the schema, values round-trip — exact; (c) id contract: ids unique, immutable, never reused, bit 63 clear — exact (single rank now; multi-rank leg lands with T21).
  - notes: Exemplar `warpx/Source/Particles/WarpXParticleContainer.H:193-199`. Per-patch vs single-container organization and concrete template arguments are the implementer's recorded decision (spec:303-307); `idcpu` bit-63 clearance assumes < 8388608 ranks.

- [ ] T08 — Device-kernel execution idiom (`[MCNX-GPU-02]`, `[MCNX-GPU-07]`: the reusable tile-loop + `amrex::ParallelFor` pattern with by-value trivially-copyable captures that all later physics kernels — gather, push, sample, deposit — must follow)
  - spec: specs/particle-container-and-gpu.md
  - tests: (a) a trivial per-packet kernel driven through the helper produces exact expected component values, verified in both CPU/MPI and GPU builds; (b) inspection gate: zero host per-particle loops, zero use of `interp.hxx`'s `amrex::AmrParticleContainer<3,2>`.
  - notes: Structural exemplar `warpx/Source/Particles/PhysicalParticleContainer.cpp:1387-1568`. Consolidation point — GEO push, PKT emission, INT events, HYD deposition all reuse this one idiom.

- [ ] T09 — `TestMCNuX` keyword-extension idiom + Schwarzschild writer (`[MCNX-VER-02]`, `[MCNX-VER-03]`: `SHARES: ADMBaseX` + `EXTENDS KEYWORD initial_data`/`initial_lapse` adding `"Schwarzschild"` + private parameters (`M`, …); `CCTK_EQUALS`-gated writers `IN ADMBaseX_InitialData` (metric/curvature) and `IN ADMBaseX_InitialGauge` (lapse) with `WRITES: …(everywhere)`; isotropic closed forms in `TestMCNuX/src/*.cxx`; writers only — no transport physics in TestMCNuX)
  - spec: specs/verification-suite-design.md
  - tests: (a) written `ψ`, `γ_ij`, `K_ij = 0`, `α`, `β^i = 0` vs the closed forms at pinned `M` — `1e-14`; (b) the extended keyword actually binds (parfile sets `ADMBaseX::initial_data = "Schwarzschild"` and the writer runs); (c) golden parity `1e-12`.
  - notes: `EXTENDS KEYWORD` has zero working precedent in the checked-out CarpetX/WeakLibInterp trees — run the one-time parse/accept/bind probe now pinned in the amended `[MCNX-VER-02]`, and on failure use its pre-authorized fallback (own selector keywords, `ADMBaseX::initial_data/initial_lapse = "none"`). Target variables `CarpetX/ADMBaseX/interface.ccl:12,16,20`.

- [ ] T10 — `TestMCNuX` hydro-profile writers (`[MCNX-VER-04]`: `SHARES: HydroBaseX` + `EXTENDS KEYWORD initial_hydro` adding `"uniform sphere"` and `"equilibration box"`; writers `IN HydroBaseX_InitialData` populating `rho`, `vel`, `temperature` (MeV), `Ye`, explicitly zeroing `eps`, `press`, `entropy`, `Bvec`, `Avecx/y/z`)
  - spec: specs/verification-suite-design.md
  - tests: (a) `ρ, T, Yₑ, vel` vs closed forms — `1e-14`; (b) zeroed variables exactly `0`; (c) golden parity `1e-12`.
  - notes: `HydroBaseX_InitialData` has an `AT initial` fallback when ODESolvers is inactive, unlike ADMBaseX (`CarpetX/HydroBaseX/schedule.ccl:15-27`). Fallback for the keyword (per amended spec): `initial_hydro` has no `"none"` value — retain `"vacuum"` with the profile writer ordered `AFTER HydroBaseX_initial_data`. The equilibration box is the fixture for `stats-equilibration-*` and the TRP benchmarks.

## Tier 3 — transport physics

- [ ] T11 — Metric gather + spatial derivatives (`[MCNX-GEO-03]`: ≥trilinear interpolation of `α, β^i, γ_ij` from ADMBaseX vertex-centered grid functions at arbitrary packet positions, derivative evaluation exact on linear fields, same-gathered-snapshot invariant for every term of one RHS evaluation)
  - spec: specs/geodesic-propagation.md
  - tests: (synthetic grid in the selftest — no TestMCNuX dependency) (a) interpolation and derivatives exact on linear fields to `1e-14`; (b) snapshot consistency: all terms of one RHS evaluation see identical gathered values — exact; (c) off-vertex and cell-corner positions covered.
  - notes: Centering provenance: `CarpetX/ADMBaseX/interface.ccl:12,16,20` (no CENTERING tag) + `CarpetX/CarpetX/src/driver.cxx:430` ("Default: vertex-centred").

- [ ] T12 — Null-geodesic RHS + `p^t` null closure (`[MCNX-GEO-01]`, `[MCNX-GEO-02]`: device-callable free functions for `dx^i/dt = γ^{ij}p_j/p^t − β^i` and `dp_i/dt = −α(∂_iα)p^t + (∂_iβ^k)p_k − ½(∂_iγ^{jk})p_jp_k/p^t`, with `p^t = √(γ^{ij}p_ip_j)/α` recomputed at every RHS evaluation and never stored as state)
  - spec: specs/geodesic-propagation.md
  - tests: (a) null identity `g_{μν}p^μp^ν = 0` after closure — `1e-14`; (b) RHS matches hand-computed analytic values at a pinned Minkowski point and a pinned Schwarzschild point — `1e-14`; (c) schema gate: `p^t` absent from the packet state/SoA schema — exact (inspection).

- [ ] T13 — ≥2nd-order integrator + push kernel (`[MCNX-GEO-04]`: RK2 (minimum compliant choice) over `Δt_push`, wired into T05's group through T08's kernel idiom over T07's container)
  - spec: specs/geodesic-propagation.md
  - tests: (a) flat-spacetime straight-line propagation exact to `1e-14` over pinned steps; (b) measured convergence order ≥2 (error drops ≥4× under `Δt/2`) on a Schwarzschild background; (c) `minkowski-freestream` benchmark golden parity `1e-12`.

- [ ] T14 — `schwarzschild-pt` conservation benchmark (`[MCNX-GEO-05]`: benchmark parfile + golden data exercising T09 + T13)
  - spec: specs/geodesic-propagation.md
  - tests: `p_t = −α²p^t + β^i p_i` conserved on the stationary background to the per-benchmark tolerance pinned in specs/verification-suite-design.md; run output golden parity `1e-12`.

- [ ] T15 — Coefficient interface + analytic gray-opacity mode (`[MCNX-VER-05]`: the `(η, κ_a, κ_s)` coefficient-triple interface plus a parameter-selected analytic source (`kappa_a0`, `kappa_s0`, `eta_scale` per species) inside thorn MCNuX)
  - spec: specs/verification-suite-design.md
  - tests: (a) returned triple equals the analytic parameter values per species/energy — `1e-14`; (b) the source selector switches cleanly and defaults are pinned; (c) interface shape is the single consumption point for INT/TRP (`[MCNX-INT-05]` provenance rule).
  - notes: Deliberately sequenced before the WeakLibInterp-backed implementation (T25): hard prerequisite for `emission-fixedseed`, `interactions-fixedseed`, `stats-emission`, `stats-beam`, `stats-equilibration-*`, `trp-alpha-limit-*`. Per the amended OPA spec, the interface emits physical units (κ in cm⁻¹); consumers convert.

- [ ] T16 — Tetrad transform, Gram–Schmidt in fixed (x,y,z) order (`[MCNX-PKT-04]` core: shared device-callable transform from `u^μ` and `g_μν` — timelike leg `u^μ`, spatial legs Gram–Schmidt in the pinned order)
  - spec: specs/packet-representation-and-sampling.md
  - tests: (a) static-flat/identity case — exact; (b) null preservation and `ν = −p_μ u^μ` recovery — `1e-14`; (c) isotropy of transformed directions — 4σ.
  - notes: Consolidation point — `[MCNX-INT-04]` scattering redraw must reuse this exact function, not reimplement it (`specs/neutrino-matter-interactions.md:136-154`).

- [ ] T17 — Emission count law + packet creation (`[MCNX-PKT-01..03]`, `[MCNX-PKT-05]`, `[MCNX-PKT-06]`: `N_p = g_s √−g ΔV Δt η_b(s)/E_p` with floor + Bernoulli remainder, one draw per `(cell,step,s,b)` with the pinned `K_cell` packing and `(q=K_cell, e=n, k=3b+s)` assignment; per-packet init — bin-center `ν`, weight `N = E_p/ν`, uniform position, creation draws in the fixed `k=0..5` order; `E_p` as a `param.ccl` parameter)
  - spec: specs/packet-representation-and-sampling.md
  - tests: (a) draw-map audit — the `(q,e,k)` tuple trace matches the pinned order exactly; (b) expected-count statistics and energy/species spectra within 4σ at `N_p = 2^20`, seed `1296518744`; (c) `g = 4` applied exactly once — νx count ratio 4× within 4σ, plus a source sweep that the factor appears at exactly one site; (d) weight `N = E_p/ν` — `1e-14`; (e) `emission-fixedseed` golden parity `1e-12` and single-rank replay determinism — exact.
  - notes: η consumed from T15's interface in physical units; PKT applies its own single conversion per the amended OPA spec.

- [ ] T18 — Interaction draw law + event competition, pure functions (`[MCNX-INT-01]`, `[MCNX-INT-02]` arithmetic: `Δt_{s,a} = −ln(1 − u_{s,a}) p^t/(κ_{s,a} ν)` with κ converted to code units (`κ[code] = κ[cm⁻¹] × 1.476625038e5`), competition over `min(Δt_s, Δt_a)`, cell-exit time and `Δt_rem`)
  - spec: specs/neutrino-matter-interactions.md
  - tests: (a) Δt matches the closed form for pinned `(κ, ν, p^t, u)` — `1e-14`; (b) zero-opacity channel yields `+∞` with no NaN — exact; (c) tie resolves to scattering — exact; (d) both draws always consumed even when one channel is inactive — draw-count check, exact.
  - notes: Dependency-light — unit-testable against hand-fed coefficients before the event loop exists. INT-side conversion is the single application point per the amended OPA spec.

- [ ] T19 — Interaction event loop: episodes, absorption, scattering, RNG map (`[MCNX-INT-02..06]`: episode (re)start at creation/step start/post-scatter/cell crossing with frozen episode kinematics; discrete absorption — packet removed entirely, no continuous weight decay, content handed to the HYD ledger, id retired; elastic scattering — isotropic fluid-frame redraw at fixed `ν` via T16; coefficients consumed only from T15's interface (primed values when TRP relabeling is active); event counter `e` incremented per episode with the fixed `k=0..3` draw order via T03)
  - spec: specs/neutrino-matter-interactions.md
  - tests: (a) beam attenuation vs `exp(−κ_a s)` — 4σ; (b) collision statistics — 4σ; (c) scattering isotropy — 4σ; (d) elasticity: `ν` unchanged and null identity preserved — `1e-14`; (e) discreteness: absorbed packets vanish, no weight decay anywhere — exact; (f) draw-map audit — exact; (g) `interactions-fixedseed` golden parity `1e-12`.

- [ ] T20 — Exchange ledger deposition + conservation closure (`[MCNX-HYD-02]`, `[MCNX-HYD-05]`: contributor kernel(s) accumulating per-cell `ΔP^t, ΔP_i, ΔL` from creation/absorption/scattering into `G^t = −ΔP^t/(ΔVΔt)` etc., undensitized/coordinate-volume normalized, `(G^t up, G_i down)`, `+=` into T06's variables)
  - spec: specs/hydro-coupling-source-terms.md
  - tests: (a) ledger closure: per step and per channel (energy, momentum, lepton number), grid total equals packet-side audit total to relative `1e-13`; (b) sign fixture: absorption raises fluid energy / lowers neutrino energy — exact sign check; (c) `interactions-fixedseed` and `cadence-lag` golden parity `1e-12`.

- [ ] T21 — Redistribution and after-step ownership invariant (`[MCNX-GPU-05]`, `[MCNX-GPU-06]`: `Redistribute`/`RedistributeLocal(N)` after the step's last position-changing operation and after every regrid, with `N` derived from a proven per-step motion bound, not assumed; scheduled inside T05's group)
  - spec: specs/particle-container-and-gpu.md
  - tests: (a) `ownership-multibox` benchmark — after every step each packet resides on the rank/level/grid/tile owning its position, multi-box and multi-rank — exact; (b) bounded-motion audit: measured max per-step displacement `< N` cells — assertion; (c) regrid hook exercised; (d) id contract holds across ranks (bit 63 clear, no reuse) — exact.

- [ ] T22 — Opacity relabeling + α control + α→1 limit (`[MCNX-TRP-02..05]` with `[MCNX-TRP-01]` baseline selector: pure post-map `(η, κ_a, κ_s, α) → (η′ = αη, κ_a′ = ακ_a, κ_s′ = κ_s + (1−α)κ_a)` in the binding 2021 convention (`arXiv:2103.16588v2` Eq. 43 — not the flipped 2020-letter α); per-cell control `α = min(1, ξ/(κ_a Δt_c))` with `ξ` and a fixed-α override as new `param.ccl` parameters; scheme selector giving the explicit baseline when deselected)
  - spec: specs/trapped-regime-treatment.md
  - tests: (a) the two exact relabeling invariants — exact/`1e-14`; (b) convention tripwire: a pinned numeric case that fails under the 2020 flipped convention — exact; (c) control-enforcement diagnostic `κ_a′Δt_c ≤ ξ` — exact; (d) `trp-alpha-limit-explicit` vs `trp-alpha-limit-relabeled` at α=1 produce bitwise-identical event sequences and outputs — exact.
  - notes: Scope deferrals recorded in the spec: the β-dependent Eq. 50 stability bound is deliberately not enforced by the baseline control; `α ∈ (0,1]` paramcheck rejection is a corpus assumption, not a source requirement.

- [ ] T23 — Diffusion displacement sampler + B(f_free)/θ₂ momentum model (`[MCNX-TRP-08]`, `[MCNX-TRP-10]`: two-branch displacement draw with the closed-form `F` and numerical inversion, the no-scatter atom branch; the 21-value `B(f_free)` table with linear interpolation and the closed-form `cosθ₂` draw)
  - spec: specs/trapped-regime-treatment.md
  - tests: (a) inversion accuracy ≤ `1e-12` in `r̃`; (b) sampled displacement moments vs the closed-form law — 4σ (`stats-diffusion`); (c) `B` interpolation vs the 21 pinned values — `1e-14`; (d) atom branch pins `cosθ₂ = 1` at `f_free = 1` — exact.
  - notes: Both pieces are self-contained math, unit-testable before the event loop is wired.

- [ ] T24 — Diffusion regime entry, two-leg position update, five-draw RNG map (`[MCNX-TRP-06]`, `[MCNX-TRP-07]`, `[MCNX-TRP-09]`: regime criterion evaluated after every performed scattering (with active `κ_s′`), advection leg + free-streaming leg position update via T13's pusher and T16's tetrad, fixed five-draw `k=0..4` map over T03)
  - spec: specs/trapped-regime-treatment.md
  - tests: (a) criterion fires exactly at the pinned condition — exact; (b) draw-map audit — exact; (c) `stats-equilibration-*` on T10's equilibration box — 4σ; (d) `trp-overlap-consistency` benchmark — golden parity `1e-12`.

- [ ] T25 — WeakLibInterp-backed opacity/EOS layer (`[MCNX-OPA-01]`, `[MCNX-OPA-02]`, `[MCNX-OPA-04]`, `[MCNX-OPA-05]`: call sites for `wli::EosInterpolateSingleVariable3DPoint`, `EmAbInterpolateSingleVariable4DPoint` (caller pre-logs E/ρ/T, `Ye` raw), `IsoOffset` + `IsoInterpolateSingleVariable5DPoint` (iMom=0); table residency via `wli::ResidentTable<D>::upload`/`TableView<D>` capture owned by a startup/parameter-recovery routine; Kirchhoff closure with `μ_ν` from the `Electron/Proton/Neutron Chemical Potential` datasets; νx mapping `κ_a = η = 0`, `κ_s = ½`-mean of the electron-type values; MCNuX must not reimplement interpolation; temperature converted MeV→Kelvin at every entry point via T02)
  - spec: specs/opacity-eos-evaluation.md
  - tests: (a) kernel parity: MCNuX wrapper vs direct WeakLibInterp call at pinned points — exact/`1e-14`; (b) unit-boundary check: temperature-argument call-site sweep proving `×1.160451812e10` is applied everywhere, plus a numeric check that a MeV input reaches the library as Kelvin — `1e-14`; (c) Kirchhoff closure vs closed form at pinned `(ρ,T,Yₑ,ν)` — `1e-12`; (d) νx: `κ_a = η = 0` exact, `κ_s(νx) = ½(κ_s(νe)+κ_s(ν̄e))` — `1e-14`; (e) fed through T15's interface, `interactions-fixedseed` still golden-parity clean.
  - notes: Per the amended spec, this layer emits physical units and applies no geometrized factor — INT/PKT convert at their points of use. Key WeakLibInterp loci: `src/eos/wli_eos.H:54-80`, `src/opacity/wli_opacity_emab_iso.H:75-117,228-232,248-290`, `src/core/wli_table.H:60-115`, fixture `specs/fixtures/wl-EOS-SFHo-15-25-50.h5ls`.

- [ ] T26 — Table-range enforcement + EOS-inversion error protocol (`[MCNX-OPA-06]`, `[MCNX-OPA-07]`: MCNuX-side clamping of every query coordinate into axis range before any WeakLibInterp call, transparency floor at `ρ < ρ_min`, clamp-diagnostic counters; wrapper around `ComputeTemperatureWith_*` consuming `EosInversionResult{T,Error}` with hard failure on `Error≠0` in test mode)
  - spec: specs/opacity-eos-evaluation.md
  - tests: (a) out-of-range query returns the clamped-coordinate reference value, never an extrapolation — exact; (b) `ρ < ρ_min` → `κ = η = 0` — exact; (c) clamp counters increment once per clamped axis — exact; (d) non-positive log input never yields NaN — NaN-equality check, exact; (e) `Error≠0` aborts in test mode — error-code check, exact.
  - notes: WeakLibInterp clamps only the bracket index and leaves the delta unclamped — it silently extrapolates and silently produces NaN on non-positive log input (`WeakLibInterp/src/core/wli_interp.H:65-87`); enforcement is MCNuX's responsibility per the corpus boundary policy. The clamp utility must wrap all five entry-point families, not just EmAb/Iso.

- [ ] T27 — Optional `Tmunu` contributor, default off (`[MCNX-HYD-06]`: parameter-gated routine `IN TmunuBaseX_AddToTmunu` computing per-cell `q_{μν}` and the 8-point cell-to-vertex average into `eTtt/eTti/eTij`, plus the conditional `configuration.ccl` capability on TmunuBaseX)
  - spec: specs/hydro-coupling-source-terms.md
  - tests: (a) closed-form `q_{μν}` for a pinned packet population — `1e-12`; (b) 8-point average exact on a linear field — `1e-14`; (c) default-off gate: with the parameter off, no TmunuBaseX dependency is required and `eTtt/eTti/eTij` are untouched — exact.
  - notes: Keep default off — the double-counting protocol against a partner applying four-force sources is an unresolved corpus question the spec itself defers (see Decisions). Precedents `CarpetX/TmunuBaseX/schedule.ccl:25-33`, `interface.ccl:11-15`.

## Tier 4 — suite completion and delivery gates

- [ ] T28 — `stats-*` benchmark set + statistical reduction harness (`[MCNX-VER-06]`, `[MCNX-VER-07]`, `[MCNX-RNG-07]`, `[MCNX-RNG-08]`: `stats-emission`, `stats-beam`, `stats-diffusion`, `stats-equilibration-*` under `MCNuX/test/`, each reducing a stochastic run to a small set of golden numeric estimates)
  - spec: specs/verification-suite-design.md
  - tests: each estimator within 4σ at `N_p = 2^20`, primary seed `1296518744`, plus a rerun at secondary seed `20260730`; the reduced numbers compared to committed golden data at `1e-12`.
  - notes: Every parfile must satisfy the amended `[MCNX-VER-06]` activation floor and include a non-vacuity golden quantity so an unbound schedule fails the `1e-12` diff instead of passing silently.

- [ ] T29 — GPU build mode + delivery gates (`[MCNX-BLD-03]` mode leg, `[MCNX-BLD-06]`, `[MCNX-CTX-02]`: an OPTIONAL GPU/CUDA capability path and a GPU build invocation alongside `agent_scripts/build.sh`'s CPU config; the non-harness gates)
  - spec: specs/build-and-integration.md
  - tests: (a) clean build in both CPU/MPI and a GPU backend mode against the pinned ThornList, with `unit-selftest` passing in both; (b) activation-split audit: no non-test parfile lists `TestMCNuX` in `ActiveThorns` — exact sweep (vacuous today, zero `.par` files exist); (c) configure without CarpetX fails — non-zero exit gate.

## Decisions & non-blocking spec items

- Resolved this planning run by spec amendment (all merged, `validate_specs.sh` exit 0): (1) foreign-thorn `READS:`/`WRITES:` convention pinned in `specs/carpetx-thorn-integration.md` `[MCNX-CTX-04]`; (2) OPA coefficient units pinned consumer-side in `specs/opacity-eos-evaluation.md`; (3) benchmark-parfile activation floor and `EXTENDS KEYWORD` probe/fallback pinned in `specs/verification-suite-design.md`. The `standing.par` citation flagged stale by research was re-checked and resolves — no action.
- Which thorn hosts the runtime `unit-selftest` battery: `[MCNX-BLD-03]` says "the corpus"; this plan assumes thorn MCNuX with the parfile under `MCNuX/test/` per `[MCNX-BLD-05]`. Candidate `specs/build-and-integration.md` amendment if a build iteration hits ambiguity; likewise whether non-harness CI gates are in delivery scope.
- Tmunu double-counting: no protocol exists for enabling `[MCNX-HYD-06]` while a GRMHD partner also applies four-force sources; the spec defers this to a future corpus change — T27 stays parameter-gated default-off, and enabling it in production with an evolving spacetime is out of scope until the HYD spec is amended.
- TRP deferrals (informational, spec-recorded): β-dependent Eq. 50 stability bound not enforced by the baseline α control; `α ∈ (0,1]` paramcheck rejection is a corpus assumption.
- T07 implementer decisions to record in-code: per-patch vs single container organization, concrete `ParticleContainerPureSoA` template arguments; id headroom assumes < 8388608 ranks.
- `MCNuX/src/stub.cxx` carries no license header despite the GPLv3 `LICENSE`; no spec mandates one — delivery-convention call, not a build increment.
- `specs/carpetx-thorn-integration.md`'s `~line` citations are self-admittedly approximate; only `driver.hxx:596` was spot-verified.

## Discovered since last plan
