# TODO — MCNuX build plan

Greenfield. The spec corpus (`specs/README.md` + 12 coded specs) is complete, self-consistent, and coverage-closed — no specs need authoring. Executable code is a stub: `MCNuX/src/stub.cxx` (two precision static_asserts) plus empty `.ccl` shells; `MCNuX/configuration.ccl` already satisfies [MCNX-BLD-02]/[MCNX-CTX-02] (`REQUIRES CarpetX`/`REQUIRES WeakLibInterp`). `MCNuX/test/` does not exist yet. `agent_scripts/{build,test}.sh` + sandbox tooling are complete and tolerate the current 0-test state. Everything below is implementation + test artifacts against existing specs, in priority order.

## Standing facts every item inherits

- Acceptance is entirely programmatic: exact/bitwise (RNG KATs, enumerations, α=1 parity), machine ~1e-14, golden 1e-12 (Cactus harness ABSTOL=RELTOL=1e-12), relaxed 1e-10, 4σ statistical at `N_p=1048576` seed `1296518744` (secondary `20260730`), conservation 1e-13/step.
- Shared utilities are implemented ONCE and reused: constants/units/species header (T1); Philox `u(S,q,e,k)` (T2); tetrad transform in fixed (x,y,z) Gram–Schmidt order (T13, reused by scattering [MCNX-INT-04] and diffusion [MCNX-TRP-09/10]); 4σ z-score/estimator machinery (T26, built with first stats benchmark); 3×3 inverse-metric + null-closure helper (T13).
- `MCNuX/test/` + `test.ccl` scaffolding is created ONCE by T5, not per-task. New `.cxx` files must be added to `MCNuX/src/make.code.defn` (currently `SRCS = stub.cxx`). The first shared header fixes the repo's header-location convention; apply uniformly afterward.
- Binding anti-patterns: never use AMReX's stateful RNG (`amrex/Src/Base/AMReX_Random.H`, [MCNX-RNG-05]); never adopt CarpetX's AoS/host-loop interpolator (`CarpetX/CarpetX/src/interp.hxx:38`, [MCNX-GPU-07]/[MCNX-CTX-07]). Grid access via the `ghext` relative-include idiom (precedent `CarpetX/ODESolvers/src/solve.hxx:4-5`, [MCNX-CTX-01]).
- The Cactus regression harness (flesh) resolves only inside the agent sandbox (`agent_scripts/sandbox/setup.sh`); harness-dependent tests (T5 onward) run there, not in the analysis environment.
- Consistency sweep (`specs/conventions-and-units.md:132`) forbids duplicated constant literals — a standing review gate once T1 lands.

## Tier 1 — Foundations (pure functions, guards, first test home)

- [x] T1: Constants/units/species header — constexpr five defining constants, nine derived geometrized conversion factors, MeV↔Kelvin `1.160451812e10`, species enum `{νe=0,ν̄e=1,νx=2}` with degeneracies `(1,1,4)` and lepton numbers `(+1,-1,0)`; no runtime unit parameters ([MCNX-CNV-02..06])
  - spec: specs/conventions-and-units.md
  - tests: compile-time/selftest recomputation of each pinned 10-digit factor from defining constants to machine tier (~1e-14 rel); round-trip conversion identities; species-table exact equality; species axis later cross-checked vs `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls` (folds into T9). Fold all into T5's `unit-selftest` when it lands.
  - notes: DONE — `MCNuX/src/mcnux_units.hxx` (header-only, `namespace MCNuX`, no cctk/AMReX includes), 47 static_asserts in three groups, included from `stub.cxx` so every build re-verifies. Header convention now fixed: flat `MCNuX/src/mcnux_*.hxx`, `namespace MCNuX`. Tolerance tiers live in `MCNuX::detail`: `rtol_pinned = 1e-9` (vs truncated 10-digit spec literals — `1e-14` vs pinned values is arithmetically impossible, deviations up to 2.1e-10; [MCNX-CNV-04] allows extra digits), `rtol_machine = 1e-14` (full-precision identities/round-trips). T5 must reuse these, not re-derive. Consistency sweep clean: each pinned literal appears exactly once (in its static_assert anchor). Naming is directional (`*_code_to_cgs`; opacity reverse direction documented in-header).
- [ ] T2: Philox4x32-10 RNG primitive — pure device-callable `u(S,q,e,k)`: round function + constants ([MCNX-RNG-01]), counter/key packing ([MCNX-RNG-03]), `[0,1)` double via `(u64 >> 11) * 2^-53` ([MCNX-RNG-04]); header wired into make.code.defn
  - spec: specs/rng-and-statistical-acceptance.md
  - tests: three KAT vectors bit-for-bit exact ([MCNX-RNG-02], spec `:61-69`); purity/independence evidence — repeated calls, permuted batch order, host-vs-device bitwise agreement ([MCNX-RNG-06]); fold into `unit-selftest` (T5). [MCNX-RNG-05] no-stateful-RNG is an inspection gate on all later physics tasks.
  - notes: hard prerequisite of T14–T16/T18/T22; no exemplar exists in the reference stack — implement from spec text alone (`:51-92`).
- [ ] T3: Paramcheck guard on subcycling — `MCNuX/param.ccl`: `SHARES: CarpetX` + `USES BOOLEAN use_subcycling` (idiom `CarpetX/ODESolvers/param.ccl:23-27`); `CCTK_PARAMCHECK` routine aborting with explanatory error when `yes` ([MCNX-CTX-06])
  - spec: specs/carpetx-thorn-integration.md
  - tests: exact — run with `use_subcycling=yes` aborts at paramcheck with the message; `=no` proceeds.
  - notes: small, independently testable, unblocks schedule work; spec anchor `:161-166`.
- [ ] T4: Transport cadence schedule group — MCNuX-owned group `AT evol AFTER ODESolvers_Solve` in `MCNuX/schedule.ccl` with an initial observable routine (step counter/diagnostic), `OPTIONS: level` (or `global`), complete `READS:`/`WRITES:`; driver-include wiring to reach `ghext` ([MCNX-CTX-01/03/04/05])
  - spec: specs/carpetx-thorn-integration.md
  - tests: builds and runs under the pinned ThornList (`agent_scripts/sandbox/mcnux.th`); cadence observable fires exactly once per coarsest-level Δt after `ODESolvers_Solve` (full `cadence-lag` golden benchmark deferred to T20). Missing READS/WRITES → CarpetX runtime error is the enforcement backpressure.
  - notes: anchors `specs/carpetx-thorn-integration.md:133-159`, `CarpetX/CarpetX/src/driver.hxx:77,596`.
- [ ] T5: First benchmark scaffolding — create `MCNuX/test/test.ccl` (`NPROCS 1`, precedent `CarpetX/TestODESolvers/test/test.ccl:1-16`), `run_selftests` parameter, selftest battery routine writing pass/fail TSV, `unit-selftest.par` + golden dir ([MCNX-VER-01], [MCNX-BLD-05])
  - spec: specs/verification-suite-design.md
  - tests: `agent_scripts/test.sh` discovers ≥1 test and passes at harness defaults ABSTOL=RELTOL=1e-12; battery includes T1 factor/round-trip checks and T2 KAT/closed-form checks as exact rows.
  - notes: consumes T1–T4; closes the "0 tests discovered" state and gives every later pure-function task a home; also `specs/build-and-integration.md:155-160`. Harness runnable only in sandbox.
- [ ] T25a: Run `bash specs/tools/validate_specs.sh` in the sandbox once to confirm the corpus passes
  - spec: specs/build-and-integration.md
  - tests: validator exit 0 is the test (exact).
  - notes: cheap, do early; never mechanically executed by research agents.

## Tier 2 — Independent fixtures & particle container

- [ ] T6: TestMCNuX Schwarzschild initial data — `SHARES: ADMBaseX` + `EXTENDS KEYWORD initial_data`/`initial_lapse` adding `"Schwarzschild"`; writers `IN ADMBaseX_InitialData`/`ADMBaseX_InitialGauge` with full `WRITES: (everywhere)`; isotropic `ψ=1+M/(2r)`, `γ_ij=ψ⁴δ_ij`, `K_ij=0`, `α=(1−M/(2r))/(1+M/(2r))`, `β^i=0` ([MCNX-VER-02/03])
  - spec: specs/verification-suite-design.md
  - tests: smoke parfile + golden TSV comparing grid output against the closed forms (machine/golden tier).
  - notes: fully unblocked — depends only on stable base thorns; anchors `specs/verification-suite-design.md:139-181`, `CarpetX/ADMBaseX/schedule.ccl:12,16,68-104`, `CarpetX/ADMBaseX/interface.ccl:12-20`; shift needs no writer (default `"zero"`).
- [ ] T7: TestMCNuX hydro profiles — `SHARES: HydroBaseX` + `EXTENDS KEYWORD initial_hydro` adding `"uniform sphere"`, `"equilibration box"`; writers `IN HydroBaseX_InitialData` for `(ρ,T,Yₑ)` profiles, `vel=0`, other fields zeroed, temperature in MeV; parameters for `M`, sphere center/radius, `(ρ₀,T₀,Yₑ₀)`, ambient values ([MCNX-VER-04])
  - spec: specs/verification-suite-design.md
  - tests: smoke parfile + golden comparing profile fields to formulas (golden tier).
  - notes: independent build+test cycle from T6 (separate keyword family/writer); anchors `specs/verification-suite-design.md:182-201`, `CarpetX/HydroBaseX/interface.ccl:11-29`, `CarpetX/HydroBaseX/schedule.ccl:17,59-75`.
- [ ] T8: Pure-SoA packet container + AMReX precision gate — container deriving `amrex::ParticleContainerPureSoA` with the six-component schema (position ×3 built-in, momentum ×3 real, weight real, species int, event-counter int, id via `idcpu`); one exemplar `amrex::ParallelFor` tile kernel establishing the device-kernel pattern ([MCNX-GPU-01/02/03]); complete [MCNX-BLD-03] with `static_assert(sizeof(amrex::Real)==8)` / `sizeof(amrex::ParticleReal)==8`
  - spec: specs/particle-container-and-gpu.md
  - tests: compile-time layout pin — static_asserts that the type derives `ParticleContainerPureSoA` and tile-view is `std::is_trivially_copyable` (spec `:236-240`); CPU/MPI build green via `agent_scripts/build.sh` (GPU-backend build deferred to T25d).
  - notes: pattern `warpx/Source/Particles/WarpXParticleContainer.H:193`; schema at spec `:127-154`, kernel model `:156-171`; construction needs `AmrCore` via `ghext` (T4); decide the precision gate's permanent home — `stub.cxx` likely retires (currently only double/CCTK_REAL halves exist at `MCNuX/src/stub.cxx:7-8`).

## Tier 3 — Opacity/EOS layer

- [ ] T9: WeakLibInterp call-boundary wrapper — thin wrappers over the five `_Point` families with pinned argument/units/log10/index contract; MeV→Kelvin (via T1 factor) at EVERY temperature argument including NES/Pair; table-range awareness hooks ([MCNX-OPA-01/02/03])
  - spec: specs/opacity-eos-evaluation.md
  - tests: kernel-parity selftests vs direct WeakLibInterp calls (exact/machine tier); unit-boundary check — Kelvin conversion applied exactly once per argument; species-axis grounding vs committed `.h5ls` fixture shapes ([MCNX-OPA-08]). Fold into `unit-selftest`.
  - notes: WeakLibInterp is consume-only, never modified; entry points verified at `WeakLibInterp/src/eos/wli_eos.H:54`, `wli_eos_inversion.H:265-419`, `wli_opacity_emab_iso.H:76,229,249`, `wli_opacity_nes_pair.H:66,193,235`, `wli_opacity_brem.H:102-130`; spec anchor `:174-220`.
- [ ] T10: Baseline coefficient assembly + νx mapping — `κ_a(s,E)`, `κ_s(s,E)`, `η(s,E)` for νe/ν̄e from EmAb+Iso with Kirchhoff closure using EOS chemical potentials at (ρ,T,Yₑ); νx: `κ_a=0`, `η=0`, `κ_s` = ½-mean of the two electron-type Iso values ([MCNX-OPA-04/05])
  - spec: specs/opacity-eos-evaluation.md
  - tests: Kirchhoff-closure identity checks; νx zero/half-mean identities (machine tier); selftest rows.
  - notes: depends on T9; NES/Pair/Brem thermal-emission assembly is explicitly out of scope — do not attempt (spec `:340-349`).
- [ ] T11: Table-range enforcement — transparency floor (ρ<ρ_min ⇒ all coefficients zero, no table touch), component-wise clamping elsewhere, per-axis clamp diagnostic counters, no-NaN guarantee ([MCNX-OPA-06]); `EosInversionResult{T,Error}` hard-fail-in-test-mode protocol ([MCNX-OPA-07]) may be stubbed since baseline transport receives T in MeV externally
  - spec: specs/opacity-eos-evaluation.md
  - tests: range-policy selftests — out-of-range probes return clamped/zero values with correct counter increments, never NaN (exact).
  - notes: WeakLibInterp extrapolates permissively, so enforcement is entirely MCNuX's (README boundary policy); spec anchor `:174-272`.
- [ ] T12: Analytic/gray opacity mode + source-agnostic dispatch — runtime-parameter-selected coefficient interface (table vs analytic); analytic-mode formulas pinned at `specs/verification-suite-design.md:202-233` ([MCNX-VER-05])
  - spec: specs/verification-suite-design.md
  - tests: analytic-mode coefficient values vs closed forms (machine tier); dispatch selection exactness.
  - notes: enabling fixture for all fixed-seed/statistical transport benchmarks (T14–T16, T18) without table dependence; only needs the interface shape, so it may precede T10/T11 under sequencing pressure.

## Tier 4 — Transport core

- [ ] T13: Tetrad construction + null-closure helpers — (a) 3×3 analytic inverse metric `γ^{ij}` + `p^t = √(γ^{ij}p_ip_j)/α`, recomputed not cached ([MCNX-GEO-02]); (b) orthonormal tetrad with timelike leg `u^μ`, Gram–Schmidt in FIXED coordinate order (x,y,z) — determinism-critical pin — plus isotropic fluid-frame direction draw and transform to coordinate frame ([MCNX-PKT-04])
  - spec: specs/packet-representation-and-sampling.md
  - tests: unit selftests on synthetic metric/velocity inputs — tetrad orthonormality, null-vector norm, ν-recovery identity, flat-static limit (machine tier).
  - notes: shared by INT scattering ([MCNX-INT-04]) and TRP diffusion ([MCNX-TRP-09/10]); implement exactly once; anchors `specs/packet-representation-and-sampling.md:119-140`, `specs/geodesic-propagation.md:87-97`.
- [ ] T26: Statistical-acceptance shared machinery — shared 4σ estimator/z-score reduction used by `stats-emission`/`stats-beam`/`stats-scatterbox`/`stats-diffusion` at `N_p=2^20`, seeds `1296518744`/`20260730` ([MCNX-RNG-07], [MCNX-VER-07])
  - spec: specs/rng-and-statistical-acceptance.md
  - tests: exercised by every stats benchmark; contract at `specs/verification-suite-design.md:255-279`.
  - notes: build WITH the first stats benchmark (T14) and reuse — do not duplicate per benchmark; listed separately so consolidation is explicit.
- [ ] T14: Emission sampling: count law + packet creation — packet state schema onto T8 container ([MCNX-PKT-01]); emission count `N_p = g_s √−g ΔV Δt η_b(s)/E_p` with floor/Bernoulli remainder, g=4 applied exactly once at creation ([MCNX-PKT-02]); bin-center energy, weight `N=E_p/ν`, uniform position/time draws ([MCNX-PKT-03]); pinned `(q,e,k)` creation draw map with `K_cell` bit packing ([MCNX-PKT-05]); id uniqueness via `idcpu`, bit 63 never set ([MCNX-GPU-04])
  - spec: specs/packet-representation-and-sampling.md
  - tests: `emission-fixedseed` benchmark — exact per-packet event/state sequence at pinned seed, single rank, golden parity 1e-12; `stats-emission` — count/energy/species estimators within 4σ at `N_p=2^20`, seed `1296518744`; draw-map audit trace (exact); id-contract test (bit 63, uniqueness).
  - notes: uses T12 analytic η; direction draw consumes T13; anchors spec `:90-162`, `specs/particle-container-and-gpu.md:180-188`.
- [ ] T15: Geodesic push: gather + RHS + integrator — trilinear-minimum interpolation of vertex-centered ADMBaseX `α, β^i, γ_ij` + order-≥2 derivatives as trivially-copyable device functor; RHS `dx^i/dt`, `dp_i/dt`; integrator free but ≥2nd order ([MCNX-GEO-01/03/04])
  - spec: specs/geodesic-propagation.md
  - tests: linear-field gather exactness (machine tier, synthetic `Array4`s); flat-spacetime free-stream exactness (machine tier); `minkowski-freestream` benchmark (needs no TestMCNuX code).
  - notes: centering confirmed `CarpetX/ADMBaseX/interface.ccl:12,16,20`; gather exemplar `warpx/Source/Particles/Gather/FieldGather.H:348`; uses T13 null closure at every RHS; do NOT template on CarpetX's interpolator ([MCNX-GPU-07]); spec anchor `:73-115`.
- [ ] T16: `schwarzschild-pt` benchmark — parfile + golden data: packet orbits on T6's Schwarzschild background, `p_t` drift ≤1e-5 at `h=M/16`, `T=100M` ([MCNX-GEO-05])
  - spec: specs/geodesic-propagation.md
  - tests: the benchmark itself (explicit numeric tolerance pinned per benchmark) + golden parity 1e-12.
  - notes: depends on T6 + T15; separate build+test cycle from both; anchors spec `:117-132`, `specs/verification-suite-design.md:241-242`.

## Tier 5 — Coupling & interactions

- [ ] T17: Source-term grid variables + zero-then-add schedule — `MCNuX/interface.ccl`: cell-centered `MCNuX::rad_force{rf_t,rf_x,rf_y,rf_z}`, `MCNuX::lepton_source{lep_src}` (`CENTERING={ccc}`, checkpointed); zero-then-add protocol mirroring `CarpetX/TmunuBaseX/schedule.ccl:1-34` (zero routine + named extension group, [MCNX-HYD-03]); deposition/ledger accumulation protocol ([MCNX-HYD-01/02]) provable with a synthetic event list before real physics
  - spec: specs/hydro-coupling-source-terms.md
  - tests: zero-then-add ordering observable (variables zeroed at one schedule point, accumulated, only then read); synthetic-event deposition unit test — deposited totals equal injected event sums to floating-point-reduction tolerance; READS/WRITES completeness enforced by CarpetX at runtime.
  - notes: spec anchors `:88-97,124-131`.
- [ ] T18a: Interaction draw + episode competition + RNG map — `Δt_{s,a} = −ln(1−u_{s,a}) p^t/(κ_{s,a} ν)`; per-cell episode boundaries (creation/step-start, post-scatter, cell-crossing — normative deviation from papers, golden-determining); shorter-of firing; zero-opacity `+∞` sentinel; tie-break to scattering; pinned `(k=0..3)` draw order, undrawn-remainder discard ([MCNX-INT-01/02/06])
  - spec: specs/neutrino-matter-interactions.md
  - tests: draw-map audit trace (exact); unit tests with synthetic `p^t`/cell-exit inputs; `interactions-fixedseed` golden parity.
  - notes: spec anchors `:83-126,165-179`; coefficients exclusively via T12/T10 interface, cell-centered ([MCNX-INT-05]).
- [ ] T18b: Absorption — full packet removal (no weight decay), full-content transfer to T17 ledger, id retirement ([MCNX-INT-03])
  - spec: specs/neutrino-matter-interactions.md
  - tests: `stats-beam` — beam attenuation e^{−κΔ} within 4σ; discreteness check (weights never decay); ledger transfer sums exact per event.
- [ ] T18c: Elastic scattering — isotropic fluid-frame redraw at fixed ν via T13 tetrad; momentum-change into ledger ([MCNX-INT-04])
  - spec: specs/neutrino-matter-interactions.md
  - tests: `stats-scatterbox` — isotropy + elasticity (ν preserved) within 4σ; fixed-seed event-sequence parity.
- [ ] T19: Conservation-ledger closure check — per-step audit: event-side sums vs grid totals for energy/momentum/lepton number, relative 1e-13/step; test-configuration-only allowed; escape tallies (per-species energy/number at outer boundary, README boundary policy) belong to this ledger surface ([MCNX-HYD-05])
  - spec: specs/hydro-coupling-source-terms.md
  - tests: closure to 1e-13/step in fixed-seed benchmarks (T14/T18 configurations).
  - notes: depends on T17 + T18; spec anchors `:177-188,295-296`.
- [ ] T20: `cadence-lag` benchmark — golden benchmark demonstrating source terms advance exactly once per coarsest Δt with the one-step lag semantics ([MCNX-CTX-05] verification)
  - spec: specs/verification-suite-design.md
  - tests: golden parity 1e-12.
  - notes: depends on T4 + T14/T17; spec anchor `:238-247`.

## Tier 6 — Trapped regime, GPU invariants, delivery gates

- [ ] T21a: TRP relabeling pure map — post-map `η′=αη, κ_a′=ακ_a, κ_s′=κ_s+(1−α)κ_a` with two exact invariants to ~1e-14; binding 2021 α convention (arXiv:2103.16588v2), guard against flipped 2020-letter convention ([MCNX-TRP-01..04] map portion)
  - spec: specs/trapped-regime-treatment.md
  - tests: invariant/monotonicity selftests (machine tier) — pure function, testable immediately after T5.
  - notes: spec anchor `:95-143`; β-dependent Eq. 50 bound enforcement is an explicit spec open question (`:450-455`) — NOT a deliverable.
- [ ] T21b: α selection + bitwise-limit benchmarks — per-cell `α = min(1, ξ/(κ_a Δt_c))`, fixed-α override parameter, `κ_a′Δt_c ≤ ξ` per-cell diagnostic; `trp-alpha-limit-explicit`/`trp-alpha-limit-relabeled` benchmarks ([MCNX-TRP-05])
  - spec: specs/trapped-regime-treatment.md
  - tests: α=1 run reproduces explicit run BITWISE at fixed seed, single rank (spec `:173-181`; `specs/verification-suite-design.md:247,375`).
  - notes: bitwise-limit benchmark needs T18; α-selection anchor `:145-171`.
- [ ] T22a: TRP diffusion kinematics pure functions — truncated-Gaussian displacement branch selector + antiderivative `F` inversion (≤1e-12 in r̃), causal cap `r_d ≤ Δt′_rem` ([MCNX-TRP-08]); two-leg position update ([MCNX-TRP-09]); correlated momentum draw with 21-point `B_i` table, atom-branch `cosθ₂=1` at `f_free=1` ([MCNX-TRP-10])
  - spec: specs/trapped-regime-treatment.md
  - tests: inversion accuracy 1e-12; table/branch pins exact; unit tests on synthetic inputs via T13 tetrad.
  - notes: spec anchors `:219-321`.
- [ ] T22b: TRP regime entry + episode dispatch + 5-draw map — `κ_s′Δt′_rem > τ_diff` evaluated after EVERY performed scattering (documented deviation, spec `:182-217,507-522`), one diffusion episode per (packet, step), no absorption channel inside; fixed five-draw order ([MCNX-TRP-06/07])
  - spec: specs/trapped-regime-treatment.md
  - tests: `stats-diffusion` (4σ); `trp-overlap-consistency`; diffusion-leg draw-map audit in `interactions-fixedseed`.
  - notes: depends on T18.
- [ ] T23: GPU redistribution + ownership invariants — local `Redistribute` gated on bounded motion; after-step ownership invariant: every packet on the rank/level/grid/tile owning its position ([MCNX-GPU-05/06])
  - spec: specs/particle-container-and-gpu.md
  - tests: ownership-invariant test after each step; bounded-motion audit; single-rank determinism check.
  - notes: depends on T4 step loop + T15 motion; spec anchor `:190-209`.
- [ ] T24: Optional Tmunu contribution — parameter-gated (default off) routine `IN TmunuBaseX_AddToTmunu` (`CarpetX/TmunuBaseX/schedule.ccl:31-33`), MC stress-energy estimator (arXiv:1708.08452v2 Eq. 16 as restated in spec) with exact 8-point cell-to-vertex average `eT(v) += Σ_{d∈{0,1}³} q(v−d)/8` ([MCNX-HYD-06])
  - spec: specs/hydro-coupling-source-terms.md
  - tests: deposition-identity unit test for the 8-point average (machine tier); gated-off default leaves eT* untouched (exact).
  - notes: deferrable; unresolved double-counting question when enabled with four-force sources on evolving spacetime — spec amendment required before that combined path; do not enable both in any committed parfile (spec `:190-222,330-336`).
- [ ] T25b: Configure-time negative-path gate — script proving that removing the CarpetX/WeakLibInterp provider fails configuration ([MCNX-BLD-02] verification)
  - spec: specs/build-and-integration.md
  - tests: gate script's own pass/fail (exact); script location convention is implementation freedom.
- [ ] T25c: Activation-split audit sweep — `TestMCNuX` appears in ActiveThorns only within `MCNuX/test/*.par` ([MCNX-BLD-06])
  - spec: specs/build-and-integration.md
  - tests: sweep script pass/fail (exact).
  - notes: relevant once parfiles exist (after T5/T14).
- [ ] T25d: GPU-backend build of the pinned ThornList ([MCNX-BLD-03/04])
  - spec: specs/build-and-integration.md
  - tests: GPU-backend configure+build green; precision gate static_asserts hold with AMReX headers.
  - notes: currently only a CPU/MPI cfg exists in the sandbox.

## Decisions & non-blocking spec items

- INHERITS vs bare READS: whether `READS:` on `ADMBaseX::`/`HydroBaseX::` groups requires `INHERITS:` in `interface.ccl` under CarpetX flesh — unresolved; T4/T15/T17 implementer determines empirically at first build.
- PARAMCHECK vs schedule-tree ordering: whether an unresolvable `AFTER ODESolvers_Solve` (when `use_subcycling=yes`) errors before the T3 paramcheck abort fires — check at T4.
- Header/gate permanent home once `stub.cxx` retires — T1/T8 implementer decides; convention then binds.
- `.gitignore` build products: likely unnecessary (build lives in mounted `$CACTUSX` tree per `agent_scripts/sandbox/README.md:39`); confirm at first sandbox build.
- Validator path quirk: `specs/README.md:103` says `tools/validate_specs.sh`; actual script is `specs/tools/validate_specs.sh`. Cosmetic; fix only via a spec-edit task if the validator itself flags it.
- Explicitly spec-deferred (NOT deliverables): packet-count control beyond `E_p` (`specs/packet-representation-and-sampling.md:245-248`); β-dependent Eq. 50 stability bound (`specs/trapped-regime-treatment.md:450-455`); NES/Pair/Brem νx thermal-emission assembly (`specs/opacity-eos-evaluation.md:340-349`); HYD-06 simultaneous-use double-counting (`specs/hydro-coupling-source-terms.md:330-336`).

## Discovered since last plan

- [ ] Spec clarification: `specs/conventions-and-units.md:129` demands `1e-14` relative agreement against the *pinned 10-digit* factor values, which is arithmetically impossible (truncated literals deviate from full-precision recomputation by up to 2.1e-10; worst: emissivity 2.06e-10). T1 resolved it per [MCNX-CNV-04] (extra digits allowed): `rtol_pinned=1e-9` vs pinned literals, `rtol_machine=1e-14` for full-precision identities — rationale commented in `MCNuX/src/mcnux_units.hxx`. Spec-edit task: amend `:129` to state the two-tier tolerance so T5's `unit-selftest` doesn't inherit an unsatisfiable criterion. Scope includes authoring the amendment via spec-author (no new spec file; edit in place).
  - spec: specs/conventions-and-units.md
  - tests: `bash specs/tools/validate_specs.sh` exit 0 after the edit (exact).
  - notes: discovered 2026-08-14 during T1; low urgency — code already carries the rationale.
