# TODO — MCNuX build plan

Foundations are in: T1–T5 + T25a landed (units header, Philox RNG header, subcycling paramcheck guard + gate, transport cadence group, `unit-selftest` battery with 36 rows, spec-corpus validator green). The spec corpus (README + 12 coded specs) remains complete and coverage-closed — no specs need authoring; the only spec defects are cosmetic text fixes (T28). Codebase is clean: no drift, no duplicated constant literals, `TestMCNuX` is an empty compiling shell, `MCNuX_TransportStep` exists but holds only the cadence observer. Everything below is implementation + test artifacts against existing specs, in priority order.

## Standing facts every item inherits

- Acceptance is entirely programmatic: exact/bitwise (RNG KATs, enumerations, α=1 parity), machine ~1e-14, golden 1e-12 (Cactus harness ABSTOL=RELTOL=1e-12), relaxed 1e-10, 4σ statistical at `N_p=1048576` seed `1296518744` (secondary `20260730`), conservation 1e-13/step.
- Shared utilities are implemented ONCE and reused: constants/units/species = `MCNuX/src/mcnux_units.hxx` (done, reuse `species_lepton_number`/`degeneracy` — never re-author); Philox `u(S,q,e,k)` = `MCNuX/src/mcnux_rng.hxx` (done); tetrad transform in fixed (x,y,z) Gram–Schmidt order + 3×3 inverse-metric/null-closure helper (T13, reused by scattering [MCNX-INT-04] and diffusion [MCNX-TRP-09/10]); 4σ z-score/estimator machinery (T26, built with first stats benchmark T14). Any new physical constant (e.g. OPA's `hc = 1.239841984e-10 MeV cm`, `specs/opacity-eos-evaluation.md:338-379`) lands in `mcnux_units.hxx` exactly once.
- Selftest battery protocol: rows append at the END of `MCNuX/src/mcnux_selftest.cxx`, bump `SIZE=` in `MCNuX/interface.ccl` (currently 36), regenerate the golden TSV from a sandbox harness run — never hand-edit; routine hard-errors on table/SIZE mismatch. Header convention: flat `MCNuX/src/mcnux_*.hxx`, `namespace MCNuX`, compile-time selftests aggregate through one `#include` in `stub.cxx`.
- Parfile facts every parfile inherits: `Cactus::presync_mode = "mixed-error"` (required to survive CarpetX's own paramcheck), `IO::out_fileinfo = "axis labels"`, `IO::out_dir = $parfile`, `IO::parfile_write = no`, `CarpetX::out_metadata = no`; `ActiveThorns = "CarpetX IOUtil MCNuX"` suffices for Minkowski; `CarpetX::restrict_during_sync = no` whenever subcycling is on.
- Binding anti-patterns: never use AMReX's stateful RNG (`amrex/Src/Base/AMReX_Random.H`, [MCNX-RNG-05]); never adopt CarpetX's AoS/host-loop interpolator (`CarpetX/CarpetX/src/interp.hxx:38`, [MCNX-GPU-07]/[MCNX-CTX-07]). Grid access via the `ghext` relative-include idiom proven in `mcnux_cadence.cxx` (`#include "../../../CarpetX/CarpetX/src/driver.hxx"`); `REQUIRES CarpetX` alone suffices — never add `REQUIRES AMReX`.
- Grid tallies: unordered atomic accumulation, no bitwise cross-rank/device assertions ([MCNX-RNG-08]) — review-check on every deposition task.
- The Cactus regression harness (flesh) resolves only inside the agent sandbox (`agent_scripts/sandbox/setup.sh`); harness-dependent tests run there, not in the analysis environment. Non-harness gates follow the `agent_scripts/gate_paramcheck.sh` pattern (transient output under `$CACTUSX/GATE/mcnux/`, never in the repo) and get CI wiring per the existing `ci.yml` precedent.
- Consistency sweep (`specs/conventions-and-units.md:132`) forbids duplicated constant literals — standing review gate; currently clean (`mcnux_units.hxx`/`mcnux_rng.hxx` are the sole literal sites).

## Tier 1 — Core enablers (all unblocked)

- [x] T8: Pure-SoA packet container + AMReX precision gate — container deriving `amrex::ParticleContainerPureSoA` with the six-component schema (position ×3 built-in, momentum ×3 real, weight real, species int, event-counter int, id via `idcpu`); one exemplar `amrex::ParallelFor` tile kernel establishing the device-kernel pattern ([MCNX-GPU-01/02/03]); complete [MCNX-BLD-03] with `static_assert(sizeof(amrex::Real)==8)` / `sizeof(amrex::ParticleReal)==8`
  - spec: specs/particle-container-and-gpu.md
  - tests: compile-time layout pin — static_asserts that the type derives `ParticleContainerPureSoA` and tile-view is `std::is_trivially_copyable` (spec `:236-240`); precision-gate legs as compile-time asserts plus two new selftest battery rows (golden TSV regenerated); CPU/MPI build green via `agent_scripts/build.sh` (GPU-backend build deferred to T25d).
  - notes: DONE 2026-08-19. `MCNuX/src/mcnux_packets.hxx` (schema `RIdx`/`IIdx`, `PacketContainer : ParticleContainerPureSoA<7,2>`, layout pin `:173-216`, AMReX precision gate `:225-234`) + `mcnux_packets.cxx` (`make_packet_container()` via `ghext->patchdata.at(0).amrcore.get()`, exemplar `ParallelFor` tile kernel `touch_packets()` over `getParticleTileData()`); battery rows 36/37, SIZE=38, golden TSV verified byte-identical to a fresh harness product. First AMReX-particle compile against the CPU/MPI sandbox config confirmed green; `REQUIRES CarpetX` alone puts AMReX particle headers on the include path. Precision-gate home RATIFIED: AMReX legs permanent in `mcnux_packets.hxx`, language-level legs stay in `stub.cxx`, which is NOT retired (remains the compile-time aggregation point). Container/kernel deliberately dead at runtime (no schedule.ccl call site — spec-complete, wiring is follow-on work, see Discovered). `mcnux_selftest.cxx` is now an AMReX consumer (`#include <AMReX_REAL.H>`). Default allocator (std::allocator on CPU) kept; T25d may revisit with `SetArena`. Retry hygiene: after an interrupted build, delete `$CACTUSX/configs/mcnux/build/MCNuX/` before `build.sh` — otherwise stale objects make it a no-op relink and static_asserts never rerun. Unblocks T14, T15, T23, T25d.
- [ ] T9: WeakLibInterp call-boundary wrapper — thin wrappers over the five `_Point` families with pinned argument/units/log10/index contract; MeV→Kelvin (via T1 factor `mev_to_kelvin`) at EVERY temperature argument including NES/Pair; table-range awareness hooks ([MCNX-OPA-01/02/03])
  - spec: specs/opacity-eos-evaluation.md
  - tests: kernel-parity selftests vs direct WeakLibInterp calls (exact/machine tier); unit-boundary check — Kelvin conversion applied exactly once per argument; species-axis grounding vs committed `.h5ls` fixture shapes ([MCNX-OPA-08]) — live table species count/order vs exactly two electron-type EmAb/Iso datasets, closing the deferred `specs/conventions-and-units.md:131,148` cross-check. Fold into `unit-selftest`.
  - notes: WeakLibInterp is consume-only, never modified; entry points verified at `WeakLibInterp/src/eos/wli_eos.H:54`, `wli_eos_inversion.H:265-419`, `wli_opacity_emab_iso.H:76,229,249`, `wli_opacity_nes_pair.H:66,193,235`, `wli_opacity_brem.H:102-130`; spec anchor `:174-220`. Table io/residency via `wli_io_eos.H:78-100`, `wli_io_opacity.H:64-124`, `wli_table.H:60-115`. First task in this area sets the `mcnux_opacity*`-style convention (functor vs free functions is implementation freedom, spec `:325-336`).
- [ ] T13: Tetrad construction + null-closure helpers — (a) 3×3 analytic inverse metric `γ^{ij}` + `p^t = √(γ^{ij}p_ip_j)/α`, recomputed not cached ([MCNX-GEO-02]); (b) orthonormal tetrad with timelike leg `u^μ`, Gram–Schmidt in FIXED coordinate order (x,y,z) — determinism-critical pin — plus isotropic fluid-frame direction draw and transform to coordinate frame ([MCNX-PKT-04])
  - spec: specs/packet-representation-and-sampling.md
  - tests: unit selftests on synthetic metric/velocity inputs — tetrad orthonormality, null-vector norm, ν-recovery identity, flat-static limit (machine tier); selftest battery rows + golden TSV regen.
  - notes: shared by INT scattering ([MCNX-INT-04]) and TRP diffusion ([MCNX-TRP-09/10]); implement exactly once; anchors `specs/packet-representation-and-sampling.md:119-140`, `specs/geodesic-propagation.md:87-97`.
- [ ] T17: Source-term grid variables + zero-then-add schedule — `MCNuX/interface.ccl`: cell-centered `MCNuX::rad_force{rf_t,rf_x,rf_y,rf_z}`, `MCNuX::lepton_source{lep_src}` (`CENTERING={ccc}`, checkpointed); zero-then-add protocol mirroring `CarpetX/TmunuBaseX/schedule.ccl:1-34` (zero routine + named extension group, [MCNX-HYD-03]); deposition/ledger accumulation protocol ([MCNX-HYD-01/02]) provable with a synthetic event list before real physics
  - spec: specs/hydro-coupling-source-terms.md
  - tests: zero-then-add ordering observable (variables zeroed at one schedule point, accumulated, only then read; two-step re-run shows no cross-step accumulation leak); synthetic-event deposition unit test — deposited totals equal injected event sums to floating-point-reduction tolerance, incl. single-event sign fixtures per the [MCNX-HYD-02] convention; READS/WRITES completeness enforced by CarpetX at runtime.
  - notes: spec anchors `:88-97,124-131`. Zero+contributor routines live inside `MCNuX_TransportStep`; reuse `species_lepton_number` from `mcnux_units.hxx` — do not re-author. Group split (single vs per-variable) is implementation freedom (spec `:293-296`); deposition must honor [MCNX-RNG-08]. Unblocks T18b, T19, T20, T24.
- [ ] T21a: TRP relabeling pure map — post-map `η′=αη, κ_a′=ακ_a, κ_s′=κ_s+(1−α)κ_a` with two exact invariants to ~1e-14; binding 2021 α convention (arXiv:2103.16588v2), guard against flipped 2020-letter convention ([MCNX-TRP-01..04] map portion)
  - spec: specs/trapped-regime-treatment.md
  - tests: invariant/monotonicity selftests (machine tier) — pure function, testable immediately after T5.
  - notes: spec anchor `:95-143`; β-dependent Eq. 50 bound enforcement is an explicit spec open question (`:450-455`) — NOT a deliverable. Selftest row earmark already exists in the battery's index map; SIZE bump + golden regen on landing.

## Tier 2 — Independent fixtures, gates & spec hygiene (unblocked)

- [ ] T6: TestMCNuX Schwarzschild initial data — `SHARES: ADMBaseX` + `EXTENDS KEYWORD initial_data`/`initial_lapse` adding `"Schwarzschild"`; writers `IN ADMBaseX_InitialData`/`ADMBaseX_InitialGauge` with full `WRITES: (everywhere)`; isotropic `ψ=1+M/(2r)`, `γ_ij=ψ⁴δ_ij`, `K_ij=0`, `α=(1−M/(2r))/(1+M/(2r))`, `β^i=0` ([MCNX-VER-02/03])
  - spec: specs/verification-suite-design.md
  - tests: smoke parfile + golden TSV comparing grid output against the closed forms (machine/golden tier).
  - notes: fully unblocked — depends only on stable base thorns; anchors `specs/verification-suite-design.md:139-181`, `CarpetX/ADMBaseX/schedule.ccl:12,16,68-104`, `CarpetX/ADMBaseX/interface.ccl:12-20`; shift needs no writer (default `"zero"`). First `.cxx` under `TestMCNuX/src/` — update its `SRCS`. Risk (spec `:456-461`): cross-repository `EXTENDS KEYWORD` mechanism unregressed in flesh — this task confirms it empirically.
- [ ] T7: TestMCNuX hydro profiles — `SHARES: HydroBaseX` + `EXTENDS KEYWORD initial_hydro` adding `"uniform sphere"`, `"equilibration box"`; writers `IN HydroBaseX_InitialData` for `(ρ,T,Yₑ)` profiles, `vel=0`, other fields zeroed, temperature in MeV; parameters for `M`, sphere center/radius, `(ρ₀,T₀,Yₑ₀)`, ambient values ([MCNX-VER-04])
  - spec: specs/verification-suite-design.md
  - tests: smoke parfile + golden comparing profile fields to formulas (golden tier); grid pinned so sharp-surface cells are ≤1% volume effect (spec `:456-493` assumption 4).
  - notes: independent build+test cycle from T6 (separate keyword family/writer); anchors `specs/verification-suite-design.md:182-201`, `CarpetX/HydroBaseX/interface.ccl:11-29`, `CarpetX/HydroBaseX/schedule.ccl:17,59-75`. Needed by `stats-*`, `interactions-fixedseed`, T29.
- [ ] T27a: Validity-machinery gate — gate script running an existing-machinery parfile with `CarpetX::poison_undefined_values = yes` (clean pass) plus a deliberately under-declared READS/WRITES meta-test variant asserting the driver's checksum/poison error fires (`specs/carpetx-thorn-integration.md:204-208`)
  - spec: specs/carpetx-thorn-integration.md
  - tests: the gate itself is the test — exact-criterion pass/fail on both legs.
  - notes: currently satisfied by nothing — no parfile sets `poison_undefined_values`, no meta-test exists; needs no new physics. Precedent `agent_scripts/gate_paramcheck.sh`; planner decision taken: gate-now (spec treats the bullet as independent), CI wiring per precedent.
- [ ] T27b: Scheduling-mode multi-box artifact — committed reproducible artifact (harness benchmark parfile with ≥2 boxes per rank + golden `transport_step_count` output, or a gate script asserting the invocation count) proving exactly one observer invocation per (patch,level) traversal (`specs/carpetx-thorn-integration.md:209-211`)
  - spec: specs/carpetx-thorn-integration.md
  - tests: golden `transport_step_count` TSV (or gate exact-count check).
  - notes: verified only manually during T4 (64-box rerun, 10 firings in 10 iterations) — not reproducible via `agent_scripts/test.sh` today. May later merge with T23's `ownership-multibox`, but needs no particle code and can close now.
- [ ] T25b: Configure-time negative-path gate — script building a ThornList variant omitting the CarpetX (and/or WeakLibInterp) provider and asserting configure failure ([MCNX-BLD-02] verification)
  - spec: specs/build-and-integration.md
  - tests: gate script's own pass/fail (exact); script location convention is implementation freedom.
  - notes: spec anchor `:186-189`; CI wiring per gate precedent.
- [ ] T25c: Activation-split audit sweep + layout audit — script asserting `TestMCNuX` appears in ActiveThorns only within `MCNuX/test/*.par` ([MCNX-BLD-06]) + smoke configure/build with a ThornList omitting `MCNuX/TestMCNuX`; fold in the layout-audit bullet (`specs/build-and-integration.md:203-206` — exactly two thorns, no `TestMCNuX/test/`, 13 spec files + validator)
  - spec: specs/build-and-integration.md
  - tests: sweep script pass/fail (exact).
  - notes: relevant once parfiles exist (after T5/T14). Parfiles now exist (`unit-selftest.par`); property currently holds by inspection. Planner decision taken: layout audit folds into this script rather than staying manual.
- [ ] T28: Spec-corpus text fixes — (a) amend `specs/conventions-and-units.md:129`, which demands `1e-14` relative agreement against the *pinned 10-digit* factor values — arithmetically impossible (truncated literals deviate from full-precision recomputation by up to 2.1e-10; worst: emissivity 2.06e-10); T1 resolved it per [MCNX-CNV-04] (extra digits allowed) with `rtol_pinned=1e-9` vs pinned literals and `rtol_machine=1e-14` for full-precision identities — the amendment states this two-tier tolerance so the `unit-selftest` criterion is satisfiable; (b) fix `specs/README.md:103` validator invocation (`tools/validate_specs.sh` → `specs/tools/validate_specs.sh`)
  - spec: specs/conventions-and-units.md
  - tests: `bash specs/tools/validate_specs.sh` exit 0 after the edit (exact).
  - notes: discovered 2026-08-14 during T1; low urgency — code already carries the rationale (commented in `MCNuX/src/mcnux_units.hxx`). Scope includes authoring the amendment via spec-author (no new spec file; edit in place); bundle both edits in one iteration.

## Tier 3 — Opacity/EOS layer

- [ ] T10: Baseline coefficient assembly + νx mapping — `κ_a(s,E)`, `κ_s(s,E)`, `η(s,E)` for νe/ν̄e from EmAb+Iso with Kirchhoff closure using EOS chemical potentials at (ρ,T,Yₑ); νx: `κ_a=0`, `η=0`, `κ_s` = ½-mean of the two electron-type Iso values ([MCNX-OPA-04/05])
  - spec: specs/opacity-eos-evaluation.md
  - tests: Kirchhoff-closure identity checks at pinned table points; νx zero/half-mean identities (machine tier); selftest rows.
  - notes: depends on T9; NES/Pair/Brem thermal-emission assembly is explicitly out of scope — do not attempt (spec `:340-349`). If the `hc` constant is needed it goes in `mcnux_units.hxx` once.
- [ ] T11: Table-range enforcement — transparency floor (ρ<ρ_min ⇒ all coefficients zero, no table touch), component-wise clamping elsewhere, per-axis clamp diagnostic counters, no-NaN guarantee ([MCNX-OPA-06]); `EosInversionResult{T,Error}` hard-fail-in-test-mode protocol ([MCNX-OPA-07]) may be stubbed since baseline transport receives T in MeV externally
  - spec: specs/opacity-eos-evaluation.md
  - tests: range-policy selftests — out-of-range probes return clamped/zero values with correct counter increments, never NaN (exact); inversion-protocol fixture.
  - notes: WeakLibInterp extrapolates permissively — does zero range-checking by design (`wli_eos.H:50-52`) — so enforcement is entirely MCNuX's (README boundary policy); spec anchor `:174-272`. Inversion leg stub-eligible: baseline path never inverts (spec `:261-263`).
- [ ] T12: Analytic/gray opacity mode + source-agnostic dispatch — per-species `kappa_a0`/`kappa_s0`/`eta_scale` parameters in `MCNuX/param.ccl` + Kirchhoff-form analytic coefficients; runtime-parameter-selected coefficient interface (table vs analytic); formulas pinned at `specs/verification-suite-design.md:202-233` ([MCNX-VER-05])
  - spec: specs/verification-suite-design.md
  - tests: analytic-mode coefficient values vs closed forms (machine tier); dispatch selection exactness (parameter flips source, values change accordingly).
  - notes: enabling fixture for all fixed-seed/statistical transport benchmarks (T14–T16, T18) without table dependence; only needs the interface shape, so it may precede T10/T11 under sequencing pressure. OPA requires the interface be source-agnostic (`specs/opacity-eos-evaluation.md:119-124`).

## Tier 4 — Transport core

- [ ] T15: Geodesic push: gather + RHS + integrator — trilinear-minimum interpolation of vertex-centered ADMBaseX `α, β^i, γ_ij` + order-≥2 derivatives as trivially-copyable device functor; RHS `dx^i/dt`, `dp_i/dt`; integrator free but ≥2nd order ([MCNX-GEO-01/03/04])
  - spec: specs/geodesic-propagation.md
  - tests: linear-field gather exactness (machine tier, synthetic `Array4`s); flat-spacetime free-stream exactness (machine tier); `minkowski-freestream` benchmark (needs no TestMCNuX code).
  - notes: centering confirmed `CarpetX/ADMBaseX/interface.ccl:12,16,20`; gather exemplar `warpx/Source/Particles/Gather/FieldGather.H:348`; uses T13 null closure at every RHS; do NOT template on CarpetX's interpolator ([MCNX-GPU-07]); spec anchor `:73-115`. Same gathered snapshot for every RHS term. Needs T8 + T13; first cross-thorn READS resolves the INHERITS-vs-bare-READS open question empirically.
- [ ] T16: `schwarzschild-pt` benchmark — parfile + golden data: packet orbits on T6's Schwarzschild background, `p_t` drift ≤1e-5 at `h=M/16`, `T=100M` ([MCNX-GEO-05])
  - spec: specs/geodesic-propagation.md
  - tests: the benchmark itself (explicit numeric tolerance pinned per benchmark) + convergence order ≥2 across Δt, Δt/2, Δt/4 + golden parity 1e-12.
  - notes: depends on T6 + T15; separate build+test cycle from both; anchors spec `:117-132`, `specs/verification-suite-design.md:241-242`. The 1e-5 drift bound is a design bound, re-pinnable at golden-creation time per spec `:456-493`.
- [ ] T26: Statistical-acceptance shared machinery — shared 4σ estimator/z-score reduction used by `stats-emission`/`stats-beam`/`stats-scatterbox`/`stats-diffusion`/`stats-equilibration` at `N_p=2^20`, seeds `1296518744`/`20260730` ([MCNX-RNG-07], [MCNX-VER-07])
  - spec: specs/rng-and-statistical-acceptance.md
  - tests: exercised by every stats benchmark; contract at `specs/verification-suite-design.md:255-279`.
  - notes: build WITH the first stats benchmark (T14) and reuse — do not duplicate per benchmark; listed separately so consolidation is explicit. Per-check estimate/expected/σ/z written as grid scalars/norm outputs; σ formula comes from each consuming benchmark's leaf spec; API shape free (`specs/rng-and-statistical-acceptance.md:125-130`).
- [ ] T14: Emission sampling: count law + packet creation — packet state schema onto T8 container ([MCNX-PKT-01]); emission count `N_p = g_s √−g ΔV Δt η_b(s)/E_p` with floor/Bernoulli remainder, g=4 applied exactly once at creation via `MCNuX::degeneracy` ([MCNX-PKT-02]); bin-center energy, weight `N=E_p/ν`, uniform position/time draws ([MCNX-PKT-03]); pinned `(q,e,k)` creation draw map (count draw, then k=0..5 per packet) with `K_cell` bit packing ([MCNX-PKT-05]); id uniqueness via `idcpu`, bit 63 never set ([MCNX-GPU-04])
  - spec: specs/packet-representation-and-sampling.md
  - tests: `emission-fixedseed` benchmark — exact per-packet event/state sequence at pinned seed, single rank, golden parity 1e-12; `stats-emission` — count/energy/species estimators within 4σ at `N_p=2^20`, seed `1296518744`, archived z-values then golden-frozen at 1e-12 (no seed-shopping on a false 4σ failure); draw-map audit trace (exact); id-contract test (bit 63, uniqueness).
  - notes: uses T12 analytic η; direction draw consumes T13; anchors spec `:90-162`, `specs/particle-container-and-gpu.md:180-188`. Builds T26's machinery with it.

## Tier 5 — Coupling & interactions

- [ ] T18a: Interaction draw + episode competition + RNG map — `Δt_{s,a} = −ln(1−u_{s,a}) p^t/(κ_{s,a} ν)`; per-cell episode boundaries (creation/step-start, post-scatter, cell-crossing — normative deviation from papers, golden-determining); shorter-of firing vs cell exit and Δt_rem; zero-opacity `+∞` sentinel, no NaN; tie-break to scattering; pinned `(k=0..3)` draw order, undrawn-remainder discard, e incremented per episode ([MCNX-INT-01/02/06])
  - spec: specs/neutrino-matter-interactions.md
  - tests: draw-map audit trace (exact); unit tests with synthetic `p^t`/cell-exit inputs; NaN-free zero-opacity fixture; `interactions-fixedseed` golden parity; `stats-beam` via T26 machinery.
  - notes: spec anchors `:83-126,165-179`; coefficients exclusively via T12/T10 interface, cell-centered ([MCNX-INT-05]) — never re-derived. Preserve spec-normative deviations verbatim (per-cell episode redraw, frozen episode kinematics, tie→scattering, spec `:249-271`).
- [ ] T18b: Absorption — full packet removal (no weight decay), full-content transfer to T17 ledger, id retirement (never reused) ([MCNX-INT-03])
  - spec: specs/neutrino-matter-interactions.md
  - tests: `stats-beam` — beam attenuation e^{−κΔ} within 4σ; discreteness check (weights never decay); ledger transfer sums exact per event.
  - notes: part of the `interactions-fixedseed` golden benchmark (`specs/verification-suite-design.md:244`, needs T7 hydro data).
- [ ] T18c: Elastic scattering — isotropic fluid-frame redraw at fixed ν via T13 tetrad; momentum-change into ledger ([MCNX-INT-04])
  - spec: specs/neutrino-matter-interactions.md
  - tests: `stats-scatterbox` — isotropy + elasticity (ν preserved, machine tier) within 4σ; fixed-seed event-sequence parity.
  - notes: reuse T13's tetrad — never reimplement (spec `:136-154`).
- [ ] T19: Conservation-ledger closure check — per-step audit: event-side sums vs grid totals for energy/momentum/lepton number, relative 1e-13/step; test-configuration-only allowed; escape tallies (per-species energy/number at outer boundary, README boundary policy) belong to this ledger surface ([MCNX-HYD-05])
  - spec: specs/hydro-coupling-source-terms.md
  - tests: closure to 1e-13/step in fixed-seed benchmarks (T14/T18 configurations), asserted as a golden pass flag.
  - notes: depends on T17 + T18; spec anchors `:177-188,295-296`. Synthetic-event arithmetic validation can land with/after T17; real closure needs real producers. Always-on vs test-only audit computation is an implementer decision (spec `:295-296`).
- [ ] T20: `cadence-lag` benchmark — golden benchmark demonstrating source terms advance exactly once per coarsest Δt with the one-step lag semantics ([MCNX-CTX-05] verification): step-function fluid perturbation at iteration n first affects source terms at n+1
  - spec: specs/verification-suite-design.md
  - tests: golden parity 1e-12 (counter + source-term tally checksum TSV).
  - notes: depends on T4 + T14/T17; spec anchor `:238-247`. Counter half already exists (`transport_step_count` in `MCNuX/interface.ccl`).

## Tier 6 — Trapped regime, GPU invariants, delivery gates

- [ ] T21b: α selection + TRP parameter surface + bitwise-limit benchmarks — per-cell `α = min(1, ξ/(κ_a Δt_c))`, fixed-α override parameter, `κ_a′Δt_c ≤ ξ` per-cell diagnostic; `MCNuX/param.ccl` entries for ξ (default 1), τ_diff (default 3), scheme selector (names free, semantics binding, spec `:434-436` — param.ccl currently has zero TRP entries); `trp-alpha-limit-explicit`/`trp-alpha-limit-relabeled` benchmarks ([MCNX-TRP-04/05])
  - spec: specs/trapped-regime-treatment.md
  - tests: selection-rule selftest rows (exact); α=1 run reproduces explicit run BITWISE at fixed seed, single rank (spec `:173-181`; `specs/verification-suite-design.md:247,375`), golden 1e-12.
  - notes: bitwise-limit benchmark needs T18; α-selection anchor `:145-171`.
- [ ] T22a: TRP diffusion kinematics pure functions (standalone math half) — truncated-Gaussian displacement branch selector + antiderivative `F` inversion (≤1e-12 in r̃), causal cap `r_d ≤ Δt′_rem` ([MCNX-TRP-08]); correlated momentum θ₂ draw with 21-point `B_i` table + linear interpolation, atom-branch `cosθ₂=1` at `f_free=1` ([MCNX-TRP-10])
  - spec: specs/trapped-regime-treatment.md
  - tests: selftest rows — F/F⁻¹ round-trip ≤1e-12; branch-selector fixtures; `B_i` interpolation fixtures; atom-branch pin (exact).
  - notes: spec anchors `:219-321`. Split this plan run: pure scalar math is unblocked and testable now via battery rows; the integration leg (two-leg position update, [MCNX-TRP-09]) moved to T22b with the rest of the blocked dispatch work.
- [ ] T22b: TRP regime entry + episode dispatch + position/momentum integration + 5-draw map — `κ_s′Δt′_rem > τ_diff` evaluated after EVERY performed scattering (documented deviation, spec `:182-217,507-522`), one diffusion episode per (packet, step), no absorption channel inside; two-leg position update (advection leg + free-streaming leg via T15's pusher, [MCNX-TRP-09]); tetrad-transformed momentum outcome via T13; fixed five-draw order ([MCNX-TRP-06/07])
  - spec: specs/trapped-regime-treatment.md
  - tests: `stats-diffusion` (three-clause, 4σ via T26); `trp-overlap-consistency`; scheme identities; diffusion-leg draw-map audit in `interactions-fixedseed` (spec `:361-406`).
  - notes: depends on T18. Also needs T13, T15, T22a. Preserve the spec's four logged deviations verbatim (`:443-557`).
- [ ] T23: GPU redistribution + ownership invariants — local `Redistribute`/`RedistributeLocal(N)` inside `MCNuX_TransportStep` after position-changing ops, gated on bounded motion (N from a proven motion bound); after-step ownership invariant: every packet on the rank/level/grid/tile owning its position ([MCNX-GPU-05/06])
  - spec: specs/particle-container-and-gpu.md
  - tests: `ownership-multibox` harness benchmark (≥2 boxes, ownership after step, golden); bounded-motion audit; single-rank determinism check (bitwise-identical populations across two fixed-seed runs — rides on T14's fixedseed benchmark).
  - notes: depends on T4 step loop + T15 motion; spec anchor `:190-209`; `specs/verification-suite-design.md:245`.
- [ ] T24: Optional Tmunu contribution — parameter-gated (default off) routine `IN TmunuBaseX_AddToTmunu` (`CarpetX/TmunuBaseX/schedule.ccl:31-33`), MC stress-energy estimator (arXiv:1708.08452v2 Eq. 16 as restated in spec) with exact 8-point cell-to-vertex average `eT(v) += Σ_{d∈{0,1}³} q(v−d)/8` ([MCNX-HYD-06])
  - spec: specs/hydro-coupling-source-terms.md
  - tests: deposition-identity unit test for the 8-point average (machine tier) + sign fixtures; gated-off default leaves eT* untouched (exact).
  - notes: deferrable; unresolved double-counting question when enabled with four-force sources on evolving spacetime — spec amendment required before that combined path; do not enable both in any committed parfile (spec `:190-222,330-336`). Never mark this and the four-force coupling simultaneously "done" without addressing the open question.
- [ ] T25d: GPU-backend build of the pinned ThornList + device RNG purity leg ([MCNX-BLD-03/04], [MCNX-GPU-02], [MCNX-RNG-06] device leg)
  - spec: specs/build-and-integration.md
  - tests: GPU-backend configure+build green; precision gate static_asserts hold with AMReX headers; per-packet-operator device-kernel inspection; device-vs-host bitwise equality of `MCNuX::u` over the `mcnux_rng.hxx` `detail::purity_batch` tuples (exact).
  - notes: currently only a CPU/MPI cfg exists in the sandbox. CI compile-only cuda/rocm matrix already exists (`.github/workflows/ci.yml:173-247`); missing legs are the sandbox GPU harness run and the deferred host-vs-device RNG-purity check (deferral recorded in `mcnux_rng.hxx:32-40`; discovered 2026-08-14 during T2 — the CPU/MPI sandbox cannot exercise this leg, so T2 shipped with all host-side purity legs asserted at compile time). GPU-vs-golden norm tolerances explicitly undecided pending first GPU run (spec-recorded open).
- [ ] T29: `stats-equilibration-{explicit,relabeled}` benchmarks — equilibration-box configuration driving the full emission/absorption/scattering chain to thermal equilibrium; arbiter of OPA's stimulated-absorption and chemical-potential conventions (`specs/opacity-eos-evaluation.md:338-379`) ([MCNX-VER-07])
  - spec: specs/verification-suite-design.md
  - tests: 4σ-then-golden protocol via T26 machinery at `N_p=2^20`, seed `1296518744`; each variant documents its own σ formula from its governing leaf spec (`:255-279`).
  - notes: new this plan run — no prior TODO owner. Needs T7 (equilibration box), T10 (table coefficients), the T14/T18 chain, and T21a/T21b for the relabeled variant. Reuse T26 — no new statistical machinery.

## Decisions & non-blocking spec items

- INHERITS vs bare READS: whether `READS:` on `ADMBaseX::`/`HydroBaseX::` groups requires `INHERITS:` in `interface.ccl` under CarpetX flesh — still open; T15/T17 implementer determines empirically at first cross-thorn READS.
- RESOLVED at T4 — PARAMCHECK vs schedule-tree ordering: the flesh silently ignores an unresolvable `AFTER ODESolvers_Solve` (no schedule-resolution error precedes the T3 paramcheck abort; verified both with ODESolvers inactive and with `use_subcycling=yes` where only `ODESolvers_Solve_Subcycling` exists). Caveat: with the AFTER target absent the group's placement in `[CCTK_EVOL]` is arbitrary (observed *before* `ODESolvers_Solve_Subcycling`) — T3's guard is therefore load-bearing for [MCNX-CTX-05], not merely defensive; any future subcycling extension must revisit the ordering clause, not just drop the guard.
- RESOLVED at T8 — precision-gate permanent home: AMReX legs live in `MCNuX/src/mcnux_packets.hxx` (re-verified by every AMReX-consuming TU); language-level legs (`sizeof(double)`, `sizeof(CCTK_REAL)`) stay in `stub.cxx`, which is NOT retired — it remains the compile-time selftest aggregation point for every shared `mcnux_*.hxx`. Convention now binds.
- RESOLVED — `.gitignore` build products unnecessary: multiple sandbox builds (T3–T5) completed with the repo staying clean (build lives in mounted `$CACTUSX` tree per `agent_scripts/sandbox/README.md:39`).
- Gate CI wiring: CLAUDE.md calls gates "sandbox only" but `gate_paramcheck.sh` already runs in CI (`ci.yml:163-164`) — precedent decided: new gates (T27a, T25b, T25c) get CI wiring too.
- Scripted no-duplicate-literal sweep joining `validate_specs.sh` vs staying a manual review gate — tooling scope, not blocking; revisit if a violation ever slips through review.
- T4 recovery caveat (not exercised): if CarpetX doesn't mark recovered scalars valid under mixed-error, schedule the init routine also `AT CCTK_POST_RECOVER_VARIABLES`.
- Explicitly spec-deferred (NOT deliverables): packet-count control beyond `E_p` (`specs/packet-representation-and-sampling.md:245-248`); β-dependent Eq. 50 stability bound (`specs/trapped-regime-treatment.md:450-455`); NES/Pair/Brem νx thermal-emission assembly (`specs/opacity-eos-evaluation.md:340-349`); HYD-06 simultaneous-use double-counting (`specs/hydro-coupling-source-terms.md:330-336`).

## Discovered since last plan

- [ ] Container runtime wiring — give T8's `PacketContainer` a durable owner and a `make_packet_container()` call site in `MCNuX/schedule.ccl` (at/after `CCTK_BASEGRID`, once `MakeNewGrids` has run; timing/single-patch constraints recorded at `MCNuX/src/mcnux_packets.cxx:39-64`); currently the container and exemplar kernel are deliberately dead at runtime (spec-legal for T8, whose layout verification is compile-time-only). Natural landing spot: fold into T14's packet-creation work or schedule as its own small precursor — planner decides.
  - spec: specs/particle-container-and-gpu.md
  - tests: harness-visible evidence the container is constructed once per run without abort (smoke leg on an existing parfile suffices); full runtime behavior remains T14/T23 scope.
  - notes: discovered 2026-08-19 during T8 close-out. Related facts: T8's other follow-ons (packet creation + id contract, ownership invariants, GPU legs) are already owned by T14/T23/T25d — no new items needed for those.
