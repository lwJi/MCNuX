# TODO — MCNuX build plan

Fresh plan (prior TODO.md was deliberately discarded at c187327). Current state: every pure-function layer — RNG, units, tetrad, opacity/coefficients/table-range, geodesic math, ledger arithmetic, TRP relabel map, emission/interaction/fluid-gather pure functions — is complete, static_assert-covered, and mirrored in the 128-row selftest battery with committed goldens. BLD/CNV/CTX wiring and TestMCNuX are fully conformant. What is missing is the entire runtime physics pipeline: no scheduled routine creates packets from the emission law, none samples/applies interaction episodes, no real event feeds the source-term ledger, boundary-exiting packets are silently dropped, and ~17 of the ~23 spec-pinned benchmarks do not exist. Every task below is greenfield wiring/test-authoring on finished building blocks; none is a bug fix.

## Standing facts every item inherits

- **Fluid-gather reuse (C1):** the gather/tile-walk helpers now live in the shared header `mcnux_gather.hxx` (hoisted from `mcnux_geodesic.cxx` by emit-core): `require_driver`, `MetricGroups`/`metric_groups()`, `HydroGroups`/`hydro_groups()`, `make_gather`/`make_fluid_gather` (int-box + PacketIter overloads), `for_each_packet_tile[_raw]`, shared `fill_packet_diag`. Reuse from there, never duplicate.
- **Atomic deposition idiom (C2):** the [MCNX-GPU-02] idiom is built ONCE in `mcnux_deposit.hxx` (`SourceGroups`/`source_groups()`, `SourceViews`, `deposit_delta()` = `source_from_delta` + five `Gpu::Atomic::AddNoRet`), first production consumer `mcnux_emission.cxx`. `episode-driver` and `absorb-remove` reuse it with `scattering_delta`/`absorption_delta` (`mcnux_srcterms.hxx:62-87`).
- **Selftest growth:** the index→name row table in `mcnux_selftest.cxx` is APPEND-ONLY; new checks go in per-domain `mcnux_selftest_<domain>.cxx` files. The battery is exactly full (rows 0..SIZE−1; SIZE=140 as of emit-bins) — appending rows means bumping `interface.ccl` SIZE to the new count in the same change and re-capturing the `unit-selftest` golden TSV.
- **Coefficient access:** always via `RangedTableCoefficients`/the analytic dispatch in `mcnux_coefficients.hxx`, never bare `BaselineTableCoefficients`.
- **Old task-id breadcrumbs in code comments** (`T21b` in `mcnux_trp.hxx:36-41`, `T23` in `mcnux_geodesic.cxx:401-434`, `T25d` in `mcnux_rng.hxx:35-40`, plus old T14a/T18a/T33 commit prefixes) refer to the discarded pre-reset numbering, NOT this plan. Mapping: T21b→`trp-alpha`, T23→`escape-tallies`, T25d→documented RNG host-vs-device GPU deferral (tracked, not planned — sandbox has no GPU path; GPU compile checks are CI-only by design).
- Benchmarks whose acceptance is subjective/statistical use the repo's sanctioned 4σ pattern (`mcnux_stats.hxx`: pinned N_p=2^20, seed 1296518744, bound 4.0, non-strict ≤); sigma must be analytically propagated/documented, never computed by `mcnux_stats.hxx`.

## Tier 1 — Runtime transport event chain

- [x] [ledger-closure] Conservation-ledger closure diagnostic ([MCNX-HYD-05]) — buildable now, no dependencies
  - spec: specs/hydro-coupling-source-terms.md
  - tests: per-step five-channel (P^t, P_x, P_y, P_z, L) balance: event-side audit sum vs grid total Σ_cells source·ΔV·Δt at relative 1e-13 (ε₀ floor for all-zero) per `hydro-coupling-source-terms.md:177-188`; wire into `source-zero-add` (synthetic event list exactly known) and extend its golden with the closure observable (re-capture); selftest rows for the audit arithmetic on the synthetic fixture list
  - notes: DONE. Landed: `LedgerAudit` accumulator + `ledger_closure()` verdict (rtol 1e-13, documented ε₀=1e-30) and the shared `synthfix` synthetic event list in `mcnux_srcterms.hxx`; `MCNuX_LedgerClosure` routine in `mcnux_srcterms.cxx` (event-side audit from pre-negation `LedgerDelta`s, `amrex::ReduceSum` grid reduction with nghost=0, `CCTK_VERROR` on failing channel), scheduled twice — `AT initial AFTER MCNuX_ZeroSourceTerms` (defined it-0 output + runtime ε₀ branch) and `IN MCNuX_TransportStep AFTER MCNuX_AddToSourceTerms`; new `mcnux_ledger_diag` ARRAY SIZE=5 observable, `test_ledger_closure` param, selftest rows 128–133, `source-zero-add` + `unit-selftest` goldens re-captured. Reuse pointers for later consumers: event paths (emit-core/episode-driver/absorb-remove) must feed the audit from `LedgerDelta` pre-negation, never from `SourceContribution` (self-cancels); closure multiply must use the same ΔV·Δt the deposit divided by. CAVEAT: the grid reduction is rank-local and unigrid-only (commented, per `fill_packet_diag` precedent) — `interactions-fixedseed`/`ownership-multibox` need a rank reduction + level masking.
- [x] [emit-bins] Emission runtime pure helpers: energy-bin grid + uniform→physical mapping
  - spec: specs/packet-representation-and-sampling.md
  - tests: new constexpr static_asserts + appended selftest rows (exact-tier): bin-center/edge identities; `CreationUniforms` mapping endpoints (u=0 → cell lo / t_n, u→1 limit) and linearity
  - notes: DONE. Landed in `mcnux_emission.hxx`: bin-grid section (`num_energy_bins`, `EnergyBin`, `energy_bin(double const*, int, int)`, `bin_width`, `bin_center_energy(EnergyBin)` overload) + mapping section (`map_uniform_to_interval`, `creation_position(u, cell_lo, dx)`, `creation_time(u, t_n, dt)`); selftest rows 134..139 in new `mcnux_selftest_emission_bins.cxx`, SIZE 134→140, `unit-selftest` golden re-captured. Open questions SETTLED for emit-core: (1) bin edges are the tabulated energy nodes themselves — N nodes → N−1 bins, bin b = [E[b], E[b+1]], η interpolated at the arithmetic center (documented as an implementation decision, no spec pin); (2) `bin_integrated_eta` cadence is exactly once per (species, bin, cell, step) in the count-law path. Emit-core still owns: log10→linear conversion of table nodes (`BaselineTableViews.emab_LogEs` are log10; 10^x not constexpr), table wiring via `RangedTableCoefficients`, analytic-source binning, and the creation-loop call sites. Mapping consumes only the six existing `CreationUniforms` doubles (zero extra RNG draws); `creation_position` takes the `(cell_lo, dx)` shape matching `CellFluidGather`.
- [x] [emit-core] Runtime emission/creation loop (`mcnux_emission.cxx`, scheduled IN MCNuX_TransportStep)
  - spec: specs/packet-representation-and-sampling.md
  - tests: build + existing harness stays green; minimal fixed-seed smoke parfile capturing `mcnux_packet_diag` for a tiny analytic-opacity box (first-capture golden)
  - notes: DONE. Landed: `MCNuX_Emission` in new `mcnux_emission.cxx` (count pass → `Scan::ExclusiveSum` → contiguous `ParticleType::NextID` id reservation → tile resize → fill pass; WarpX append idiom, no prior in-repo precedent), scheduled `if (enable_emission)` IN MCNuX_AddToSourceTerms; the group is now `AFTER MCNuX_TransportObserver` (an AFTER on the routine alone would NOT bind — observer isn't a sibling inside the group), pinning the 0-based step index n = `*transport_step_count`−1. C1 hoisted to `mcnux_gather.hxx`, C2 built in `mcnux_deposit.hxx` (see standing facts). Id scheme: `SetParticleIDandCPU(NextID-reserved id, MyProc())`; per-packet RNG key is NOT raw idcpu but `packet_rng_key(id,cpu)=(cpu<<39)|id` (`mcnux_emission.hxx`, static_assert-pinned bit-63-clear + cell-key disjointness) — later event paths must key on this, not idcpu. Analytic fixture energy grid pinned: {5,15,25} MeV → 2 bins (in `mcnux_emission.cxx`; superseded per-table when table-residency lands). Count-law units: dV→cm³, dt→s, η stays Kirchhoff cgs; ν converted MeV→code at the tetrad-transform call site. `emission_step_audit()` (`LedgerAudit`, pre-negation deltas, same ΔV·Δt as deposit) is filled per step but NOT yet consumed by `MCNuX_LedgerClosure` (delta below). Paramcheck guards: `enable_emission` × table source, × `test_synthetic_packets` (fixture ids 1..8 bypass NextID). New params: `enable_emission`, `test_emission`, `rng_seed` (default = stats primary 1296518744). Smoke test `emit-smoke` (goldens committed): `eta_scale[0]=5.5e-54`, `kappa_a0[0]=1.0` (required — analytic η = eta_scale·kappa_a0·Kirchhoff, kappa_a0 defaults to 0) → exactly 7 packets/2 steps, both bins, ≈4% margins to the nearest count-changing scales. Emission diag reuses shared `fill_packet_diag` (SIZE=8 cap unchanged).
- [ ] [exit-time] Cell-exit-time pure helper (prerequisite for episode competition, [MCNX-INT-02])
  - spec: specs/neutrino-matter-interactions.md
  - tests: static_asserts + selftest rows — exact-tier axis-aligned flat cases; machine-tier agreement with an RK4-stepped boundary crossing
  - notes: time to next cell boundary given cell geometry + current position/momentum from the push; no home yet (candidate: alongside `geodesic_step` in `mcnux_geodesic.hxx` or in `mcnux_interactions.hxx`).
- [ ] [episode-driver] Episode driver + scattering event application (scheduled IN MCNuX_TransportStep)
  - spec: specs/neutrino-matter-interactions.md
  - tests: selftest rows for the scatter-redraw pure composition (ν preserved exactly, null closure machine-tier, flat-limit isotropy); fixed-seed scatter-box parfile + golden (per-packet event log / `mcnux_packet_diag`) verifying elasticity + determinism (`neutrino-matter-interactions.md:197-223`)
  - notes: depends on exit-time (emit-core for a realistic population; can also run on the synthetic fixture). Wires [MCNX-INT-01/02/04/05/06] live: per packet gather fluid (C1), κ_a/κ_s via RangedTableCoefficients, ν via `fluid_frame_energy`, cell-exit time, `episode_uniforms`/`interaction_time`/`compete_episode` (`mcnux_interactions.hxx:122-220`), increment `IntIdx::event_counter` at episode start per [MCNX-INT-06]; on Scatter: `scatter_uniforms`, redraw p^μ at fixed ν via `fluid_frame_direction`/`build_tetrad`/`fluid_momentum_to_coordinate`/`lower_spatial_momentum` (`mcnux_tetrad.hxx:244-320`), deposit `scattering_delta` (C2), begin new episode. Preserve the three pinned deviations: per-cell episodes, frozen episode kinematics, tie-to-scatter (`neutrino-matter-interactions.md:249-271`). Scattering-only first — all packets survive, no removal machinery yet. If too large, split the scatter-redraw pure composition into its own selftest-only increment first.
- [ ] [absorb-remove] Absorption event application + packet removal/retirement
  - spec: specs/neutrino-matter-interactions.md
  - tests: fixed-seed beam-attenuation-style parfile + golden: packet count decreases discretely, absorbed energy/lepton number appears in `rad_force`/`lepton_source` (discreteness + collision-statistics observables, `neutrino-matter-interactions.md:197-223`); id non-reuse assertion over the run
  - notes: depends on episode-driver. On `EpisodeEnd::Absorb`: deposit full packet content via `absorption_delta` (`mcnux_srcterms.hxx:72-77`, C2), remove packet with id retired per [MCNX-INT-03] + [MCNX-GPU-04]. No removal API is used anywhere in the repo today — confirm the retirement contract in `specs/particle-container-and-gpu.md` before assuming AMReX swap-and-pop suffices.

## Tier 2 — Event-chain benchmarks & population management

- [ ] [emit-fixedseed] `emission-fixedseed` benchmark + draw-map audit
  - spec: specs/verification-suite-design.md
  - tests: golden parity at 1e-12 (deterministic): pinned tiny emission step with per-packet TSV (species, bin, ν, N, position, time, direction); exact-tier trace of (q,e,k) tuples consumed (`packet-representation-and-sampling.md:203-205`)
  - notes: depends on emit-core. Closes the [MCNX-VER-06] row at `verification-suite-design.md:243` and [MCNX-PKT-05/06] deterministic verification; the statistical leg is merged into stats-writer (both need the same reduction writer).
- [ ] [schwarzschild-pt] `schwarzschild-pt` benchmark — buildable now, zero blockers
  - spec: specs/verification-suite-design.md
  - tests: harness golden parity + drift/convergence thresholds: pinned trajectory family at isotropic r∈[6M,10M] on TestMCNuX "Schwarzschild" data, T=100M, h=M/16, per-packet `p_t` via `mcnux_packet_diag`, drift ≤1e-5, plus Δt/Δt/2/Δt/4 convergence legs (order ≥2)
  - notes: golden-parity anchor for [MCNX-GEO-01/02/04/05] (`verification-suite-design.md:242,360-364`). Needs a Schwarzschild-appropriate variant of `PacketFixture` (`mcnux_geodesic.cxx:313-392`; current fixture is Minkowski-small-box-tuned) and its own domain — `schwarzschild-id.par`'s box (r≈[3.46,6.93]) does not cover [6M,10M]. Spec anticipates re-pinning h/T after first capture (`verification-suite-design.md:467-473`).
- [ ] [escape-tallies] Escape tallies + bounded-redistribution audit ([MCNX-GPU-05/06], GEO escape clause)
  - spec: specs/particle-container-and-gpu.md
  - tests: freestream-style parfile where a known packet subset exits: golden tallies equal the exact expected energy/number; audit observable in every push benchmark
  - notes: (a) detect packets crossing the outer boundary during/before `Redistribute()` (`mcnux_geodesic.cxx:401-434`, in-code "T23"), accumulate per-species energy+number escape tallies (`specs/README.md:20`; `geodesic-propagation.md:66-69`), remove explicitly instead of AMReX silent drop; decide tally storage surface (likely HYD-owned grid scalars — open question). (b) derive per-step motion bound N (max coord speed × dt / min cell width + margin), switch to `RedistributeLocal(N)` where valid, record derived-N vs measured max displacement (audit mode, `particle-container-and-gpu.md:250-254`). (a) and (b) are separable build+test cycles — split if either grows.
- [ ] [stats-writer] Statistical-reduction writer + `stats-emission` benchmark ([MCNX-PKT-06])
  - spec: specs/rng-and-statistical-acceptance.md
  - tests: the 4σ acceptance itself (non-deterministic backpressure — sanctioned pattern) on realized packet counts/energy/isotropy per `verification-suite-design.md:264-270`; programmatic given pinned seeds
  - notes: depends on emit-core. Scheduled routine computing (estimate, expected, sigma, z) via `MCNuX::zscore` (`mcnux_stats.hxx:77-95`) and writing TSV/norm output — no stats-* benchmark can exist without it. Writer + first benchmark are one task (writer is untestable alone).
- [ ] [interactions-fixedseed] `interactions-fixedseed` benchmark + event-driven ledger closure
  - spec: specs/verification-suite-design.md
  - tests: golden parity 1e-12: pinned-seed run with real creation/absorption/scattering events, per-packet event log TSV golden; ledger closes against the packet-side audit sum (ledger-closure diagnostic)
  - notes: depends on emit-core, episode-driver, absorb-remove, ledger-closure. Closes [MCNX-VER-06] row `verification-suite-design.md:244`, upgrades [MCNX-HYD-02] verification (matrix row :390). Retire/gate `MCNuX_SyntheticDeposit` (`schedule.ccl:216-226`) as the sole contributor.

## Tier 3 — Trapped-regime treatment

- [ ] [trp-alpha] Trapped regime stage 1: α selection rule + parameters ([MCNX-TRP-04], wiring for TRP-01/02)
  - spec: specs/trapped-regime-treatment.md
  - tests: selftest rows for the selection rule (exact endpoints: α=1 when ξ large; proportionality); paramcheck gate for invalid α; wiring verified by trp-alpha-parity
  - notes: depends on emit-core, episode-driver. Per-cell α = min(1, ξ/(κ_a Δt_c)) (needs cell light-crossing time Δt_c), diagnostic κ_a′Δt_c ≤ ξ, new `param.ccl` entries: ξ, fixed-α override, scheme-selector KEYWORD, τ_diff (for trp-diffusion), α∈(0,1] paramcheck (`mcnux_trp.hxx:36-41` names this "T21b"). Substitute primed coefficients transparently at the emission (η_b input) and event-sampling (κ inputs) call sites, gated by the selector — `relabel()` (`mcnux_trp.hxx:60-64`) currently has zero call sites. [MCNX-TRP-04] Eq.-50 β-dependent bound: spec defers enforcement ("choose ξ conservatively") — continue the documented deferral.
- [ ] [trp-alpha-parity] `trp-alpha-limit-{explicit,relabeled}` bitwise-parity benchmark ([MCNX-TRP-05])
  - spec: specs/trapped-regime-treatment.md
  - tests: paired-parfile golden parity: relabeled-at-α=1 run vs explicit run, identical per-packet event sequences, ABSTOL=RELTOL=1e-12 (`verification-suite-design.md:247`)
  - notes: depends on trp-alpha. The algebraic identity is already selftested; this is the run-level check.
- [ ] [trp-diffusion] Trapped regime stage 2: diffusion approximation ([MCNX-TRP-06..10]) + its benchmarks
  - spec: specs/trapped-regime-treatment.md
  - tests: exact/machine-tier selftest rows (static_asserts + appended rows) per pure piece; then 4σ statistical benchmarks `stats-diffusion`, `stats-equilibration-{explicit,relabeled}`, and `trp-overlap-consistency` via the stats-writer
  - notes: depends on episode-driver, trp-alpha, stats-writer. Four separable build+test cycles: (a) entry criterion κ_s′Δt′_rem > τ_diff evaluated after every performed scattering (`trapped-regime-treatment.md:192-217`); (b) two-branch displacement law — erf-based F inversion, p_ns, r̃_max + pinned five-draw RNG map (:298-339); (c) two-leg position update (advection leg reuses `mcnux_fluid.hxx` Valencia lift, free-stream leg reuses `geodesic_step`); (d) correlated momentum model — 21-point B_i table + linear interp + θ₂ closed-form draw (NOT uniform in cosθ; `fluid_frame_direction` is not reusable here). Equilibration benchmarks need the unused `"equilibration box"` TestMCNuX profile (`hydro_profiles.cxx:39`).

## Tier 4 — Remaining verification, integration & extensions

- [ ] [stats-interactions] Remaining stats benchmarks: `stats-beam`, `stats-scatterbox`
  - spec: specs/verification-suite-design.md
  - tests: 4σ acceptance via the stats-writer at pinned seeds — attenuation law; scattering isotropy (`verification-suite-design.md:264-270`; `neutrino-matter-interactions.md:197-223`)
  - notes: depends on episode-driver, absorb-remove, stats-writer.
- [ ] [cadence-lag] `cadence-lag` benchmark ([MCNX-CTX-05])
  - spec: specs/carpetx-thorn-integration.md
  - tests: harness golden parity 1e-12: real multi-stage ODESolvers RK4 with a trivial evolved state; golden-check `mcnux_transport_diag::transport_step_count` (never output in any current parfile) — one increment per CCTK_EVOL iteration, not per RK stage; injected fluid step-change at iteration n first affects source reads at n+1; checkpoint/recovery leg (`carpetx-thorn-integration.md:196-203`; `verification-suite-design.md:246`)
  - notes: reuses `source-zero-add` scaffolding; step-change injection fixture may need building — check overlap with ownership-multibox's restart leg.
- [ ] [ownership-multibox] `ownership-multibox` benchmark + [MCNX-HYD-03] meta-tests
  - spec: specs/verification-suite-design.md
  - tests: (a) ≥2-box-per-rank parfile with packets crossing box boundaries; per-invocation counter proving one invocation per (patch, level) for container-walking routines + `OK()`-style placement assertion after each step; (b) contributor-split meta-test (two artificial contributors ≡ single-contributor run to FP-reduction tolerance) + restart/recovery test reproducing partner-visible `rad_force`/`lepton_source` (`hydro-coupling-source-terms.md:248-254`)
  - notes: (a) closes [MCNX-CTX-03] AND [MCNX-GPU-05/06] verification in one artifact — the benchmark is double-cited (`verification-suite-design.md:245,405`), schedule it once. (b) is a separable second cycle — split into two tasks if both grow.
- [ ] [table-residency] Table-residency layer (real HDF5 → `BaselineTableViews`)
  - spec: specs/opacity-eos-evaluation.md
  - tests: parfile-level benchmark evaluating `RangedTableCoefficients` against a real table load, golden TSV; clamp-counter observables
  - notes: load production `wl-EOS-SFHo-15-25-50.h5` / `wl-Op-...-{EmAb,Iso}.h5` into device-resident `BaselineTableViews` (`mcnux_table_coeffs.hxx:63-65` flags "a residency layer later"; ownership straddles OPA and GPU specs — confirm in `specs/particle-container-and-gpu.md` before scoping). Enables table-sourced (non-analytic) runs of the emission loop and episode driver.
- [ ] [gates] Non-harness gates: configure-without-CarpetX; validity-machinery meta-test
  - spec: specs/build-and-integration.md
  - tests: the gate scripts' exact-criterion checks themselves
  - notes: (a) the one explicitly-required missing gate: configure-without-CarpetX failure ([MCNX-CTX-02]/[MCNX-BLD-02], `verification-suite-design.md:249-253`) — new `agent_scripts/` gate script. (b) deliberately-under-declared-routine meta-test proving poison/presync enforcement fires ([MCNX-CTX-04], `carpetx-thorn-integration.md:204-208`) — home (test vs CI gate) unresolved by specs (`carpetx-thorn-integration.md:406` vs `verification-suite-design.md:251-252`), decide when scoping. Other matrix `gate:` rows may stay inspection-only per Implementation freedom.
- [ ] [tmunu] [MCNX-HYD-06] Tmunu contract (optional, default-off, diagnostic-grade)
  - spec: specs/hydro-coupling-source-terms.md
  - tests: selftest closed-form rows (one-packet, uniform field, linear field, vertex-total identity) per `hydro-coupling-source-terms.md:263-269`
  - notes: low priority. q_{μν} MC stress-energy estimator + parameter-gated contributor IN TmunuBaseX_AddToTmunu with 8-point cell-to-vertex average; add TmunuBaseX to `interface.ccl` INHERITS + `configuration.ccl` (already in ThornList). Spec's own feedback-subtraction open question (`hydro-coupling-source-terms.md:330-336`) means this lands as diagnostic-grade only.

## Decisions & non-blocking spec items

- **Spec typo (spec-repair, not a build increment):** primary-seed hex gloss `0x4D434E58` at `specs/rng-and-statistical-acceptance.md:104` is wrong; code correctly uses decimal 1296518744 = 0x4D474E58 and documents the discrepancy (`mcnux_stats.hxx:17-22`). Fix the spec gloss and run `specs/tools/validate_specs.sh`.
- **Cosmetic comment check:** `.h5ls` dataset-order citation at `mcnux_opacity.hxx:226-228` vs the (correct, static_assert-pinned) νe=0/ν̄e=1 slots — verify fixture order only if touching that file.
- **Packet-id retirement contract:** absorb-remove must confirm `specs/particle-container-and-gpu.md`'s concrete retirement semantics before assuming AMReX swap-and-pop suffices.
- **Escape-tally storage:** values likely land in HYD-owned grid variables while the invariant is GPU-owned ([MCNX-GPU-05]); escape-tallies must pick the surface.
- **[MCNX-RNG-06] host-vs-device bitwise leg:** documented deferral to a GPU build (`mcnux_rng.hxx:35-40`, in-code "T25d") — tracked, not an oversight; sandbox has no GPU path (GPU compile checks are CI-only by design).

## Discovered since last plan

- [ ] [gpu-annotations] Host/device annotations for the emission call chain (CI GPU builds)
  - spec: specs/particle-container-and-gpu.md
  - tests: CI GPU compile check (CUDA/HIP) stays the verifier — no sandbox GPU path
  - notes: found by emit-core. `build_tetrad`/`fluid_frame_direction` (`mcnux_tetrad.hxx`), `analytic_coefficients`/`evaluate_coefficients` (`mcnux_coefficients.hxx`), and `RangedTableCoefficients::operator()` lack host/device annotations; fine on the CPU sandbox build, but the emission count/fill kernels call them from device lambdas — a CUDA/HIP build needs `MCNUX_HOST_DEVICE` (or equivalent) on that chain.
- [ ] [emission-level-masking] Emission double-counts on refined meshes (unigrid assumption)
  - spec: specs/carpetx-thorn-integration.md
  - tests: covered by the ownership-multibox / multi-level benchmark family once masking lands
  - notes: found by emit-core. `MCNuX_Emission` walks every refinement level; with real mesh refinement, coarse cells covered by finer levels would double-emit. Same rank-local/unigrid caveat already recorded for the ledger-closure grid reduction — fix both with the same level-masking machinery.
- [ ] [ledger-emission-audit] Wire `emission_step_audit()` into `MCNuX_LedgerClosure`
  - spec: specs/hydro-coupling-source-terms.md
  - tests: extend the ledger-closure observable/golden to close against real emission events (precursor to interactions-fixedseed)
  - notes: found by emit-core. The audit accumulator is declared in `mcnux_srcterms.hxx`, defined and filled per step in `mcnux_emission.cxx` (pre-negation `LedgerDelta`s, same ΔV·Δt as the deposit), but `MCNuX_LedgerClosure` still audits only the `synthfix` synthetic list.
