# Opacity and EOS evaluation (domain leaf: the WeakLibInterp consumption contract)

> Domain leaf spec. Self-contained: an agent can implement MCNuX's evaluation of weak-interaction opacities and EOS quantities through the WeakLibInterp library from this file alone, referencing `conventions-and-units.md` for the unit systems, conversion factors, and the binding species enumeration, and `rng-and-statistical-acceptance.md` only where a check reuses the pinned seeds (evaluation itself consumes no random draws). The tolerance tiers named here (machine `~1e-14`, golden/parity `1e-12`, relaxed `1e-10`, exact) are defined in `README.md` ("Global correctness contract"), which is the canonical arbiter on any conflict. WeakLibInterp itself is consumed, never modified; its own spec corpus governs its internals.

## Purpose & scope

This spec defines how MCNuX obtains microphysics from the WeakLibInterp interpolation library: the five `_Point` entry-point families and their caller-side argument conventions, the range-enforcement obligation that WeakLibInterp's permissive extrapolation policy places on MCNuX, the EOS-inversion integer error-code protocol, the table residency (read-once/upload-once) contract, the binding mapping from the corpus's three effective neutrino species to the datasets actually present in the production tables (closing the species-axis question left open in `conventions-and-units.md`), and the assembly of the transport coefficients — absorption opacity κ_a(s, E), elastic-scattering opacity κ_s(s, E), and emissivity spectrum η_s(E) — that the neutrino–matter interactions and packet-sampling specs consume.

In scope:

- The five entry-point families (EOS evaluate, EOS inversion, EmAb/Iso opacity, NES/Pair kernels, Brem kernel) with argument conventions, units, and return conventions.
- Range enforcement: the reference physical validity ranges and the pinned clamp/zero policy MCNuX applies **before** every table call.
- EOS-inversion error-code handling (codes 0/01/02/03/10/11/13, `T = 0` on failure) and the pinned consumer response.
- Table residency: read-once, upload-once, bitwise-stable evaluation for the run.
- The species-axis mapping {νe, ν̄e, νx} → table datasets, confirmed against the committed table-structure snapshots.
- The baseline transport-coefficient assembly: κ_a from EmAb, κ_s from Iso, η_s from κ_a via detailed balance (Kirchhoff), with the equilibrium chemical potentials from the EOS.

Out of scope:

- WeakLibInterp's interpolation internals, storage layout, and parity-vs-Fortran contract (owned by WeakLibInterp's own `specs/`; consumed as a black box here).
- How the assembled coefficients are used to sample events (see `neutrino-matter-interactions.md`) or emission (see `packet-representation-and-sampling.md`); the low-T `T^6` emissivity extrapolation is owned by the interactions spec and applied downstream of the assembly here.
- The trapped-regime relabeling of η, κ_a, κ_s (the trapped-regime spec; it transforms the coefficients this spec produces).
- Unit-system definitions and conversion factors (restated below only as consumed; `conventions-and-units.md` is binding).
- How grid fluid data reaches the evaluation call sites (the CarpetX-integration spec) and how evaluation is batched over cells/packets (implementation freedom).

## Source of truth

- **WeakLibInterp public headers** — the authoritative entry points and argument conventions (all compute entry points are header-inline `AMREX_GPU_HOST_DEVICE` functions over flat `double` arrays):
  - `WeakLibInterp/src/eos/wli_eos.H` — `EosInterpolateSingleVariable3DPoint` / `EosInterpolateDifferentiateSingleVariable3DPoint`: 3D (ρ, T, Yₑ) evaluate; takes **raw physical** ρ and T (the routine applies `log10` internally); Yₑ linear.
  - `WeakLibInterp/src/eos/wli_eos_inversion.H` — `ComputeTemperatureWith_{DEY,DPY,DSY}_{Guess,NoGuess}` returning `EosInversionResult{T, Error}` (and `_NoError` variants returning bare `T`): recover T from (ρ, X, Yₑ), X ∈ {E, P, S}; raw physical inputs; integer error code; `T = 0` on any failure.
  - `WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` — `EmAbInterpolateSingleVariable4DPoint` (4D: E, ρ, T, Yₑ) and `IsoInterpolateSingleVariable5DPoint` (5D: E, moment, ρ, T, Yₑ; integer moment slice): E, ρ, T **already `log10`'d by the caller**; Yₑ linear.
  - `WeakLibInterp/src/opacity/wli_opacity_nes_pair.H` — `NESPairInterpolateSingleVariable2D2DAlignedPoint` plus the symmetry-fill composers `NESDetailedBalanceFillPoint` (upper triangle × `exp((E − E′)/T)`, E and T **physical MeV** in the factor) and `PairCrossingSymmetryFillPoint` (exact index/component relabeling, no factor): T and η_e = μₑ/T already `log10`'d; energy indices integer; **T indexed in MeV** in these tables.
  - `WeakLibInterp/src/opacity/wli_opacity_brem.H` — `BremInterpolateSingleVariable2D2DAlignedSummedPoint`: ρ, T already `log10`'d; summed over 3 effective densities with weights `Alpha = [1, 1, 28/3]`.
  - `WeakLibInterp/src/core/wli_table.H` — the residency pieces: `ResidentTable<D>` (device buffer, filled once via `upload()`) and `TableView<D>` (trivially copyable pointer+extents handle captured by value in kernels).
  - `WeakLibInterp/src/io/wli_io_eos.cpp`, `WeakLibInterp/src/io/wli_io_opacity.H` — root-rank HDF5 read + MPI broadcast (`read_eos_table`, `read_emab_table`, `read_scat_iso_table`, `read_scat_nes_table`, `read_scat_pair_table`, `read_scat_brem_table`) and the committed on-disk schema constants.
- **WeakLibInterp's own specs** — the boundary policy and tolerance provenance this spec builds on: `WeakLibInterp/specs/fortran-parity-and-tolerances.md` (tolerance tiers `1e-12` parity / `1e-10` relaxed / `~1e-14` machine; the **permissive** boundary policy: bracket index clamped, delta not — out-of-range inputs extrapolate linearly, non-positive log arguments produce NaN, **no error is ever raised**), `WeakLibInterp/specs/opacity-emab-iso.md`, `WeakLibInterp/specs/opacity-nes-pair.md`, `WeakLibInterp/specs/opacity-brem.md`, `WeakLibInterp/specs/eos-interpolation.md`, `WeakLibInterp/specs/eos-inversion.md` (the error-code protocol restated below), `WeakLibInterp/specs/table-format-and-io.md`, `WeakLibInterp/specs/amrex-device-interface.md`.
- **Committed production-table structure snapshots** — the authoritative species-axis evidence (Fortran shapes in parentheses): `WeakLibInterp/specs/fixtures/wl-EOS-SFHo-15-25-50.h5ls`, `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls` (`/EmAb/Electron Neutrino` and `/EmAb/Electron Antineutrino`, each (40, 185, 81, 30), `nOpacities = 2`), `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-Iso.h5ls` (`/Scat_Iso_Kernels/Electron Neutrino` and `…/Electron Antineutrino`, each (40, 2, 185, 81, 30), `nOpacities = 2`, `nMoments = 2`), `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-NES.h5ls` and `…-Pair.h5ls` (`Kernels` (40, 40, 4, 81, 120): no species axis; 4 kernel components), `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-Brem.h5ls` (`S_sigma` (40, 40, 1, 185, 81): no species axis).
- **Cactus delivery** — `WeakLibInterp/cactus/thorns/WeakLibInterp/configuration.ccl`: the library is delivered as an ExternalLibraries-style capability; the consumer thorn's entire link obligation is `REQUIRES WeakLibInterp` in its own `configuration.ccl` (the build contract is owned by the build-and-integration spec).
- **Detailed-balance emissivity construction** — derived and pinned in this spec (the production tables store opacity, not emissivity, so the construction is MCNuX's); provenance: Kirchhoff's law for a stimulated-absorption-corrected opacity (the weaklib EmAb datasets are corrected absorption opacities — the legacy on-disk group name is `/EmAb_CorrectedAbsorption`, per `WeakLibInterp/src/io/wli_io_opacity.H`), standard neutrino-transport usage as in Bruenn 1985 (ApJS 58:771) and NuLib-based codes (Foucart 2018 [arXiv:1708.08452](https://arxiv.org/abs/1708.08452) §2).
- `conventions-and-units.md` — the caller-side `log10` table, per-channel temperature units (Kelvin for EOS/EmAb/Iso/Brem, MeV for NES/Pair), the conversion factors, and the binding species enumeration restated below.

## Inputs & outputs

### The five entry-point families (consumption summary)

| Family | Entry points (header) | Table axes | Caller passes | Returns |
|---|---|---|---|---|
| EOS evaluate | `EosInterpolateSingleVariable3DPoint` (+`Differentiate`) — `WeakLibInterp/src/eos/wli_eos.H` | 3D (ρ, T, Yₑ) | **raw physical** ρ [g/cm³], T [K]; Yₑ linear | value in the dependent variable's table units (or value + 3 partials) |
| EOS inversion | `ComputeTemperatureWith_{DEY,DPY,DSY}_{Guess,NoGuess}` — `WeakLibInterp/src/eos/wli_eos_inversion.H` | 3D, X ∈ {E, P, S} | raw physical ρ, X, Yₑ (+ optional T guess [K]) | `EosInversionResult{T [K], Error}`; `T = 0` on failure |
| EmAb / Iso opacity | `EmAbInterpolateSingleVariable4DPoint`, `IsoInterpolateSingleVariable5DPoint` — `WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` | 4D (E, ρ, T, Yₑ) / 5D (E, moment, ρ, T, Yₑ) | E [MeV], ρ [g/cm³], T [K] **already `log10`'d**; Yₑ linear; moment integer | opacity value [cm⁻¹] per species (per dataset units; see Open questions) |
| NES / Pair kernels | `NESPairInterpolateSingleVariable2D2DAlignedPoint`, `NESDetailedBalanceFillPoint`, `PairCrossingSymmetryFillPoint` — `WeakLibInterp/src/opacity/wli_opacity_nes_pair.H` | 5D (E′, E, kernel, T, η_e) | T [**MeV**], η_e = μₑ/T already `log10`'d; (iE′, iE, kernel) integer | kernel value in table units |
| Brem kernel | `BremInterpolateSingleVariable2D2DAlignedSummedPoint` — `WeakLibInterp/src/opacity/wli_opacity_brem.H` | 5D (E′, E, moment, ρ, T) | ρ (3 effective densities) [g/cm³], T [K] already `log10`'d; `Alpha = [1, 1, 28/3]` | kernel value in table units |

Consequences restated from `conventions-and-units.md` (binding at every call site): fluid temperature held in MeV is converted to Kelvin (× `1.160451812e10`) before EOS/EmAb/Iso/Brem calls but **not** before NES/Pair calls; the `log10` is applied by the caller exactly where the table says so — passing raw values to a pre-`log10` channel (or vice versa) is a correctness bug, not a tolerance issue. The electron degeneracy η_e = μₑ/T (dimensionless; μₑ from the EOS `Electron Chemical Potential` dataset, T in MeV) is computed by MCNuX before NES/Pair calls.

### Reference physical validity ranges and the clamp/zero policy (binding)

The production tables' validity ranges (from `WeakLibInterp/specs/fortran-parity-and-tolerances.md`, restated in `conventions-and-units.md`):

| Axis | Valid range |
|---|---|
| ρ | [`1.66054e3`, `3.16409e15`] g/cm³ |
| T | [`1.16045e9`, `1.83919e12`] K (equivalently [0.1, 158.5] MeV) |
| Yₑ | [0.01, 0.6] |
| E (opacity energy axes) | the loaded table's `/EnergyGrid/Values` span (production: 40 log-spaced groups) |

WeakLibInterp **never** enforces these: its boundary policy is permissive extrapolation (no error, no clamp; non-positive log arguments produce silently propagating NaN). Range enforcement is therefore **MCNuX's responsibility**, with this pinned policy applied before every opacity/EOS call:

- **ρ below range** → the cell is *interaction-free*: κ_a = κ_s = η = 0 for all species and energies; no table call is made (near-vacuum matter is transparent; extrapolated table values there are meaningless).
- **ρ, T, Yₑ, E above range; T, Yₑ, E below range** → the input is clamped to the nearest range edge before the call, and a per-run diagnostic counter of clamped evaluations is maintained (observable).
- **Non-finite input (NaN/Inf in ρ, T, Yₑ, or E)** → hard error (abort with diagnostics); NaN never reaches a `_Point` call and never propagates into transport.

(The low-temperature emissivity treatment for T < `T_low` = 0.5 MeV sits *on top* of this policy and is owned by `neutrino-matter-interactions.md`; opacities at low T use the clamped table value.)

### EOS-inversion error-code protocol (restated, binding)

The complete code set of the inversion family (provenance `WeakLibInterp/specs/eos-inversion.md`; codes are produced in this priority order — uninitialized before NaN, NaN before bounds, ρ before X before Yₑ):

| Code | Meaning |
|---|---|
| 0 | returned successfully; recovered T is valid |
| 01 | first argument ρ outside table bounds |
| 02 | second argument X (E, P, or S) outside table bounds |
| 03 | third argument Yₑ outside table bounds |
| 10 | EOS inversion not initialized |
| 11 | NaN in argument(s) |
| 13 | unable to find any root |

On any non-zero code the routine returns `T = 0`; for a `_NoError` caller, `T = 0` is the **only** failure signal. Codes 04–09 and 12 are never produced.

Pinned MCNuX response: every inversion call checks the code (or, if a `_NoError` variant is used, tests `T = 0`); a recovered T from a failing call **never** enters any physics formula. On failure, the cell uses the HydroBaseX `temperature` value if it is finite and in-range (after the clamp policy above); otherwise the cell is treated as interaction-free for the step and a per-run diagnostic counter of inversion failures is maintained. In verification fixtures (family B), an unexpected non-zero code is a hard test failure.

### Species-axis mapping (binding; closes the `conventions-and-units.md` open question)

The committed snapshots confirm the production EmAb and Iso tables carry **exactly two** species datasets — `Electron Neutrino` and `Electron Antineutrino` (`nOpacities = 2`) — and the NES/Pair/Brem kernel tables carry **no species axis** at all (kernel components only). The three-effective-species transport contract of `conventions-and-units.md` ({νe, ν̄e, νx}, νx degeneracy g = 4) therefore stands, with this binding dataset mapping:

| Species s | EmAb (κ_a) | Iso (κ_s) |
|---|---|---|
| 0 (νe) | `Electron Neutrino` dataset | `Electron Neutrino` dataset |
| 1 (ν̄e) | `Electron Antineutrino` dataset | `Electron Antineutrino` dataset |
| 2 (νx) | none — κ_a(νx, E) = 0 and η_νx(E) = 0 in the baseline (no charged-current channel for heavy leptons in these tables) | arithmetic mean of the two datasets: κ_s(νx, E) = ½ [κ_s(νe, E) + κ_s(ν̄e, E)] (νx lumps two neutrino and two antineutrino species; the datasets differ only by correction terms such as weak magnetism) |

The νx degeneracy g = 4 is folded in at packet creation per `conventions-and-units.md` and is **not** applied anywhere in this spec's per-species coefficients. NES/Pair/Brem evaluation contracts are pinned above (they are part of the five-family consumption surface), but their assembly into transport coefficients — in particular thermal pair emissivity for νx — is **deferred** (see Open questions): in the baseline interaction set, νx packets are neither emitted nor absorbed by tabulated physics, and νx transport is exercised with synthetic emissivities in tests.

### Baseline transport-coefficient assembly (normative)

Per cell (fluid state ρ, T, Yₑ after the clamp policy) and species s, on the loaded table's energy grid E ∈ {E_1 … E_B}:

**1. Absorption opacity.** κ_a(s, E) — the EmAb evaluation at (E, ρ, T, Yₑ) of the mapped dataset (zero for νx), in cm⁻¹. The stored quantity is the stimulated-absorption-corrected opacity.

**2. Elastic-scattering opacity.** κ_s(s, E) — the Iso evaluation at (E, moment = 0, ρ, T, Yₑ) of the mapped dataset(s) (zeroth Legendre moment = the transport-relevant total elastic opacity; the first-moment slice is not used in the baseline), in cm⁻¹.

**3. Emissivity spectrum via detailed balance (Kirchhoff).** Because κ_a is the corrected opacity, the fluid-frame spectral energy emissivity is

```text
η_s(E) = c κ_a(s, E) · (4π E³ / (hc)³) · f_eq(E; T, μ_ν,s) ,
f_eq(E; T, μ_ν,s) = 1 / ( exp( (E − μ_ν,s)/T ) + 1 )
```

in `MeV cm⁻³ s⁻¹ MeV⁻¹`, with E, T, μ in MeV and the binding constants c = `2.99792458e10` cm/s (defining, `conventions-and-units.md`) and hc = `1.23984198e-10` MeV·cm (exact by SI definition — h, c, and the MeV are exact; the derivation discipline of MCNX-CNV-01 applies to it). The equilibrium neutrino chemical potentials, from the EOS dependent variables `Electron Chemical Potential`, `Proton Chemical Potential`, `Neutron Chemical Potential`:

```text
μ_ν,0 = μ_e + μ_p − μ_n   (νe),    μ_ν,1 = −μ_ν,0   (ν̄e),    μ_ν,2 = 0   (νx)
```

**4. Group emissivities.** The per-group emissivities the packet sampler consumes (`packet-representation-and-sampling.md`) are η_{s,b} = ∫ over group b of η_s(E) dE, with group edges taken from the loaded table's energy-grid geometry; the quadrature rule within a group is implementation freedom (a one-point group-center × width rule is acceptable) but must be deterministic and fixed for a run.

### Table residency (observable contract)

- Each required table (EOS + the opacity channels the configuration enables) is read **once per run** — root-rank HDF5 read, then MPI broadcast (the WeakLibInterp `read_*` functions provide exactly this) — and uploaded **once** to device memory (`ResidentTable<D>::upload`, one host-to-device copy); kernels receive the trivially copyable `TableView<D>` by value.
- No table file I/O occurs after initialization; no per-step re-upload occurs absent a table-changing event (there is none in the current parameter surface).
- Evaluation is a pure function of its inputs for the whole run: the same (channel, species, E, ρ, T, Yₑ) query returns the bitwise-identical `double` at every step, on every rank, matching the single-rank CPU reproducibility scheme of `rng-and-statistical-acceptance.md`.

## Correctness requirements

- **[MCNX-OPA-01] Argument-convention parity (golden `1e-12`).** For every consumed channel, MCNuX's evaluation path (unit conversion → clamp policy → `log10` where required → `_Point` call → assembly) reproduces a direct reference call of the corresponding WeakLibInterp `_Point` function on identical physical inputs at `rtol = 1e-12`, `atol = 1e-30` (WeakLibInterp's own parity tier). In particular: T converted to Kelvin for EOS/EmAb/Iso/Brem and left in MeV for NES/Pair; E, ρ, T pre-`log10`'d exactly where the family table above says so and nowhere else; Yₑ passed linear everywhere; the Iso moment and NES/Pair energy/kernel indices passed as integers. Reference: the family table above; provenance the cited WeakLibInterp headers.
- **[MCNX-OPA-02] Range enforcement (exact).** Zero out-of-range and zero non-finite values ever reach a `_Point` call: below-range ρ yields the interaction-free result with no call; other out-of-range inputs are clamped to the range edge (a clamped call's result equals the edge evaluation bitwise); non-finite inputs abort. The clamped-evaluation and inversion-failure diagnostic counters are observable. Reference: the clamp/zero policy above; provenance the permissive boundary policy of `WeakLibInterp/specs/fortran-parity-and-tolerances.md`.
- **[MCNX-OPA-03] Inversion error protocol (exact).** MCNuX's inversion handling produces/consumes exactly the code set {0, 01, 02, 03, 10, 11, 13}; a `T` from a call with non-zero code (equivalently `T = 0` from a `_NoError` variant) never enters any physics formula; the pinned fallback (HydroBaseX temperature if valid, else interaction-free + counter) is taken. Reference: the protocol table above; provenance `WeakLibInterp/specs/eos-inversion.md`.
- **[MCNX-OPA-04] Species-axis mapping (exact + machine).** The dataset mapping table is honored exactly: νe/ν̄e read their named datasets; κ_a(νx) = 0 and η_νx = 0 identically in the baseline; κ_s(νx, E) equals the arithmetic mean of the two datasets' values at `rtol = 1e-14`; no factor of g = 4 appears in any per-species coefficient (it enters once, at νx packet creation, per `conventions-and-units.md`). Reference: the mapping table above; evidence the committed `.h5ls` snapshots (`nOpacities = 2`).
- **[MCNX-OPA-05] Detailed-balance emissivity identity (machine).** At every evaluated (s, E, cell): η_s(E) / (c κ_a(s, E)) equals `(4π E³/(hc)³) f_eq(E; T, μ_ν,s)` at `rtol = 1e-14` (with the μ_ν,s construction above), wherever κ_a > 0; and μ_ν,1 = −μ_ν,0 exactly. This is the algebraic assembly identity; its physical consequence (radiation equilibrates to f_eq) is exercised by the interactions spec's thermalization benchmark. Reference: equation 3 above.
- **[MCNX-OPA-06] Residency and evaluation stability (exact).** One read and one device upload per table per run (instrumented); a fixed probe set of queries evaluated at step 0 and at a later step k returns bitwise-identical results; the same probe set evaluated on 1 and 2 ranks returns bitwise-identical results (tables are broadcast, evaluation is pure). Reference: the residency contract above.

## Verification

- **MCNX-OPA-01**: a family-B unit check builds small **synthetic** tables (known analytic content, e.g. exactly log-linear dependence) through the WeakLibInterp host structures, then evaluates a fixed probe lattice (≥ `1e3` points per channel spanning the interior) through both the MCNuX path and direct `_Point` reference calls, asserting `rtol = 1e-12`, `atol = 1e-30` per point. A production-table leg (same probes on `wl-EOS-SFHo-15-25-50.h5` and companions, gated by an environment variable pointing at the table set) reports SKIP loudly when tables are absent — never a silent pass.
- **MCNX-OPA-02**: the same harness with probe points outside every range edge (each axis, both sides) and NaN probes: asserts the interaction-free result for below-range ρ, bitwise edge-equality for clamped axes, the abort on non-finite input, and the counter values (exact integers).
- **MCNX-OPA-03**: a unit check drives the inversion family across in-bounds, each out-of-bounds axis, NaN, and a no-root synthetic table, asserting the exact code per case, `T = 0` on every failure, and the pinned fallback selection logic (exact).
- **MCNX-OPA-04**: a unit check on synthetic two-dataset tables with distinct known content: asserts νe/ν̄e dataset selection (values distinguishable by construction), κ_a(νx) ≡ 0, η_νx ≡ 0, and the νx Iso mean at `rtol = 1e-14`.
- **MCNX-OPA-05**: a unit check sweeps (s, E, T, μ) over a probe lattice and asserts the identity at `rtol = 1e-14`, including the μ_ν,1 = −μ_ν,0 sign.
- **MCNX-OPA-06**: an instrumented 2-rank fixture run asserts read/upload counts and the bitwise stability comparisons (exact).
- All of the above appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `WeakLibInterp/` paths (headers, specs, and `.h5ls` fixtures), and contains the required claim strings (the caller-side `log10` convention, all seven inversion error codes as table cells, the validity ranges, the residency claims).

## Implementation freedom

- **Evaluation granularity and caching** — per-cell precomputation of (κ_a, κ_s, η) spectra once per transport step vs. on-demand per-packet evaluation, batching over cells/species/groups, and any host/device split — provided every observable value satisfies MCNX-OPA-01/05 and the residency contract.
- The in-group quadrature rule for η_{s,b} (one-point, multi-point) — provided it is deterministic, fixed for a run, and the interactions spec's equilibration benchmark passes.
- Which optional channels (NES/Pair/Brem) are loaded when the configuration does not consume them, and how table paths are parameterized.
- How the clamp policy is implemented (branch, min/max, masked) and how the diagnostic counters are represented — provided counts are exact and observable.
- Whether inversion is used at all (a configuration that trusts HydroBaseX `temperature` directly need never invert) — the protocol of MCNX-OPA-03 binds only actual inversion calls.
- Data structures wrapping `TableView`/`ResidentTable`, and any prefetching or memoization — provided bitwise evaluation stability (MCNX-OPA-06) holds.

## Open questions / assumptions

- **Dataset units are taken as cm⁻¹ (assumption, confirm at first load).** The EmAb and Iso zeroth-moment datasets are assumed to store opacities in cm⁻¹ (their `/…/Units` metadata datasets are authoritative). The implementation must log the Units strings at first table load; if they differ from cm⁻¹, only the single conversion constant applied at the assembly boundary changes — every requirement above is invariant. (WeakLibInterp itself defers units to table metadata.)
- **Chemical-potential rest-mass convention (assumption, confirm at first load).** μ_ν,0 = μ_e + μ_p − μ_n assumes the table's `Proton Chemical Potential` and `Neutron Chemical Potential` are stored in a consistent rest-mass convention such that the combination is the physical νe equilibrium chemical potential. This must be confirmed once against the loaded table (sanity check: in hot NSE matter at Yₑ ≈ 0.1, μ_ν,0 should be large and negative-to-moderate in the expected sense; a rest-mass-offset error shifts it by ≈ 1.293 MeV, the neutron–proton mass difference). If the convention differs, the μ_ν,0 formula gains the documented constant offset here — a spec change, not silent code drift.
- **νx pair-process emissivity is deferred (the corpus's largest open physics item).** With no charged-current channel, baseline νx packets are never emitted or absorbed from tabulated physics; production νx emission requires assembling thermal pair emissivity/absorptivity from the NES/Pair/Brem kernels (their evaluation contracts are already pinned above). That assembly, its blocking-factor treatment, and its verification benchmarks are a future spec change here; until then νx transport machinery is exercised with synthetic emissivities in tests.
- **νx Iso mean is an approximation (assumption, non-blocking).** The ½(νe + ν̄e) elastic opacity for νx averages the weak-magnetism asymmetry over the lumped 2ν + 2ν̄ content; the residual error is far below the corpus's transport fidelity. If a production table set ever provides four species datasets, the mapping table gains explicit νx/ν̄x rows (the three-effective-species transport contract stands regardless).
- **First-moment Iso slice unused (assumption).** The baseline elastic-scattering treatment (isotropic re-emission in the fluid frame) uses only the zeroth moment; anisotropic elastic scattering using the first moment would be a spec change in `neutrino-matter-interactions.md` consuming the moment = 1 slice already reachable through the pinned Iso entry point.
