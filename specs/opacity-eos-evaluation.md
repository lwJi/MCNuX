# Opacity and EOS evaluation (WeakLibInterp consumption contract)

> Leaf spec. Self-contained: an agent can implement MCNuX's opacity/EOS evaluation layer
> — the complete WeakLibInterp consumption surface, the baseline transport-coefficient
> assembly, and the νx dataset mapping — from this file alone. It restates only the
> conventions it uses and references [conventions-and-units](./conventions-and-units.md)
> (units, constants, species) and
> [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (statistical
> acceptance) for the rest. `README.md` is canonical if any restated convention here
> conflicts with it.

## Purpose & scope

This spec pins how MCNuX evaluates weak-interaction opacities and EOS quantities through
the WeakLibInterp library:

- The evaluation contract for **all five** WeakLibInterp entry-point families — EOS
  interpolation, EOS inversion, EmAb/Iso, NES/Pair, Brem — covering argument
  conventions, units, log10 ownership, and index conventions, so the consumption surface
  is complete and durable even where the baseline does not yet consume a family.
- The **baseline transport-coefficient assembly**: per-species absorption opacity
  `κ_a(s, E)`, elastic scattering opacity `κ_s(s, E)`, and spectral emissivity
  `η(s, E)` from the EmAb and Iso tables only (Kirchhoff closure for η).
- The **νx dataset mapping** against the two-species production tables (exact tier).
- The **table-range enforcement policy** (WeakLibInterp extrapolates permissively;
  range enforcement is MCNuX's responsibility, per `README.md`).

Out of scope:

- How the coefficients are consumed — event sampling and outcomes
  ([neutrino-matter-interactions](./neutrino-matter-interactions.md)); the
  trapped-regime relabeling η′/κ_a′/κ_s′
  ([trapped-regime-treatment](./trapped-regime-treatment.md)).
- Where the degeneracy factor g = 4 is applied — **exactly once, at packet creation**
  ([packet-representation-and-sampling](./packet-representation-and-sampling.md));
  it appears **never in per-species transport coefficients** (this spec owns that
  prohibition's coefficient side).
- Unit definitions and pinned constants ([conventions-and-units](./conventions-and-units.md));
  this spec applies them at the call boundary.
- Device-kernel execution idioms and table residency mechanics
  ([particle-container-and-gpu](./particle-container-and-gpu.md)); grid access and
  scheduling ([carpetx-thorn-integration](./carpetx-thorn-integration.md)).
- The analytic (gray) opacity mode's benchmark values and the concrete verification
  problems ([verification-suite-design](./verification-suite-design.md)); this spec pins
  only that the coefficient interface below is source-agnostic (table or analytic).

## Source of truth

WeakLibInterp is consumed, never modified. Its own spec corpus is the governing contract
for kernel semantics; on interpretation conflicts its precedence rule (*table > thornado
> weaklib*, `WeakLibInterp/specs/README.md`) applies. Reference artifacts:

- `WeakLibInterp/src/eos/wli_eos.H` — EOS 3D point evaluate/differentiate
  (`EosInterpolateSingleVariable3DPoint`); the library log10s ρ and T internally, Yₑ
  linear.
- `WeakLibInterp/src/eos/wli_eos_inversion.H` — temperature recovery
  (`ComputeTemperatureWith_{DEY,DPY,DSY}_{NoGuess,Guess}` →
  `EosInversionResult{T, Error}`; any error ⇒ `T = 0`).
- `WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` — EmAb 4D and Iso 5D point
  kernels (`EmAbInterpolateSingleVariable4DPoint`, Iso + `IsoOffset` 2D offset
  selector); all pre-log10'd coordinates.
- `WeakLibInterp/src/opacity/wli_opacity_nes_pair.H` — the channel-neutral NES/Pair
  bilinear kernel (`NESPairInterpolateSingleVariable2D2DAlignedPoint`) plus the
  `NESDetailedBalanceFillPoint` / `PairCrossingSymmetryFillPoint` symmetry fills.
- `WeakLibInterp/src/opacity/wli_opacity_brem.H` — Brem per-density kernel and
  `BremInterpolateSingleVariable2D2DAlignedSummedPoint` over 3 effective densities with
  weights `Alpha = [1, 1, 28/3]`.
- `WeakLibInterp/src/core/wli_interp.H`, `WeakLibInterp/src/core/wli_index.H`,
  `WeakLibInterp/src/core/wli_table.H` — bracket/delta (clamped index, unclamped
  delta ⇒ edge-cell extrapolation), column-major `flat_index` (first index fastest),
  universal recovery `10**stored − OS`, and the trivially-copyable `TableView`.
- `WeakLibInterp/specs/opacity-nes-pair.md` — **the NES/Pair temperature axis is
  Kelvin** (documented interpretation-conflict resolution; the MeV-labeled dead-code
  comment upstream in weaklib is overridden).
- `WeakLibInterp/test/test_parallelfor_wrappers.cpp` — the device consumption pattern
  (host-built `TableView` captured by value into one `amrex::ParallelFor`; bit-identity
  between host loop and device launcher).
- Committed structural snapshots (cited in place — MCNuX vendors no fixture copies),
  pinned by `WeakLibInterp/specs/fixtures/tables.provenance`:
  `WeakLibInterp/specs/fixtures/wl-EOS-SFHo-15-25-50.h5ls`,
  `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls`,
  `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-Iso.h5ls`,
  `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-NES.h5ls`,
  `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-Pair.h5ls`,
  `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-Brem.h5ls` — structural
  snapshots of the production tables `wl-EOS-SFHo-15-25-50.h5`,
  `wl-Op-SFHo-15-25-50-E40-EmAb.h5`, `wl-Op-SFHo-15-25-50-E40-Iso.h5`,
  `wl-Op-SFHo-15-25-50-E40-NES.h5`, `wl-Op-SFHo-15-25-50-E40-Pair.h5`,
  `wl-Op-SFHo-15-25-50-E40-Brem.h5`.
- Thornado production precedent for the νx mapping (informational provenance recorded
  from the design review; thornado is deliberately **outside** this corpus's five
  validator-checked reference repositories, so these pointers are not validator-resolved
  citations): heavy-lepton species receive exactly zero tabulated EmAb opacity (the
  `iS <= iNuE_Bar` guard with `opEC = Zero` otherwise, thornado
  `Modules/Opacities/NeutrinoOpacitiesComputationModule.F90`, ~lines 887–904); the two
  Iso datasets are reused via a chirality map keyed on lepton-number sign
  (`iC = nChirals - (LeptonNumber(iS)+1)/nChirals`, thornado
  `Modules/TwoMoment/TwoMoment_NeutrinoMatterSolverModule.F90`, ~line 3419); equilibrium
  distributions use μ_ν = 0 for heavy species (same opacities module, ~lines 195–199).

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — defines the roles of η, κ_a, κ_s consumed by the MC event sampling (its Eq. 18 emission and Eq. 24 interaction-time draws) |

## Inputs & outputs

### Coefficient interface (what the rest of MCNuX sees)

The observable product of this spec is the per-species transport-coefficient triple,
evaluated at a fluid state and fluid-frame neutrino energy:

| Symbol | Meaning | Units (microphysics system) |
|---|---|---|
| κ_a(s, E) | absorption opacity, species s, fluid-frame energy E | cm⁻¹ |
| κ_s(s, E) | elastic (isoenergetic) scattering opacity | cm⁻¹ |
| η(s, E) | spectral energy emissivity per unit volume (per single species, g = 1) | MeV cm⁻³ s⁻¹ MeV⁻¹ |

Inputs: fluid state `(ρ, T, Yₑ)` with ρ in g/cm³, **T in MeV** (MCNuX-internal), Yₑ
dimensionless; species `s ∈ {νe, ν̄e, νx}` (indices 0, 1, 2); energy E in MeV. The
interface is **source-agnostic**: a runtime parameter selects the production weaklib
tables (this spec) or the analytic per-species η/κ_a/κ_s formulas (a first-class MCNuX
capability; formulas and benchmarks pinned in
[verification-suite-design](./verification-suite-design.md)). Everything below the
interface in this spec concerns the table source. Conversion of κ and η to geometrized
code units uses exactly the pinned factors of
[conventions-and-units](./conventions-and-units.md).

### The five WeakLibInterp entry-point families (complete consumption surface)

All kernels are scalar, allocation-free, `double`-typed pure functions over a
`TableView` (flat column-major storage, **first index fastest**; `h5ls` prints the
reverse, C-order shape). Everything is stored `log10(value + OS)` and recovered
`10**stored − OS`; the offset `OS` is per-dataset (1D `Offsets` for EOS/EmAb, 2D
`Offsets[nOpacities, nMoments]` for every scattering kernel).

| Family | Entry points | Coordinate arguments (in order) | Units | log10 owner |
|---|---|---|---|---|
| EOS evaluate | `EosInterpolateSingleVariable3DPoint` | (ρ, T, Yₑ) **raw** | g/cm³, **Kelvin**, — | **library** (logs ρ, T internally; Yₑ linear) |
| EOS inversion | `ComputeTemperatureWith_{DEY,DPY,DSY}_{NoGuess,Guess}` | (ρ, X, Yₑ) raw, X = ε/P/s | g/cm³, X-dependent, — | library |
| EmAb | `EmAbInterpolateSingleVariable4DPoint` | (LogE, LogD, LogT, Yₑ) | log10 of MeV, g/cm³, **Kelvin**; Yₑ raw | **caller (MCNuX)** |
| Iso | Iso point kernel + `IsoOffset(iOp, iMom)` | (LogE, LogD, LogT, Yₑ) + moment index `iMom` | as EmAb; `iMom` raw int | caller |
| NES/Pair | `NESPairInterpolateSingleVariable2D2DAlignedPoint` (+ `NESDetailedBalanceFillPoint`, `PairCrossingSymmetryFillPoint`) | (LogT, LogEta) + raw ints (iEp, iE, kernel) | log10 of **Kelvin**, dimensionless η_e = μ_e/(k_B T) | caller |
| Brem | `BremInterpolateSingleVariable2D2DAlignedSummedPoint` | (LogD[3] effective densities, LogT) + raw ints (iEp, iE) | log10 of g/cm³, **Kelvin** | caller |

Index conventions (binding): all integer indices are zero-based at this interface;
energy indices `iEp, iE` address the 40-point tabulated energy grid directly (the 2D
kernels interpolate only the thermodynamic plane — energies and moment index are
discrete and **never interpolated**); Brem's thermodynamic plane is (ρ, T) with ρ
**before** T (the transpose of NES/Pair's (T, η_e) plane); Brem sums 3 effective
densities with fixed weights `Alpha = [1, 1, 28/3]` applied to **recovered** (not
log-space) values.

Fortran-order table shapes, pinned by the committed snapshots: EOS
`(nD, nT, nY) = (185, 81, 30)` per dependent variable; EmAb `(nE, nD, nT, nY) =
(40, 185, 81, 30)` per species; Iso `(nE, nMom, nD, nT, nY) = (40, 2, 185, 81, 30)` per
species; NES/Pair `(nEp, nE, nMom, nT, nEta) = (40, 40, 4, 81, 120)`; Brem
`(nEp, nE, nMom, nD, nT) = (40, 40, 1, 185, 81)`.

### Boundary unit conversions (applied by MCNuX at every call)

- Temperature: `T[K] = T[MeV] × 1.160451812e10` (exact factor, pinned in
  [conventions-and-units](./conventions-and-units.md)) before **every** WeakLibInterp
  temperature argument — all channels, **including NES/Pair**, whose temperature axis is
  Kelvin per `WeakLibInterp/specs/opacity-nes-pair.md`. Assuming MeV there is a
  documented, easy-to-make error.
- Electron degeneracy for NES/Pair: `η_e = μ_e / T` with both in MeV (equivalently
  `μ_e/(k_B T)` with T in Kelvin and `k_B = 8.617333262e-11 MeV/K`); `μ_e` read from the
  EOS table's `Electron Chemical Potential` dataset (MeV).
- Energies in MeV, ρ in g/cm³, Yₑ dimensionless. Log10 is taken of the value in exactly
  these units.

## Correctness requirements

- **[MCNX-OPA-01] Complete, exact consumption surface.** MCNuX consumes the
  WeakLibInterp `_Point` kernels of the five families exactly as tabulated above:
  argument order, zero-based integer indices, discrete (never-interpolated) energy and
  moment axes, per-dataset offsets, and Brem's `Alpha = [1, 1, 28/3]`
  recovered-value weighting. MCNuX calls the WeakLibInterp kernels; it does **not**
  reimplement the interpolation. For identical inputs and table data, every value MCNuX
  uses equals the corresponding WeakLibInterp kernel output bitwise (exact tier).

- **[MCNX-OPA-02] Kelvin at every temperature argument.** Every temperature passed to
  any WeakLibInterp entry point — EOS evaluate, inversion, EmAb, Iso, NES/Pair, Brem —
  is `T[MeV] × 1.160451812e10` Kelvin. No channel receives MeV temperature.

- **[MCNX-OPA-03] Log10 ownership.** The EOS evaluate/inversion families receive raw
  (un-logged) ρ, T, X, Yₑ; the library owns their log-space handling. All four opacity
  families receive coordinates **already log10'd by MCNuX** in the units of this spec
  (Yₑ always raw/linear, never logged). MCNuX must guarantee no non-positive value
  reaches a log10 (WeakLibInterp propagates silent NaN); the range policy of
  [MCNX-OPA-06] is what guarantees it.

- **[MCNX-OPA-04] Baseline assembly (EmAb + Iso only).** The baseline transport
  coefficients for the electron-type species are assembled from exactly two families:

  ```text
  κ_a(s, E)  =  EmAb(s; E, ρ, T, Yₑ)                        s ∈ {νe, ν̄e}     [cm⁻¹]
  κ_s(s, E)  =  Iso(s; E, ρ, T, Yₑ; iMom = 0)               s ∈ {νe, ν̄e}     [cm⁻¹]
  η(s, E)    =  c · κ_a(s, E) · f_eq(E; s) · 4π E³/(hc)³                       [MeV cm⁻³ s⁻¹ MeV⁻¹]
  ```

  with the Kirchhoff closure (restated, normative; thornado-precedent form) using the
  equilibrium Fermi–Dirac occupancy

  ```text
  f_eq(E; s) = 1 / ( exp( (E − μ_ν(s)) / T ) + 1 )          (E, μ_ν, T all in MeV)
  μ_ν(νe) = μ_e + μ_p − μ_n ,   μ_ν(ν̄e) = −μ_ν(νe) ,   μ_ν(νx) = 0
  ```

  where μ_e, μ_p, μ_n are the EOS table's `Electron Chemical Potential`,
  `Proton Chemical Potential`, and `Neutron Chemical Potential` datasets (MeV) evaluated
  at the same (ρ, T, Yₑ), and the constants are `c = 2.99792458e10 cm/s` and
  `hc = 1.239841984e-10 MeV cm` (see Open questions for the hc pin). The Iso moment
  index is `iMom = 0` (the l = 0 Legendre coefficient — the total elastic rate,
  matching the isotropic outgoing redraw of
  [neutrino-matter-interactions](./neutrino-matter-interactions.md)); the l = 1 slice
  is not consumed by the baseline. The bin-integrated emissivity consumed by emission
  sampling is `η_b(s) = η(s, E_b) ΔE_b` with `E_b` the bin center and `ΔE_b` the bin
  width of the tabulated energy grid. NES, Pair, and Brem are **not** assembled into
  baseline rates (structured open question below).

- **[MCNX-OPA-05] νx mapping (exact tier).** The production EmAb/Iso tables carry
  exactly two species datasets (`Electron Neutrino`, `Electron Antineutrino`); νx is
  untabulated. The binding mapping is:

  ```text
  κ_a(νx, E) = 0        (from tabulated physics)
  η(νx, E) = 0          (from tabulated physics)
  κ_s(νx, E) = ½ [κ_s(νe, E) + κ_s(ν̄e, E)]
  ```

  The ½-mean is the exact lump-average of the production chirality map (two Iso
  datasets keyed on lepton-number sign) over the 2ν + 2ν̄ heavy species that νx
  aggregates. The degeneracy factor g = 4 enters **exactly once, at packet creation**
  ([packet-representation-and-sampling](./packet-representation-and-sampling.md)) and
  **never in per-species transport coefficients** — the mapping above is per-species
  (g = 1) by construction. Tabulated νx emission is identically zero until the
  NES/Pair/Brem open question is taken up; νx transport is exercised in tests via the
  analytic opacity mode ([verification-suite-design](./verification-suite-design.md)).

- **[MCNX-OPA-06] Table-range enforcement (MCNuX's responsibility).** WeakLibInterp
  clamps only the bracket index, never the interpolation delta — out-of-range queries
  extrapolate linearly from the edge cell, silently. MCNuX therefore enforces ranges
  before every evaluation, with this binding policy:

  - **Transparency floor:** if `ρ < ρ_min(table)` (the table's lowest tabulated
    density), the cell is treated as transparent for the table source:
    `κ_a = κ_s = η = 0` for all species and energies. No table evaluation occurs.
  - **Clamping elsewhere:** for `ρ ≥ ρ_min`, the query coordinates (ρ, T, Yₑ, E, η_e as
    applicable) are clamped into the tabulated axis ranges (component-wise, in the same
    log/linear space the axis is stored in) before evaluation. No evaluated coordinate
    ever lies outside its table axis range; no NaN from a WeakLibInterp kernel may
    propagate into transport.
  - **Diagnostics:** the implementation maintains per-axis clamp counters (a
    diagnostic, not a physics output) so silent saturation is observable.

- **[MCNX-OPA-07] EOS-inversion error protocol.** Any use of the inversion family must
  consume the full `EosInversionResult{T, Error}`: `Error ≠ 0` implies the returned
  `T = 0` is not a temperature and must not be used; in test configurations an
  inversion error on the baseline path is a hard failure. (The baseline transport path
  does not require inversion — temperature arrives in MeV from the GRMHD partner per
  [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) — but the contract is
  pinned so any future consumer inherits it.)

- **[MCNX-OPA-08] Species-axis grounding.** The species enumeration's table grounding
  is the committed snapshots: `wl-Op-SFHo-15-25-50-E40-EmAb.h5ls` and
  `wl-Op-SFHo-15-25-50-E40-Iso.h5ls` list exactly the two electron-type datasets and no
  heavy-lepton dataset; NES/Pair carry a single species-less `Kernels` dataset
  (nMoments = 4) plus the 120-point `EtaGrid`; Brem carries the single `S_sigma`
  dataset. If a future table revision changes any species axis, [MCNX-OPA-05] and the
  species enumeration of [conventions-and-units](./conventions-and-units.md) must be
  re-confirmed before that table is adopted (see Open questions).

### Conventions restated (the subset this leaf uses)

- Microphysics units: E and T in MeV internally, ρ in g/cm³, Yₑ dimensionless;
  temperature converted MeV → Kelvin (`× 1.160451812e10`, exact) at every WeakLibInterp
  temperature argument; `k_B = 8.617333262e-11 MeV/K`.
- Species `{νe, ν̄e, νx}` with indices (0, 1, 2), degeneracies (1, 1, 4), lepton numbers
  (+1, −1, 0); g = 4 never appears in this spec's outputs.
- Opacity/emissivity → geometrized code-unit conversion uses exactly the pinned factor
  set of [conventions-and-units](./conventions-and-units.md) (κ[code] = κ[cm⁻¹] ×
  1.476625038e5).
- All reals are IEEE-754 binary64 (`double`), matching WeakLibInterp's pinned value
  type.

## Verification

- **Kernel-parity (exact tier).** For a pinned set of in-range sample points spanning
  each family, the value MCNuX uses equals the directly-invoked WeakLibInterp `_Point`
  kernel output bitwise. Host and device evaluations agree bitwise for identical inputs
  (the WeakLibInterp `test_parallelfor_wrappers.cpp` pattern).
- **Unit-boundary checks (machine tier, ~1e-14 relative).** The MeV→Kelvin conversion
  at each call site reproduces `T[MeV] × 1.160451812e10`; `η_e` computed both as
  `μ_e[MeV]/T[MeV]` and `μ_e/(k_B T[K])` agrees to machine tier.
- **νx mapping identities (exact / machine tier).** `κ_a(νx, E) = 0` and `η(νx, E) = 0`
  exactly; `κ_s(νx, E)` equals the ½-mean of the two evaluated electron-type values to
  machine tier (`1e-14` relative) for every sample point.
- **Kirchhoff consistency (machine tier + statistical).** At any in-range state,
  `η(s, E) / κ_a(s, E) = c f_eq(E; s) 4πE³/(hc)³` holds to machine tier by
  construction; the equilibration-box benchmark
  ([verification-suite-design](./verification-suite-design.md)) verifies that the
  assembled coefficients drive the packet spectrum to `f_eq` within the 4σ bar of
  [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md).
- **Range-policy unit checks (exact).** Queries with `ρ < ρ_min` return exactly zero
  coefficients and touch no table; queries just outside any other axis bound return
  exactly the value at the clamped coordinate; no NaN is produced for any input in the
  policy's domain; clamp counters increment as specified.
- **Inversion protocol (exact).** For pinned failing inputs, `Error ≠ 0` and `T = 0`
  are surfaced, and the baseline test-mode behavior (hard failure) triggers.

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in
order; contains the required literals `× 1.160451812e10`, `including NES/Pair`,
`κ_a(νx, E) = 0`, `η(νx, E) = 0`, `κ_s(νx, E) = ½ [κ_s(νe, E) + κ_s(ν̄e, E)]`, `28/3`,
and `never in per-species transport coefficients`; that its cited `WeakLibInterp/` and
`CarpetX/` paths resolve; that the structural claims above (two electron-type species
datasets in EmAb/Iso, the EOS chemical-potential dataset, the NES/Pair `EtaGrid` and
`Kernels` shapes, Brem's `S_sigma`) are anchored verbatim in the committed `.h5ls`
snapshots under the resolved WeakLibInterp root; and that its `[MCNX-OPA-NN]` ids are
declared here and covered by the [verification-suite-design](./verification-suite-design.md)
matrix.

## Implementation freedom

- Table residency and upload strategy (when/where `ResidentTable`/`TableView` uploads
  happen, arena choice), caching of assembled coefficients (per cell, per bin, per
  step, or none), batching and kernel fusion — provided observable values equal the
  pure-function contract of [MCNX-OPA-01].
- Where the MeV→Kelvin conversion and the log10s are computed (per call, per cell,
  hoisted) — observable results must be as if applied at every call.
- The internal representation of the coefficient triple and of the source-selection
  parameter (names, enum layout), and how the analytic mode is dispatched.
- Precomputing `f_eq`/chemical-potential lookups vs evaluating on the fly.
- How clamp diagnostics are aggregated and reported.

## Open questions / assumptions

- **NES/Pair/Brem rate assembly (structured open question, deliberately not resolved
  here).** The consumption contract for the NES, Pair, and Brem kernels is pinned above
  so the surface is complete, but their assembly into inelastic scattering / pair /
  bremsstrahlung **rates** — and hence νx thermal emission — is not part of the
  baseline (the baseline event set is absorption + elastic scattering,
  [neutrino-matter-interactions](./neutrino-matter-interactions.md)). A future
  inelastic extension adds the assembly equations and their verification problems to
  this spec and to [verification-suite-design](./verification-suite-design.md) without
  changing the entry-point contract. Until then, tabulated νx emission is identically
  zero and νx coverage runs through the analytic opacity mode.
- **Species-axis confirmation against the production tables (open question, carried
  per the species-enumeration decision).** The νx mapping assumes the production
  EmAb/Iso tables carry exactly the two electron-type species datasets pinned by the
  committed `.h5ls` snapshots (`nOpacities = 2`). The validator re-greps the snapshots
  on every run; a table revision with a different species axis must fail loudly here
  before adoption.
- **Assumption: hc pin extends the constant set.** The Kirchhoff normalization needs
  the Planck constant, which the pinned defining-constant set of
  [conventions-and-units](./conventions-and-units.md) does not carry. This spec pins
  `h = 6.62607015e-34 J s` (SI-2019 exact) ⇒ `hc = 1.239841984e-10 MeV cm` (10
  digits, derived with the pinned `e` and `c`; recomputation must agree to machine
  tier). This is an SI-exact extension, not a competing constant set.
- **Assumption: stimulated-absorption convention.** The tabulated EmAb value is
  consumed as the absorption opacity `κ_a` satisfying the Kirchhoff closure
  `η = c κ_a f_eq 4πE³/(hc)³` exactly (thornado's production convention: emissivity is
  computed from the tabulated opacity times the equilibrium occupancy). Whether the
  table's value includes a stimulated-absorption correction is thereby absorbed into
  the closure; the equilibration-box benchmark is the arbiter — if it fails its 4σ
  bar systematically, this assumption is the first suspect and its resolution updates
  this spec.
- **Assumption: chemical-potential convention.** `μ_ν(νe) = μ_e + μ_p − μ_n` uses the
  table's chemical-potential datasets in the table's own rest-mass convention
  (thornado-precedent combination). Any rest-mass-offset inconsistency between the
  datasets would surface in the same equilibration benchmark.
- **Assumption: transparency floor.** Treating `ρ < ρ_min(table)` as vacuum
  (coefficients exactly zero) mirrors the production consumer's density guard and
  avoids edge-cell extrapolation into near-vacuum, where log-space extrapolation is
  least trustworthy. Benchmarks that need opacity in low-density regions use the
  analytic mode.
