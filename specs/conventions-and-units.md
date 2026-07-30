# Conventions and units (cross-cutting contract)

> Cross-cutting spec. Every other spec restates, in its own "Conventions restated" block, only the subset of these conventions it actually uses and references this file for the rest. An agent can implement MCNuX's unit handling, constant set, and species enumeration from this file alone. `README.md` is canonical if anything here conflicts with it.

## Purpose & scope

This spec pins the corpus-wide conventions every MCNuX spec and every MCNuX golden file depends on:

- Metric signature, index conventions, and the 3+1 notation all restated equations use.
- The two-unit system (geometrized for spacetime/transport, MeV + cgs for microphysics), the pinned defining-constant set, and the derived conversion factors — with every derivation shown and the derived numbers binding.
- The MeV → Kelvin conversion applied at every WeakLibInterp temperature argument.
- The neutrino species enumeration, degeneracy factors, and lepton numbers.

Out of scope:

- Random-number and statistical-acceptance conventions ([rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)).
- The per-channel WeakLibInterp argument conventions, log10 ownership, and table-range policy ([opacity-eos-evaluation](./opacity-eos-evaluation.md)); this spec pins only the *units* crossing that boundary.
- Packet state and where the degeneracy factor is applied in sampling ([packet-representation-and-sampling](./packet-representation-and-sampling.md)); this spec pins only the enumeration and the exactly-once rule.
- The MeV requirement on the GRMHD partner's `HydroBaseX::temperature` as an *interface obligation* ([hydro-coupling-source-terms](./hydro-coupling-source-terms.md)); this spec pins the unit convention it rests on.

## Source of truth

The defining constants are pinned from published standards (restated below as binding values, per the README's restate-and-pin rule):

- `c` and the SI-2019 exact constants (`k_B`, `e`): SI definition (exact by definition).
- `G`: CODATA 2018.
- `GM_sun`: IAU 2015 Resolution B3 nominal solar mass parameter (exact by convention).

Reference-stack provenance for the conventions this spec must interoperate with:

- `CarpetX/HydroBaseX/README:18-32` — the stack's geometrized `c = G = 1`, `M = M_sun` convention and its (older) constant set.
- `CarpetX/CarpetX/src/io_openpmd.cxx:162-176` — CarpetX's CODATA-2018 constant set used for openPMD unit metadata.
- `CarpetX/HydroBaseX/interface.ccl:19` — `temperature` declared with **no unit annotation** (see Open questions).
- `WeakLibInterp/specs/opacity-nes-pair.md:84` — the NES/Pair temperature axis is **Kelvin** (documented interpretation-conflict resolution; do not be misled by the MeV-labeled dead-code comment upstream of WeakLibInterp).
- `WeakLibInterp/specs/fixtures/wl-EOS-SFHo-15-25-50.h5ls`, `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls` — committed structural snapshots of the production tables (cited in place; MCNuX vendors no copies): thermodynamic axes in cgs + Kelvin, and exactly two species datasets (`Electron Neutrino`, `Electron Antineutrino`) in the EmAb group.

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — precedent for lumping the four heavy-lepton species into one effective νx species in GR-MC transport |
| arXiv:1806.02349v1 | Foucart et al. 2018, PRD 98, 063007 — cited here only for the α symbol-collision warning (its α is a spectral-fit index) |

## Inputs & outputs

This spec defines conventions, not a callable surface. The tables below are the binding values.

### Metric, indices, and 3+1 notation

- **Signature** `(−,+,+,+)`. Greek indices `μ, ν, …` run over spacetime `0–3`; Latin indices `i, j, k` over space `1–3`. Einstein summation. Partial derivatives written `∂_μ`.
- **3+1 decomposition**: lapse `α`, shift `β^i`, spatial metric `γ_ij`, with

  ```text
  ds² = −α² dt² + γ_ij (dx^i + β^i dt)(dx^j + β^j dt)
  ```

  `γ = det(γ_ij)`; the unit normal to the slices is `n^μ = (1/α, −β^i/α)`.
- **Momenta**: `p^μ` is the packet's coordinate-frame four-momentum; indices are lowered/raised with `g_μν`/`g^μν`. The fluid-frame (comoving) neutrino energy is `ν = −p_μ u^μ` with `u^μ` the fluid four-velocity.
- **Symbol collision warning — three α's.** In geometric context, unqualified `α` is the lapse. The trapped-regime relabeling fraction ([trapped-regime-treatment](./trapped-regime-treatment.md)) and the spectral-fit index (arXiv:1806.02349) also use `α`; every spec that uses a non-lapse `α` restates its meaning locally.

### The two unit systems

1. **Geometrized units** `G = c = M_sun = 1` for spacetime and transport: coordinates, metric quantities, packet positions and four-momenta, time steps, opacities as used by the transport kernels, and source terms on the grid.
2. **MeV + cgs** for microphysics: neutrino energies in MeV, temperature in **MeV** (held in MeV everywhere inside MCNuX), rest-mass density ρ in g/cm³, electron fraction Yₑ dimensionless.

Conversions between the two systems use exactly the pinned factor set below. There are no runtime unit parameters and no configurable constant set — golden data must not be configuration-dependent.

### Pinned defining constants (exact inputs of every derivation)

| Constant | Pinned value | Status |
|---|---|---|
| c | 2.99792458e10 cm/s | SI exact |
| k_B | 1.380649e-16 erg/K | SI-2019 exact |
| e | 1.602176634e-19 C (⇒ 1 MeV = 1.602176634e-6 erg, exact) | SI-2019 exact |
| G | 6.67430e-8 cm³ g⁻¹ s⁻² | CODATA 2018 |
| GM_sun | 1.3271244e26 cm³ s⁻² | IAU 2015 B3 nominal, exact by convention |

### Derived conversion factors (binding values; derivations shown)

One geometrized code unit of each quantity equals the cgs value listed. Derived values are quoted to 10 significant digits and are **binding as written**; the derivations let any agent recompute them.

| Quantity | Derivation | 1 code unit = |
|---|---|---|
| length | L = GM_sun/c² | 1.476625038e5 cm |
| time | t = L/c = GM_sun/c³ | 4.925490948e-6 s |
| mass | M_sun = GM_sun/G | 1.988409871e33 g |
| energy | M_sun c² = GM_sun c²/G | 1.787093669e54 erg = 1.115416135e60 MeV |
| mass density | M_sun/L³ | 6.175828479e17 g/cm³ |
| energy density, pressure | M_sun c²/L³ | 5.550557829e38 erg/cm³ |
| opacity (inverse length) | 1/L | 6.772199944e-6 cm⁻¹ (κ[code] = κ[cm⁻¹] × 1.476625038e5) |
| rate (inverse time) | 1/t | 2.030254467e5 s⁻¹ |
| emissivity (energy/volume/time) | M_sun c²/(L³ t) | 1.126904483e44 erg cm⁻³ s⁻¹ |

Uncertainty note: the length and time factors rest only on the exact `GM_sun` and `c`; the mass, density, pressure, and emissivity factors inherit G's ~2e-5 relative uncertainty — irrelevant at the accuracy level where these factors act. What correctness requires is that exactly this one pinned, documented factor set is used consistently throughout MCNuX.

### Temperature conversions at the WeakLibInterp boundary

- MeV → Kelvin:

  ```text
  T[K] = T[MeV] × 1.160451812e10
  ```

  Derivation (exact under the SI-2019 `k_B` and `e`): 1 MeV = 1e6 × 1.602176634e-19 J = 1.602176634e-13 J; the factor is (1.602176634e-13 J) / (1.380649e-23 J/K) = 1.1604518121550083e10 K, pinned to 10 digits as `1.160451812e10`.
- Boltzmann constant in microphysics units: `k_B = 8.617333262e-11 MeV/K` (the exact reciprocal of the factor above), for any Boltzmann factor written with T in Kelvin.
- This conversion is applied before **every** WeakLibInterp temperature argument — all channels, **including NES/Pair**, whose temperature axis is Kelvin per `WeakLibInterp/specs/opacity-nes-pair.md:84`. Assuming MeV there is a documented, easy-to-make error. Energies passed to WeakLibInterp are in MeV, ρ in g/cm³, Yₑ dimensionless; which arguments arrive pre-log10'd is owned by [opacity-eos-evaluation](./opacity-eos-evaluation.md).

### Neutrino species enumeration

Three effective species, in this binding order:

| Index | Species | Represents | Degeneracy g | Lepton number ℓ |
|---|---|---|---|---|
| 0 | νe | electron neutrino | 1 | +1 |
| 1 | ν̄e | electron antineutrino | 1 | −1 |
| 2 | νx | νμ, ν̄μ, ντ, ν̄τ (lumped) | 4 | 0 |

νx carries degeneracy factor g = 4 in energy tallies and contributes zero net lepton number (lumping precedent: arXiv:1708.08452).

## Correctness requirements

- **[MCNX-CNV-01] Notation.** All restated equations in the corpus, and all MCNuX diagnostics that expose geometric quantities, use the signature, index, and 3+1 conventions above: signature `(−,+,+,+)`, lapse `α`, shift `β^i`, spatial metric `γ_ij`, coordinate-frame `p^μ`, fluid-frame energy `ν = −p_μ u^μ`.
- **[MCNX-CNV-02] Two-unit system.** Spacetime/transport quantities are geometrized `G = c = M_sun = 1`; microphysics quantities are MeV + cgs with temperature held in MeV internally. No runtime parameter selects units or constants.
- **[MCNX-CNV-03] Pinned constant set.** Every conversion anywhere in MCNuX derives from exactly the five pinned defining constants above; no other constant set (including any carried elsewhere in the reference stack — see Open questions) is used.
- **[MCNX-CNV-04] Derived factors binding.** The derived conversion factors table is normative as written (10 significant digits). Recomputing any factor from the defining constants must reproduce the quoted value to the machine tier (relative `1e-14`); implementations may carry more digits but must agree with the pinned values to that tier.
- **[MCNX-CNV-05] Kelvin at the table boundary.** Every temperature argument of every WeakLibInterp entry point — EOS, inversion, EmAb/Iso, NES/Pair, Brem — receives `T[MeV] × 1.160451812e10` Kelvin. No WeakLibInterp channel receives MeV temperature.
- **[MCNX-CNV-06] Species enumeration.** The species set is exactly `{νe, ν̄e, νx}` with the indices, degeneracies, and lepton numbers of the table above. The factor g = 4 for νx enters **exactly once, at packet creation** (emission sampling weights, [packet-representation-and-sampling](./packet-representation-and-sampling.md)) — never in per-species transport coefficients ([opacity-eos-evaluation](./opacity-eos-evaluation.md) owns the νx coefficient mapping).

## Verification

- **Factor recomputation (machine tier, `~1e-14` relative).** For each row of the derived-factor table, recompute the factor from the five defining constants in IEEE-754 double arithmetic and compare against the pinned 10-digit value: `|recomputed − pinned| ≤ 1e-14·|pinned| + 1e-30`. Same check for `1.160451812e10` (MeV → K) and `8.617333262e-11` (k_B in MeV/K), which must also be mutual reciprocals to the machine tier.
- **Round-trip identities (machine tier).** cgs → code → cgs is the identity to `1e-14` relative for every factor; `T[MeV] → T[K] → T[MeV]` likewise.
- **Species-table cross-check.** The committed WeakLibInterp snapshots pin the production-table species axes: `WeakLibInterp/specs/fixtures/wl-Op-SFHo-15-25-50-E40-EmAb.h5ls` contains exactly the two datasets `Electron Neutrino` and `Electron Antineutrino` and no heavy-lepton dataset — consistent with νx being untabulated and handled by the mapping in [opacity-eos-evaluation](./opacity-eos-evaluation.md). Any change to that species axis invalidates this enumeration's table grounding (see Open questions).
- **Consistency sweep.** No MCNuX source file may define an independent physical-constant value: verification includes a source sweep confirming every constant literal traces to the single pinned set (implementation location is free; duplication is not).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order; contains the pinned literals `2.99792458e10`, `1.380649e-16`, `6.67430e-8`, `1.3271244e26`, `1.160451812e10`, and `g = 4`; cites reference paths that resolve under the CarpetX and WeakLibInterp roots; and declares its `[MCNX-CNV-NN]` ids in the Correctness requirements section, each covered by the [verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Where in the code conversions are applied (at ingest, at the call boundary, cached, folded into other factors), how constants are represented (constexpr, header, generated), and naming — provided every observable value is consistent with the single pinned factor set to the machine tier.
- Carrying more than 10 significant digits internally (e.g. computing `L` from the defining constants at full double precision) — the pinned 10-digit values are the agreement anchor, not a mandated truncation.
- Whether temperature is converted to Kelvin once per cell per step or per WeakLibInterp call — observable results must be as if every call received the converted value.

## Open questions / assumptions

- **Three-way constant discrepancy in the reference stack (recorded, resolved by pinning).** The stack carries three mutually inconsistent constant sets: `CarpetX/HydroBaseX/README:20-26` (CODATA-2006-era: G = 6.67428e-11 m³ kg⁻¹ s⁻², M_sun = 1.98892e30 kg), CarpetX's openPMD metadata (`CarpetX/CarpetX/src/io_openpmd.cxx:162-176`, CODATA 2018: G = 6.67430e-11, M_sun = 1.98847e30 kg), and a third, FIL-lineage set carried by downstream production GRMHD codes (not cited here — this corpus is GRMHD-agnostic and cites no production GRMHD code). MCNuX adopts none of them wholesale: the pinned set above (SI-2019 exact + CODATA-2018 G + IAU-2015 GM_sun) is normative for MCNuX. Consequence: MCNuX's mass unit (1.988409871e33 g) differs from HydroBaseX's README value in the 4th digit; the discrepancy is far below any physical effect but matters for golden data, which is why the set is pinned here.
- **`HydroBaseX::temperature` has no declared unit** (`CarpetX/HydroBaseX/interface.ccl:19`). MCNuX states as an explicit interface requirement that any GRMHD partner populates it in **MeV**; the requirement is owned by [hydro-coupling-source-terms](./hydro-coupling-source-terms.md). Assumption recorded: HydroBaseX's silence is treated as delegating the unit choice to the coupled codes.
- **Species-axis confirmation against the production tables (open question, owned by [opacity-eos-evaluation](./opacity-eos-evaluation.md)).** The enumeration above assumes the production weaklib opacity tables carry exactly the two electron-type species datasets pinned by the committed `.h5ls` snapshots. That spec carries the structured confirmation question and the νx coefficient mapping that depends on it.
- **Assumption: exactness lineage.** `GM_sun` is exact by IAU convention and `c`, `k_B`, `e` by SI definition; only `G`-derived factors carry measurement uncertainty. If CODATA revises G, MCNuX does **not** track it — the pinned set is frozen; a re-pin would be a deliberate, golden-data-invalidating corpus change.
