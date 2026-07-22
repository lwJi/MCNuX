# Conventions and units (cross-cutting notation/units contract)

> Cross-cutting spec. Every other MCNuX spec references this file for the shared notation, unit systems, conversion constants, species enumeration, and packet nomenclature, and restates only the subset it actually uses. The corpus-wide numeric tolerance tiers (machine `~1e-14`, golden `1e-12`, relaxed `1e-10`, exact), the 4σ statistical acceptance bar, and the conservation tolerances are defined in `README.md` ("Global correctness contract"); every tier named in this file resolves there. If a leaf spec and this file ever disagree, `README.md` is the canonical arbiter; absent that, this file governs notation and units.

## Purpose & scope

This spec defines the conventions every MCNuX spec (and the eventual implementation) inherits: the metric signature and index conventions, the 3+1 variable naming and its mapping to the CarpetX base thorns (ADMBaseX, HydroBaseX, TmunuBaseX) with grid centerings, the two unit systems (geometrized for spacetime/transport, cgs+MeV for microphysics) with their binding defining constants and derived conversion factors, the WeakLibInterp caller-side argument conventions (which inputs arrive already `log10`'d, per-channel temperature units), the binding neutrino species enumeration, and the packet component nomenclature all other specs use.

In scope:

- Metric signature, index conventions, and 3+1 variable symbols mapped to concrete CarpetX grid variables and centerings.
- The defining physical constants (binding), the derived geometrized↔cgs/MeV conversion factors, and where unit conversion happens.
- The WeakLibInterp call-boundary conventions: which arguments the caller must pre-`log10`, per-channel temperature units, reference table ranges.
- The binding species enumeration {νe, ν̄e, νx} with degeneracy and lepton-number assignments.
- Packet component names, symbols, and units (nomenclature reservation).

Out of scope:

- Packet component types, valid ranges, and storage layout (owned by the packet-representation and particle-container specs).
- The physics equations that use this notation (owned by the domain leaf specs).
- WeakLibInterp's boundary/extrapolation policy and the range-enforcement obligation it implies (owned by the opacity/EOS evaluation spec; the unit facts needed to honor it are stated here).
- The RNG and statistical-acceptance contract (see `rng-and-statistical-acceptance.md`).

## Source of truth

- **CarpetX base-thorn variable declarations** — the authoritative names, groupings, and centerings of the metric, fluid, and stress-energy grid variables:
  - `CarpetX/ADMBaseX/interface.ccl` — groups `metric` (gxx, gxy, gxz, gyy, gyz, gzz), `curv` (kxx…kzz), `lapse` (alp), `shift` (betax, betay, betaz), and their time derivatives; all vertex-centered (`VVV`, the CarpetX default when `CENTERING` is omitted).
  - `CarpetX/HydroBaseX/interface.ccl` — cell-centered (`CCC`) scalars rho, eps, press, temperature, entropy, Ye and 3-vectors vel, Bvec.
  - `CarpetX/TmunuBaseX/interface.ccl` — vertex-centered eTtt, eTti, eTij.
- **WeakLibInterp call-boundary conventions** — the authoritative argument conventions of the interpolation entry points MCNuX consumes:
  - `WeakLibInterp/src/eos/wli_eos.H` — EOS evaluate takes **raw physical** ρ and T (the routine applies `log10` internally); Yₑ linear.
  - `WeakLibInterp/src/eos/wli_eos_inversion.H` — EOS inversion recovers T from (ρ, X, Yₑ) with X one of {specific internal energy E, pressure P, entropy S}; raw physical inputs; integer error codes on failure.
  - `WeakLibInterp/src/opacity/wli_opacity_emab_iso.H` — EmAb/Iso opacity take E, ρ, T **already `log10`'d by the caller**; Yₑ linear.
  - `WeakLibInterp/src/opacity/wli_opacity_nes_pair.H` — NES/Pair kernels take T and η = μₑ/(k_B T) already `log10`'d, where μₑ is the electron chemical potential (an EOS dependent variable) and η the resulting dimensionless electron degeneracy parameter; **T indexed in MeV** in these tables.
  - `WeakLibInterp/src/opacity/wli_opacity_brem.H` — Brem kernel takes ρ and T already `log10`'d.
  - `WeakLibInterp/specs/fortran-parity-and-tolerances.md` — the reference physical table ranges restated below.
- **Defining physical constants** (external provenance, restated as binding values below): the exact SI defining constants (CODATA/SI-2019: c, k_B, e), the IAU 2015 Resolution B3 nominal solar mass parameter (GM)_☉, and the CODATA 2018 Newtonian constant G.
- **Species-set provenance**: the three-effective-species set {νe, ν̄e, νx} with heavy-lepton degeneracy 4 follows Miller, Ryan & Dolence 2019 (nubhlight), [arXiv:1903.09273](https://arxiv.org/abs/1903.09273); the enumeration below is normative regardless.

## Inputs & outputs

This spec does not define a callable surface; it defines the conventions all MCNuX surfaces obey.

### Metric signature and index conventions

- Metric signature **(−,+,+,+)**.
- Greek indices μ, ν, … run over spacetime components 0–3; Latin indices i, j, k, … run over spatial components 1–3. Einstein summation applies.
- Component 0 is the coordinate time direction t. Covariant components carry lower indices (p_μ), contravariant upper (p^μ); indices are raised/lowered with the four-metric g_μν.
- The 3+1 (ADM) line element is `ds² = −α² dt² + γ_ij (dx^i + β^i dt)(dx^j + β^j dt)`, with lapse α, shift β^i, and spatial metric γ_ij; γ^{ij} is the inverse of γ_ij (not the spatial part of g^{μν} raised with the four-metric); K_ij is the extrinsic curvature. √−g = α√γ with γ = det(γ_ij).

### 3+1 variables mapped to CarpetX grid variables

| Quantity | Symbol | CarpetX variable(s) | Thorn / group | Centering |
|---|---|---|---|---|
| Lapse | α | `alp` | ADMBaseX::lapse | vertex (VVV) |
| Shift | β^i | `betax, betay, betaz` | ADMBaseX::shift | vertex |
| Spatial metric | γ_ij | `gxx, gxy, gxz, gyy, gyz, gzz` | ADMBaseX::metric | vertex |
| Extrinsic curvature | K_ij | `kxx, kxy, kxz, kyy, kyz, kzz` | ADMBaseX::curv | vertex |
| Rest-mass density | ρ | `rho` | HydroBaseX::rho | cell (CCC) |
| Fluid velocity | v^i | `velx, vely, velz` | HydroBaseX::vel | cell |
| Specific internal energy | ε_int | `eps` | HydroBaseX::eps | cell |
| Pressure | P | `press` | HydroBaseX::press | cell |
| Temperature | T | `temperature` | HydroBaseX::temperature | cell |
| Electron fraction | Yₑ | `Ye` | HydroBaseX::Ye | cell |
| Stress-energy | T_μν | `eTtt`, `eTti…`, `eTij…` | TmunuBaseX | vertex |

MCNuX-owned source-term grid variables are **cell-centered** (matching HydroBaseX's `CCC`); the metric is vertex-centered, so any code path evaluating both at the same point owes an explicit centering bridge (the obligation is pinned in the CarpetX-integration spec; the centering facts here are the binding inputs to it).

### Unit systems and defining constants

Two unit systems, with a fixed boundary between them:

1. **Geometrized (code) units** — `G = c = M_sun = 1` — for spacetime geometry, grid coordinates, coordinate time, and packet transport.
2. **Microphysics units** — ρ in `g/cm^3`, temperatures and neutrino energies in `MeV` (with ħ = c = 1 for neutrino kinematics) — the units in which MCNuX *holds* fluid and radiation quantities for all EOS/opacity/rate work. Individual table axes may demand other units at the call itself (notably Kelvin temperature axes; see the per-channel table below, which governs).

The **defining constants** (binding; the first three are exact by definition, the fourth exact-by-convention, the fifth the CODATA 2018 measured value):

| Constant | Binding value | Status |
|---|---|---|
| Speed of light c | `2.99792458e10 cm/s` | exact (SI) |
| 1 MeV in erg | `1.602176634e-6 erg` | exact (SI) |
| Boltzmann constant k_B | `8.617333262e-11 MeV/K` | exact (SI-derived) |
| Solar mass parameter (GM)_☉ | `1.3271244e26 cm^3/s^2` | exact by convention (IAU 2015 B3 nominal) |
| Newtonian constant G | `6.67430e-8 cm^3 g^-1 s^-2` | CODATA 2018 |

The **derived conversion factors** (normative rule: compute them from the defining constants; the decimals quoted here are informative, given to exactly 7 significant figures — except the temperature factor, which is an exact ratio of two exact constants and is quoted to 10 significant figures):

| Code unit | Expression | Value |
|---|---|---|
| Length | (GM)_☉ / c² | `1.476625e5 cm` |
| Time | (GM)_☉ / c³ | `4.925491e-6 s` |
| Mass | (GM)_☉ / G | `1.988410e33 g` |
| Density | mass / length³ | `6.175828e17 g/cm^3` |
| Energy | mass · c² | `1.787094e54 erg` = `1.115416e60 MeV` |
| Temperature (K per MeV) | 1 / k_B | `1.160451812e10 K/MeV` (exact) |

Unit conversion happens **at the call boundary**: geometrized quantities are converted to microphysics units immediately before a WeakLibInterp table call or a rate evaluation, and rates/exchanges are converted back to geometrized units before deposition onto the grid. No spec or implementation may introduce a third unit system for any externally observable quantity.

### WeakLibInterp call-boundary conventions

The caller-side argument conventions, per channel (the header files cited under Source of truth are authoritative):

| Channel family | Table axes | Caller passes | Temperature unit in table |
|---|---|---|---|
| EOS evaluate | (ρ, T, Yₑ) | **raw physical** ρ, T; Yₑ linear | Kelvin |
| EOS inversion | (ρ, X, Yₑ), X ∈ {E, P, S} (here E = specific internal energy, not neutrino energy) | raw physical; integer error codes on failure | Kelvin |
| EmAb / Iso opacity | (E, ρ, T, Yₑ) [+ moment] | E, ρ, T **already `log10`'d**; Yₑ linear | Kelvin |
| NES / Pair kernels | (E′, E, kernel, T, η) | T, η already `log10`'d; energy indices integer | **MeV** |
| Brem kernel | (E′, E, moment, ρ, T) | ρ, T already `log10`'d | Kelvin |

Axis symbols used above: **E, E′** — neutrino energies, always in **MeV** on every opacity-table axis (the caller applies `log10` to the MeV value where the table above says so; on the NES/Pair/Brem axes the pre-collapsed energy indices are integers). **X** — the EOS dependent variable being inverted for temperature: specific internal energy E, pressure P, or entropy S. **moment** — the integer Legendre-moment (angular-expansion) index of a scattering kernel, used directly as a table index and never interpolated; its cardinality is a property of the loaded table (2 moments in the production Iso tables, 1 in Brem; NES/Pair carry 4 kernel components addressed the same way). **η = μₑ/(k_B T)** — the dimensionless electron degeneracy parameter, with μₑ the electron chemical potential supplied by the EOS.

Consequences every consumer spec inherits:

- Fluid temperature held in MeV must be converted to Kelvin (multiply by `1.160451812e10`) before EOS/EmAb/Iso/Brem calls, but **not** before NES/Pair calls.
- The `log10` is applied by the caller exactly where the table above says so; passing raw values to a pre-`log10` channel (or vice versa) is a correctness bug, not a tolerance issue.
- Reference physical validity ranges of the production tables (from `WeakLibInterp/specs/fortran-parity-and-tolerances.md`): ρ ∈ [`1.66054e3`, `3.16409e15`] g/cm³; T ∈ [`1.16045e9`, `1.83919e12`] K (equivalently [0.1, 158.5] MeV); Yₑ ∈ [0.01, 0.6]. WeakLibInterp extrapolates permissively outside these ranges (no error, no clamp), so range enforcement is MCNuX's responsibility — the enforcement requirement is owned by the opacity/EOS evaluation spec; the ranges here are the binding numbers it uses.

### Neutrino species enumeration (binding)

Exactly three effective species, with fixed integer labels:

| Index | Species | Symbol | Degeneracy g | Lepton number ℓ |
|---|---|---|---|---|
| 0 | electron neutrino | νe | 1 | +1 |
| 1 | electron antineutrino | ν̄e | 1 | −1 |
| 2 | heavy-lepton (effective) | νx | 4 | 0 |

νx represents νμ, ν̄μ, ντ, and ν̄τ collectively. The degeneracy factor `g = 4` is folded in **once, at packet creation**: heavy-lepton emissivities are multiplied by g when νx packets are sampled, so a νx packet's particle-count weight N already counts physical neutrinos summed over all four heavy-lepton species. No downstream tally, source term, or diagnostic may apply any additional factor of g (that would double-count). Lepton-number accounting uses ℓ from this table: νe carries +1, ν̄e carries −1, νx carries 0 net lepton number.

### Packet component nomenclature (names reserved for all specs)

| Symbol | Name | Meaning | Units |
|---|---|---|---|
| x^i | position | coordinate position of the packet | code (geometrized) length |
| p_t, p_i | four-momentum | covariant coordinate-frame components of the packet four-momentum p_μ | MeV |
| ε | fluid-frame energy | ε = −p_μ u^μ with u^μ the fluid four-velocity (derived, not stored independently) | MeV |
| N | weight | number of physical neutrinos the packet represents (g already folded in for νx) | dimensionless |
| s | species | integer species label per the enumeration above | — |
| (q, e) | RNG identity | unique 64-bit packet id q and event counter e (see `rng-and-statistical-acceptance.md`) | — |

Carrying p_μ in MeV against code-unit coordinates is consistent because the geodesic equation is invariant under a constant rescaling of p_μ; the conversion factor `1 MeV = 1/1.115416e60` code energy units enters only where packet four-momenta source geometrized grid quantities (stress-energy and source-term deposition). Types, valid ranges, and storage layout for these components are owned by the packet-representation and particle-container specs; the names, index placement (covariant p_μ), and units above are binding corpus-wide.

## Correctness requirements

- **[MCNX-CNV-01] Constants are derived, not transcribed.** Every physical constant and unit-conversion factor used by the implementation must be computed from the five defining constants in this file (or be one of them verbatim) and must agree with the value so derived to the machine-precision tier (`rtol ~1e-14`, defined in `README.md`). The informative decimals quoted in the derived-factor table are given to exactly 7 significant figures (the temperature factor to 10, being an exact ratio of exact constants) and must agree with the implementation's derived values to `1e-6` relative. Reference: the defining-constants table above; provenance SI-2019 exact constants, IAU 2015 B3, CODATA 2018.
- **[MCNX-CNV-02] Species enumeration is exact.** The implementation's species labels, degeneracy factors, and lepton numbers must match the enumeration table exactly (integer equality): indices {0, 1, 2} ↔ {νe, ν̄e, νx}, g = {1, 1, 4}, ℓ = {+1, −1, 0}; `g = 4` applied exactly once, at νx packet creation, and never again downstream. Reference: the binding enumeration above; provenance nubhlight [arXiv:1903.09273](https://arxiv.org/abs/1903.09273).
- **Notation is single-valued.** Any spec or implementation artifact using a symbol from this file (α, β^i, γ_ij, K_ij, p_μ, ε, N, s, q, e) must use it with the meaning, index placement, and units fixed here. Signature is (−,+,+,+) everywhere; no spec may locally flip it.
- **Unit conversion at the boundary.** Conversions between geometrized and microphysics units happen at the WeakLibInterp call boundary and the grid-deposition boundary as stated above; the temperature-unit facts (Kelvin for EOS/EmAb/Iso/Brem, MeV for NES/Pair) and the caller-side `log10` table are hard requirements on every call site.

## Verification

- **MCNX-CNV-01**: a unit check recomputes every conversion factor the implementation exposes from the five defining constants and asserts agreement at `rtol = 1e-14`, `atol = 1e-30`; it additionally asserts the quoted decimals of this file against the derived values at `rtol = 1e-6`. Runs host-only, no tables, no grid.
- **MCNX-CNV-02**: a unit check asserts the species table (indices, g, ℓ) by integer equality against the implementation's species metadata, and asserts on a synthetic emission/tally round trip that the νx energy tally carries exactly one factor of g = 4 relative to the single-species rate (deterministic bookkeeping form, relaxed tier `1e-10`).
- Both checks appear as rows of the coverage matrix in `verification-suite-design.md` (closure enforced by `specs/tools/validate_specs.sh`).
- Mechanical: `bash specs/tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order, names concrete numeric tolerances, resolves its cited `CarpetX/` and `WeakLibInterp/` paths, and contains the required claim strings (species degeneracy, temperature conversion, `log10` convention).

## Implementation freedom

- The internal unit system for *storage* (an implementation may hold packet energies geometrized internally, or fluid data in cgs) — provided every externally observable quantity and every WeakLibInterp call honors the boundary conventions above at the stated tolerances.
- Constant naming, code organization, and whether conversion factors are precomputed or evaluated on the fly.
- The numeric representation of species labels in memory (enum, int component, …) — provided the observable integer values match the enumeration.
- Any additional *internal* diagnostics in any units, provided they are not part of a specified interface or golden output.

## Open questions / assumptions

- **HydroBaseX variable units (assumption, to confirm at integration).** `CarpetX/HydroBaseX/interface.ccl` declares names and centerings but not units; producing GRMHD codes conventionally supply `rho` in geometrized code units (multiply by `6.175828e17` to get g/cm³), `temperature` in MeV, and `Ye` dimensionless. MCNuX assumes exactly that. If the actual producer differs, only the conversion applied at the HydroBaseX read boundary changes; the microphysics-side units in this file are unaffected. Confirm against the producing thorn's documentation when coupling is first exercised.
- **Species axis of the production opacity tables (deferred to the opacity/EOS spec).** The enumeration {νe, ν̄e, νx} must be confirmed against the species axis actually present in the production weaklib opacity tables (pinned by the committed snapshots under WeakLibInterp's `specs/fixtures/` directory) when the opacity/EOS evaluation spec is drafted; if the tables distinguish four species (νx split), the enumeration here gains an explicit table-to-species mapping but the three-effective-species transport contract stands.
- **G's uncertainty (non-blocking).** The geometrized↔cgs *mass* conversion inherits the CODATA uncertainty of G (relative `2.2e-5`); this is irrelevant to internal consistency because MCNX-CNV-01 pins all factors to the same binding G value, and it only matters if outputs are compared against literature values computed with a different G.
