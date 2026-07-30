# MCNuX specifications

This `specs/` set is the durable contract for **MCNuX**, a GPU Monte Carlo neutrino transport code implemented as a Cactus thorn on the CarpetX driver, with transport packets carried by the AMReX particle module and weak-interaction opacities and EOS quantities evaluated through the WeakLibInterp library. It is written to drive a **Ralph loop**: a long-running loop in which a *fresh* agent re-reads the specs on every iteration with no memory of prior runs. Each spec pins down *what "correct" means and how to verify it* while leaving *how to implement it* open.

Open any single spec file in isolation: it is self-contained. The two cross-cutting specs ([conventions-and-units](./conventions-and-units.md), [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md)) carry the shared unit/notation and RNG/statistical contracts; every leaf spec restates only the conventions it uses and references those two for the rest. **This README is the canonical source if any spec ever conflicts with it.**

## Global correctness contract (inherited by every spec)

- **Restated equations are normative (equation provenance).** Every governing equation is written out in full in the spec that uses it, in the corpus's single notation (fixed by [conventions-and-units](./conventions-and-units.md)), with its citation — bare arXiv id plus equation number — carried as provenance only. On any conflict between spec text and paper, **the spec text is normative**. Cite-only equations (paper + number, agent reads the paper) are forbidden. Standing assumption, carried by every physics spec: equations were extracted from the ar5iv renderings of the pinned arXiv versions and the restated forms are normative; a spot-check of a typeset PDF that finds a discrepancy updates the spec, not the citation.
- **Citation form.** Inline citations are bare arXiv ids with equation numbers (e.g. `arXiv:2103.16588 Eq. 48`) — never author-year labels ("Foucart 2018 Eq. 16" is genuinely ambiguous between two registered papers that both have a meaningful Eq. 16). Each spec's Source of truth section carries a registry table of exactly the papers that spec cites, pinning the arXiv **version** once per paper; the master registry below covers every published source cited anywhere in the corpus. Every equation number in a spec refers to that spec's pinned version. Non-arXiv sources are registered by DOI.
- **Interpretation-conflict precedence.** If a leaf spec conflicts with this README, the README governs. If spec text conflicts with a cited paper, the spec text governs (above). When *reference sources* disagree among themselves on interpretation — units, conventions, live-vs-dead code paths — the resolution is pinned in the affected spec and the conflict is logged in that spec's "Open questions / assumptions" section (worked example: the WeakLibInterp NES/Pair temperature axis is **Kelvin**, per `WeakLibInterp/specs/opacity-nes-pair.md`, overriding a misleading MeV comment upstream).
- **Tolerance tiers for deterministic checks** (`|got − expected| ≤ max(abstol, reltol·|expected|)` unless a spec states otherwise):
  - **Machine tier** — `~1e-14` (a few ULP) — closed-form exactness checks: recomputing pinned conversion factors from the defining constants, node identities, exact algebraic invariants.
  - **Golden parity** — `1e-12` — fixed-seed, single-rank deterministic benchmark outputs against committed golden data; this matches the Cactus regression harness defaults `ABSTOL=RELTOL=1e-12` (`flesh/lib/sbin/RunTestUtils.pl`).
  - **Relaxed** — `1e-10` — checks that legitimately accumulate more floating-point error (long integrations, interpolation-of-interpolation); a spec may invoke this tier only with a stated rationale.
  - **Exact** — bitwise / integer / NaN-equality — RNG known-answer vectors, per-packet event sequences at fixed seed on a single rank, enumerations, error codes.
- **Statistical acceptance bar: 4σ at pinned packet counts and pinned seeds.** Every Monte Carlo statistical check compares an estimator against its expected value with tolerance four standard errors, evaluated at the packet count and RNG seed pinned by [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) (default packet count `N_p = 1048576 = 2^20`, primary seed `1296518744`). That spec is the global statistical contract every stochastic domain spec cites.
- **Conservation tolerances.** The exchange ledger must balance identically up to floating-point reduction: for each conserved channel (energy, momentum, lepton number), the sum of what packets lose per step equals what the source-term grid variables gain, to relative `1e-13` per step. Analytic conservation benchmarks (e.g. Schwarzschild `p_t` conservation) carry explicit numeric tolerances pinned per benchmark in [verification-suite-design](./verification-suite-design.md).
- **Units conventions.** Two unit systems, pinned in [conventions-and-units](./conventions-and-units.md): geometrized `G = c = M_sun = 1` for spacetime and transport; MeV + cgs for microphysics, with temperature held in **MeV** internally and converted to **Kelvin** (× `1.160451812e10`, exact) at every WeakLibInterp temperature argument. Exactly one pinned, documented conversion-factor set is used throughout MCNuX; there are no runtime unit parameters.
- **Boundary policy.** A packet whose position leaves the outer domain boundary is removed from transport and accumulated into escape tallies (per-species energy and number); no default reflection or periodicity. WeakLibInterp table lookups extrapolate permissively out of range, so **table-range enforcement is MCNuX's responsibility**; the enforcement policy is pinned in [opacity-eos-evaluation](./opacity-eos-evaluation.md). After each transport step every packet resides on the rank/level/grid/tile owning its position ([particle-container-and-gpu](./particle-container-and-gpu.md)).
- **Reproducibility.** All physics random draws come from the pinned stateless counter-based generator of [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md); per-packet event sequences are exactly reproducible at fixed seed in single-rank CPU test mode. There is **no** bitwise cross-run/cross-rank/cross-device requirement: grid tallies deposited by unordered atomics are reproducible only to floating-point-reduction tolerance.

## Master paper registry

Every published source cited anywhere in this corpus, version-pinned. Per-spec registries must agree with this table. Authors/year/title are an informational gloss only; the pinned id is the citation.

| Pinned id | Informational gloss |
|---|---|
| arXiv:1708.08452v2 | Foucart 2018, MNRAS 475, 4186 — the GR-MC methods paper; its Eq. 16 is the MC stress-energy estimator |
| arXiv:1806.02349v1 | Foucart, Duez, Kidder, Nguyen, Pfeiffer, Scheel 2018, PRD 98, 063007 — closure-evaluation companion; its Eq. 16 is a spectral fit, *not* the estimator (registered chiefly to disambiguate) |
| arXiv:2008.08089v3 | Foucart et al. 2020, ApJL 902, L27 — the letter introducing the trapped-regime relabeling, with the **flipped** α convention (its α is this corpus's 1−α) |
| arXiv:2103.16588v2 | Foucart et al. 2021, ApJ 920, 82 — SpEC MC implementation paper; the **binding** α convention for the opacity relabeling and its stability bound |
| arXiv:1903.09273v1 | Miller, Ryan, Dolence 2019 — νbhlight |
| arXiv:1505.05119v1 | Ryan, Dolence, Gammie 2015, ApJ 807, 31 — bhlight |
| arXiv:1507.03606v2 | Richers et al. 2015, ApJ 813, 38 — Sedonu |
| arXiv:1706.06187v2 | Richers et al. 2017, ApJ 847, 133 — MC vs discrete-ordinates comparison |
| 10.1145/2063384.2063405 | Salmon, Moraes, Dror, Shaw 2011, SC'11 — "Parallel random numbers: as easy as 1, 2, 3" (Random123; defines Philox4x32-10) |

## Spec index

The corpus is exactly this README plus the 12 coded specs below — nothing else. Each coded spec owns the `[MCNX-XXX-NN]` requirement ids bearing its code: every id *declared* (in a `## Correctness requirements` section) uses the declaring spec's own code, so an id's prefix always tells a fresh agent which spec owns the requirement.

| Code | Spec | Description |
|---|---|---|
| CNV | [conventions-and-units](./conventions-and-units.md) | Cross-cutting: metric signature and index conventions, the two-unit system with pinned defining constants and derived conversion factors, MeV↔Kelvin at the WeakLibInterp boundary, the three-species enumeration with g = 4 for νx. |
| RNG | [rng-and-statistical-acceptance](./rng-and-statistical-acceptance.md) | Cross-cutting: Philox4x32-10 pinned (counter/key packing, [0,1) double construction, known-answer vectors), single-rank exact reproducibility, and the global 4σ statistical acceptance bar at pinned packet counts and seeds. |
| PKT | [packet-representation-and-sampling](./packet-representation-and-sampling.md) | Domain leaf: packet state (position, coordinate-frame four-momentum, weight, species, RNG identity), emission sampling of energies, angles, and species. |
| GEO | [geodesic-propagation](./geodesic-propagation.md) | Domain leaf: geodesic push on the dynamical 3+1 metric gathered from ADMBaseX grid functions. |
| INT | [neutrino-matter-interactions](./neutrino-matter-interactions.md) | Domain leaf: absorption and elastic scattering — optical-depth event sampling and outcomes. |
| TRP | [trapped-regime-treatment](./trapped-regime-treatment.md) | Domain leaf: optically-thick scheme — explicit baseline plus the opacity-relabeling extension with its stability bound (binding 2021 α convention). |
| OPA | [opacity-eos-evaluation](./opacity-eos-evaluation.md) | Domain leaf: the WeakLibInterp evaluation contract across all five entry-point families, the baseline EmAb + Iso coefficient assembly, and the νx dataset mapping. |
| HYD | [hydro-coupling-source-terms](./hydro-coupling-source-terms.md) | Domain leaf: energy/momentum/lepton-number source-term grid variables, the zero-then-add handoff protocol to a GRMHD partner, and the optional Tmunu contribution. |
| GPU | [particle-container-and-gpu](./particle-container-and-gpu.md) | Technical leaf: pure-SoA packet container, device-kernel execution model, redistribution invariants. |
| CTX | [carpetx-thorn-integration](./carpetx-thorn-integration.md) | Technical leaf: thorn/schedule/parameter/grid-access contracts against the CarpetX driver and base thorns. |
| VER | [verification-suite-design](./verification-suite-design.md) | Technical leaf: benchmark and acceptance-test design, the `TestMCNuX` sibling test thorn, and the coverage matrix closing over every declared `[MCNX-*]` id. |
| BLD | [build-and-integration](./build-and-integration.md) | Delivery: build and dependency contract, thorn assembly, the pinned test ThornList, golden-data location. |

The validator (`tools/validate_specs.sh`) enforces that this README links **exactly** these files — no orphan spec files on disk, no missing links, no broken links — and that the [verification-suite-design](./verification-suite-design.md) coverage matrix references every `[MCNX-XXX-NN]` id declared anywhere (the closure check). Coverage grows by registry append in `tools/validate_specs.sh`.

## Required structure for every spec (the section contract)

Every coded spec contains these seven `## ` sections, in this order:

```text
## Purpose & scope → ## Source of truth → ## Inputs & outputs →
## Correctness requirements → ## Verification → ## Implementation freedom →
## Open questions / assumptions
```

and instantiates this skeleton (this README is the one file with a different anatomy):

````markdown
# <Spec title>

> Leaf spec. Self-contained: an agent can implement <X> from this file alone. It restates
> only the conventions it uses and references <cross-cutting specs> for the rest.
> `README.md` is canonical if any restated convention here conflicts with it.

## Purpose & scope
<covers / does NOT cover, naming the sibling spec that owns each exclusion>

## Source of truth
<prose + repo file:line references>
| Pinned id | Informational gloss |
|---|---|
| arXiv:2103.16588v2 | Foucart et al. 2021, "Implementation of Monte-Carlo transport…" |

## Inputs & outputs
## Correctness requirements
- [MCNX-XXX-NN] <requirement, with restated equation + provenance>
### Conventions restated (the subset this leaf uses)
## Verification
### Mechanical (validator)
<what the validator checks about this file>
## Implementation freedom
## Open questions / assumptions
````

Cross-cutting and technical specs use the same skeleton with tier-appropriate blockquote wording. Requirement ids have the fixed form `[MCNX-XXX-NN]` (`XXX` = the owning spec's code from the index table, `NN` = two digits); they are declared only inside the owner's `## Correctness requirements` section, and every declared id appears in the [verification-suite-design](./verification-suite-design.md) coverage matrix.

## Validating the spec set

`bash tools/validate_specs.sh` runs the corpus linter: 7-section order per spec; every cited reference-repo path resolves; per-spec required-string claims; citation-registry closure (inline ids ⊂ per-spec registry ⊂ master registry, versions agreeing); `[MCNX-XXX-NN]` ownership, declaration, and coverage closure; README link and exact-set closure; and structural claims anchored in WeakLibInterp's committed table snapshots (`WeakLibInterp/specs/fixtures/tables.provenance`, cited in place — MCNuX vendors no fixture copies). Exit 0 iff every check passes.

Cited repository paths are written with a repo-name prefix and resolve against five reference roots, each overridable by environment variable (a stale override self-heals via a per-repo probe file; defaults assume the sibling-checkout layout):

| Prefix | Env override | Probe file |
|---|---|---|
| `amrex/` | `MCNX_AMREX_ROOT` | `amrex/Src/Particle/AMReX_ParticleTile.H` |
| `warpx/` | `MCNX_WARPX_ROOT` | `warpx/Source/Particles/WarpXParticleContainer.H` |
| `CarpetX/` | `MCNX_CARPETX_ROOT` | `CarpetX/CarpetX/src/driver.hxx` |
| `flesh/` | `MCNX_FLESH_ROOT` | `flesh/lib/sbin/RunTestUtils.pl` |
| `WeakLibInterp/` | `MCNX_WEAKLIB_INTERP_ROOT` | `WeakLibInterp/specs/fixtures/tables.provenance` |

These five roots cover **every** repository path cited anywhere in the corpus — in particular, no production GRMHD code is cited anywhere; all coupling contracts are stated against the CarpetX base thorns plus explicit interface requirements on whatever GRMHD partner couples to MCNuX.
