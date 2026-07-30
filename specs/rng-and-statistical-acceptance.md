# RNG and statistical acceptance (cross-cutting contract)

> Cross-cutting spec. Every stochastic domain spec cites this file for random draws and statistical pass/fail criteria; leaf specs restate only the draw conventions they use. An agent can implement the MCNuX random-number generator and the statistical acceptance harness from this file alone. `README.md` is canonical if anything here conflicts with it.

## Purpose & scope

This spec pins:

- The random-number generator, as a **pure function**: algorithm (Philox4x32-10), constants, round function, the exact packing of (seed, packet id, event counter, draw index) into counter/key words, and the exact output-word → [0, 1) IEEE-754 double construction.
- Known-answer test (KAT) vectors as exact-equality correctness fixtures.
- The reproducibility contract: single-rank exact reproducibility of per-packet event sequences; explicitly **no** bitwise cross-run/cross-rank/cross-device requirement.
- The global statistical acceptance bar (4σ) with its pinned packet counts and pinned seeds.
- The prohibition on AMReX's stateful device RNG for physics draws, with rationale.

Out of scope:

- *Which* physics events consume draws, and each event's draw order and count — every consuming domain spec ([packet-representation-and-sampling](./packet-representation-and-sampling.md), [neutrino-matter-interactions](./neutrino-matter-interactions.md), [trapped-regime-treatment](./trapped-regime-treatment.md)) documents its own (event, draw) consumption against the primitive defined here.
- Concrete statistical benchmark problems and their expected values ([verification-suite-design](./verification-suite-design.md)).
- The packet id's storage and lifetime in the particle container ([particle-container-and-gpu](./particle-container-and-gpu.md)); this spec consumes the id as an opaque 64-bit value.

## Source of truth

The generator is **pinned, not left free**. Rationale, recorded per the design decision: fixed-seed golden data is a function of the generator, so leaving the family free would let a fresh iteration legitimately swap generators and invalidate every golden file; and no counter-based generator exists anywhere in the reference stack to serve as an exemplar (AMReX's device RNG is a stateful per-thread pool; νbhlight uses per-thread Mersenne Twister). The algorithm below is restated in full from the Random123 publication and its reference implementation (`philox.h` and the published `kat_vectors` file of the Random123 distribution); the restated constants, round function, and KAT vectors are normative as written.

The **forbidden** alternative, cited as provenance for the prohibition: AMReX's stateful per-thread device RNG — `amrex/Src/Base/AMReX_Random.H` (device-callable `Random(RandomEngine const&)` overloads over persistent `curand`/`hiprand` per-thread states) and `amrex/Src/Base/AMReX_RandomEngine.H`. Its values depend on which hardware thread executes which packet, i.e. on execution order — exactly what this contract excludes for physics draws.

| Pinned id | Informational gloss |
|---|---|
| 10.1145/2063384.2063405 | Salmon, Moraes, Dror, Shaw 2011, SC'11 — "Parallel random numbers: as easy as 1, 2, 3" (Random123; defines Philox4x32-10 and publishes its KAT vectors) |

## Inputs & outputs

The MCNuX RNG is the pure function

```text
u : (S, q, e, k) → r ∈ [0, 1)
```

| Argument | Type / range | Meaning |
|---|---|---|
| S | uint64 | run seed (runtime parameter; pinned values below for golden/statistical runs) |
| q | uint64 | packet id — the packet's AMReX 64-bit `idcpu` identity (opaque here; uniqueness contract in [particle-container-and-gpu](./particle-container-and-gpu.md)) |
| e | uint32 | event counter — indexes the physics events of packet q's life (0 = creation/emission sampling; incremented per event as the consuming domain spec pins) |
| k | uint32 | draw index within event e — consumed consecutively from 0 in the order the consuming spec documents |
| r | double (IEEE-754 binary64) | uniform variate in [0, 1), 53 significant bits, error-free construction |

The output is a deterministic function of the four arguments only: independent of execution order, thread assignment, rank count, device, and batching. Distinct `(S, q, e, k)` tuples yield statistically independent variates (the Philox counter-mode guarantee). All non-uniform distributions (exponential via −ln(r), isotropic angles, spectral sampling, …) are constructed from these uniforms by the consuming spec's stated transformation and draw order.

## Correctness requirements

- **[MCNX-RNG-01] Algorithm pinned: Philox4x32-10.** The generator block function is Philox4x32-10 exactly: state = four 32-bit counter words `(c0, c1, c2, c3)` and two 32-bit key words `(k0, k1)`; **10 rounds**; multiplication constants `M0 = 0xD2511F53`, `M1 = 0xCD9E8D57`; Weyl key-schedule constants `W0 = 0x9E3779B9`, `W1 = 0xBB67AE85`. One round maps `(c0, c1, c2, c3)` with key `(k0, k1)` to

  ```text
  hi0·2³² + lo0 = M0 · c0        (full 32×32 → 64-bit product)
  hi1·2³² + lo1 = M1 · c2
  (c0, c1, c2, c3) ← (hi1 ^ c1 ^ k0,  lo1,  hi0 ^ c3 ^ k1,  lo0)
  ```

  Round 1 uses the key as given; before each subsequent round the key is bumped: `k0 ← k0 + W0`, `k1 ← k1 + W1` (mod 2³²). The output block `(x0, x1, x2, x3)` is the counter after the 10th round.

- **[MCNX-RNG-02] Known-answer vectors, exact equality.** The implementation must reproduce the published Random123 KAT vectors bit-for-bit (counter and key words given little-endian word order, values in hex):

  | Counter (c0 c1 c2 c3) | Key (k0 k1) | Output (x0 x1 x2 x3) |
  |---|---|---|
  | 00000000 00000000 00000000 00000000 | 00000000 00000000 | 6627e8d5 e169c58d bc57ac4c 9b00dbd8 |
  | ffffffff ffffffff ffffffff ffffffff | ffffffff ffffffff | 408f276d 41c83b0e a20bc7c6 6d5451fd |
  | 243f6a88 85a308d3 13198a2e 03707344 | a4093822 299f31d0 | d16cfe09 94fdcceb 5001e420 24126ea1 |

  These are exact-equality fixtures: any mismatch in any word is a failure. If the round-function prose above and these fixtures ever disagree, the fixtures govern and this spec must be corrected.

- **[MCNX-RNG-03] Counter/key packing, exact.** The tuple `(S, q, e, k)` maps to Philox inputs as

  ```text
  b  = k >> 1                      (block index: draw pair within event)
  c0 = b
  c1 = e
  c2 = q & 0xFFFFFFFF              (packet id, low 32 bits)
  c3 = q >> 32                     (packet id, high 32 bits)
  k0 = S & 0xFFFFFFFF              (seed, low 32 bits)
  k1 = S >> 32                     (seed, high 32 bits)
  ```

  Each Philox block yields two doubles (next requirement); draw `k` uses word pair `j = k & 1` of block `b`.

- **[MCNX-RNG-04] [0, 1) double construction, exact.** From output block `(x0, x1, x2, x3)` and word-pair index `j ∈ {0, 1}`:

  ```text
  u64 = x_{2j}  +  2³² · x_{2j+1}      (unsigned 64-bit)
  r   = (u64 >> 11) × 2^-53
  ```

  i.e. the top 53 bits of the 64-bit word scaled by `2^-53`. The integer-to-double scaling is error-free (every value of `u64 >> 11` is exactly representable in binary64), yields exactly 53 significant bits, and gives `r ∈ [0, 1)`: `r = 0.0` iff `u64 < 2^11`, and the maximum value is `1 − 2^-53 < 1`. No other construction (multiplication by `2^-64`, float intermediates, rejection) is permitted.

- **[MCNX-RNG-05] All physics draws are this pure function; stateful RNG forbidden.** Every random draw that affects physics — emission sampling, event-time sampling, scattering redraws, any stochastic decision on packets — is an evaluation of `u(S, q, e, k)` with a documented `(e, k)` assignment in the consuming spec. AMReX's stateful per-thread RNG (`amrex/Src/Base/AMReX_Random.H`, `amrex/Src/Base/AMReX_RandomEngine.H`, including `ParallelForRNG`) is **forbidden for physics draws**. Rationale (binding, recorded): its output depends on thread scheduling, so it cannot produce execution-order-independent per-packet streams, and it would make golden data irreproducible.

- **[MCNX-RNG-06] Single-rank exact reproducibility.** In single-rank CPU test mode at fixed seed, the per-packet event sequence — every event's type, drawn values, and resulting packet state updates, in each packet's own event order — is exactly reproducible run-to-run (bitwise; exact tier). This is what makes fixed-seed golden outputs work with the Cactus regression harness at its `ABSTOL=RELTOL=1e-12` defaults (`flesh/lib/sbin/RunTestUtils.pl`).

- **[MCNX-RNG-07] Statistical acceptance bar: 4σ at pinned counts and seeds.** Every statistical Monte Carlo check in the corpus passes iff

  ```text
  |estimate − expected| ≤ 4σ
  ```

  where σ is the standard error of the estimator (analytically propagated, or the documented sample estimate) at the **pinned packet count** and **pinned seed**. Defaults, binding unless a benchmark in [verification-suite-design](./verification-suite-design.md) explicitly pins different values: packet count `N_p = 1048576` (2^20); primary seed `S = 1296518744` (0x4D434E58); secondary seed `S = 20260730` for seed-robustness cross-checks. Golden statistical data is always generated with the primary seed. A statistical check may not be re-run with fresh seeds to pass ("seed shopping"); at pinned seed the check is deterministic.

- **[MCNX-RNG-08] No bitwise cross-run/cross-rank/cross-device requirement.** Grid-deposited tallies (source terms, moments) use unordered atomic accumulation; their values are reproducible across rank counts, devices, and runs only to floating-point-reduction tolerance, and no MCNuX check may demand bitwise equality of grid tallies across those axes (that would prescribe deposition ordering — an internal). Cross-configuration comparisons of grid quantities use norm-based tolerances stated by the comparing spec.

### Conventions restated (the subset this spec uses)

- All reals are IEEE-754 binary64 (`double`); Philox words are unsigned 32-bit with wraparound arithmetic mod 2³²; `S` and `q` are unsigned 64-bit.
- Neutrino physics quantities drawn from these variates carry the units of [conventions-and-units](./conventions-and-units.md); this spec's output `r` is dimensionless.

## Verification

- **KAT test (exact).** Evaluate the Philox4x32-10 block function on the three KAT vectors of [MCNX-RNG-02]; require bit-exact equality of all four output words for each. Run on every supported backend (host and device); the function is identical everywhere.
- **Double-construction closed forms (exact).** `u64 = 0 → r = 0.0`; `u64 = 2^11 − 1 → r = 0.0`; `u64 = 2^11 → r = 2^-53`; `u64 = 2^64 − 1 → r = 1 − 2^-53` (the largest double below 1.0); `r < 1.0` for all inputs; `r` is nondecreasing in `u64`.
- **Purity/independence (exact).** Evaluating `u(S, q, e, k)` for a fixed tuple must give bitwise-identical results across: repeated calls, permuted evaluation order over a batch of tuples, host vs device, and any batching strategy.
- **Uniformity at the pinned bar (statistical, 4σ).** With `N = 1048576` draws at the primary seed `1296518744` over a documented tuple sweep: sample mean within 4σ of 1/2 (σ = 1/√(12N)), sample second moment within 4σ of 1/3. These are acceptance smoke checks at the pinned bar, not a test-suite replacement for the published generator's statistical quality.
- **Reproducibility harness (exact / golden parity 1e-12).** A fixed-seed single-rank benchmark run twice produces identical per-packet event logs (exact) and identical golden TSV/norm outputs within the harness tolerance `1e-12` ([verification-suite-design](./verification-suite-design.md) pins the concrete benchmarks).

### Mechanical (validator)

`bash tools/validate_specs.sh` asserts this file carries the 7 mandated sections in order; contains the literals `Philox4x32-10`, `0xD2511F53`, `0xCD9E8D57`, `0x9E3779B9`, `0xBB67AE85`, the KAT words `6627e8d5`, `408f276d`, `d16cfe09`, the bar `4σ`, the pinned seed `1296518744`, and the DOI `10.1145/2063384.2063405`; that the DOI appears in this file's Source of truth registry and the README master registry; that its cited `amrex/` and `flesh/` paths resolve; and that its `[MCNX-RNG-NN]` ids are declared here and covered by the [verification-suite-design](./verification-suite-design.md) matrix.

## Implementation freedom

- Kernel fusion, batching, caching, and vectorization of draws; generating both doubles of a block at once and caching the second; recomputing vs storing; where the generator runs (host or device) — provided every observable drawn value equals the pure function `u(S, q, e, k)`.
- Internal, non-physics uses of randomness (e.g. randomized diagnostics timers) are unconstrained, provided they cannot affect physics results or golden outputs.
- The API shape (functor, free function, class), word-level micro-optimizations (e.g. `mulhilo` via 64-bit multiply vs intrinsics), and loop structure — the KAT vectors and purity checks are the arbiters.
- How the event counter `e` is stored/advanced per packet (it is packet state; [packet-representation-and-sampling](./packet-representation-and-sampling.md) owns the component, this spec owns its meaning).

## Open questions / assumptions

- **Assumption: no counter-based exemplar in the reference stack (verified).** Neither AMReX, WarpX, CarpetX, the flesh, nor WeakLibInterp contains a counter-based RNG usable as a production exemplar; this is why the algorithm is pinned from the published Random123 source rather than mirrored from stack code. The restated constants and KAT vectors above were extracted from the Random123 reference implementation and are normative as written.
- **Assumption: draw budget fits the packing.** The packing gives each (packet, event) pair 2³² blocks = 2³³ doubles, and each packet 2³² events; no plausible MCNuX event consumes anywhere near that. If a future scheme exceeded it, the packing — and with it all golden data — would need a deliberate re-pin.
- **Assumption: 4σ flake budget.** At 4σ, a single statistical check false-fails with probability ≈ 6.3e-5; across the corpus's ~dozens of statistical checks the aggregate false-failure probability stays below ~1e-2 per full-suite run. Checks are deterministic at pinned seed, so a false failure is a fixed property of (benchmark, seed) discovered once, not a flaky test; the remedy is documented re-pinning of that benchmark's expected band, never seed shopping.
