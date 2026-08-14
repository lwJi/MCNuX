#ifndef MCNUX_RNG_HXX
#define MCNUX_RNG_HXX

#include <cstdint>

// The MCNuX random-number generator: the pure function u(S, q, e, k) -> [0, 1).
//
// Governing spec: specs/rng-and-statistical-acceptance.md, requirements
// [MCNX-RNG-01..06].
//
// [MCNX-RNG-01] pins the algorithm to Philox4x32-10 (Random123,
// doi:10.1145/2063384.2063405) with its multiplication and Weyl constants;
// [MCNX-RNG-02] pins three published known-answer vectors as exact-equality
// fixtures; [MCNX-RNG-03] pins the packing of (seed, packet id, event counter,
// draw index) into counter/key words; [MCNX-RNG-04] pins the output-word ->
// binary64 construction. The static_asserts at the bottom of this header *are*
// the test suite for those four requirements: any translation unit that
// includes this header re-verifies them at build time.
//
// [MCNX-RNG-05]: this file deliberately does not include, name, or call
// AMReX's stateful per-thread device RNG (AMReX_Random.H, AMReX_RandomEngine.H,
// ParallelForRNG). Physics draws must be evaluations of the pure function
// below, because a stateful engine's output depends on thread scheduling and
// would make golden data irreproducible.
//
// Device callability: every function here is `constexpr` and touches only
// fixed-width unsigned integers and one double multiply -- no allocation, no
// globals, no library calls. AMReX's nvcc toolchain compiles with
// `--expt-relaxed-constexpr` (amrex/Tools/GNUMake/comps/nvcc.mak), so plain
// constexpr functions are callable from device code; annotating them with
// AMREX_GPU_HOST_DEVICE would require an AMReX include and is deliberately
// avoided so that host code, device code, and standalone unit tests can all
// include this header (same convention as mcnux_units.hxx).
//
// Deferral, recorded: the host-vs-device bitwise leg of [MCNX-RNG-06]'s purity
// check cannot be exercised while the sandbox build is CPU/MPI-only; it is
// deferred to the GPU build (T25d). The repeated-call, permuted-order, and
// batching-independence legs are asserted below at compile time, which -- being
// constant-folded by the compiler's own arithmetic -- is the strongest evidence
// available in a CPU-only configuration.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes.

namespace MCNuX {

// One Philox4x32 output block: four 32-bit words in little-endian word order.
struct PhiloxBlock {
  std::uint32_t x0;
  std::uint32_t x1;
  std::uint32_t x2;
  std::uint32_t x3;
};

namespace detail {

// ---------------------------------------------------------------------------
// Pinned Philox4x32-10 constants  [MCNX-RNG-01]
// ---------------------------------------------------------------------------
// rng-and-statistical-acceptance.md:51. These four literals appear in MCNuX
// exactly once -- here.

inline constexpr std::uint32_t philox_M0 = 0xD2511F53u; // multiplier for c0
inline constexpr std::uint32_t philox_M1 = 0xCD9E8D57u; // multiplier for c2
inline constexpr std::uint32_t philox_W0 = 0x9E3779B9u; // Weyl bump for k0
inline constexpr std::uint32_t philox_W1 = 0xBB67AE85u; // Weyl bump for k1

inline constexpr int philox_rounds = 10;

// Half of a full 32x32 -> 64 bit product.
struct HiLo {
  std::uint32_t hi;
  std::uint32_t lo;
};

// Full 32x32 -> 64 bit unsigned product, split into high and low words. Done
// by widening to uint64_t rather than by an intrinsic: portable, legal in a
// constant expression, and exact (a 64-bit product of 32-bit operands never
// overflows). Word-level micro-optimization is explicitly free
// (rng-and-statistical-acceptance.md:129); the KAT vectors are the arbiter.
constexpr HiLo mulhilo(std::uint32_t a, std::uint32_t b) noexcept {
  const std::uint64_t p =
      static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b);
  return HiLo{static_cast<std::uint32_t>(p >> 32),
              static_cast<std::uint32_t>(p & 0xFFFFFFFFu)};
}

// One Philox4x32 round with key (k0, k1)
// (rng-and-statistical-acceptance.md:53-57):
//   hi0:lo0 = M0 * c0,  hi1:lo1 = M1 * c2
//   (c0, c1, c2, c3) <- (hi1 ^ c1 ^ k0,  lo1,  hi0 ^ c3 ^ k1,  lo0)
// All arithmetic is unsigned, i.e. wraparound mod 2^32 as the spec requires
// (rng-and-statistical-acceptance.md:110); signed types would be undefined
// behaviour on overflow and unusable in a constant expression.
constexpr PhiloxBlock philox_round(const PhiloxBlock &c, std::uint32_t k0,
                                   std::uint32_t k1) noexcept {
  const HiLo p0 = mulhilo(philox_M0, c.x0);
  const HiLo p1 = mulhilo(philox_M1, c.x2);
  return PhiloxBlock{static_cast<std::uint32_t>(p1.hi ^ c.x1 ^ k0), p1.lo,
                     static_cast<std::uint32_t>(p0.hi ^ c.x3 ^ k1), p0.lo};
}

// The Philox inputs produced by the packing of [MCNX-RNG-03].
struct PhiloxInput {
  std::uint32_t c0;
  std::uint32_t c1;
  std::uint32_t c2;
  std::uint32_t c3;
  std::uint32_t k0;
  std::uint32_t k1;
};

constexpr std::uint32_t lo32(std::uint64_t v) noexcept {
  return static_cast<std::uint32_t>(v & 0xFFFFFFFFu);
}

constexpr std::uint32_t hi32(std::uint64_t v) noexcept {
  return static_cast<std::uint32_t>(v >> 32);
}

} // namespace detail

// ---------------------------------------------------------------------------
// Block function  [MCNX-RNG-01]
// ---------------------------------------------------------------------------
// Philox4x32-10: ten rounds over counter (c0, c1, c2, c3) with key (k0, k1).
// Round 1 uses the key as given; the key is bumped by the Weyl constants
// *before* each of rounds 2..10 (nine bumps in total, never after the last
// round). The returned block is the counter state after round 10.
constexpr PhiloxBlock philox4x32_10(std::uint32_t c0, std::uint32_t c1,
                                    std::uint32_t c2, std::uint32_t c3,
                                    std::uint32_t k0,
                                    std::uint32_t k1) noexcept {
  PhiloxBlock c{c0, c1, c2, c3};
  for (int r = 0; r < detail::philox_rounds; ++r) {
    if (r > 0) {
      k0 = static_cast<std::uint32_t>(k0 + detail::philox_W0);
      k1 = static_cast<std::uint32_t>(k1 + detail::philox_W1);
    }
    c = detail::philox_round(c, k0, k1);
  }
  return c;
}

// ---------------------------------------------------------------------------
// Counter/key packing  [MCNX-RNG-03]
// ---------------------------------------------------------------------------
// (S, q, e, k) -> Philox inputs (rng-and-statistical-acceptance.md:71-83).
// The block index is b = k >> 1: each block yields two doubles, so draws
// k = 2b and k = 2b+1 share block b and are distinguished by the word-pair
// index j = k & 1.
constexpr detail::PhiloxInput philox_input(std::uint64_t S, std::uint64_t q,
                                           std::uint32_t e,
                                           std::uint32_t k) noexcept {
  const std::uint32_t b = k >> 1;
  return detail::PhiloxInput{b,
                             e,
                             detail::lo32(q),
                             detail::hi32(q),
                             detail::lo32(S),
                             detail::hi32(S)};
}

// The Philox block backing draw index k of event e of packet q under seed S.
constexpr PhiloxBlock philox_block(std::uint64_t S, std::uint64_t q,
                                   std::uint32_t e, std::uint32_t k) noexcept {
  const detail::PhiloxInput in = philox_input(S, q, e, k);
  return philox4x32_10(in.c0, in.c1, in.c2, in.c3, in.k0, in.k1);
}

// ---------------------------------------------------------------------------
// [0, 1) double construction  [MCNX-RNG-04]
// ---------------------------------------------------------------------------
// r = (u64 >> 11) * 2^-53 (rng-and-statistical-acceptance.md:85-92). The
// integer-to-double conversion is error-free: u64 >> 11 < 2^53, so every value
// is exactly representable in binary64, and scaling by a power of two is exact.
// No other construction is permitted -- in particular no multiplication by
// 2^-64, no float intermediate, and no rejection loop.
constexpr double uniform_from_u64(std::uint64_t u64) noexcept {
  return static_cast<double>(u64 >> 11) * 0x1p-53;
}

// Word pair j of a block, assembled little-endian:
//   u64 = x_{2j} + 2^32 * x_{2j+1}.
constexpr std::uint64_t block_word_pair(const PhiloxBlock &x,
                                        std::uint32_t j) noexcept {
  const std::uint32_t lo = (j == 0u) ? x.x0 : x.x2;
  const std::uint32_t hi = (j == 0u) ? x.x1 : x.x3;
  return static_cast<std::uint64_t>(lo) +
         (static_cast<std::uint64_t>(hi) << 32);
}

// ---------------------------------------------------------------------------
// The MCNuX random-number generator  [MCNX-RNG-05, MCNX-RNG-06]
// ---------------------------------------------------------------------------
// u(S, q, e, k) -> r in [0, 1). Deterministic in its four arguments alone:
// independent of execution order, thread assignment, rank count, device, and
// batching. Every physics draw in MCNuX is an evaluation of this function with
// an (e, k) assignment documented by the consuming domain spec.
//
//   S  run seed
//   q  packet id (the packet's AMReX 64-bit idcpu identity; opaque here)
//   e  event counter within packet q's life (0 = creation/emission sampling)
//   k  draw index within event e, consumed consecutively from 0
constexpr double u(std::uint64_t S, std::uint64_t q, std::uint32_t e,
                   std::uint32_t k) noexcept {
  return uniform_from_u64(block_word_pair(philox_block(S, q, e, k), k & 1u));
}

// ---------------------------------------------------------------------------
// Compile-time verification (rng-and-statistical-acceptance.md:113-119)
// ---------------------------------------------------------------------------

namespace detail {

// (a) Known-answer vectors  [MCNX-RNG-02]
//     rng-and-statistical-acceptance.md:61-69. Exact-equality fixtures: any
//     mismatch in any word is a failure. Counter, key, and output are listed in
//     little-endian *word* order; each hex token below is a plain uint32 value
//     (no byte swapping). These fixtures govern over the round-function prose.

constexpr bool kat_matches(const PhiloxBlock &got, std::uint32_t x0,
                           std::uint32_t x1, std::uint32_t x2,
                           std::uint32_t x3) noexcept {
  return got.x0 == x0 && got.x1 == x1 && got.x2 == x2 && got.x3 == x3;
}

} // namespace detail

static_assert(detail::kat_matches(philox4x32_10(0x00000000u, 0x00000000u,
                                                0x00000000u, 0x00000000u,
                                                0x00000000u, 0x00000000u),
                                  0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu,
                                  0x9b00dbd8u),
              "Philox4x32-10 KAT 1 (zero counter, zero key) failed");
static_assert(detail::kat_matches(philox4x32_10(0xffffffffu, 0xffffffffu,
                                                0xffffffffu, 0xffffffffu,
                                                0xffffffffu, 0xffffffffu),
                                  0x408f276du, 0x41c83b0eu, 0xa20bc7c6u,
                                  0x6d5451fdu),
              "Philox4x32-10 KAT 2 (all-ones counter, all-ones key) failed");
static_assert(detail::kat_matches(philox4x32_10(0x243f6a88u, 0x85a308d3u,
                                                0x13198a2eu, 0x03707344u,
                                                0xa4093822u, 0x299f31d0u),
                                  0xd16cfe09u, 0x94fdccebu, 0x5001e420u,
                                  0x24126ea1u),
              "Philox4x32-10 KAT 3 (pi-digits counter and key) failed");

// A wrong number of rounds, or a key bumped after the last round instead of
// before rounds 2..10, is the most likely transcription error and would break
// the KATs above; these two asserts pin the constants themselves so that a
// mistyped literal cannot be masked by a compensating error elsewhere.
static_assert(detail::philox_rounds == 10, "Philox4x32-10 is exactly 10 rounds");
static_assert(detail::philox_M0 == 0xD2511F53u &&
                  detail::philox_M1 == 0xCD9E8D57u &&
                  detail::philox_W0 == 0x9E3779B9u &&
                  detail::philox_W1 == 0xBB67AE85u,
              "Philox4x32-10 constants disagree with the pinned values");

// (b) [0, 1) double construction, closed forms  [MCNX-RNG-04]
//     rng-and-statistical-acceptance.md:116. Exact equality -- no tolerance:
//     the construction is error-free by design, so any tolerance here would
//     hide a defect.

static_assert(uniform_from_u64(0u) == 0.0, "u64 = 0 must give exactly 0.0");
static_assert(uniform_from_u64((1ull << 11) - 1u) == 0.0,
              "u64 = 2^11 - 1 must still give exactly 0.0");
static_assert(uniform_from_u64(1ull << 11) == 0x1p-53,
              "u64 = 2^11 must give exactly 2^-53");
static_assert(uniform_from_u64(~0ull) == 1.0 - 0x1p-53,
              "u64 = 2^64 - 1 must give the largest double below 1.0");
static_assert(uniform_from_u64(~0ull) < 1.0,
              "the maximum variate must be strictly below 1.0");
static_assert(1.0 - 0x1p-53 < 1.0 && 1.0 - 0x1p-54 == 1.0,
              "1 - 2^-53 must be the largest double below 1.0");
static_assert(uniform_from_u64(3ull << 11) == 3.0 * 0x1p-53,
              "the scaling must be exact for small multiples of 2^-53");
static_assert(uniform_from_u64(1ull << 63) == 0.5,
              "the top bit alone must give exactly 1/2");
static_assert(uniform_from_u64((1ull << 63) + (1ull << 62)) == 0.75,
              "u64 = 3*2^62 must give exactly 3/4");

namespace detail {

// r is nondecreasing in u64, and always in [0, 1). Probed across the whole
// 64-bit range including both ends and the 2^11 truncation boundary.
constexpr bool uniform_is_monotone_and_in_range() noexcept {
  const std::uint64_t vs[] = {0ull,
                              1ull,
                              (1ull << 11) - 1ull,
                              1ull << 11,
                              (1ull << 11) + 1ull,
                              1ull << 20,
                              1ull << 32,
                              0x0123456789abcdefull,
                              1ull << 62,
                              1ull << 63,
                              (1ull << 63) + (1ull << 62),
                              ~0ull - 1ull,
                              ~0ull};
  const int n = static_cast<int>(sizeof(vs) / sizeof(vs[0]));
  for (int i = 0; i < n; ++i) {
    const double r = uniform_from_u64(vs[i]);
    if (!(r >= 0.0) || !(r < 1.0))
      return false;
    if (i > 0 && !(uniform_from_u64(vs[i - 1]) <= r))
      return false;
  }
  return true;
}

} // namespace detail

static_assert(detail::uniform_is_monotone_and_in_range(),
              "r must lie in [0, 1) and be nondecreasing in u64");

namespace detail {

// (c) Purity / permutation independence  [MCNX-RNG-06]
//     rng-and-statistical-acceptance.md:117. Doubles are compared with `==`:
//     the construction is a nonnegative integer times a power of two, so it can
//     produce neither NaN nor -0.0, and `==` therefore *is* bitwise equality
//     for these values.

struct DrawKey {
  std::uint64_t S;
  std::uint64_t q;
  std::uint32_t e;
  std::uint32_t k;
};

inline constexpr int purity_batch_size = 8;

// A batch spanning the pinned seeds, a 64-bit packet id straddling the word
// split, several events, and both parities of the draw index.
constexpr DrawKey purity_batch(int i) noexcept {
  const DrawKey ts[purity_batch_size] = {
      {1296518744ull, 0ull, 0u, 0u},
      {1296518744ull, 0ull, 0u, 1u},
      {1296518744ull, 0ull, 0u, 2u},
      {1296518744ull, 0ull, 1u, 0u},
      {1296518744ull, 0x00000001FFFFFFFFull, 7u, 5u},
      {20260730ull, 0x00000001FFFFFFFFull, 7u, 5u},
      {0xFEDCBA9876543210ull, 0x0123456789ABCDEFull, 4294967295u, 4294967294u},
      {0ull, 0ull, 0u, 0u}};
  return ts[i];
}

// Every leg of the purity requirement that a CPU-only build can exercise:
// repeated calls, permuted evaluation order, reverse order, and a "batched"
// evaluation that computes whole blocks and extracts both doubles at once
// (the caching strategy that rng-and-statistical-acceptance.md:127 permits).
constexpr bool purity_holds() noexcept {
  double fwd[purity_batch_size] = {};
  for (int i = 0; i < purity_batch_size; ++i) {
    const DrawKey t = purity_batch(i);
    fwd[i] = u(t.S, t.q, t.e, t.k);
  }

  // Repeated evaluation of the same tuple is bitwise identical.
  for (int i = 0; i < purity_batch_size; ++i) {
    const DrawKey t = purity_batch(i);
    if (!(u(t.S, t.q, t.e, t.k) == fwd[i]))
      return false;
    if (!(u(t.S, t.q, t.e, t.k) == u(t.S, t.q, t.e, t.k)))
      return false;
  }

  // Evaluation in a permuted order yields the same value per tuple, so the
  // result depends on the tuple alone and not on when it was drawn.
  const int perm[purity_batch_size] = {5, 2, 7, 0, 3, 6, 1, 4};
  double permuted[purity_batch_size] = {};
  for (int i = 0; i < purity_batch_size; ++i) {
    const int j = perm[i];
    const DrawKey t = purity_batch(j);
    permuted[j] = u(t.S, t.q, t.e, t.k);
  }
  for (int i = 0; i < purity_batch_size; ++i)
    if (!(permuted[i] == fwd[i]))
      return false;

  // Reverse order, likewise.
  for (int i = purity_batch_size - 1; i >= 0; --i) {
    const DrawKey t = purity_batch(i);
    if (!(u(t.S, t.q, t.e, t.k) == fwd[i]))
      return false;
  }

  // Batching independence: taking both doubles out of one shared block gives
  // exactly what draw-at-a-time evaluation gives.
  for (int i = 0; i < purity_batch_size; ++i) {
    const DrawKey t = purity_batch(i);
    const PhiloxBlock x = philox_block(t.S, t.q, t.e, t.k);
    if (!(uniform_from_u64(block_word_pair(x, t.k & 1u)) == fwd[i]))
      return false;
  }

  // Every variate stays in range.
  for (int i = 0; i < purity_batch_size; ++i)
    if (!(fwd[i] >= 0.0) || !(fwd[i] < 1.0))
      return false;

  // Distinct tuples give distinct variates. Not a statistical claim -- just
  // evidence that the packing does not collapse different (S, q, e, k) onto
  // one stream (e.g. by dropping e, or by ignoring the high word of q).
  for (int i = 0; i < purity_batch_size; ++i)
    for (int j = i + 1; j < purity_batch_size; ++j)
      if (fwd[i] == fwd[j])
        return false;

  return true;
}

// Consecutive draws k = 2b and k = 2b+1 come from the *same* Philox block b,
// as word pairs j = 0 and j = 1 respectively; draws in different blocks use
// different counters. This pins the [MCNX-RNG-03] pairing that the purity
// checks above would otherwise not distinguish from "one block per draw".
constexpr bool pairing_holds() noexcept {
  const std::uint64_t S = 1296518744ull;
  const std::uint64_t q = 0x0123456789ABCDEFull;
  const std::uint32_t e = 3u;
  for (std::uint32_t b = 0u; b < 4u; ++b) {
    const std::uint32_t k0 = 2u * b;
    const std::uint32_t k1 = 2u * b + 1u;
    const PhiloxBlock x = philox_block(S, q, e, k0);

    // Both draws of the pair are backed by one and the same block.
    const PhiloxBlock y = philox_block(S, q, e, k1);
    if (!(x.x0 == y.x0 && x.x1 == y.x1 && x.x2 == y.x2 && x.x3 == y.x3))
      return false;

    // ... and that block's counter is exactly the packing of [MCNX-RNG-03].
    const PhiloxInput in = philox_input(S, q, e, k1);
    if (!(in.c0 == b && in.c1 == e && in.c2 == lo32(q) && in.c3 == hi32(q) &&
          in.k0 == lo32(S) && in.k1 == hi32(S)))
      return false;

    // The even draw takes word pair 0, the odd draw word pair 1.
    if (!(u(S, q, e, k0) == uniform_from_u64(block_word_pair(x, 0u))))
      return false;
    if (!(u(S, q, e, k1) == uniform_from_u64(block_word_pair(x, 1u))))
      return false;

    // The two draws of a pair differ, and so do the blocks of adjacent pairs.
    if (u(S, q, e, k0) == u(S, q, e, k1))
      return false;
    const PhiloxBlock z = philox_block(S, q, e, k1 + 1u);
    if (x.x0 == z.x0 && x.x1 == z.x1 && x.x2 == z.x2 && x.x3 == z.x3)
      return false;
  }
  return true;
}

} // namespace detail

static_assert(detail::purity_holds(),
              "u(S, q, e, k) must be a pure function of its arguments: "
              "identical under repeated calls, permuted order, and batching");
static_assert(detail::pairing_holds(),
              "draws k = 2b and k = 2b+1 must share Philox block b as word "
              "pairs 0 and 1");

} // namespace MCNuX

#endif // MCNUX_RNG_HXX
