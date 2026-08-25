#ifndef MCNUX_EMISSION_HXX
#define MCNUX_EMISSION_HXX

// Emission sampling — the pure functions of
// specs/packet-representation-and-sampling.md (count law, per-packet
// initialization quantities, cell key, creation RNG draw map). No scheduling,
// no cell loop, no container writes: the runtime creation loop that assembles
// these pieces (sqrt(-g) from the metric gather, the bin grid, the packet
// buffer) is a separate task.
//
// Governing spec requirements:
//   * [MCNX-PKT-02] — the emission count law (provenance arXiv:1708.08452
//     Eq. 18, restated normatively in the spec):
//         N_p = g_s sqrt(-g) dV dt eta_b(s) / E_p ,  g_s = (1, 1, 4)
//     Exactly floor(N_p) packets are created, plus one more with probability
//     N_p - floor(N_p) (one uniform draw per (cell, step, s, b)). The nu_x
//     degeneracy g = 4 enters the corpus EXACTLY ONCE, here, through
//     MCNuX::degeneracy(Species) (mcnux_units.hxx pins creation as the sole
//     locus) — never in per-species transport coefficients and never again
//     downstream (the packet weight N already counts all lumped states).
//     The Bernoulli convention is `u_draw < frac` with u in [0, 1), which
//     realizes probability exactly N_p - floor(N_p): u == 0 fires whenever
//     frac > 0, and an exactly-integer N_p can never round up.
//   * [MCNX-PKT-03] — per-packet initialization: fluid-frame energy
//     nu = 0.5 (E_{b-1} + E_b) (the bin center, EXACT — energies are never
//     drawn continuously in the baseline) and weight N = E_p/nu (each packet
//     carries fluid-frame energy E_p). This assigned nu is a DIFFERENT
//     quantity from mcnux_interactions.hxx's fluid_frame_energy() (which
//     computes -p.u from an already-built momentum); do not conflate them.
//   * [MCNX-PKT-05] — the creation RNG consumption map. All draws are
//     evaluations of the pure u(S, q, e, k) of mcnux_rng.hxx; AMReX's
//     stateful device RNG is forbidden for them. Cell-level count draw:
//     q = K_cell, e = n (transport-step index from 0), k = 3b + s.
//     Per-packet creation draws at event e = 0 of the new packet's own id:
//     k = 0, 1, 2 -> position x, y, z; k = 3 -> creation time; k = 4 -> u1
//     (cos theta); k = 5 -> u2 (phi), consumed in exactly this order and no
//     others. The cell key packing
//         K_cell = 2^63 + patch*2^56 + level*2^48
//                  + (k_c - k_0)*2^32 + (j_c - j_0)*2^16 + (i_c - i_0)
//     relies on the spec's capacity assumption (patch < 128, level < 256,
//     per-dimension level extents < 65536) and on bit 63 disambiguating cell
//     keys from packet idcpu values, which never have bit 63 set
//     ([MCNX-GPU-04] of specs/particle-container-and-gpu.md: AMReX packs
//     40-bit id + 24-bit cpu; bit 63 would require cpu >= 2^23 ranks).
//
// Bin-grid boundary (deliberate scope cut): every function here is
// BIN-AGNOSTIC — the count law consumes a caller-supplied scalar eta_b.
// specs/opacity-eos-evaluation.md [MCNX-OPA-04] pins the producer formula
// eta_b(s) = eta(s, E_b) * dE_b with E_b on the tabulated energy grid;
// bin_integrated_eta below is the sanctioned bridge from the SPECTRAL eta of
// the frozen coefficient interface (mcnux_coefficients.hxx, Coefficients.eta
// semantics unchanged) to the bin-integrated eta_b this header consumes. The
// concrete grid (dE_b from the tabulated energy nodes, the edge-bin width
// convention, analytic-source binning) belongs to the runtime creation loop.
//
// RNG discipline (mcnux_tetrad.hxx precedent): the wrappers below do ONLY
// (q, e, k) bookkeeping; direction math is the existing
// fluid_frame_direction(u1, u2) of mcnux_tetrad.hxx taking already-drawn
// doubles.
//
// Naming note: the pinned creation draw-map constants carry a `creation`
// infix (draw_k_creation_phi = 5, not draw_k_phi) because
// mcnux_interactions.hxx already owns draw_k_phi = 3 in this namespace for
// the DIFFERENT per-episode (e >= 1) draw map; creation is e = 0 with its own
// k space (emission k = 4/5 vs scattering k = 2/3).
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// C++, constexpr throughout (the math here is rational), so host code, device
// code, and standalone unit tests can all include it. The static_asserts at
// the bottom are the compile-time half of the T14a test surface; the runtime
// half lives in the `unit-selftest` battery (mcnux_selftest_emission.cxx,
// rows 118..127).

#include "mcnux_rng.hxx"   // u(S, q, e, k)
#include "mcnux_units.hxx" // Species, degeneracy, species_index

#include <cassert>
#include <cstdint>

namespace MCNuX {

// ---------------------------------------------------------------------------
// Bin-integrated emissivity bridge  [MCNX-OPA-04] producer formula
// ---------------------------------------------------------------------------

// eta_b = eta(s, E_b) * dE_b: converts the SPECTRAL emissivity of the
// coefficient interface (Coefficients.eta, MeV cm^-3 s^-1 MeV^-1 in the
// Kirchhoff form's units) to the bin-integrated eta_b the count law consumes.
// The bin width dE_b comes from the caller's grid (deferred to the creation
// loop); this helper is the only sanctioned bridge between the two.
constexpr double bin_integrated_eta(double eta_spectral, double dE_b) noexcept {
  return eta_spectral * dE_b;
}

// ---------------------------------------------------------------------------
// Emission count law  [MCNX-PKT-02]
// ---------------------------------------------------------------------------

// N_p = g_s sqrt(-g) dV dt eta_b / E_p. The degeneracy g_s = (1, 1, 4)
// enters here EXACTLY ONCE in the corpus, via degeneracy(Species). The caller
// supplies sqrt(-g) = alpha sqrt(gamma) as one scalar (assembled by the
// creation loop from the metric gather) and evaluates the product in one
// documented unit system — the dimensionless result is unit-system
// independent to machine tier.
constexpr double packet_count(Species s, double sqrt_neg_g, double dV,
                              double dt, double eta_b, double E_p) noexcept {
  return static_cast<double>(degeneracy(s)) * sqrt_neg_g * dV * dt * eta_b /
         E_p;
}

// floor(N_p) packets plus one more iff u_draw < N_p - floor(N_p). With
// u_draw in [0, 1) the `<` comparison realizes the extra packet with
// probability exactly N_p - floor(N_p): the uniform grid {j * 2^-53} has
// exactly frac * 2^53 points below frac. u_draw == 0.0 (a reachable exact
// RNG output) fires whenever frac > 0; an exactly-integer N_p (frac == 0)
// never rounds up. Requires N_p >= 0 (guaranteed by construction: every
// count-law factor is non-negative), under which int64 truncation IS floor.
constexpr std::int64_t emission_count(double N_p, double u_draw) noexcept {
  assert(N_p >= 0.0);
  const std::int64_t base = static_cast<std::int64_t>(N_p);
  const double frac = N_p - static_cast<double>(base);
  return base + (u_draw < frac ? 1 : 0);
}

// ---------------------------------------------------------------------------
// Per-packet initialization quantities  [MCNX-PKT-03]
// ---------------------------------------------------------------------------

// Fluid-frame energy of every packet created in bin [E_lo, E_hi]: the bin
// center, exact — baseline energies are never drawn continuously.
constexpr double bin_center_energy(double E_lo, double E_hi) noexcept {
  return 0.5 * (E_lo + E_hi);
}

// Packet weight N = E_p/nu (> 0): each packet carries fluid-frame energy
// E_p; for nu_x, N counts all four lumped heavy-lepton states (g = 4 already
// entered through the count law, never again here).
constexpr double packet_weight(double E_p, double nu) noexcept {
  return E_p / nu;
}

// ---------------------------------------------------------------------------
// Cell key  [MCNX-PKT-05]
// ---------------------------------------------------------------------------

// The spec's capacity assumption, checked before packing: patch in [0, 128),
// level in [0, 256), each cell offset from the level's domain lower corner in
// [0, 65536). Exceeding any of these requires a deliberate,
// golden-data-invalidating re-pin of the packing — never a silent overflow.
constexpr bool cell_key_capacity_ok(int patch, int level, int ic, int jc,
                                    int kc, int i0, int j0, int k0) noexcept {
  const long long di = static_cast<long long>(ic) - i0;
  const long long dj = static_cast<long long>(jc) - j0;
  const long long dk = static_cast<long long>(kc) - k0;
  return patch >= 0 && patch < 128 && level >= 0 && level < 256 && di >= 0 &&
         di < 65536 && dj >= 0 && dj < 65536 && dk >= 0 && dk < 65536;
}

// K_cell = 2^63 + patch*2^56 + level*2^48 + (k_c-k_0)*2^32 + (j_c-j_0)*2^16
//          + (i_c-i_0), packed with shifts (mcnux_rng.hxx philox_input
// precedent). Under the capacity guard the six fields occupy disjoint bit
// ranges — flag bit 63 | patch bits 56..62 | level bits 48..55 | dk bits
// 32..47 | dj bits 16..31 | di bits 0..15 — so the ORs below equal the
// spec's sum exactly. Bit 63 is always set, disambiguating cell keys from
// packet idcpu values ([MCNX-GPU-04]: real ids never have bit 63 set).
constexpr std::uint64_t cell_key(int patch, int level, int ic, int jc, int kc,
                                 int i0, int j0, int k0) noexcept {
  assert(cell_key_capacity_ok(patch, level, ic, jc, kc, i0, j0, k0));
  const std::uint64_t di = static_cast<std::uint64_t>(ic - i0);
  const std::uint64_t dj = static_cast<std::uint64_t>(jc - j0);
  const std::uint64_t dk = static_cast<std::uint64_t>(kc - k0);
  return (std::uint64_t{1} << 63) |
         (static_cast<std::uint64_t>(patch) << 56) |
         (static_cast<std::uint64_t>(level) << 48) | (dk << 32) | (dj << 16) |
         di;
}

// ---------------------------------------------------------------------------
// Creation RNG draw map  [MCNX-PKT-05]
// ---------------------------------------------------------------------------

// The pinned per-packet creation draw indices, consumed in exactly this
// order at event e = 0 of the new packet's own id. (`creation` infix: see
// the naming note in the header comment — interactions' per-episode map owns
// the bare draw_k_phi name at a different (e, k) assignment.)
inline constexpr std::uint32_t draw_k_creation_pos_x = 0u;   // position x
inline constexpr std::uint32_t draw_k_creation_pos_y = 1u;   // position y
inline constexpr std::uint32_t draw_k_creation_pos_z = 2u;   // position z
inline constexpr std::uint32_t draw_k_creation_time = 3u;    // creation time
inline constexpr std::uint32_t draw_k_creation_costheta = 4u; // u1 (cos theta)
inline constexpr std::uint32_t draw_k_creation_phi = 5u;      // u2 (phi)

// The six per-packet creation uniforms, k = 0..5 in the pinned order;
// creation is e = 0 and consumes no other draws. Bookkeeping only: the
// position/time/direction math consumes these as plain doubles
// (fluid_frame_direction(u1, u2) of mcnux_tetrad.hxx for the direction pair).
struct CreationUniforms {
  double u_pos_x, u_pos_y, u_pos_z, u_time, u1, u2;
};

constexpr CreationUniforms creation_uniforms(std::uint64_t S,
                                             std::uint64_t q) noexcept {
  return {u(S, q, 0u, draw_k_creation_pos_x),
          u(S, q, 0u, draw_k_creation_pos_y),
          u(S, q, 0u, draw_k_creation_pos_z),
          u(S, q, 0u, draw_k_creation_time),
          u(S, q, 0u, draw_k_creation_costheta),
          u(S, q, 0u, draw_k_creation_phi)};
}

// The cell-level count draw (the Bernoulli of [MCNX-PKT-02]): q = K_cell,
// e = n (the transport-step index, counted from 0), k = 3b + s. NOTE: unlike
// the per-packet map above, this k is COMPUTED from the bin index b and the
// species index s — it is not a fixed constant.
constexpr double cell_count_uniform(std::uint64_t S, std::uint64_t K_cell,
                                    std::uint32_t n, int b,
                                    Species s) noexcept {
  return u(S, K_cell, n,
           static_cast<std::uint32_t>(3 * b + species_index(s)));
}

// ---------------------------------------------------------------------------
// Compile-time verification (the constexpr half of the T14a test surface;
// the runtime half pins the same fixtures as golden rows in the
// `unit-selftest` battery, mcnux_selftest_emission.cxx rows 118..127)
// ---------------------------------------------------------------------------

namespace detail {

// Count law: g enters exactly once and is visible as an exact factor of 4
// (nu_x) vs 1 (nu_e, nu_e_bar) at identical inputs; the g = 1 species agree
// bitwise.
constexpr bool count_law_degeneracy_holds() noexcept {
  const double sqrt_neg_g = 1.2, dV = 0.001, dt = 0.05, eta_b = 3.75,
               E_p = 10.0;
  const double n_e = packet_count(Species::NuE, sqrt_neg_g, dV, dt, eta_b, E_p);
  const double n_eb =
      packet_count(Species::NuEBar, sqrt_neg_g, dV, dt, eta_b, E_p);
  const double n_x = packet_count(Species::NuX, sqrt_neg_g, dV, dt, eta_b, E_p);
  return n_eb == n_e && n_x == 4.0 * n_e && n_e > 0.0;
}

// Floor/Bernoulli arithmetic at exact binary fixtures: both branches of the
// `u_draw < frac` comparison, the u == 0 edge, and the exactly-integer case.
constexpr bool emission_count_table_holds() noexcept {
  return emission_count(2.75, 0.5) == 3 &&  // u below frac -> +1
         emission_count(2.75, 0.75) == 2 && // u == frac: strict < -> no +1
         emission_count(2.75, 0.9) == 2 &&  // u above frac -> no +1
         emission_count(2.75, 0.0) == 3 &&  // u == 0 fires when frac > 0
         emission_count(3.0, 0.0) == 3 &&   // integer N_p never rounds up
         emission_count(0.25, 0.0) == 1 &&  // pure-Bernoulli cell, fires
         emission_count(0.25, 0.25) == 0 && // pure-Bernoulli cell, no fire
         emission_count(0.0, 0.0) == 0;     // zero emissivity -> no packets
}

// The pinned cell-key fixture: patch 3, level 5, offsets (7, 11, 13) ->
// hand-composed 0x83 | 05 | 000D | 000B | 0007 byte layout.
constexpr bool cell_key_fixture_holds() noexcept {
  return cell_key(3, 5, 17, 21, 33, 10, 10, 20) == 0x8305000D000B0007ull &&
         cell_key(0, 0, 0, 0, 0, 0, 0, 0) == 0x8000000000000000ull;
}

// Bit 63 is set on every packable key (disambiguation from packet idcpu).
constexpr bool cell_key_flag_bit_holds() noexcept {
  return (cell_key(127, 255, 65535, 65535, 65535, 0, 0, 0) >> 63) == 1u &&
         (cell_key(0, 0, 0, 0, 0, 0, 0, 0) >> 63) == 1u;
}

// The capacity predicate accepts the corners and rejects each overflow and
// a negative offset.
constexpr bool cell_key_capacity_table_holds() noexcept {
  return cell_key_capacity_ok(0, 0, 0, 0, 0, 0, 0, 0) &&
         cell_key_capacity_ok(127, 255, 65535, 65535, 65535, 0, 0, 0) &&
         !cell_key_capacity_ok(128, 0, 0, 0, 0, 0, 0, 0) &&
         !cell_key_capacity_ok(-1, 0, 0, 0, 0, 0, 0, 0) &&
         !cell_key_capacity_ok(0, 256, 0, 0, 0, 0, 0, 0) &&
         !cell_key_capacity_ok(0, 0, 65536, 0, 0, 0, 0, 0) &&
         !cell_key_capacity_ok(0, 0, 0, 65536, 0, 0, 0, 0) &&
         !cell_key_capacity_ok(0, 0, 0, 0, 65536, 0, 0, 0) &&
         !cell_key_capacity_ok(0, 0, 5, 5, 5, 6, 0, 0); // ic < i0
}

// The wrappers are exactly the pinned direct u(S, q, e, k) calls
// (mcnux_interactions.hxx draw_map_wrappers_hold precedent).
constexpr bool creation_draw_map_wrappers_hold() noexcept {
  const std::uint64_t S = 1296518744ull;
  const std::uint64_t q = 0x0123456789ABCDEFull;
  const CreationUniforms cu = creation_uniforms(S, q);
  if (!(cu.u_pos_x == u(S, q, 0u, 0u) && cu.u_pos_y == u(S, q, 0u, 1u) &&
        cu.u_pos_z == u(S, q, 0u, 2u) && cu.u_time == u(S, q, 0u, 3u) &&
        cu.u1 == u(S, q, 0u, 4u) && cu.u2 == u(S, q, 0u, 5u)))
    return false;
  const std::uint64_t K = cell_key(3, 5, 17, 21, 33, 10, 10, 20);
  for (int b = 0; b < 2; ++b)
    for (int s = 0; s < NUM_SPECIES; ++s)
      if (cell_count_uniform(S, K, 7u, b, static_cast<Species>(s)) !=
          u(S, K, 7u, static_cast<std::uint32_t>(3 * b + s)))
        return false;
  return true;
}

} // namespace detail

static_assert(draw_k_creation_pos_x == 0u && draw_k_creation_pos_y == 1u &&
                  draw_k_creation_pos_z == 2u && draw_k_creation_time == 3u &&
                  draw_k_creation_costheta == 4u && draw_k_creation_phi == 5u,
              "[MCNX-PKT-05] creation draw-index assignment is pinned: "
              "k = 0..2 position x/y/z, k = 3 time, k = 4 cos theta, "
              "k = 5 phi");
static_assert(detail::count_law_degeneracy_holds(),
              "[MCNX-PKT-02] the count law N_p = g_s sqrt(-g) dV dt eta_b/E_p "
              "must carry g = (1, 1, 4) exactly (nu_x = 4x nu_e at identical "
              "inputs; this is g's single point of entry)");
static_assert(detail::emission_count_table_holds(),
              "[MCNX-PKT-02] emission_count must be floor(N_p) plus one iff "
              "u_draw < N_p - floor(N_p) (strict <, u in [0, 1))");
static_assert(bin_center_energy(10.5, 12.5) == 11.5,
              "[MCNX-PKT-03] packet energy is the bin center "
              "0.5 (E_lo + E_hi), exact");
static_assert(packet_weight(10.0, 2.5) == 4.0,
              "[MCNX-PKT-03] packet weight is N = E_p/nu");
static_assert(bin_integrated_eta(2.5, 4.0) == 10.0,
              "[MCNX-OPA-04] eta_b = eta(s, E_b) * dE_b (bin-integration "
              "bridge from the spectral coefficient interface)");
static_assert(detail::cell_key_fixture_holds(),
              "[MCNX-PKT-05] K_cell packing must match the hand-composed "
              "bit-exact fixture 0x8305000D000B0007");
static_assert(detail::cell_key_flag_bit_holds(),
              "[MCNX-PKT-05]/[MCNX-GPU-04] every cell key has bit 63 set "
              "(disambiguation from packet idcpu values, which never do)");
static_assert(detail::cell_key_capacity_table_holds(),
              "[MCNX-PKT-05] cell-key capacity guard: patch < 128, "
              "level < 256, offsets in [0, 65536), no negatives");
static_assert(detail::creation_draw_map_wrappers_hold(),
              "[MCNX-PKT-05] creation_uniforms/cell_count_uniform must be "
              "exactly the pinned u(S, q, e, k) evaluations (k = 0..5 order; "
              "count draw k = 3b + s)");

} // namespace MCNuX

#endif // MCNUX_EMISSION_HXX
