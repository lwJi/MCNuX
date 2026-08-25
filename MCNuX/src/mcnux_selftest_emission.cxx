// ---------------------------------------------------------------------------
// Rows 118..127 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the emission count law, per-packet initialization
// quantities, cell-key packing, and creation RNG draw map of
// specs/packet-representation-and-sampling.md [MCNX-PKT-02/03/05] (T14a),
// through the mcnux_emission.hxx pure functions. Deterministic pinned
// fixtures only (the `unit-selftest` parfile is a single-shot itlast=0 run):
// the floor/Bernoulli branches, the exact g = (1, 1, 4) degeneracy factor at
// its single point of entry, the bin-center/weight/bin-integration
// arithmetic, a hand-computed bit-exact K_cell literal plus flag-bit and
// capacity-predicate rows, and the bitwise k = 0..5 creation draw-map trace
// against direct u(S, q, 0, k) calls. No scheduling, no cell loop, no
// containers.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_emission.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_units.hxx"

#include <cstdint>

namespace MCNuX {

namespace {

// Pinned draw-map fixture identity (the mcnux_selftest_interactions.cxx
// precedent: any pinned (S, q) works — the direct-u comparisons are exact).
constexpr std::uint64_t fix_S = 1296518744ull;
constexpr std::uint64_t fix_q = 42ull;

// Pinned count-law fixture (exact binary factors so the products carry no
// rounding surprises across compilers).
constexpr double fix_sqrt_neg_g = 1.25;
constexpr double fix_dV = 0.5;
constexpr double fix_dt = 0.25;
constexpr double fix_eta_b = 44.0;
constexpr double fix_E_p = 2.5;
// N_p(NuE) = 1 * 1.25 * 0.5 * 0.25 * 44.0 / 2.5 = 2.75, exactly.

} // namespace

void append_emission_rows(Battery &b) {
  // --- Row 118: [MCNX-PKT-02] count-law floor/Bernoulli arithmetic at the
  // pinned fixture N_p = 2.75: exact floor, both Bernoulli branches around
  // the remainder 0.75 (u below -> +1, u at/above -> +0; strict <), the
  // u == 0 edge (fires iff frac > 0), and the exactly-integer case. ---
  {
    const double N_p = packet_count(Species::NuE, fix_sqrt_neg_g, fix_dV,
                                    fix_dt, fix_eta_b, fix_E_p);
    const bool ok = N_p == 2.75 &&                  //
                    emission_count(N_p, 0.5) == 3 &&  // u < frac
                    emission_count(N_p, 0.75) == 2 && // u == frac: strict <
                    emission_count(N_p, 0.9) == 2 &&  // u > frac
                    emission_count(N_p, 0.0) == 3 &&  // u == 0, frac > 0
                    emission_count(3.0, 0.0) == 3 &&  // integer N_p
                    emission_count(0.25, 0.0) == 1 && // pure Bernoulli fires
                    emission_count(0.25, 0.25) == 0 && // and does not
                    emission_count(0.0, 0.0) == 0;
    b.add_boolean("emission.count.floor_bernoulli", ok);
  }

  // --- Row 119: [MCNX-PKT-02] the degeneracy g = (1, 1, 4) at its single
  // point of entry: nu_x count is exactly 4x the nu_e count at identical
  // inputs (the exact-tier factor the analytic nu_x benchmark verifies
  // statistically), and nu_e_bar agrees with nu_e bitwise. ---
  {
    const double n_e = packet_count(Species::NuE, fix_sqrt_neg_g, fix_dV,
                                    fix_dt, fix_eta_b, fix_E_p);
    const double n_eb = packet_count(Species::NuEBar, fix_sqrt_neg_g, fix_dV,
                                     fix_dt, fix_eta_b, fix_E_p);
    const double n_x = packet_count(Species::NuX, fix_sqrt_neg_g, fix_dV,
                                    fix_dt, fix_eta_b, fix_E_p);
    b.add_boolean("emission.count.g4_single_entry",
                  n_eb == n_e && n_x == 4.0 * n_e && n_e == 2.75);
  }

  // --- Row 120: [MCNX-PKT-03] bin-center packet energy nu = 0.5 (E_lo +
  // E_hi), exact (baseline energies are never drawn continuously). ---
  b.add_exact("emission.init.bin_center", bin_center_energy(10.5, 12.5), 11.5);

  // --- Row 121: [MCNX-PKT-03] packet weight N = E_p/nu (each packet carries
  // fluid-frame energy E_p), exact at a binary fixture. ---
  b.add_exact("emission.init.weight", packet_weight(fix_E_p, 0.625), 4.0);

  // --- Row 122: [MCNX-OPA-04] the bin-integration bridge eta_b =
  // eta(s, E_b) * dE_b from the spectral coefficient interface, exact. ---
  b.add_exact("emission.eta.bin_integral", bin_integrated_eta(2.5, 4.0), 10.0);

  // --- Row 123: [MCNX-PKT-05] K_cell packing, bit-exact against the
  // hand-computed literal: patch 3, level 5, cell (17, 21, 33) on a level
  // with lower corner (10, 10, 20) -> offsets (7, 11, 13) ->
  // 2^63 + 3*2^56 + 5*2^48 + 13*2^32 + 11*2^16 + 7 = 0x8305000D000B0007. ---
  {
    const std::uint64_t got = cell_key(3, 5, 17, 21, 33, 10, 10, 20);
    b.add_boolean("emission.key.bit_exact",
                  got == 0x8305000D000B0007ull &&
                      cell_key(0, 0, 0, 0, 0, 0, 0, 0) ==
                          0x8000000000000000ull);
  }

  // --- Row 124: [MCNX-PKT-05]/[MCNX-GPU-04] the 2^63 flag bit is set on
  // every packable key, including the all-max corner — cell keys are
  // disambiguated from packet idcpu values, which never have bit 63 set. ---
  {
    const bool ok =
        (cell_key(0, 0, 0, 0, 0, 0, 0, 0) >> 63) == 1u &&
        (cell_key(127, 255, 65535, 65535, 65535, 0, 0, 0) >> 63) == 1u &&
        (cell_key(3, 5, 17, 21, 33, 10, 10, 20) >> 63) == 1u;
    b.add_boolean("emission.key.flag_bit", ok);
  }

  // --- Row 125: [MCNX-PKT-05] capacity guard: the predicate accepts both
  // in-range corners and rejects each single-field overflow (patch 128,
  // level 256, any offset 65536), a negative patch, and a cell below the
  // level's lower corner. ---
  {
    const bool ok = cell_key_capacity_ok(0, 0, 0, 0, 0, 0, 0, 0) &&
                    cell_key_capacity_ok(127, 255, 65535, 65535, 65535, 0, 0,
                                         0) &&
                    !cell_key_capacity_ok(128, 0, 0, 0, 0, 0, 0, 0) &&
                    !cell_key_capacity_ok(-1, 0, 0, 0, 0, 0, 0, 0) &&
                    !cell_key_capacity_ok(0, 256, 0, 0, 0, 0, 0, 0) &&
                    !cell_key_capacity_ok(0, 0, 65536, 0, 0, 0, 0, 0) &&
                    !cell_key_capacity_ok(0, 0, 0, 65536, 0, 0, 0, 0) &&
                    !cell_key_capacity_ok(0, 0, 0, 0, 65536, 0, 0, 0) &&
                    !cell_key_capacity_ok(0, 0, 5, 5, 5, 6, 0, 0);
    b.add_boolean("emission.key.capacity", ok);
  }

  // --- Row 126: [MCNX-PKT-05] creation draw map, bitwise trace: the six
  // creation_uniforms fields must equal the direct u(S, q, 0, k) calls for
  // k = 0..5 in exactly the pinned order (positions x, y, z; time;
  // cos theta; phi), and the pinned constants carry that order. Creation is
  // e = 0 and consumes no other draws (the int.rng.episode_draw_map
  // precedent for e >= 1). ---
  {
    const CreationUniforms cu = creation_uniforms(fix_S, fix_q);
    const double got[6] = {cu.u_pos_x, cu.u_pos_y, cu.u_pos_z,
                           cu.u_time,  cu.u1,      cu.u2};
    bool ok = draw_k_creation_pos_x == 0u && draw_k_creation_pos_y == 1u &&
              draw_k_creation_pos_z == 2u && draw_k_creation_time == 3u &&
              draw_k_creation_costheta == 4u && draw_k_creation_phi == 5u;
    for (std::uint32_t k = 0u; ok && k < 6u; ++k)
      ok = got[k] == u(fix_S, fix_q, 0u, k);
    b.add_boolean("emission.rng.creation_draw_map", ok);
  }

  // --- Row 127: [MCNX-PKT-05] cell-level count draw: q = K_cell, e = n,
  // k = 3b + s, bitwise against direct u calls over bins b = 0, 1 and all
  // three species (k computed from (b, s), not a fixed constant). ---
  {
    const std::uint64_t K = cell_key(3, 5, 17, 21, 33, 10, 10, 20);
    const std::uint32_t n = 7u;
    bool ok = true;
    for (int bin = 0; ok && bin < 2; ++bin)
      for (int s = 0; ok && s < NUM_SPECIES; ++s)
        ok = cell_count_uniform(fix_S, K, n, bin, static_cast<Species>(s)) ==
             u(fix_S, K, n, static_cast<std::uint32_t>(3 * bin + s));
    b.add_boolean("emission.rng.cell_count_draw", ok);
  }
}

} // namespace MCNuX
