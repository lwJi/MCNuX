// ---------------------------------------------------------------------------
// Rows 156..161 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the [MCNX-GPU-05] escape-tally arithmetic of
// specs/particle-container-and-gpu.md through the mcnux_escape.hxx pure
// functions: the species-major slot map, the half-open outside_domain
// predicate (interior inside, exactly-on-lower-face inside,
// exactly-on-upper-face escaped, beyond-one-axis escaped), and the
// EscapeTally::accumulate arithmetic on a hand-computed two-species fixture.
// All inputs are exact binary fractions, so every row is exact (no
// tolerance).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_escape.hxx"

namespace MCNuX {

void append_escape_rows(Battery &b) {
  // Row 156 — the species-major slot map: slot = 3 * species + channel over
  // the full NUM_SPECIES x {count, number, energy} grid, plus the exact
  // slot count (the declared SIZE of MCNuX::mcnux_escape_diag).
  bool slots_ok = escape_num_slots == 9;
  for (int s = 0; s < NUM_SPECIES; ++s)
    for (int c = 0; c < escape_num_channels; ++c)
      slots_ok = slots_ok && escape_slot(s, c) == 3 * s + c;
  slots_ok = slots_ok &&
             escape_slot(Species::NuX, escape_channel_energy) ==
                 escape_num_slots - 1;
  b.add_boolean("escape.slot.species_major", slots_ok);

  // Rows 157..160 — the half-open escape predicate on the [-1, 1]^3
  // benchmark domain (exact binary bounds): interior points are inside; a
  // packet exactly ON the lower domain face is inside; exactly ON the upper
  // face it has escaped (the half-open convention that pins WHICH iteration
  // an exact face landing escapes at); beyond a single axis suffices.
  const double lo[3] = {-1.0, -1.0, -1.0};
  const double hi[3] = {1.0, 1.0, 1.0};
  const double interior[3] = {0.125, -0.25, 0.5};
  const double on_lo[3] = {-1.0, 0.25, -0.5};
  const double on_hi[3] = {1.0, 0.25, -0.5};
  const double one_axis[3] = {0.0, 1.25, 0.0};
  b.add_boolean("escape.domain.interior_inside",
                !outside_domain(interior, lo, hi));
  b.add_boolean("escape.domain.lo_face_inside", !outside_domain(on_lo, lo, hi));
  b.add_boolean("escape.domain.hi_face_escaped",
                outside_domain(on_hi, lo, hi));
  b.add_boolean("escape.domain.one_axis_beyond",
                outside_domain(one_axis, lo, hi));

  // Row 161 — EscapeTally::accumulate on the hand fixture: two nu_e
  // escapes (N, p^t) = (1.5, 2.0) and (0.5, 4.0), one nu_x escape
  // (2.0, 0.25); expected slots nu_e {2, 2, 5}, nu_e_bar all zero, nu_x
  // {1, 2, 0.5} — every value an exact binary fraction.
  EscapeTally t{};
  t.accumulate(species_index(Species::NuE), 1.5, 2.0);
  t.accumulate(species_index(Species::NuE), 0.5, 4.0);
  t.accumulate(species_index(Species::NuX), 2.0, 0.25);
  b.add_boolean("escape.tally.accumulate",
                t.v[0] == 2.0 && t.v[1] == 2.0 && t.v[2] == 5.0 &&
                    t.v[3] == 0.0 && t.v[4] == 0.0 && t.v[5] == 0.0 &&
                    t.v[6] == 1.0 && t.v[7] == 2.0 && t.v[8] == 0.5);
}

} // namespace MCNuX
