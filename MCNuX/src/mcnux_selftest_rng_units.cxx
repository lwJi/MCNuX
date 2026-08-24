// ---------------------------------------------------------------------------
// Rows 0..37 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx): the RNG contract of
// specs/rng-and-statistical-acceptance.md [MCNX-RNG-01..06], the unit-system
// contract of specs/conventions-and-units.md [MCNX-CNV-02..04, -06], and the
// binary64 precision gate of specs/build-and-integration.md [MCNX-BLD-03].
// Every check is a call into the shared headers; no constant, fixture, or
// tolerance is restated here.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_particles.hxx"
#include "mcnux_rng.hxx"
#include "mcnux_units.hxx"

#include <cctk.h>

namespace MCNuX {

void append_rng_units_rows(Battery &b) {
  // Rows 0..2 — Philox4x32-10 known-answer vectors, exact word equality.
  // rng-and-statistical-acceptance.md [MCNX-RNG-01], [MCNX-RNG-02].
  const char *const kat_names[detail::num_kat_vectors] = {
      "rng.philox.kat1", "rng.philox.kat2", "rng.philox.kat3"};
  for (int i = 0; i < detail::num_kat_vectors; ++i)
    b.add_boolean(kat_names[i], detail::kat_holds(i));

  // Rows 3..6 — the [0, 1) construction's closed forms, exact.
  // [MCNX-RNG-04].
  const char *const cf_names[detail::num_uniform_closed_forms] = {
      "rng.uniform.u64_zero", "rng.uniform.u64_2p11_minus_1",
      "rng.uniform.u64_2p11", "rng.uniform.u64_max"};
  for (int i = 0; i < detail::num_uniform_closed_forms; ++i) {
    const detail::UniformClosedForm f = detail::uniform_closed_form(i);
    b.add_exact(cf_names[i], uniform_from_u64(f.u64), f.r);
  }

  // Row 7 — u(S, q, e, k) is a pure function of its arguments: repeated,
  // permuted, reversed, and batched evaluation all agree bitwise.
  // [MCNX-RNG-05], [MCNX-RNG-06].
  b.add_boolean("rng.purity", detail::purity_holds());

  // Row 8 — draws k = 2b and k = 2b+1 share Philox block b as word pairs 0
  // and 1, over the counter/key packing of [MCNX-RNG-03].
  b.add_boolean("rng.pairing", detail::pairing_holds());

  // Rows 9..20 — each derived conversion factor, recomputed from the five
  // defining constants, against the spec's pinned table anchor.
  // conventions-and-units.md [MCNX-CNV-03], [MCNX-CNV-04].
  b.add_pinned("units.factor.length_code_to_cgs", length_code_to_cgs,
               pinned::length_code_to_cgs);
  b.add_pinned("units.factor.time_code_to_cgs", time_code_to_cgs,
               pinned::time_code_to_cgs);
  b.add_pinned("units.factor.mass_code_to_cgs", mass_code_to_cgs,
               pinned::mass_code_to_cgs);
  b.add_pinned("units.factor.energy_code_to_cgs", energy_code_to_cgs,
               pinned::energy_code_to_cgs);
  b.add_pinned("units.factor.energy_code_to_MeV", energy_code_to_MeV,
               pinned::energy_code_to_MeV);
  b.add_pinned("units.factor.mass_density_code_to_cgs",
               mass_density_code_to_cgs, pinned::mass_density_code_to_cgs);
  b.add_pinned("units.factor.energy_density_code_to_cgs",
               energy_density_code_to_cgs, pinned::energy_density_code_to_cgs);
  b.add_pinned("units.factor.opacity_code_to_cgs", opacity_code_to_cgs,
               pinned::opacity_code_to_cgs);
  b.add_pinned("units.factor.rate_code_to_cgs", rate_code_to_cgs,
               pinned::rate_code_to_cgs);
  b.add_pinned("units.factor.emissivity_code_to_cgs", emissivity_code_to_cgs,
               pinned::emissivity_code_to_cgs);
  b.add_pinned("units.factor.mev_to_kelvin", mev_to_kelvin,
               pinned::mev_to_kelvin);
  b.add_pinned("units.factor.k_B_MeV_per_K", k_B_MeV_per_K,
               pinned::k_B_MeV_per_K);

  // Rows 21..32 — cgs -> code -> cgs (and back) is the identity to the machine
  // tier for the same factors, probed over several decades. [MCNX-CNV-02].
  b.add_boolean("units.roundtrip.length_code_to_cgs",
                detail::round_trips(length_code_to_cgs));
  b.add_boolean("units.roundtrip.time_code_to_cgs",
                detail::round_trips(time_code_to_cgs));
  b.add_boolean("units.roundtrip.mass_code_to_cgs",
                detail::round_trips(mass_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_code_to_cgs",
                detail::round_trips(energy_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_code_to_MeV",
                detail::round_trips(energy_code_to_MeV));
  b.add_boolean("units.roundtrip.mass_density_code_to_cgs",
                detail::round_trips(mass_density_code_to_cgs));
  b.add_boolean("units.roundtrip.energy_density_code_to_cgs",
                detail::round_trips(energy_density_code_to_cgs));
  b.add_boolean("units.roundtrip.opacity_code_to_cgs",
                detail::round_trips(opacity_code_to_cgs));
  b.add_boolean("units.roundtrip.rate_code_to_cgs",
                detail::round_trips(rate_code_to_cgs));
  b.add_boolean("units.roundtrip.emissivity_code_to_cgs",
                detail::round_trips(emissivity_code_to_cgs));
  b.add_boolean("units.roundtrip.mev_to_kelvin",
                detail::round_trips(mev_to_kelvin));
  b.add_boolean("units.roundtrip.k_B_MeV_per_K",
                detail::round_trips(k_B_MeV_per_K));

  // Row 33 — the species enumeration, degeneracies (including nu_x's g = 4),
  // and lepton numbers, exact. [MCNX-CNV-06].
  b.add_boolean("units.species_table", detail::species_table_holds());

  // Rows 34..35 — the binary64 precision gate of
  // specs/build-and-integration.md [MCNX-BLD-03], measured at runtime rather
  // than only asserted in stub.cxx.
  b.add_exact("precision.sizeof_double", double(sizeof(double)), 8.0);
  b.add_exact("precision.sizeof_cctk_real", double(sizeof(CCTK_REAL)), 8.0);

  // Rows 36..37 — the AMReX legs of the same gate, per the coverage matrix of
  // specs/verification-suite-design.md (both sizeof(amrex::Real) and
  // sizeof(amrex::ParticleReal)). Compile-time asserted in
  // mcnux_particles.hxx; measured here as golden rows.
  b.add_exact("precision.sizeof_amrex_real", double(sizeof(amrex::Real)), 8.0);
  b.add_exact("precision.sizeof_amrex_particle_real",
              double(sizeof(amrex::ParticleReal)), 8.0);
}

} // namespace MCNuX
