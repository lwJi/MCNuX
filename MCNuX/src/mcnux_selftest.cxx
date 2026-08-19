// The MCNuX unit self-test battery — the runtime half of the `unit-selftest`
// benchmark of specs/verification-suite-design.md [MCNX-VER-01],
// [MCNX-VER-06].
//
// The battery re-runs, at runtime, the checks that the shared headers already
// assert at compile time, and writes them out as a table of rows so that the
// Cactus regression harness can diff them against committed golden data
// (MCNuX/test/unit-selftest/, produced by CarpetX's TSV writer at
// ABSTOL=RELTOL=1e-12). The spec's binding shape is "one row per named check
// (measured error, pass flag)"; here a row is one index of the 1d grid array
// MCNuX::mcnux_selftest, carrying `selftest_error` and `selftest_pass`.
//
// Why both compile-time and runtime? The static_asserts prove the predicates
// hold for the *compiler's* constant evaluator; these rows additionally pin
// the *measured* errors as golden numbers, so that a change which keeps a
// check passing but moves its error (a re-derived factor, a different
// evaluation order) is still caught at 1e-12.
//
// ---------------------------------------------------------------------------
// Row table (index -> name).  APPEND-ONLY.
// ---------------------------------------------------------------------------
// Golden data is keyed by row index, so rows are never reordered, renamed in
// meaning, or removed: later tasks append their checks at the end and
// regenerate the golden file (and bump SIZE in interface.ccl). The full
// battery the spec eventually requires (tetrad/null identities, Tmunu closed
// forms, nu_x mapping, relabeling identities, sign/range/inversion fixtures)
// is built up this way.
//
//   0..2    rng.philox.kat{1,2,3}                 [MCNX-RNG-01], [MCNX-RNG-02]
//   3..6    rng.uniform.{u64_zero, u64_2p11_minus_1,
//                        u64_2p11, u64_max}       [MCNX-RNG-04]
//   7       rng.purity                            [MCNX-RNG-05], [MCNX-RNG-06]
//   8       rng.pairing                           [MCNX-RNG-03]
//   9..20   units.factor.<name>                   [MCNX-CNV-03], [MCNX-CNV-04]
//  21..32   units.roundtrip.<name>                [MCNX-CNV-02]
//  33       units.species_table                   [MCNX-CNV-06]
//  34..35   precision.sizeof_{double, cctk_real}  [MCNX-BLD-03]
//  36..37   precision.sizeof_amrex_{real,
//                                   particlereal} [MCNX-BLD-03]
//
// Tiers, per specs/README.md and the tolerance discussion in mcnux_units.hxx:
//   * exact   — pass requires a measured error of identically 0 (KATs, the
//               [0, 1) closed forms, purity/pairing, the species table, the
//               precision rows).
//   * machine — the round-trip rows; MCNuX::detail::round_trips judges at
//               detail::rtol_machine = 1e-14 internally, so the row is
//               boolean and its measured error is 0 or 1.
//   * pinned  — the factor-recomputation rows, judged at
//               detail::rtol_pinned = 1e-9. Their measured errors are ~5e-10
//               by construction: the spec's table values carry 10 significant
//               digits, so agreement cannot be asserted at the machine tier
//               (see the rtol_pinned comment in mcnux_units.hxx).

#include "mcnux_rng.hxx"
#include "mcnux_units.hxx"

#include <AMReX_REAL.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>

#include <cstddef>

namespace MCNuX {

namespace {

// One row of the battery.
struct CheckRow {
  const char *name;
  double error; // measured error, in the row's own tier
  bool pass;
};

// Relative error |got/want - 1|, with an absolute fallback where the anchor is
// zero. This is the reported quantity for the `pinned` tier; the boolean tiers
// report 0 (pass) or 1 (fail), and the exact double comparisons report the
// absolute difference (which is 0 exactly when the check passes).
constexpr double rel_error(double got, double want) noexcept {
  const double d = detail::cabs(got - want);
  return want == 0.0 ? d : d / detail::cabs(want);
}

// Ordered, fixed-capacity accumulator; `capacity` only bounds the builder, the
// authoritative length is the number of rows actually appended, which is
// cross-checked against the interface.ccl array size below.
class Battery {
public:
  static constexpr int capacity = 256;

  void add_boolean(const char *name, bool ok) { add(name, ok ? 0.0 : 1.0, ok); }

  // Exact double comparison: no tolerance at all, and the measured error is
  // the absolute difference so a near miss is still visible in the golden data.
  void add_exact(const char *name, double got, double want) {
    add(name, detail::cabs(got - want), got == want);
  }

  // Factor recomputation against a pinned 10-significant-digit table anchor.
  void add_pinned(const char *name, double got, double want) {
    add(name, rel_error(got, want),
        detail::approx_eq(got, want, detail::rtol_pinned));
  }

  int size() const { return n_; }
  const CheckRow &operator[](int i) const { return rows_[i]; }

private:
  void add(const char *name, double error, bool pass) {
    if (n_ >= capacity)
      CCTK_VERROR("MCNuX self-test battery exceeded its builder capacity of "
                  "%d rows while adding \"%s\"",
                  capacity, name);
    rows_[n_].name = name;
    rows_[n_].error = error;
    rows_[n_].pass = pass;
    ++n_;
  }

  CheckRow rows_[capacity] = {};
  int n_ = 0;
};

// ---------------------------------------------------------------------------
// The battery itself. Every check is a call into the shared headers; no
// constant, fixture, or tolerance is restated here.
// ---------------------------------------------------------------------------
void run_battery(Battery &b) {
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

  // Rows 34..37 — the binary64 precision gate of
  // specs/build-and-integration.md [MCNX-BLD-03], measured at runtime rather
  // than only asserted at compile time (rows 34..35 mirror the language-level
  // legs in stub.cxx, rows 36..37 the AMReX legs in mcnux_packets.hxx).
  // amrex::Real and amrex::ParticleReal are gated by independent build macros,
  // so neither row implies the other.
  b.add_exact("precision.sizeof_double", double(sizeof(double)), 8.0);
  b.add_exact("precision.sizeof_cctk_real", double(sizeof(CCTK_REAL)), 8.0);
  b.add_exact("precision.sizeof_amrex_real", double(sizeof(amrex::Real)), 8.0);
  b.add_exact("precision.sizeof_amrex_particlereal",
              double(sizeof(amrex::ParticleReal)), 8.0);
}

// Size of the MCNuX::mcnux_selftest array as declared in interface.ccl.
int selftest_array_size() {
  const int gi = CCTK_GroupIndex("MCNuX::mcnux_selftest");
  if (gi < 0)
    CCTK_VERROR("MCNuX::mcnux_selftest is not a known group");
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("MCNuX::mcnux_selftest must be a 1-dimensional grid array");
  return int(*sizes[0]);
}

} // namespace

extern "C" void MCNuX_SelfTest(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_SelfTest;
  DECLARE_CCTK_PARAMETERS;

  Battery battery;
  run_battery(battery);

  const int nrows = battery.size();
  const int nslots = selftest_array_size();
  if (nrows != nslots)
    CCTK_VERROR("MCNuX self-test battery has %d rows, but "
                "MCNuX::mcnux_selftest is declared with SIZE=%d in "
                "MCNuX/interface.ccl; the two must agree (and the golden data "
                "in MCNuX/test/unit-selftest/ must be regenerated whenever "
                "rows are appended)",
                nrows, nslots);

  int nfailed = 0;
  for (int i = 0; i < nrows; ++i) {
    const CheckRow &row = battery[i];
    selftest_error[i] = row.error;
    selftest_pass[i] = row.pass ? 1.0 : 0.0;

    // The TSV table is numeric only (it is diffed column-wise by the harness),
    // so the index -> name mapping is echoed to the log instead.
    CCTK_VINFO("MCNuX selftest row %3d  %-40s  error = %.17g  %s", i, row.name,
               row.error, row.pass ? "PASS" : "FAIL");
    if (!row.pass)
      ++nfailed;
  }

  if (nfailed > 0)
    CCTK_VWARN(CCTK_WARN_ALERT,
               "MCNuX unit self-test battery: %d of %d checks FAILED (see the "
               "per-row lines above; the golden diff of the `unit-selftest` "
               "benchmark is the scoring mechanism)",
               nfailed, nrows);
  else
    CCTK_VINFO("MCNuX unit self-test battery: all %d checks passed", nrows);
}

} // namespace MCNuX
