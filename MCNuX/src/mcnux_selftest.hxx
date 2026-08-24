// Shared machinery of the MCNuX unit self-test battery
// (specs/verification-suite-design.md [MCNX-VER-01], [MCNX-VER-06]).
//
// The battery is split across translation units so that appending checks
// touches only a small file: this header carries the row accumulator and the
// declarations of the per-domain row appenders; each appender lives in its own
// mcnux_selftest_<domain>.cxx; and mcnux_selftest.cxx owns the authoritative
// APPEND-ONLY index -> name row table, the ordered appender call sequence
// (which IS the row ordering), and the scheduled routine. Unlike the other
// mcnux_*.hxx headers this one carries no compile-time selftests — it is
// host-only test scaffolding — but it is still #include'd from stub.cxx so
// the "every shared header aggregates through stub.cxx" rule stays uniform.

#ifndef MCNUX_SELFTEST_HXX
#define MCNUX_SELFTEST_HXX

#include "mcnux_units.hxx"

#include <cctk.h>

#include <cmath>
#include <cstddef>

namespace MCNuX {

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
constexpr double selftest_rel_error(double got, double want) noexcept {
  const double d = detail::cabs(got - want);
  return want == 0.0 ? d : d / detail::cabs(want);
}

// Ordered, fixed-capacity accumulator; `capacity` only bounds the builder, the
// authoritative length is the number of rows actually appended, which is
// cross-checked against the interface.ccl array size in mcnux_selftest.cxx.
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
    add(name, selftest_rel_error(got, want),
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

// A bounded, finite, smooth log-stored synthetic table value keyed on the
// flat column-major offset (the WeakLibInterp
// test_parallelfor_wrappers.cpp generator) — the shared fixture generator of
// the opacity, table-coefficient, and table-range appenders.
inline double synth_tval(std::size_t k) {
  return 0.35 + 0.4 * std::sin(0.6 * double(k) + 0.2);
}

// The per-domain row appenders, one per mcnux_selftest_<domain>.cxx, in the
// battery's row order. The ordered call sequence in mcnux_selftest.cxx's
// run_battery is authoritative; the row ranges here are a courtesy gloss of
// the index -> name table in that file.
void append_rng_units_rows(Battery &b);      // rows 0..37
void append_opacity_rows(Battery &b);        // rows 38..52
void append_tetrad_rows(Battery &b);         // rows 53..59
void append_srcterm_rows(Battery &b);        // rows 60..65
void append_trp_rows(Battery &b);            // rows 66..70
void append_coefficient_rows(Battery &b);    // rows 71..78
void append_geodesic_rows(Battery &b);       // rows 79..85
void append_baseline_table_rows(Battery &b); // rows 86..91
void append_table_range_rows(Battery &b);    // rows 92..99
void append_stats_rows(Battery &b);          // rows 100..105

} // namespace MCNuX

#endif // MCNUX_SELFTEST_HXX
