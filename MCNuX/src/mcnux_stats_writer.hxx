#ifndef MCNUX_STATS_WRITER_HXX
#define MCNUX_STATS_WRITER_HXX

// Shared row-writer scaffold of the stats-* verdict tables ([MCNX-VER-07],
// specs/verification-suite-design.md:255-263: every statistical check writes
// the four fields estimate / expected / sigma / z): a dumb view over the
// four column arrays of one mcnux_stats*_diag grid-array group, with the
// deterministic zero fill and the per-row ZScore mirror that every stats
// writer repeats (factored from MCNuX_StatsEmission, mcnux_emission.cxx;
// second consumer MCNuX_StatsBeam, third consumer MCNuX_StatsScatterbox,
// both mcnux_interactions.cxx).
//
// Deliberately CCTK-free (plain double* columns — CCTK_REAL is binary64 by
// the [MCNX-BLD-03] gate in stub.cxx), so the scaffold carries compile-time
// fixtures like every other shared header. What stays in the callers, by
// design:
//   * the nrows-vs-declared-SIZE cross-check (CCTK_VERROR, group lookup);
//   * the per-row CCTK_VINFO log line (the index -> meaning mapping is
//     owned by, and documented in, each writer).
//
// The z values written here are golden numbers ([MCNX-VER-07]): acceptance
// |z| <= 4 is judged once at capture, thereafter the archived values are
// diffed at 1e-12.

#include "mcnux_stats.hxx" // ZScore

#include <type_traits>

namespace MCNuX {

// View over the four columns of one stats verdict table. The pointers are
// the scheduled routine's CCTK_ARGUMENTS column arrays; nrows is the
// declared group SIZE (cross-checked by the caller before constructing the
// view).
struct StatsDiagView {
  double *estimate;
  double *expected;
  double *sigma;
  double *z;
  int nrows;

  // Deterministic zero fill of all four columns — the defined iteration-0 /
  // pre-compute baseline of every stats writer.
  constexpr void zero_fill() const noexcept {
    for (int r = 0; r < nrows; ++r) {
      estimate[r] = 0.0;
      expected[r] = 0.0;
      sigma[r] = 0.0;
      z[r] = 0.0;
    }
  }

  // Mirror one ZScore (built by the frozen MCNuX::zscore(), never
  // reimplemented) into row r. Pure pass-through: no field is recomputed.
  constexpr void put(int r, const ZScore &zs) const noexcept {
    estimate[r] = zs.estimate;
    expected[r] = zs.expected;
    sigma[r] = zs.sigma;
    z[r] = zs.z;
  }
};

// ---------------------------------------------------------------------------
// Compile-time verification (exact binary-fraction fixtures)
// ---------------------------------------------------------------------------

static_assert(std::is_trivially_copyable<StatsDiagView>::value,
              "StatsDiagView must be a dumb, trivially copyable view");

namespace detail {

// zero_fill clears every slot of every column; put writes exactly the four
// ZScore fields of exactly its row (neighbors untouched). The zscore fixture
// (2.5, 1.0, 0.5) -> z == 3.0 is the exact one of mcnux_stats.hxx.
constexpr bool stats_writer_holds() noexcept {
  double est[3] = {1.0, 2.0, 3.0};
  double exp_[3] = {4.0, 5.0, 6.0};
  double sig[3] = {7.0, 8.0, 9.0};
  double zz[3] = {10.0, 11.0, 12.0};
  const StatsDiagView v{est, exp_, sig, zz, 3};

  v.zero_fill();
  for (int r = 0; r < 3; ++r)
    if (!(est[r] == 0.0 && exp_[r] == 0.0 && sig[r] == 0.0 && zz[r] == 0.0))
      return false;

  v.put(1, zscore(2.5, 1.0, 0.5));
  return est[1] == 2.5 && exp_[1] == 1.0 && sig[1] == 0.5 && zz[1] == 3.0 &&
         est[0] == 0.0 && est[2] == 0.0 && exp_[0] == 0.0 && exp_[2] == 0.0 &&
         sig[0] == 0.0 && sig[2] == 0.0 && zz[0] == 0.0 && zz[2] == 0.0;
}

} // namespace detail

static_assert(detail::stats_writer_holds(),
              "[MCNX-VER-07] the stats-writer scaffold must zero-fill all "
              "four columns and mirror a ZScore into exactly its row");

} // namespace MCNuX

#endif // MCNUX_STATS_WRITER_HXX
