// Source-term grid-variable routines of specs/hydro-coupling-source-terms.md:
// the zero point of the zero-then-add protocol [MCNX-HYD-03] and the
// parameter-gated synthetic-event contributor proving the [MCNX-HYD-02]
// ledger path (the `source-zero-add` benchmark).
//
// Both routines walk the cell-centered (CENTERING={ccc}) groups
// MCNuX::rad_force and MCNuX::lepton_source with the CarpetX loop machinery
// (<1,1,1> loop templates — the HydroBaseX idiom for ccc variables), over the
// `everywhere` region matching their WRITES: clauses in MCNuX/schedule.ccl.
// The ledger arithmetic itself lives in mcnux_srcterms.hxx and is evaluated
// on the host before the loops; the device lambdas only accumulate the
// precomputed contributions ([MCNX-RNG-08] governs real multi-contributor
// deposition — the synthetic list writes one contribution per cell inside a
// grid loop, so plain += at the point is deterministic).

#include "mcnux_srcterms.hxx"
#include "mcnux_units.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace MCNuX {
using namespace Loop;

// The zero point of [MCNX-HYD-03]: one routine body serving both schedule
// bindings (AT initial for the variables' initial data, and at the start of
// MCNuX_TransportStep each step). Every source variable is set to zero over
// the `everywhere` region; every write after this point in the step is a
// contributor's += (no MCNuX routine assigns after the zero point).
extern "C" void MCNuX_ZeroSourceTerms(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MCNuX_ZeroSourceTerms;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      rf_t(p.I) = 0;
                                      rf_x(p.I) = 0;
                                      rf_y(p.I) = 0;
                                      rf_z(p.I) = 0;
                                      lep_src(p.I) = 0;
                                    });
}

namespace {

// Fixture cells of the synthetic event list, identified by their cell-center
// coordinates on the source-zero-add benchmark grid (8^3 cells over the
// CarpetX default domain [-1, 1]^3, dx = 0.25). Each lies on one of CarpetX's
// TSV output lines (transverse index 4, coordinate 0.125), so every deposited
// value is golden-checked:
//   emission   at cell (2,4,4), center (-0.375, 0.125, 0.125) — on the x line
//   absorption at cell (4,5,4), center ( 0.125, 0.375, 0.125) — on the y line
//   scattering at cell (4,4,6), center ( 0.125, 0.125, 0.625) — on the z line
// The cells are matched by coordinate with a quarter-cell tolerance (cell
// centers are a full dx apart), so the match is exact-by-construction and
// independent of index conventions.
// The cell objects themselves are defined as function-local constexpr values
// inside MCNuX_SyntheticDeposit (captured by value into the device lambda):
// nvcc only lets device code use namespace-scope constexpr variables of
// scalar type, so struct-typed fixtures at namespace scope are "undefined in
// device code" under CUDA (HIP/clang accepts them, which hides the defect).
struct FixtureCell {
  double x, y, z;
};

} // namespace

// The synthetic-event contributor ([MCNX-HYD-03]'s contributor idiom,
// parameter-gated by MCNuX::test_synthetic_deposit): a deterministic,
// iteration-scaled event list — the three [MCNX-HYD-02] sign fixtures, one
// per fixture cell — pushed through the shared ledger arithmetic with fixed
// synthetic (dP^t, dP_i, dL, dV, dt), so every grid value is analytically
// predictable: value(iteration n) = n * value(iteration 1). A missing or
// misordered zero point would double-count across steps and break exactly
// this linearity, which the golden diff catches at 1e-12.
extern "C" void MCNuX_SyntheticDeposit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MCNuX_SyntheticDeposit;
  DECLARE_CCTK_PARAMETERS;

  // Fixed synthetic denominators (exact binary fractions; dV dt = 0.5).
  // These are NOT the grid's dV or the run's dt — the synthetic list proves
  // the ledger arithmetic, not a physical deposition.
  constexpr double dV = 2.0;
  constexpr double dt = 0.25;

  // Iteration scale: packets-worth N grows linearly with cctk_iteration, so
  // consecutive steps deposit different, analytically-known values.
  const double scale = static_cast<double>(cctk_iteration);

  // Fixture cells (see the FixtureCell comment above for the coordinates).
  constexpr FixtureCell emission_cell{-0.375, 0.125, 0.125};
  constexpr FixtureCell absorption_cell{0.125, 0.375, 0.125};
  constexpr FixtureCell scattering_cell{0.125, 0.125, 0.625};

  // The three sign fixtures of hydro-coupling-source-terms.md:260-262,
  // evaluated on the host through mcnux_srcterms.hxx:
  //   nu_e emission   -> S_l < 0 and G^t < 0 in exactly the emission cell;
  //   nu_e absorption -> both signs reversed in the absorption cell;
  //   elastic scattering -> momentum components only in the scattering cell.
  const SourceContribution emis = source_from_delta(
      emission_delta(Species::NuE, 3.0 * scale, 2.5, 0.5, -1.25, 0.75), dV,
      dt);
  const SourceContribution absn = source_from_delta(
      absorption_delta(Species::NuE, 1.5 * scale, 3.0, -0.5, 1.0, 0.25), dV,
      dt);
  const SourceContribution scat = source_from_delta(
      scattering_delta(2.0 * scale, 1.5, 1.0, 0.25, -0.5, 1.5, -0.25, 0.75,
                       0.5),
      dV, dt);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto is_cell = [&](const FixtureCell &c) {
          using std::fabs;
          return fabs(p.x - c.x) < 0.25 * p.dx &&
                 fabs(p.y - c.y) < 0.25 * p.dy && fabs(p.z - c.z) < 0.25 * p.dz;
        };
        const auto accumulate = [&](const SourceContribution &s) {
          // Contributor idiom: read-modify-write += only ([MCNX-HYD-03]).
          rf_t(p.I) += s.Gt;
          rf_x(p.I) += s.Gx;
          rf_y(p.I) += s.Gy;
          rf_z(p.I) += s.Gz;
          lep_src(p.I) += s.Sl;
        };
        if (is_cell(emission_cell))
          accumulate(emis);
        else if (is_cell(absorption_cell))
          accumulate(absn);
        else if (is_cell(scattering_cell))
          accumulate(scat);
      });
}

} // namespace MCNuX
