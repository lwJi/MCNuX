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

// TODO: Don't include files from other thorns; create a proper interface
//
// The same relative include as mcnux_cadence.cxx (see the rationale there):
// CarpetX exposes no public capability header for `ghext`, which the
// conservation-ledger closure diagnostic below walks for its grid-side
// reduction ([MCNX-CTX-01]).
#include "../../../CarpetX/CarpetX/src/driver.hxx"

#include "mcnux_srcterms.hxx"
#include "mcnux_units.hxx"

#include <loop_device.hxx>

#include <AMReX_FabArrayUtility.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>

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

  // Iteration scale: packets-worth N grows linearly with cctk_iteration, so
  // consecutive steps deposit different, analytically-known values.
  const double scale = static_cast<double>(cctk_iteration);

  // Fixture cells (see the FixtureCell comment above for the coordinates).
  constexpr FixtureCell emission_cell{-0.375, 0.125, 0.125};
  constexpr FixtureCell absorption_cell{0.125, 0.375, 0.125};
  constexpr FixtureCell scattering_cell{0.125, 0.125, 0.625};

  // The three sign fixtures of hydro-coupling-source-terms.md:260-262 — the
  // shared synthfix event list of mcnux_srcterms.hxx (one source of truth
  // with the [MCNX-HYD-05] closure audit and the selftest rows), evaluated on
  // the host with the fixed synthetic denominators synthfix::dV, synthfix::dt
  // (dV dt = 0.5 — NOT the grid's dV or the run's dt; the synthetic list
  // proves the ledger arithmetic, not a physical deposition):
  //   nu_e emission   -> S_l < 0 and G^t < 0 in exactly the emission cell;
  //   nu_e absorption -> both signs reversed in the absorption cell;
  //   elastic scattering -> momentum components only in the scattering cell.
  const SourceContribution emis = source_from_delta(
      synthfix::emission(scale), synthfix::dV, synthfix::dt);
  const SourceContribution absn = source_from_delta(
      synthfix::absorption(scale), synthfix::dV, synthfix::dt);
  const SourceContribution scat = source_from_delta(
      synthfix::scattering(scale), synthfix::dV, synthfix::dt);

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

// ---------------------------------------------------------------------------
// Conservation-ledger closure diagnostic  (MCNuX::mcnux_ledger_diag)
// ---------------------------------------------------------------------------

namespace {

// Sum one component of one group's MultiFab over every patch and level,
// valid cells only (nghost = 0 — ghost cells would double-count the shared
// faces). Rank-local sum: every MCNuX test runs NPROCS 1 (the fill_packet_diag
// single-rank note of mcnux_geodesic.cxx); a multi-rank closure would need a
// reduction over ranks. Single-level assumption likewise: the source-zero-add
// benchmark is unigrid — with mesh refinement, coarse cells covered by finer
// levels would double-count and the per-cell dV would differ per level.
double sum_group_component(const int gi, const int comp) {
  if (!CarpetX::ghext)
    CCTK_VERROR("MCNuX reached the ledger-closure diagnostic, but the CarpetX "
                "driver singleton `ghext` is null "
                "(specs/carpetx-thorn-integration.md [MCNX-CTX-01])");

  constexpr int tl = 0;
  double total = 0.0;
  for (const auto &patchdata : CarpetX::ghext->patchdata) {
    for (const auto &leveldata : patchdata.leveldata) {
      const amrex::MultiFab &mf = *leveldata.groupdata.at(gi)->mfab.at(tl);
      total += double(amrex::ReduceSum(
          mf, amrex::IntVect(0),
          [=] AMREX_GPU_DEVICE(const amrex::Box &bx,
                               const amrex::Array4<const amrex::Real> &arr)
              -> amrex::Real {
            amrex::Real s = 0.0;
            amrex::Loop(bx,
                        [&](int i, int j, int k) { s += arr(i, j, k, comp); });
            return s;
          }));
    }
  }
  return total;
}

} // namespace

// The conservation-ledger closure diagnostic of [MCNX-HYD-05]
// (hydro-coupling-source-terms.md:177-188), parameter-gated by
// MCNuX::test_ledger_closure and scheduled after all contributors of
// MCNuX_AddToSourceTerms have run: per channel X in {P^t, P_x, P_y, P_z, L},
//
//   | Sum_cells (source value) * dV * dt  +  dX_packets |
//       <=  ledger_rtol * max(Sum_events |dX|, ledger_eps0)
//
// with the grid side reduced from the deposited MCNuX::rad_force /
// MCNuX::lepton_source variables and the event side independently rebuilt
// from the shared synthfix LedgerDelta list (pre-negation packet deltas —
// never the deposited SourceContribution, which would trivially self-cancel).
// The multiply uses the fixture's synthfix::dV * synthfix::dt = 0.5 (the
// deposit's own normalization), not the grid's cell size or the run's time
// step. The verdict table is mirrored into MCNuX::mcnux_ledger_diag (one row
// per channel) for golden output; any failing channel additionally aborts
// the run so a regression is loud even without a golden diff.
//
// The AT initial leg audits the all-zero state (zeroed variables, scale-0
// event list): both sides are exactly zero and the verdict exercises the
// documented eps0 floor.
extern "C" void MCNuX_LedgerClosure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_LedgerClosure;
  DECLARE_CCTK_PARAMETERS;

  // Event side: the synthetic list's LedgerDelta values at this iteration's
  // scale (MCNuX_SyntheticDeposit above deposits exactly this list with
  // scale = cctk_iteration).
  const double scale = static_cast<double>(cctk_iteration);
  LedgerAudit audit{};
  audit.accumulate(synthfix::emission(scale));
  audit.accumulate(synthfix::absorption(scale));
  audit.accumulate(synthfix::scattering(scale));

  // Grid side: Sum_cells (source value) * dV * dt per channel, with the
  // fixture's dV dt.
  const int gi_rad = CCTK_GroupIndex("MCNuX::rad_force");
  const int gi_lep = CCTK_GroupIndex("MCNuX::lepton_source");
  if (gi_rad < 0 || gi_lep < 0)
    CCTK_VERROR("MCNuX ledger closure needs the MCNuX::rad_force and "
                "MCNuX::lepton_source groups");

  const double dVdt = synthfix::dV * synthfix::dt;
  const double grid_total[5] = {
      sum_group_component(gi_rad, 0) * dVdt, // rf_t    -> P^t
      sum_group_component(gi_rad, 1) * dVdt, // rf_x    -> P_x
      sum_group_component(gi_rad, 2) * dVdt, // rf_y    -> P_y
      sum_group_component(gi_rad, 3) * dVdt, // rf_z    -> P_z
      sum_group_component(gi_lep, 0) * dVdt, // lep_src -> L
  };
  const double net[5] = {audit.net.dPt, audit.net.dPx, audit.net.dPy,
                         audit.net.dPz, audit.net.dL};
  const double gross[5] = {audit.gross.dPt, audit.gross.dPx, audit.gross.dPy,
                           audit.gross.dPz, audit.gross.dL};

  constexpr const char *channel[5] = {"P^t", "P_x", "P_y", "P_z", "L"};
  int nfailed = 0;
  for (int c = 0; c < 5; ++c) {
    const LedgerClosure lc = ledger_closure(grid_total[c], net[c], gross[c]);
    lg_grid_total[c] = grid_total[c];
    lg_event_net[c] = net[c];
    lg_gross[c] = gross[c];
    lg_residual[c] = lc.residual;
    lg_pass[c] = lc.pass ? 1.0 : 0.0;
    CCTK_VINFO("MCNuX ledger closure %-3s: grid = %.17g, net = %.17g, "
               "gross = %.17g, residual = %.17g  %s",
               channel[c], grid_total[c], net[c], gross[c], lc.residual,
               lc.pass ? "PASS" : "FAIL");
    if (!lc.pass)
      ++nfailed;
  }

  if (nfailed > 0)
    CCTK_VERROR("MCNuX conservation-ledger closure FAILED on %d of 5 channels "
                "at iteration %d ([MCNX-HYD-05], "
                "specs/hydro-coupling-source-terms.md:177-188; see the "
                "per-channel lines above)",
                nfailed, cctk_iteration);
}

} // namespace MCNuX
