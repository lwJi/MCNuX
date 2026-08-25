// ---------------------------------------------------------------------------
// Rows 134..139 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the emission energy-bin grid (nodes-as-edges
// convention: bin b = [E[b], E[b+1]] on the tabulated linear-MeV nodes,
// dE_b = E[b+1] - E[b], [MCNX-OPA-04]) and the uniform -> physical creation
// mapping (u = 0 -> cell lo / t_n exactly, affine lo + u * width,
// [MCNX-PKT-03]) of specs/packet-representation-and-sampling.md, through the
// mcnux_emission.hxx pure functions. Deterministic pinned fixtures only (the
// `unit-selftest` parfile is a single-shot itlast=0 run): a small non-uniform
// constexpr node array for the grid rows, exact binary constants for the
// mapping rows. Exact tier throughout (add_exact/add_boolean); no
// scheduling, no cell loop, no RNG draws (the mapping consumes already-drawn
// doubles).
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_emission.hxx"

namespace MCNuX {

namespace {

// Pinned non-uniform linear-MeV node fixture (the same shape as the
// energy_bin_grid_holds constexpr fixture in mcnux_emission.hxx): 4 nodes ->
// 3 bins, edges at the nodes.
constexpr int fix_nE = 4;
constexpr double fix_E[fix_nE] = {1.0, 3.0, 8.0, 20.0};

// Pinned mapping fixtures — exact binary values so every product/sum below
// is exact.
constexpr double fix_lo = 3.5;    // interval lower endpoint
constexpr double fix_w = 0.25;    // interval width
constexpr double fix_cell_lo = -2.0, fix_dx = 0.125; // cell bounds shape
constexpr double fix_t_n = 10.0, fix_dt = 0.5;       // step bounds

} // namespace

void append_emission_bins_rows(Battery &b) {
  // --- Row 134: [MCNX-OPA-04] nodes-as-edges: the E_lo/E_hi edges of the
  // middle bin equal the fixture nodes bitwise (and so for the outer bins —
  // no invented outer half-bins). ---
  {
    const EnergyBin b0 = energy_bin(fix_E, fix_nE, 0);
    const EnergyBin b1 = energy_bin(fix_E, fix_nE, 1);
    const EnergyBin b2 = energy_bin(fix_E, fix_nE, 2);
    const bool ok = b1.E_lo == fix_E[1] && b1.E_hi == fix_E[2] &&
                    b0.E_lo == fix_E[0] && b0.E_hi == fix_E[1] &&
                    b2.E_lo == fix_E[2] && b2.E_hi == fix_E[3];
    b.add_boolean("emission.bins.edges_are_nodes", ok);
  }

  // --- Row 135: [MCNX-OPA-04] bin count: N nodes -> N - 1 bins, exact. ---
  b.add_exact("emission.bins.count",
              static_cast<double>(num_energy_bins(fix_nE)), 3.0);

  // --- Row 136: [MCNX-PKT-03] bin-center energy of a grid-produced bin:
  // nu = 0.5 (E[b] + E[b+1]) bitwise (the EnergyBin overload delegates to
  // the pinned scalar formula). ---
  b.add_exact("emission.bins.center",
              bin_center_energy(energy_bin(fix_E, fix_nE, 1)),
              0.5 * (fix_E[1] + fix_E[2]));

  // --- Row 137: [MCNX-OPA-04] width/eta bridge: bin_width feeds
  // bin_integrated_eta consistently, eta_b = eta * (E[b+1] - E[b])
  // bitwise. ---
  b.add_exact("emission.bins.width_eta_bridge",
              bin_integrated_eta(2.5, bin_width(energy_bin(fix_E, fix_nE, 1))),
              2.5 * (fix_E[2] - fix_E[1]));

  // --- Row 138: [MCNX-PKT-03] u = 0 lower endpoints, exact: the scalar
  // workhorse hits lo, creation_position hits the cell's lower corner, and
  // creation_time hits t_n (half-open [t_n, t_n + dt): u = 0 is a reachable
  // exact RNG output). ---
  {
    const bool ok =
        map_uniform_to_interval(0.0, fix_lo, fix_w) == fix_lo &&
        creation_position(0.0, fix_cell_lo, fix_dx) == fix_cell_lo &&
        creation_time(0.0, fix_t_n, fix_dt) == fix_t_n;
    b.add_boolean("emission.map.zero_endpoints", ok);
  }

  // --- Row 139: [MCNX-PKT-03] affine linearity: the u = 0.5 midpoint is
  // lo + 0.5 * width bitwise (exact binary fixture). ---
  b.add_exact("emission.map.linearity",
              map_uniform_to_interval(0.5, fix_lo, fix_w),
              fix_lo + 0.5 * fix_w);
}

} // namespace MCNuX
