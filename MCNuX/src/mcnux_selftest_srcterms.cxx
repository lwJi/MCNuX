// ---------------------------------------------------------------------------
// Rows 60..65 and 128..133 of the MCNuX unit self-test battery (index -> name
// table in mcnux_selftest.cxx) — the [MCNX-HYD-02] source-term ledger
// arithmetic and the [MCNX-HYD-05] conservation-ledger closure audit of
// specs/hydro-coupling-source-terms.md. Pure-function fixtures through the
// mcnux_srcterms.hxx helpers on the header's own compile-time fixture sets
// (exact binary-fraction inputs, dV dt = 0.5), so every row is exact: the
// worked signs of the spec (nu_e emission -> S_l < 0, G^t < 0; nu_e
// absorption reverses both; elastic scattering changes momentum components
// only), the exact fixture values, and the closure audit (net/gross
// accumulation, exact grid-vs-event cancellation, verdict and eps0-floor
// branches), all judged with no tolerance.
// ---------------------------------------------------------------------------

#include "mcnux_selftest.hxx"

#include "mcnux_srcterms.hxx"

namespace MCNuX {

void append_srcterm_rows(Battery &b) {
  constexpr SourceContribution emis = detail::srcfix_emission;
  constexpr SourceContribution absn = detail::srcfix_absorption;
  constexpr SourceContribution scat = detail::srcfix_scattering;

  // Rows 60..61 — nu_e emission: the fluid loses a lepton and energy
  // (S_l < 0, G^t < 0), at the exact fixture values -6 and -15.
  b.add_boolean("srcterm.emission.lepton_sign",
                emis.Sl < 0.0 && emis.Sl == -6.0);
  b.add_boolean("srcterm.emission.energy_sign",
                emis.Gt < 0.0 && emis.Gt == -15.0);

  // Rows 62..63 — nu_e absorption reverses both signs and is the exact
  // negation of emission, component-wise.
  b.add_boolean("srcterm.absorption.lepton_sign",
                absn.Sl > 0.0 && absn.Sl == -emis.Sl);
  b.add_boolean("srcterm.absorption.energy_sign",
                absn.Gt > 0.0 && absn.Gt == -emis.Gt && absn.Gx == -emis.Gx &&
                    absn.Gy == -emis.Gy && absn.Gz == -emis.Gz);

  // Rows 64..65 — elastic scattering changes momentum components only:
  // S_l and G^t are exactly zero, the momentum transfer is the exact,
  // nonzero fixture value.
  b.add_exact("srcterm.scattering.lepton_zero", scat.Sl, 0.0);
  b.add_boolean("srcterm.scattering.momentum_only",
                scat.Gt == 0.0 && scat.Gx == 5.0 && scat.Gy == -2.0 &&
                    scat.Gz == -4.0);
}

// Rows 128..133 — the [MCNX-HYD-05] conservation-ledger closure audit on the
// shared synthetic-event fixture list (synthfix of mcnux_srcterms.hxx, scale
// 1: the exact event list MCNuX_SyntheticDeposit deposits at iteration 1).
// The event side is built from the pre-negation LedgerDelta values; the grid
// side is the sum of the three deposited SourceContribution values times
// dV dt (each event lands in exactly one cell) — the two must cancel
// exactly, and the closure verdict must pass on the consistent fixture,
// fail on a corrupted grid total, and judge the all-zero list against the
// documented eps0 floor. All inputs are exact binary fractions, so every
// row is exact.
void append_ledger_closure_rows(Battery &b) {
  constexpr LedgerAudit audit = synthfix::audit(1.0);
  constexpr SourceContribution e =
      source_from_delta(synthfix::emission(1.0), synthfix::dV, synthfix::dt);
  constexpr SourceContribution a =
      source_from_delta(synthfix::absorption(1.0), synthfix::dV, synthfix::dt);
  constexpr SourceContribution s =
      source_from_delta(synthfix::scattering(1.0), synthfix::dV, synthfix::dt);
  constexpr double dVdt = synthfix::dV * synthfix::dt;

  const double grid[5] = {
      (e.Gt + a.Gt + s.Gt) * dVdt, (e.Gx + a.Gx + s.Gx) * dVdt,
      (e.Gy + a.Gy + s.Gy) * dVdt, (e.Gz + a.Gz + s.Gz) * dVdt,
      (e.Sl + a.Sl + s.Sl) * dVdt};
  const double net[5] = {audit.net.dPt, audit.net.dPx, audit.net.dPy,
                         audit.net.dPz, audit.net.dL};
  const double gross[5] = {audit.gross.dPt, audit.gross.dPx, audit.gross.dPy,
                           audit.gross.dPz, audit.gross.dL};

  // Row 128 — event-side audit net, all five channels at the exact
  // hand-computed fixture values (emission (7.5, 1.5, -3.75, 2.25, 3) +
  // absorption (-4.5, 0.75, -1.5, -0.375, -1.5) + scattering
  // (0, -2.5, 1, 2, 0)).
  b.add_boolean("ledger.audit.net",
                net[0] == 3.0 && net[1] == -0.25 && net[2] == -4.25 &&
                    net[3] == 3.875 && net[4] == 1.5);

  // Row 129 — event-side audit gross activity, exact.
  b.add_boolean("ledger.audit.gross",
                gross[0] == 12.0 && gross[1] == 4.75 && gross[2] == 6.25 &&
                    gross[3] == 4.625 && gross[4] == 4.5);

  // Row 130 — the worked identity Sum_cells G * dV dt = -dX: grid total and
  // event net cancel exactly on every channel (the measured error is the
  // summed absolute residual, identically 0).
  double cancel = 0.0;
  for (int c = 0; c < 5; ++c)
    cancel += detail::cabs(grid[c] + net[c]);
  b.add_exact("ledger.closure.cancellation", cancel, 0.0);

  // Row 131 — the closure verdict passes on the consistent fixture with a
  // residual of exactly 0 on every channel.
  bool all_pass = true;
  for (int c = 0; c < 5; ++c) {
    const LedgerClosure lc = ledger_closure(grid[c], net[c], gross[c]);
    all_pass = all_pass && lc.pass && lc.residual == 0.0;
  }
  b.add_boolean("ledger.closure.pass", all_pass);

  // Row 132 — leak detection: a doubled grid total (the signature of a
  // missing or misordered zero point) fails closure on every channel (all
  // five fixture grid totals are nonzero, far beyond the 1e-13 band).
  bool all_detect = true;
  for (int c = 0; c < 5; ++c)
    all_detect =
        all_detect && !ledger_closure(2.0 * grid[c], net[c], gross[c]).pass;
  b.add_boolean("ledger.closure.detects_leak", all_detect);

  // Row 133 — the all-zero list passes via the documented eps0 floor: zero
  // gross activity flips the denominator to ledger_eps0, and the exactly-zero
  // residual passes the resulting 1e-43 bound.
  const LedgerClosure z = ledger_closure(0.0, 0.0, 0.0);
  b.add_boolean("ledger.closure.allzero_floor",
                z.pass && z.residual == 0.0 &&
                    z.bound == ledger_rtol * ledger_eps0);
}

} // namespace MCNuX
