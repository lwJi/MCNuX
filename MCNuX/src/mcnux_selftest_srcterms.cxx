// ---------------------------------------------------------------------------
// Rows 60..65 of the MCNuX unit self-test battery (index -> name table in
// mcnux_selftest.cxx) — the [MCNX-HYD-02] source-term ledger arithmetic of
// specs/hydro-coupling-source-terms.md (T17). Pure-function sign fixtures
// through the mcnux_srcterms.hxx helpers on the header's own compile-time
// fixture set (exact binary-fraction inputs, dV dt = 0.5), so every row is
// exact: the worked signs of the spec (nu_e emission -> S_l < 0, G^t < 0;
// nu_e absorption reverses both; elastic scattering changes momentum
// components only) plus the exact fixture values, judged with no tolerance.
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

} // namespace MCNuX
